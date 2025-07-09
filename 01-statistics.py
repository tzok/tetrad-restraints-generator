#!/usr/bin/env python3
"""
01-statistics.py

Skeleton script to load every *.json file in a given directory and return
a mapping from file name → dictionary with the JSON data.

Example
-------
$ python 01-statistics.py data/
"""

from __future__ import annotations

import argparse
import json
import gzip
from pathlib import Path
from typing import Dict, Any, List, Set

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from rnapolis.parser_v2 import parse_cif_atoms
from rnapolis.tertiary_v2 import Structure, Residue, calculate_torsion_angle


def build_nucleotide_map(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    """
    Build an index mapping a nucleotide's ``fullName`` → the nucleotide object.

    Parameters
    ----------
    data : Dict[str, Any]
        JSON dictionary that may contain the ``nucleotides`` list.

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Mapping from ``fullName`` to the full nucleotide entry.
        Any nucleotide without a ``fullName`` string is ignored.
    """
    index: Dict[str, Dict[str, Any]] = {}
    for nucleotide in data.get("nucleotides", []):
        full_name = nucleotide.get("fullName")
        if isinstance(full_name, str):
            index[full_name] = nucleotide
    return index


def extract_g_tetrads(
    data: Dict[str, Any], nucleotide_map: Dict[str, Dict[str, Any]]
) -> list[Dict[str, Any]]:
    """
    Extract tetrads that satisfy two criteria:

    1. Every nucleotide referenced by ``nt1``‒``nt4`` is a guanine
       (``shortName`` == ``"G"``).
    2. For every unordered pair of nucleotides inside the tetrad there exists
       a record in the global ``basePairs`` list whose ``lw`` annotation is
       either ``"cWH"`` or ``"cHW"``:

       • If ``lw == "cWH"`` we collect the pair exactly as stored
         ``(base_pair["nt1"], base_pair["nt2"])``.
       • If ``lw == "cHW"`` we collect the reversed order
         ``(base_pair["nt2"], base_pair["nt1"])``.

       Should any pair be missing or use another ``lw`` annotation the whole
       tetrad is discarded.

    Returns
    -------
    list[Dict[str, Any]]
        Each item is a dictionary with two keys:

        - ``"tetrad"`` – the original tetrad object
        - ``"pairs"``  – list with the four collected nucleotide pairs ``(nt1‒nt2, nt2‒nt3, nt3‒nt4, nt4‒nt1)``

    Parameters
    ----------
    data : Dict[str, Any]
        Original JSON content.
    nucleotide_map : Dict[str, Dict[str, Any]]
        Mapping created by :pyfunc:`build_nucleotide_map`.

    Returns
    -------
    list[Dict[str, Any]]
        List of tetrad dictionaries that satisfy the condition.
    """
    g_tetrads: list[Dict[str, Any]] = []

    # Index basePairs for quick look-up
    base_pair_lookup = {
        (bp.get("nt1"), bp.get("nt2")): bp for bp in data.get("basePairs", [])
    }

    for helix in data.get("helices", []):
        for quad in helix.get("quadruplexes", []):
            for tetrad in quad.get("tetrads", []):
                nts = [tetrad.get(f"nt{i}") for i in range(1, 5)]

                # --- condition 1: every nucleotide is a guanine ----------------
                if not all(
                    isinstance(nt, str)
                    and nt in nucleotide_map
                    and nucleotide_map[nt].get("shortName") == "G"
                    for nt in nts
                ):
                    continue

                # --- condition 2: validate all four base-pairs -----------------
                collected_pairs: list[tuple[str, str]] = []
                valid = True
                for i, j in ((0, 1), (1, 2), (2, 3), (3, 0)):
                    a, b = nts[i], nts[j]

                    bp = base_pair_lookup.get((a, b))
                    if bp is None:
                        bp = base_pair_lookup.get((b, a))

                    if bp is None:
                        valid = False
                        break

                    lw = bp.get("lw")
                    if lw == "cWH":
                        collected_pairs.append((bp["nt1"], bp["nt2"]))
                    elif lw == "cHW":
                        collected_pairs.append((bp["nt2"], bp["nt1"]))
                    else:
                        valid = False
                        break
                if not valid:
                    break

                if valid and len(collected_pairs) == 4:
                    g_tetrads.append({"tetrad": tetrad, "pairs": collected_pairs})

    return g_tetrads


def load_json_files(directory: Path) -> Dict[str, Dict[str, Any]]:
    """
    Load every ``*.json`` file in *directory* that contains at least one helix
    with a non-empty ``quadruplexes`` list.

    For every kept JSON file also build a mapping from ``fullName`` of each
    nucleotide (found in the ``nucleotides`` list) to the full nucleotide
    object.  The returned dictionary therefore maps **file name →** a new
    dictionary with two keys:

    - ``"data"`` – the raw JSON loaded from disk
    - ``"nucleotide_map"`` – ``Dict[str, Dict[str, Any]]`` created from the
      nucleotide list
    - ``"g_tetrads"`` – ``List[Dict[str, Any]]`` where each entry contains a
      G-tetrad and its collected base-pairs

    Parameters
    ----------
    directory : Path
        Folder containing JSON files.

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Mapping **file_name → loaded JSON dict**.
        If no JSON files are found the mapping will be empty.
    """
    if not directory.is_dir():
        raise NotADirectoryError(f"{directory} is not a directory")

    json_map: Dict[str, Dict[str, Any]] = {}
    for json_file in directory.glob("*.json"):
        try:
            with json_file.open("r", encoding="utf-8") as fp:
                data = json.load(fp)

            # keep only files that have at least one helix whose
            # ``quadruplexes`` field is a non-empty list
            helices = data.get("helices", [])
            if any(
                isinstance(helix.get("quadruplexes"), list)
                and len(helix["quadruplexes"]) > 0
                for helix in helices
            ):
                nucleotide_map = build_nucleotide_map(data)
                json_map[json_file.name] = {
                    "data": data,
                    "nucleotide_map": nucleotide_map,
                    "g_tetrads": extract_g_tetrads(data, nucleotide_map),
                }
        except json.JSONDecodeError as exc:
            # In a full implementation you might want to log or collect errors.
            print(f"Skipping {json_file} – invalid JSON: {exc}")

    return json_map


def collect_tetrad_residues(
    cif_path: Path, g_tetrads: List[Dict[str, Any]]
) -> List[Residue]:
    """
    Given a mmCIF.gz path and a list of G-tetrads extracted from JSON,
    return the unique residues participating in those tetrads.

    Parameters
    ----------
    cif_path : Path
        Path to the ``*.cif.gz`` structure file.
    g_tetrads : list[Dict[str, Any]]
        Output of :pyfunc:`extract_g_tetrads`.

    Returns
    -------
    List[rnapolis.tertiary_v2.Residue]
        Unique residues corresponding to the nt identifiers in the tetrads.
        If the structure file is missing or unreadable an empty list is
        returned.
    """
    if not cif_path.is_file():
        print(f"Structure file not found: {cif_path}")
        return []

    try:
        with gzip.open(cif_path, "rt", encoding="utf-8") as handle:
            atoms_df = parse_cif_atoms(handle)
    except Exception as exc:  # noqa: BLE001
        print(f"Could not parse {cif_path}: {exc}")
        return []

    structure = Structure(atoms_df)
    residue_lookup = {str(residue): residue for residue in structure.residues}

    residues: Set[Residue] = set()
    for tetrad_entry in g_tetrads:
        tetrad = tetrad_entry["tetrad"]
        for i in range(1, 5):
            nt_id = tetrad.get(f"nt{i}")
            if isinstance(nt_id, str) and nt_id in residue_lookup:
                residues.add(residue_lookup[nt_id])

    return list(residues)


def load_cache(cache_file: Path) -> Dict[str, Any]:
    """
    Load cache from *cache_file* if it exists, otherwise return empty dict.
    """
    if cache_file.is_file():
        try:
            with cache_file.open("r", encoding="utf-8") as fp:
                return json.load(fp)
        except Exception:  # noqa: BLE001
            pass
    return {}


def save_cache(cache_file: Path, cache_data: Dict[str, Any]) -> None:
    """Write *cache_data* into *cache_file* as JSON."""
    try:
        with cache_file.open("w", encoding="utf-8") as fp:
            json.dump(cache_data, fp)
    except Exception as exc:  # noqa: BLE001
        print(f"Could not write cache to {cache_file}: {exc}")


def analyse_tetrad_geometry(
    structure: Structure,
    residue_lookup: Dict[str, Residue],
    g_tetrads: List[Dict[str, Any]],
) -> tuple[list[float], list[float], list[float], list[float]]:
    """
    For the provided G-tetrads calculate geometric properties:

    1. N1–O6 inter-atomic distances for every base pair.
    2. N2–N7 inter-atomic distances for every base pair.
    3. Torsion angle between N9 atoms (ordered nt1→nt2→nt3→nt4).
    4. Torsion angle between O6 atoms (ordered nt1→nt2→nt3→nt4).

    Returns four lists with the collected values.
    """
    dist_n1_o6: list[float] = []
    dist_n2_n7: list[float] = []
    tors_n9: list[float] = []
    tors_o6: list[float] = []

    for tetrad_entry in g_tetrads:
        tetrad = tetrad_entry["tetrad"]
        nts_order = [tetrad.get(f"nt{i}") for i in range(1, 5)]

        # ----- torsion angles -------------------------------------------------
        residues_ordered = []
        for nt in nts_order:
            res = residue_lookup.get(nt)
            if res is None:
                break
            residues_ordered.append(res)
        if len(residues_ordered) == 4:
            atoms_n9 = [r.find_atom("N9") for r in residues_ordered]
            atoms_o6 = [r.find_atom("O6") for r in residues_ordered]
            if all(atoms_n9):
                tors_n9.append(
                    calculate_torsion_angle(
                        *(atom.coordinates for atom in atoms_n9)  # type: ignore[arg-type]
                    )
                )
            if all(atoms_o6):
                tors_o6.append(
                    calculate_torsion_angle(
                        *(atom.coordinates for atom in atoms_o6)  # type: ignore[arg-type]
                    )
                )

        # ----- distances for each base pair ----------------------------------
        for nt_a, nt_b in tetrad_entry["pairs"]:
            res_a = residue_lookup.get(nt_a)
            res_b = residue_lookup.get(nt_b)
            if res_a is None or res_b is None:
                continue

            atom_n1 = res_a.find_atom("N1")
            atom_o6 = res_b.find_atom("O6")
            if atom_n1 and atom_o6:
                dist_n1_o6.append(
                    np.linalg.norm(atom_n1.coordinates - atom_o6.coordinates).item()
                )

            atom_n2 = res_a.find_atom("N2")
            atom_n7 = res_b.find_atom("N7")
            if atom_n2 and atom_n7:
                dist_n2_n7.append(
                    np.linalg.norm(atom_n2.coordinates - atom_n7.coordinates).item()
                )

    return dist_n1_o6, dist_n2_n7, tors_n9, tors_o6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Process all JSON files in a directory and analyse G-tetrads "
            "using corresponding structure files."
        )
    )
    parser.add_argument(
        "json_directory",
        type=Path,
        help="Directory containing .json files",
    )
    parser.add_argument(
        "cif_directory",
        type=Path,
        help="Directory containing .cif.gz files with the same base names",
    )
    parser.add_argument(
        "--cache-file",
        type=Path,
        default=Path(".statistics_cache.json"),
        help="Path to cache file for geometry results (default: .statistics_cache.json)",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("fitted_params.json"),
        help="Path where fitted distribution parameters will be stored",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mapping = load_json_files(args.json_directory)

    # Accumulators for global histograms
    all_n1_o6: list[float] = []
    all_n2_n7: list[float] = []
    all_tors_n9: list[float] = []
    all_tors_o6: list[float] = []
    print(f"Loaded {len(mapping)} JSON files:")
    cache = load_cache(args.cache_file)

    for name, entry in mapping.items():
        print(
            f"  {name}: "
            f"{len(entry['data'].get('nucleotides', []))} nucleotides, "
            f"{len(entry['nucleotide_map'])} indexed by fullName, "
            f"{len(entry['g_tetrads'])} G-tetrads"
        )

        # -------------------- use cached geometry if available ---------------
        if name in cache:
            cached = cache[name]
            all_n1_o6.extend(cached.get("n1_o6", []))
            all_n2_n7.extend(cached.get("n2_n7", []))
            all_tors_n9.extend(cached.get("tors_n9", []))
            all_tors_o6.extend(cached.get("tors_o6", []))
            print("    Loaded geometry from cache")
            continue

        # -------------------- compute geometry --------------------------------
        if entry["g_tetrads"]:
            base_name = Path(name).stem
            cif_path = args.cif_directory / f"{base_name}.cif.gz"
            residues = collect_tetrad_residues(cif_path, entry["g_tetrads"])
            if residues:
                structure = Structure(parse_cif_atoms(gzip.open(cif_path, "rt")))
                residue_lookup = {str(r): r for r in structure.residues}

                d1, d2, t1, t2 = analyse_tetrad_geometry(
                    structure, residue_lookup, entry["g_tetrads"]
                )
                all_n1_o6.extend(d1)
                all_n2_n7.extend(d2)
                all_tors_n9.extend(t1)
                all_tors_o6.extend(t2)

                cache[name] = {
                    "n1_o6": d1,
                    "n2_n7": d2,
                    "tors_n9": t1,
                    "tors_o6": t2,
                }

                first = entry["g_tetrads"][0]
                print(f"    First tetrad: {first['tetrad']}")
                print(f"    Base pairs: {first['pairs']}")
                print(f"    Residues found in structure: {len(residues)}")
            else:
                print("    No residues found – skipping geometry")

    # ---------------------------- fit distributions --------------------------
    fitted_params: Dict[str, Any] = {}

    if all_n1_o6:
        mu, sigma = stats.norm.fit(all_n1_o6)
        fitted_params["n1_o6"] = {"mu": mu, "sigma": sigma}

    if all_n2_n7:
        mu, sigma = stats.norm.fit(all_n2_n7)
        fitted_params["n2_n7"] = {"mu": mu, "sigma": sigma}

    if all_tors_n9:
        kappa, loc, scale = stats.vonmises.fit(all_tors_n9, fscale=1)
        fitted_params["tors_n9"] = {"kappa": kappa, "loc": loc}

    if all_tors_o6:
        kappa, loc, scale = stats.vonmises.fit(all_tors_o6, fscale=1)
        fitted_params["tors_o6"] = {"kappa": kappa, "loc": loc}

    # ------------------------------ histograms ------------------------------
    if all_n1_o6:
        plt.figure()
        plt.hist(all_n1_o6, bins=30, density=True, color="steelblue", edgecolor="black")
        plt.title("N1–O6 distances")
        plt.xlabel("Distance (Å)")

        # overlay fitted normal
        params = fitted_params["n1_o6"]
        x = np.linspace(min(all_n1_o6), max(all_n1_o6), 200)
        plt.plot(x, stats.norm.pdf(x, params["mu"], params["sigma"]), "k--", lw=2)

    if all_n2_n7:
        plt.figure()
        plt.hist(all_n2_n7, bins=30, density=True, color="seagreen", edgecolor="black")
        plt.title("N2–N7 distances")
        plt.xlabel("Distance (Å)")

        params = fitted_params["n2_n7"]
        x = np.linspace(min(all_n2_n7), max(all_n2_n7), 200)
        plt.plot(x, stats.norm.pdf(x, params["mu"], params["sigma"]), "k--", lw=2)

    if all_tors_n9:
        plt.figure()
        plt.hist(
            np.degrees(all_tors_n9),
            bins=30,
            density=True,
            color="tomato",
            edgecolor="black",
        )
        plt.title("Torsion angle between N9 atoms")
        plt.xlabel("Angle (degrees)")

        params = fitted_params["tors_n9"]
        rad_range = np.linspace(min(all_tors_n9), max(all_tors_n9), 360)
        plt.plot(
            np.degrees(rad_range),
            stats.vonmises.pdf(rad_range, params["kappa"], loc=params["loc"]),
            "k--",
            lw=2,
        )

    if all_tors_o6:
        plt.figure()
        plt.hist(
            np.degrees(all_tors_o6),
            bins=30,
            density=True,
            color="orchid",
            edgecolor="black",
        )
        plt.title("Torsion angle between O6 atoms")
        plt.xlabel("Angle (degrees)")

        params = fitted_params["tors_o6"]
        rad_range = np.linspace(min(all_tors_o6), max(all_tors_o6), 360)
        plt.plot(
            np.degrees(rad_range),
            stats.vonmises.pdf(rad_range, params["kappa"], loc=params["loc"]),
            "k--",
            lw=2,
        )

    if any((all_n1_o6, all_n2_n7, all_tors_n9, all_tors_o6)):
        plt.tight_layout()
        plt.show()

    # ------------------------ save updated cache -----------------------------
    save_cache(args.cache_file, cache)

    # ---------------------- write fitted parameters --------------------------
    try:
        with args.output_json.open("w", encoding="utf-8") as fp:
            json.dump(fitted_params, fp, indent=2)
        print(f"Saved fitted parameters to {args.output_json}")
    except Exception as exc:  # noqa: BLE001
        print(f"Could not write fitted parameters: {exc}")


if __name__ == "__main__":
    main()

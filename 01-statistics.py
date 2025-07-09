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
from pathlib import Path
from typing import Dict, Any


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process all JSON files in a directory."
    )
    parser.add_argument(
        "directory",
        type=Path,
        help="Directory containing .json files",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mapping = load_json_files(args.directory)
    print(f"Loaded {len(mapping)} JSON files:")
    for name, entry in mapping.items():
        print(
            f"  {name}: "
            f"{len(entry['data'].get('nucleotides', []))} nucleotides, "
            f"{len(entry['nucleotide_map'])} indexed by fullName, "
            f"{len(entry['g_tetrads'])} G-tetrads"
        )


if __name__ == "__main__":
    main()

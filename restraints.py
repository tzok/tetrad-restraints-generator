#!/usr/bin/env python3
"""
02-restraints.py

Utility to analyse a “QRS” string – a sequence composed of dots (place-holders)
and paired letters.  Each letter must occur exactly **four** times:
    • twice in lowercase
    • twice in UPPERCASE

The script outputs a dictionary that maps every unique letter (lower-case key)
to a list with the **(index, character)** of all its occurrences in the input
string, ordered by position.

Example
-------
$ python 02-restraints.py "...q...Q...q...Q"
{'q': [(3, 'q'), (7, 'Q'), (11, 'q'), (15, 'Q')]}
"""

from __future__ import annotations

import argparse
import json
import sys
import math
import textwrap
from pathlib import Path
from typing import Dict, List, Tuple


def parse_qrs(qrs: str) -> Dict[str, List[Tuple[int, str]]]:
    """
    Build mapping *letter → list[(index, char)]* for a QRS string.

    Parameters
    ----------
    qrs : str
        Input QRS string; only ASCII letters are analysed.
        Any other character (dots ``.``, brackets, spaces, …) is ignored.

    Returns
    -------
    Dict[str, List[Tuple[int, str]]]
        Mapping described above.

    Raises
    ------
    ValueError
        If any letter does not appear exactly four times with the required
        2×lower-case + 2×upper-case pattern.
    """
    mapping: Dict[str, List[Tuple[int, str]]] = {}

    for idx, ch in enumerate(qrs):
        # Skip any placeholder that is NOT an ASCII letter (e.g. '.', brackets, spaces)
        if not ch.isalpha():
            continue
        key = ch.lower()
        mapping.setdefault(key, []).append((idx, ch))

    # Validate counts and case distribution
    for key, occurrences in mapping.items():
        if len(occurrences) != 4:
            raise ValueError(
                f"Letter '{key}' occurs {len(occurrences)} times, expected 4"
            )
        lower = sum(1 for _, c in occurrences if c.islower())
        upper = 4 - lower
        if lower != 2 or upper != 2:
            raise ValueError(
                f"Letter '{key}' must appear twice in each case "
                f"(got {lower} lower, {upper} upper)"
            )

    return mapping


def classify_letters(mapping: Dict[str, List[Tuple[int, str]]]) -> Dict[str, str]:
    """
    Classify each letter according to the pattern of upper/lower-case
    occurrences (ordered by position in the QRS string).

    Pattern → classification code
    ---------------------------------
      ULUL → O+
      LULU → O-
      ULLU → N+
      LULL → N-
      UULL → Z+
      LLUU → Z-

    Returns
    -------
    Dict[str, str]
        Mapping *letter → classification code*.

    Raises
    ------
    ValueError
        If any letter does not match one of the recognised patterns.
    """
    pattern_to_code = {
        "ULUL": "O+",
        "LULU": "O-",
        "ULLU": "N+",
        "LULL": "N-",
        "UULL": "Z+",
        "LLUU": "Z-",
    }

    class_map: Dict[str, str] = {}
    for letter, occurrences in mapping.items():
        # sort just in case (they should already be in order)
        ordered = sorted(occurrences, key=lambda t: t[0])
        pattern = "".join("U" if ch.isupper() else "L" for _, ch in ordered)
        try:
            class_map[letter] = pattern_to_code[pattern]
        except KeyError as exc:
            raise ValueError(
                f"Letter '{letter}' has unsupported pattern '{pattern}'"
            ) from exc
    return class_map


def load_fitted_params(path: Path) -> Dict[str, Dict[str, float]]:
    """Load *fitted_params.json* produced by 01-statistics.py."""
    try:
        with path.open("r", encoding="utf-8") as fp:
            return json.load(fp)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"Cannot find fitted parameters file: {path}") from exc


def build_pairs(order: str) -> List[Tuple[int, int]]:
    """Return list of 4 index pairs according to classification code."""
    patterns = {
        "O+": [(0, 1), (1, 2), (2, 3), (3, 0)],
        "O-": [(0, 3), (3, 2), (2, 1), (1, 0)],
        "N+": [(0, 1), (1, 3), (3, 2), (2, 0)],
        "N-": [(0, 2), (2, 3), (3, 1), (1, 0)],
        "Z+": [(0, 2), (2, 1), (1, 3), (3, 0)],
        "Z-": [(0, 3), (3, 1), (1, 2), (2, 0)],
    }
    return patterns[order]


def generate_restraints(qrs: str, params: Dict[str, Dict[str, float]]) -> List[str]:
    """
    Create distance and torsion restraints for *qrs* using statistics in
    *params*.

    Output format
    -------------
    Distance block:
        1;1
        idx1 N1 idx2 O6 mean minus plus
        idx1 N2 idx2 N7 mean minus plus
        … (8 lines per letter)

    Torsion block:
        1;1
        i1 N9 i2 N9 i3 N9 i4 N9 10 angle range 2
        i1 O6 i2 O6 i3 O6 i4 O6 10 angle range 2
        … (2 lines per letter)
    """
    mapping = parse_qrs(qrs)
    classes = classify_letters(mapping)

    # Extract statistics
    mu_n1_o6 = params["n1_o6"]["mu"]
    minus_n1_o6 = mu_n1_o6 - params["n1_o6"]["ci99"][0]
    plus_n1_o6 = params["n1_o6"]["ci99"][1] - mu_n1_o6

    mu_n2_n7 = params["n2_n7"]["mu"]
    minus_n2_n7 = mu_n2_n7 - params["n2_n7"]["ci99"][0]
    plus_n2_n7 = params["n2_n7"]["ci99"][1] - mu_n2_n7

    # Torsion statistics (stored in radians → convert to degrees)
    tors_n9_loc_deg = math.degrees(params["tors_n9"]["loc"])
    tors_n9_range_deg = math.degrees(params["tors_n9"]["ci99"][1]) - tors_n9_loc_deg

    tors_o6_loc_deg = math.degrees(params["tors_o6"]["loc"])
    tors_o6_range_deg = math.degrees(params["tors_o6"]["ci99"][1]) - tors_o6_loc_deg

    lines = ["1;1"]
    for letter, occurrences in mapping.items():
        ordered = sorted(occurrences, key=lambda t: t[0])
        indices = [idx + 1 for idx, _ in ordered]  # convert to 1-based
        code = classes[letter]
        for i, j in build_pairs(code):
            r_i, r_j = indices[i], indices[j]

            # N1-O6
            lines.append(
                f"{r_i:>4d} {'N1':>4} {r_j:>4d} {'O6':>4} {mu_n1_o6:7.2f} {minus_n1_o6:7.2f} {plus_n1_o6:7.2f}"
            )
            # N2-N7
            lines.append(
                f"{r_i:>4d} {'N2':>4} {r_j:>4d} {'N7':>4} {mu_n2_n7:7.2f} {minus_n2_n7:7.2f} {plus_n2_n7:7.2f}"
            )
    # ---------------------------- torsion block -----------------------------
    lines.append("1;1")
    torsion_order = {
        "O+": [0, 1, 2, 3],
        "O-": [0, 3, 2, 1],
        "N+": [0, 1, 3, 2],
        "N-": [0, 2, 3, 1],
        "Z+": [0, 2, 1, 3],
        "Z-": [0, 3, 1, 2],
    }

    for letter, occurrences in mapping.items():
        ordered = sorted(occurrences, key=lambda t: t[0])
        indices = [idx + 1 for idx, _ in ordered]
        order = torsion_order[classes[letter]]

        # Build list of torsion orders to apply by generating all four cyclic
        # rotations of the base order (i.e. shift the list progressively).
        orders_to_apply = [order[i:] + order[:i] for i in range(4)]

        for ord4 in orders_to_apply:
            i1, i2, i3, i4 = (indices[k] for k in ord4)

            # N9 torsion
            lines.append(
                f"{i1:>4d} {'N9':>4} {i2:>4d} {'N9':>4} {i3:>4d} {'N9':>4} {i4:>4d} {'N9':>4} "
                f"{10.0:7.1f} {tors_n9_loc_deg:8.2f} {tors_n9_range_deg:7.2f} {tors_n9_range_deg:7.2f}"
            )
            # O6 torsion
            lines.append(
                f"{i1:>4d} {'O6':>4} {i2:>4d} {'O6':>4} {i3:>4d} {'O6':>4} {i4:>4d} {'O6':>4} "
                f"{10.0:7.1f} {tors_o6_loc_deg:8.2f} {tors_o6_range_deg:7.2f} {tors_o6_range_deg:7.2f}"
            )

    # Split the combined list into separate distance and torsion blocks.
    # The first "1;1" starts the distance block; the second marks the
    # beginning of the torsion block.
    first_sep = lines.index("1;1", 1)
    distance_lines = lines[:first_sep]
    torsion_lines = lines[first_sep:]
    return distance_lines, torsion_lines, generate_xplor_script(qrs, params)


def _build_planar_group(indices: List[int]) -> str:
    """
    Build an XPLOR planar‐restraint *group* block for one tetrad.

    Parameters
    ----------
    indices : list[int]
        Four 1-based residue indices belonging to the tetrad.

    Returns
    -------
    str
        Formatted multiline XPLOR group block.
    """
    indices_sorted = sorted(indices)
    n9_lines = [f"        (resid {i} and name N9) or" for i in indices_sorted]
    o6_lines = [
        f"        (resid {i} and name O6){'' if idx == len(indices_sorted) - 1 else ' or'}"
        for idx, i in enumerate(indices_sorted)
    ]
    group_lines = (
        ["    group selection=("]
        + n9_lines
        + o6_lines
        + ["    )", "    weight=$planarweight", "    end"]
    )
    return "\n".join(group_lines)


def generate_xplor_script(qrs: str, params: Dict[str, Dict[str, float]]) -> str:
    """
    Produce a complete XPLOR script with planar restraints for every tetrad in
    *qrs*.  Geometry statistics *params* are currently unused but kept for
    API symmetry with *generate_restraints()*.
    """
    mapping = parse_qrs(qrs)
    groups_block = "\n".join(
        _build_planar_group([idx + 1 for idx, _ in sorted(occ, key=lambda t: t[0])])
        for _, occ in mapping.items()
    )

    template = textwrap.dedent(
        """\
parameter
    @TOPPAR:dna-rna-allatom.param
end

structure
    @$PSF_PATH$
end

evaluate ($filename="$INPUT_PDB_PATH$")
    coor @@$filename
end

noe
    nres=3000
    class all
    @$XPLOR_NOES_PATH$
end

evaluate ($planarweight=20)
    restraints planar @$XPLOR_PLANAR_PATH$
{groups}
end

noe
    ceiling=1000
    averaging  * cent
    potential  * square
    sqconstant * 1.
    sqexponent * 2
    scale      * 50.
end

restraints
    dihedral
    nassign=300
    @$XPLOR_DIHE_PATH$
    scale=1.0
end

flags
    exclude *
    include bond angl cdih impr vdw noe elec plan
end

minimize
    powell nstep=1000 nprint=100
end

evaluate ($filename="$OUTPUT_PDB_PATH$")
    write coordinates output =$filename
end

stop
"""
    )
    return template.format(groups=groups_block)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyse QRS string.")
    parser.add_argument("qrs", nargs="?", help="QRS string (dots and letters)")
    parser.add_argument(
        "--params",
        type=Path,
        default=Path("fitted_params.json"),
        help="Path to fitted_params.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    qrs_str = args.qrs
    if qrs_str is None:
        qrs_str = sys.stdin.readline().rstrip("\n")
    params = load_fitted_params(args.params)
    dist_lines, torsion_lines, xplor_script = generate_restraints(qrs_str, params)

    # Print distance restraints
    print("\n".join(dist_lines))
    print()
    # Print torsion restraints
    print("\n".join(torsion_lines))
    print()
    # Print XPLOR planar restraints
    print(xplor_script)


if __name__ == "__main__":
    main()

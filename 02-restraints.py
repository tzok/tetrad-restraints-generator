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
from pathlib import Path
from typing import Dict, List, Tuple


def parse_qrs(qrs: str) -> Dict[str, List[Tuple[int, str]]]:
    """
    Build mapping *letter → list[(index, char)]* for a QRS string.

    Parameters
    ----------
    qrs : str
        Input QRS string containing dots ``.`` and letters.

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
        if ch == ".":
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
    Create restraint lines for *qrs* using statistics in *params*.
    Returns list of lines including the mandatory header ``1;1``.
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

    lines = ["1;1"]
    for letter, occurrences in mapping.items():
        ordered = sorted(occurrences, key=lambda t: t[0])
        indices = [idx + 1 for idx, _ in ordered]  # convert to 1-based
        code = classes[letter]
        for i, j in build_pairs(code):
            r_i, r_j = indices[i], indices[j]

            # N1-O6
            lines.append(
                f"{r_i} N1 {r_j} O6 {mu_n1_o6:.3f} {minus_n1_o6:.3f} {plus_n1_o6:.3f}"
            )
            # N2-N7
            lines.append(
                f"{r_i} N2 {r_j} N7 {mu_n2_n7:.3f} {minus_n2_n7:.3f} {plus_n2_n7:.3f}"
            )
    return lines


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
    restraint_lines = generate_restraints(qrs_str, params)
    print("\n".join(restraint_lines))


if __name__ == "__main__":
    main()

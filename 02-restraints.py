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
import sys
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
            raise ValueError(f"Letter '{key}' occurs {len(occurrences)} times, expected 4")
        lower = sum(1 for _, c in occurrences if c.islower())
        upper = 4 - lower
        if lower != 2 or upper != 2:
            raise ValueError(
                f"Letter '{key}' must appear twice in each case "
                f"(got {lower} lower, {upper} upper)"
            )

    return mapping


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyse QRS string.")
    parser.add_argument(
        "qrs",
        nargs="?",
        help="QRS string consisting of dots and letters",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    qrs_str = args.qrs
    if qrs_str is None:
        qrs_str = sys.stdin.readline().rstrip("\n")
    mapping = parse_qrs(qrs_str)
    print(mapping)


if __name__ == "__main__":
    main()

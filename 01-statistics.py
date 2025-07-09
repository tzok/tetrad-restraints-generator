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


def load_json_files(directory: Path) -> Dict[str, Dict[str, Any]]:
    """
    Load every ``*.json`` file in *directory* that contains at least one helix
    with a non-empty ``quadruplexes`` list.

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
                json_map[json_file.name] = data
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
    for name, data in mapping.items():
        print(f"  {name}: {type(data).__name__} with {len(data)} top-level keys")


if __name__ == "__main__":
    main()

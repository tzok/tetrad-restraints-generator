#!/usr/bin/env python3
"""
Streamlit UI to generate distance & torsion restraints from a QRS string.

The page:
1. Loads *fitted_params.json* to pre-populate the form with recommended
   distance / torsion statistics.
2. Lets the user tweak every parameter.
3. On “Submit” it calls the generator implemented in *02-restraints.py* and
   shows the resulting block that can be copy-pasted into external software.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import streamlit as st

# Re-use the logic from the CLI tool
from 02_restraints import (
    generate_restraints,
    load_fitted_params,
)  # ruff: noqa: E402  # local import


DEFAULT_PARAM_FILE = Path("fitted_params.json")


@st.cache_data(show_spinner=False)
def _load_default_params(path: Path) -> dict:
    try:
        return load_fitted_params(path)
    except FileNotFoundError:
        st.warning(f"Parameters file not found: {path}\nUsing placeholders.")
        # minimal fallback with dummy ranges
        return {
            "n1_o6": {"mu": 3.0, "ci99": [2.5, 3.5]},
            "n2_n7": {"mu": 3.2, "ci99": [2.7, 3.7]},
            "tors_n9": {"loc": math.radians(0.0), "ci99": [math.radians(-30), math.radians(30)]},
            "tors_o6": {"loc": math.radians(0.0), "ci99": [math.radians(-30), math.radians(30)]},
        }


def ci_to_range(mu: float, ci: list[float]) -> tuple[float, float]:
    """Return (minus, plus) positive ranges from mean & CI."""
    return (mu - ci[0], ci[1] - mu)


def range_to_ci(mu: float, minus: float, plus: float) -> list[float]:
    return [mu - minus, mu + plus]


def main() -> None:
    st.title("QRS → restraint generator")

    # ---------------------------- sidebar -----------------------------------
    st.sidebar.header("Parameters")

    # load defaults
    param_file = st.sidebar.file_uploader(
        "fitted_params.json (optional)",
        type=["json"],
        label_visibility="collapsed",
    )
    if param_file is not None:
        params_raw = json.load(param_file)
    else:
        params_raw = _load_default_params(DEFAULT_PARAM_FILE)

    mu_n1o6 = params_raw["n1_o6"]["mu"]
    minus_n1o6, plus_n1o6 = ci_to_range(mu_n1o6, params_raw["n1_o6"]["ci99"])

    mu_n2n7 = params_raw["n2_n7"]["mu"]
    minus_n2n7, plus_n2n7 = ci_to_range(mu_n2n7, params_raw["n2_n7"]["ci99"])

    loc_n9_deg = math.degrees(params_raw["tors_n9"]["loc"])
    range_n9_deg = math.degrees(params_raw["tors_n9"]["ci99"][1]) - loc_n9_deg

    loc_o6_deg = math.degrees(params_raw["tors_o6"]["loc"])
    range_o6_deg = math.degrees(params_raw["tors_o6"]["ci99"][1]) - loc_o6_deg

    # Distance N1-O6
    st.sidebar.subheader("N1–O6 distance (Å)")
    mu_n1o6 = st.sidebar.number_input("Mean", value=float(mu_n1o6), step=0.01, key="n1o6_mu")
    minus_n1o6 = st.sidebar.number_input("Minus", value=float(minus_n1o6), step=0.01, key="n1o6_minus")
    plus_n1o6 = st.sidebar.number_input("Plus", value=float(plus_n1o6), step=0.01, key="n1o6_plus")

    # Distance N2-N7
    st.sidebar.subheader("N2–N7 distance (Å)")
    mu_n2n7 = st.sidebar.number_input("Mean ", value=float(mu_n2n7), step=0.01, key="n2n7_mu")
    minus_n2n7 = st.sidebar.number_input("Minus ", value=float(minus_n2n7), step=0.01, key="n2n7_minus")
    plus_n2n7 = st.sidebar.number_input("Plus ", value=float(plus_n2n7), step=0.01, key="n2n7_plus")

    # Torsion N9
    st.sidebar.subheader("N9 planarity (°)")
    loc_n9_deg = st.sidebar.number_input("Angle", value=float(loc_n9_deg), step=1.0, key="n9_loc")
    range_n9_deg = st.sidebar.number_input("± Range", value=float(range_n9_deg), step=1.0, key="n9_rng")

    # Torsion O6
    st.sidebar.subheader("O6 planarity (°)")
    loc_o6_deg = st.sidebar.number_input("Angle ", value=float(loc_o6_deg), step=1.0, key="o6_loc")
    range_o6_deg = st.sidebar.number_input("± Range ", value=float(range_o6_deg), step=1.0, key="o6_rng")

    # ------------------------------ main ------------------------------------
    qrs_input = st.text_input("QRS sequence", placeholder="...q...Q...q...Q")

    if st.button("Submit") and qrs_input:
        # Build params dictionary in the same structure as 01-statistics output
        params = {
            "n1_o6": {
                "mu": mu_n1o6,
                "ci99": range_to_ci(mu_n1o6, minus_n1o6, plus_n1o6),
            },
            "n2_n7": {
                "mu": mu_n2n7,
                "ci99": range_to_ci(mu_n2n7, minus_n2n7, plus_n2n7),
            },
            "tors_n9": {
                "loc": math.radians(loc_n9_deg),
                "ci99": [
                    math.radians(loc_n9_deg - range_n9_deg),
                    math.radians(loc_n9_deg + range_n9_deg),
                ],
            },
            "tors_o6": {
                "loc": math.radians(loc_o6_deg),
                "ci99": [
                    math.radians(loc_o6_deg - range_o6_deg),
                    math.radians(loc_o6_deg + range_o6_deg),
                ],
            },
        }

        try:
            restraints = generate_restraints(qrs_input, params)
            result_text = "\n".join(restraints)
            st.text_area("Generated restraints", value=result_text, height=300)
            st.download_button(
                "Download .rst file",
                data=result_text,
                file_name="restraints.rst",
                mime="text/plain",
            )
        except Exception as exc:  # noqa: BLE001
            st.error(f"Error: {exc}")


if __name__ == "__main__":
    main()

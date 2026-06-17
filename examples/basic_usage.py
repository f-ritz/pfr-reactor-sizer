"""
Basic usage examples for pfrsizer.

Run with:
    python -m examples.basic_usage
or (after `pip install -e .`):
    python examples/basic_usage.py
"""

import sys
from pathlib import Path

# Allow running directly from source tree without installation
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pfrsizer import Reaction, Feed, PFRConfig
from pfrsizer.solvers import solve_pfr_isothermal, solve_pfr_adiabatic
from pfrsizer.plot import plot_profiles


def main():
    print("=== pfrsizer basic usage ===\n")

    # Example 1: Isothermal first-order reaction A -> B
    print("1) Isothermal gas-phase PFR: A -> B")
    rxn1 = Reaction(
        stoichiometry={"A": -1, "B": 1},
        k0=0.05,      # 1/s
        E=0.0,
        delta_H=0.0,
        orders={"A": 1},
        name="A->B isothermal"
    )
    feed1 = Feed(F0={"A": 1.0}, T0=350.0, P0=101325 * 4)
    cfg1 = PFRConfig(target_X=0.85, max_V=60.0)

    res1 = solve_pfr_isothermal(rxn1, feed1, cfg1)
    print(f"   Required V for X=0.85: {res1.final_V:.4f} m³")
    print(f"   Final X: {res1.final_X:.4f}")
    plot_profiles(res1, show=False, save_path="example1_isothermal.png")
    print("   Saved plot: example1_isothermal.png\n")

    # Example 2: Adiabatic exothermic reaction with mole increase
    print("2) Adiabatic PFR: A -> 2B (exothermic)")
    rxn2 = Reaction(
        stoichiometry={"A": -1, "B": 2},
        k0=8.5e8,
        E=65000,         # J/mol
        delta_H=-72000,  # exothermic
        orders={"A": 1.0},
        name="A->2B adiabatic"
    )
    feed2 = Feed(
        F0={"A": 1.0, "Inert": 0.4},
        T0=330.0,
        P0=3 * 101325
    )
    Cp = {"A": 115.0, "B": 82.0, "Inert": 29.5}  # J/mol K

    cfg2 = PFRConfig(mode="adiabatic", target_X=0.80, max_V=3.5)

    res2 = solve_pfr_adiabatic(rxn2, feed2, Cp, cfg2)
    print(f"   V for X=0.80: {res2.final_V:.4f} m³")
    print(f"   Outlet T: {res2.final_T:.1f} K  (started at {feed2.T0} K)")
    print(f"   Outlet P: {res2.final_P / 101325:.3f} atm")
    plot_profiles(res2, show=False, save_path="example2_adiabatic.png")
    print("   Saved plot: example2_adiabatic.png\n")

    print("Done. Open the .png files to see the profiles.")
    print("You can also run the CLI:  pfrsizer example A_to_2B_adiabatic")


if __name__ == "__main__":
    main()

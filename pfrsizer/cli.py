"""
Command-line interface for pfrsizer.
"""

import sys
from pathlib import Path

# --- Robust support for running cli.py directly ----------------------------
def _ensure_pfrsizer_on_path():
    here = Path(__file__).resolve()
    for ancestor in [here.parent] + list(here.parents):
        if (ancestor / "pfrsizer" / "__init__.py").is_file():
            if str(ancestor) not in sys.path:
                sys.path.insert(0, str(ancestor))
            return
    sys.path.insert(0, str(here.parent.parent))

_ensure_pfrsizer_on_path()
# ---------------------------------------------------------------------------

import argparse
from typing import Optional

try:
    from rich import print as rprint
    from rich.table import Table
    from rich.console import Console
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False
    rprint = print

from pfrsizer.models import Reaction, Feed, PFRConfig
from pfrsizer.solvers import solve_pfr
from pfrsizer.plot import plot_profiles


def _print_result(result, use_rich: bool = RICH_AVAILABLE):
    if use_rich:
        console = Console()
        table = Table(title="PFR Result Summary", show_header=True, header_style="bold magenta")
        table.add_column("Quantity", style="cyan")
        table.add_column("Value", style="green")

        table.add_row("Mode", result.config.mode)
        table.add_row("Limiting species", result.limiting_species)
        table.add_row("Reactor Length L (m)", f"{result.final_z:.6g}")
        table.add_row("Reactor Volume V (m³)", f"{result.final_V:.6g}")
        table.add_row("Diameter (m)", f"{result.config.diameter if result.config else '?'}")
        table.add_row("Cross-sectional area (m²)", f"{result.config.A_c if result.config else '?':.6f}")
        table.add_row("Final Conversion X", f"{result.final_X:.5f}")
        table.add_row("Final Temperature (K)", f"{result.final_T:.2f}")
        table.add_row("Final Pressure (Pa)", f"{result.final_P:.1f} ({result.final_P/101325:.4f} atm)")
        table.add_row("Success", str(result.success))
        table.add_row("Message", result.message)

        console.print(table)

        # Outlet flows table
        flow_table = Table(title="Outlet Molar Flows (mol/s)", header_style="bold blue")
        flow_table.add_column("Species")
        flow_table.add_column("Flow rate")
        for sp, f in result.outlet_flows().items():
            flow_table.add_row(sp, f"{f:.6g}")
        console.print(flow_table)
    else:
        print("\n=== PFR Result Summary ===")
        print(f"Mode:                {result.config.mode}")
        print(f"Limiting species:    {result.limiting_species}")
        print(f"Reactor Length (m):  {result.final_z:.6g}")
        print(f"Reactor Volume (m³): {result.final_V:.6g}")
        print(f"Diameter (m):        {result.config.diameter if result.config else '?'}")
        print(f"Cross-section (m²):  {result.config.A_c if result.config else '?':.6f}")
        print(f"Final Conversion X:  {result.final_X:.5f}")
        print(f"Final T (K):         {result.final_T:.2f}")
        print(f"Final P (Pa):        {result.final_P:.1f}  ({result.final_P/101325:.4f} atm)")
        print(f"Success:             {result.success}")
        print(f"Message:             {result.message}")
        print("\nOutlet flows (mol/s):")
        for sp, f in result.outlet_flows().items():
            print(f"  {sp}: {f:.6g}")


def run_example(example: str = "A_to_B_isothermal"):
    """Run one of the built-in example cases."""
    if example == "A_to_B_isothermal":
        # Classic simple isothermal gas phase: A -> B
        rxn = Reaction(
            stoichiometry={"A": -1, "B": 1},
            k0=0.05,          # s^-1  (first order)
            E=0.0,            # isothermal, E irrelevant
            delta_H=0.0,
            orders={"A": 1},
            name="A -> B (1st order)"
        )
        feed = Feed(
            F0={"A": 1.0, "B": 0.0},
            T0=350.0,
            P0=5 * 101325,    # 5 atm
        )
        cfg = PFRConfig(
            mode="isothermal",
            target_X=0.8,
            max_L=50.0,
            diameter=0.1,
            pressure_model="constant",
        )
        result = solve_pfr(rxn, feed, config=cfg)
        print("=== Example: Isothermal first-order A -> B ===")
        _print_result(result)
        plot_profiles(result, save_path=None)
        return result

    elif example == "A_to_2B_adiabatic":
        # Exothermic gas expansion reaction A -> 2B  (volume increase + temp rise)
        rxn = Reaction(
            stoichiometry={"A": -1, "B": 2},
            k0=8.5e8,            # tuned to give reasonable rates at 330-500K
            E=65000,             # J/mol
            delta_H=-72000,      # exothermic J/mol
            orders={"A": 1},
            name="A -> 2B (exothermic)"
        )
        feed = Feed(
            F0={"A": 1.0, "B": 0.0, "Inert": 0.4},
            T0=330.0,
            P0=3 * 101325,
        )
        Cp = {"A": 115.0, "B": 82.0, "Inert": 29.5}   # J/mol/K approx values

        cfg = PFRConfig(
            mode="adiabatic",
            target_X=0.80,
            max_L=50.0,
            diameter=0.1,
            pressure_model="constant",
        )
        result = solve_pfr(rxn, feed, Cp=Cp, config=cfg)
        print("=== Example: Adiabatic A -> 2B (exothermic, volume change) ===")
        _print_result(result)
        plot_profiles(result)
        return result

    elif example == "pressure_drop":
        # Isothermal with pressure drop
        rxn = Reaction(
            stoichiometry={"A": -1, "B": 1},
            k0=0.12,
            E=0,
            delta_H=0,
            orders={"A": 1},
            name="A->B with pressure drop"
        )
        feed = Feed(F0={"A": 0.8, "B": 0.0}, T0=400.0, P0=8*101325)
        cfg = PFRConfig(
            target_X=0.9,
            max_L=100.0,
            diameter=0.1,
            pressure_model="simple_drop",
            alpha=0.08,
        )
        result = solve_pfr(rxn, feed, config=cfg)
        print("=== Example: Isothermal with pressure drop (alpha parameter) ===")
        _print_result(result)
        plot_profiles(result)
        return result

    else:
        print(f"Unknown example '{example}'. Available: A_to_B_isothermal, A_to_2B_adiabatic, pressure_drop")
        return None


def main(argv: Optional[list] = None):
    parser = argparse.ArgumentParser(
        prog="pfrsizer",
        description="PFR Reactor Sizer - design non-catalytic gas-phase plug flow reactors (isothermal & adiabatic)."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Example subcommand
    ex = subparsers.add_parser("example", help="Run a built-in example case")
    ex.add_argument(
        "name",
        nargs="?",
        default="A_to_B_isothermal",
        choices=["A_to_B_isothermal", "A_to_2B_adiabatic", "pressure_drop"],
        help="Name of the example to run"
    )
    ex.add_argument("--no-plot", action="store_true", help="Do not show plots")

    args = parser.parse_args(argv)

    if args.command == "example":
        result = run_example(args.name)
        if result is None:
            sys.exit(1)
        if not args.no_plot:
            pass

    return 0


if __name__ == "__main__":
    sys.exit(main())

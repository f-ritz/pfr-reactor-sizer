"""
Plotting utilities for PFR results.

Most profiles are plotted against reactor length z (m).
Temperature vs Conversion (T vs X) is also included (linear for adiabatic).
Volume V (m^3) is shown on a secondary axis or in the title.
"""

from typing import Optional
import matplotlib.pyplot as plt
import numpy as np

from .models import PFRResult


def plot_profiles(
    result: PFRResult,
    show: bool = True,
    save_path: Optional[str] = None,
    figsize: tuple = (11, 7),
):
    """
    Plot key profiles from a PFR simulation:
      - Conversion X vs z
      - Temperature T vs z
      - Temperature vs Conversion (T vs X)  -- new; linear for adiabatic (energy balance)
      - Pressure P vs z (if variable)
      - Molar flows vs z
      - Reaction rate vs z
    """
    z = np.array(result.z)
    if len(z) == 0:
        print("No data to plot.")
        return

    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    ax_idx = 0

    # 1. Conversion
    ax = axes[ax_idx]
    ax_idx += 1
    ax.plot(z, result.X, "b-", linewidth=2, label="Conversion X")
    ax.set_xlabel("Reactor Length z (m)")
    ax.set_ylabel("Conversion X")
    ax.set_title("Conversion Profile")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")

    # 2. Temperature vs z
    ax = axes[ax_idx]
    ax_idx += 1
    ax.plot(z, result.T, "r-", linewidth=2, label="Temperature")
    if result.feed:
        ax.axhline(result.feed.T0, color="gray", linestyle="--", label=f"T0 = {result.feed.T0:.1f} K")
        ax.legend(fontsize=8)
    ax.set_xlabel("Reactor Length z (m)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Temperature Profile")
    ax.grid(True, alpha=0.3)

    # 3. NEW: Temperature vs Conversion (T vs X)
    # For adiabatic: approximately linear (energy balance line)
    # For isothermal: flat
    ax = axes[ax_idx]
    ax_idx += 1
    ax.plot(result.X, result.T, "g-", linewidth=2)
    if result.feed:
        ax.axhline(result.feed.T0, color="gray", linestyle="--", label=f"T0 = {result.feed.T0:.1f} K")
        ax.legend(fontsize=8)
    ax.set_xlabel("X")
    ax.set_ylabel("T (K)")
    ax.set_title("Temperature vs Conversion")
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)

    # 4. Pressure (if enabled)
    if result.config.pressure_model != "constant":
        ax = axes[ax_idx]
        ax_idx += 1
        ax.plot(z, [p / 101325 for p in result.P], "g-", linewidth=2)
        ax.set_xlabel("Reactor Length z (m)")
        ax.set_ylabel("Pressure (atm)")
        ax.set_title("Pressure Profile")
        ax.grid(True, alpha=0.3)

    # 5. Molar flows (all species)
    ax = axes[ax_idx]
    ax_idx += 1
    if result.F and len(result.F) > 0:
        species = sorted(result.F[0].keys())
        for sp in species:
            Fi = [f.get(sp, 0.0) for f in result.F]
            ax.plot(z, Fi, linewidth=1.8, label=sp)
        ax.legend(loc="best", fontsize=8, ncol=2)
    ax.set_xlabel("Reactor Length z (m)")
    ax.set_ylabel("Molar flow Fi (mol/s)")
    ax.set_title("Molar Flow Profiles")
    ax.grid(True, alpha=0.3)

    # 6. Rate
    ax = axes[ax_idx]
    ax_idx += 1
    ax.plot(z, result.r, "m-", linewidth=2)
    ax.set_xlabel("z (m)")
    ax.set_ylabel("r (mol/m³/s)")
    ax.set_title("Reaction Rate")
    ax.grid(True, alpha=0.3)

    # Hide unused axes
    while ax_idx < len(axes):
        axes[ax_idx].axis("off")
        ax_idx += 1

    # Show geometry info in title
    D = result.config.diameter if result.config and result.config.diameter else "?"
    A_c = result.config.A_c if result.config else "?"
    fig.suptitle(
        f"PFR Simulation — {result.reaction.name} — {result.config.mode}\n"
        f"D = {D} m, A_c = {A_c:.6f} m², L = {result.final_z:.3f} m, V = {result.final_V:.4f} m³",
        fontsize=13
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved plot to {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_rate_and_concentrations(
    result: PFRResult,
    show: bool = True,
    save_path: Optional[str] = None,
    figsize: tuple = (10, 4),
):
    """Plot reaction rate and outlet concentrations (or profiles)."""
    z = np.array(result.z)
    if len(z) == 0:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Rate
    ax1.plot(z, result.r, "m-", linewidth=2)
    ax1.set_xlabel("z (m)")
    ax1.set_ylabel("Rate r (mol/m³/s)")
    ax1.set_title("Reaction Rate Profile")
    ax1.grid(True, alpha=0.3)

    # Final concentrations (bar)
    C_out = result.outlet_concentrations()
    species = list(C_out.keys())
    concs = [C_out[s] for s in species]
    ax2.barh(species, concs, color="teal", alpha=0.75)
    ax2.set_xlabel("Concentration (mol/m³)")
    ax2.set_title("Outlet Concentrations (approx. ideal gas)")
    ax2.grid(True, axis="x", alpha=0.3)

    fig.suptitle("PFR Rate & Concentrations", fontsize=13)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)

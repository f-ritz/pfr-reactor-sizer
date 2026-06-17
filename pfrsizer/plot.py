"""
Plotting utilities for PFR results.
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
      - Conversion X vs V
      - Temperature T vs V (if adiabatic)
      - Pressure P vs V (if variable)
      - Molar flows vs V
      - Reaction rate vs V
    """
    V = np.array(result.V)
    if len(V) == 0:
        print("No data to plot.")
        return

    n_plots = 2
    if result.config.mode == "adiabatic":
        n_plots += 1
    if result.config.pressure_model != "constant":
        n_plots += 1

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()
    ax_idx = 0

    # 1. Conversion
    ax = axes[ax_idx]
    ax_idx += 1
    ax.plot(V, result.X, "b-", linewidth=2, label="Conversion X")
    ax.set_xlabel("Reactor Volume V (m³)")
    ax.set_ylabel("Conversion X")
    ax.set_title("Conversion Profile")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")

    # 2. Temperature
    if result.config.mode == "adiabatic":
        ax = axes[ax_idx]
        ax_idx += 1
        ax.plot(V, result.T, "r-", linewidth=2, label="Temperature")
        ax.axhline(result.feed.T0, color="gray", linestyle="--", label=f"T0 = {result.feed.T0:.1f} K")
        ax.set_xlabel("Reactor Volume V (m³)")
        ax.set_ylabel("Temperature (K)")
        ax.set_title("Temperature Profile (Adiabatic)")
        ax.grid(True, alpha=0.3)
        ax.legend()

    # 3. Pressure
    if result.config.pressure_model != "constant":
        ax = axes[ax_idx]
        ax_idx += 1
        ax.plot(V, [p / 101325 for p in result.P], "g-", linewidth=2)
        ax.set_xlabel("Reactor Volume V (m³)")
        ax.set_ylabel("Pressure (atm)")
        ax.set_title("Pressure Profile")
        ax.grid(True, alpha=0.3)

    # 4. Molar flows (all species)
    ax = axes[ax_idx]
    ax_idx += 1
    species = sorted(result.F[0].keys())
    for sp in species:
        Fi = [f.get(sp, 0.0) for f in result.F]
        ax.plot(V, Fi, linewidth=1.8, label=sp)
    ax.set_xlabel("Reactor Volume V (m³)")
    ax.set_ylabel("Molar flow Fi (mol/s)")
    ax.set_title("Molar Flow Profiles")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, ncol=2)

    # Hide unused axes
    while ax_idx < len(axes):
        axes[ax_idx].axis("off")
        ax_idx += 1

    fig.suptitle(f"PFR Simulation — {result.reaction.name} — {result.config.mode}", fontsize=14)
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
    V = np.array(result.V)
    if len(V) == 0:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # Rate
    ax1.plot(V, result.r, "m-", linewidth=2)
    ax1.set_xlabel("V (m³)")
    ax1.set_ylabel("Rate r (mol/m³/s)")
    ax1.set_title("Reaction Rate Profile")
    ax1.grid(True, alpha=0.3)

    # Final concentrations (bar)
    # We compute concentrations along the reactor for one species or all
    # For simplicity show final concentrations as horizontal bars
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

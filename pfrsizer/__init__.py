"""
pfrsizer - Plug Flow Reactor Sizer

Design and analysis of non-catalytic gas-phase plug flow reactors (PFRs).
Supports isothermal and adiabatic operation, with optional pressure drop.
"""

__version__ = "0.1.0"

from .models import Reaction, Feed, PFRConfig, PFRResult
from .solvers import solve_pfr_isothermal, solve_pfr_adiabatic
from .core import ideal_gas_concentration
from .plot import plot_profiles, plot_rate_and_concentrations

__all__ = [
    "Reaction",
    "Feed",
    "PFRConfig",
    "PFRResult",
    "solve_pfr_isothermal",
    "solve_pfr_adiabatic",
    "ideal_gas_concentration",
    "plot_profiles",
    "plot_rate_and_concentrations",
]

"""
pfrsizer - Plug Flow Reactor Sizer

Design and analysis of non-catalytic gas-phase plug flow reactors (PFRs).
Supports isothermal (with computed heat duty) and adiabatic operation.
PubChem lookup for species properties + Tkinter GUI for Windows EXE distribution.
"""

__version__ = "1.0.0"

from .models import Reaction, Feed, PFRConfig, PFRResult
from .solvers import solve_pfr_isothermal, solve_pfr_adiabatic, solve_pfr
from .core import ideal_gas_concentration, inlet_volumetric_flow, compute_concentrations_from_v
from .plot import plot_profiles, plot_rate_and_concentrations
from .pubchem import lookup_species, lookup_species_list
from .nist_webbook import lookup_nist, lookup_nist_list
from .reaction_parser import parse_reaction, pretty_stoich
from .components import Component, compute_delta_H, build_components_from_pubchem

__all__ = [
    "Reaction",
    "Feed",
    "PFRConfig",
    "PFRResult",
    "solve_pfr_isothermal",
    "solve_pfr_adiabatic",
    "solve_pfr",
    "ideal_gas_concentration",
    "inlet_volumetric_flow",
    "compute_concentrations_from_v",
    "plot_profiles",
    "plot_rate_and_concentrations",
    "lookup_species",
    "lookup_species_list",
    "lookup_nist",
    "lookup_nist_list",
    "parse_reaction",
    "pretty_stoich",
    "Component",
    "compute_delta_H",
    "build_components_from_pubchem",
]

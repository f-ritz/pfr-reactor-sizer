"""
Core utility functions for PFR calculations (ideal gas, concentrations, etc.).

All equations follow Fogler's *Elements of Chemical Reaction Engineering*.
The independent variable is reactor length z (m).
"""

from __future__ import annotations
from typing import Dict
import math

from .models import Reaction, Feed, PFRConfig


R_DEFAULT = 8.314  # J / mol / K


def volumetric_flow(Ft: float, T: float, P: float, v0: float, Ft0: float, T0: float, P0: float) -> float:
    """
    Volumetric flow rate v (m^3/s) for ideal gas at current conditions.

    v = v0 * (Ft / Ft0) * (T / T0) * (P0 / P)

    Fogler eq. 4-23 (variable-density gas-phase systems).
    """
    if Ft0 <= 0 or P <= 0:
        return 0.0
    return v0 * (Ft / Ft0) * (T / T0) * (P0 / P)


def ideal_gas_concentration(Fi: float, Ft: float, T: float, P: float, R: float = R_DEFAULT) -> float:
    """
    Concentration of species i (mol/m^3) from ideal gas law.

    C_i = (Fi / Ft) * (P / (R * T))
    """
    if Ft <= 0 or T <= 0 or P <= 0:
        return 0.0
    y_i = Fi / Ft
    return y_i * (P / (R * T))


def compute_concentrations(F: Dict[str, float], T: float, P: float, R: float = R_DEFAULT) -> Dict[str, float]:
    """Compute all concentrations from current molar flows + T,P."""
    Ft = sum(max(f, 0.0) for f in F.values())
    if Ft <= 0:
        return {sp: 0.0 for sp in F}
    return {sp: ideal_gas_concentration(Fi, Ft, T, P, R) for sp, Fi in F.items()}


def compute_concentrations_from_v(F: Dict[str, float], v: float) -> Dict[str, float]:
    """
    Compute concentrations from molar flows and volumetric flow rate.

    C_i = F_i / v

    This is the most direct method — no ideal gas assumption needed if v is known.
    """
    if v <= 0:
        return {sp: 0.0 for sp in F}
    Ft = sum(max(f, 0.0) for f in F.values())
    if Ft <= 0:
        return {sp: 0.0 for sp in F}
    return {sp: max(F.get(sp, 0.0), 0.0) / v for sp in F}


def delta_n(reaction: Reaction) -> float:
    """
    Change in number of moles for the reaction (sum nu).
    Used for volume change factor.
    """
    return sum(reaction.stoichiometry.values())


def inlet_volumetric_flow(feed: Feed, R: float = R_DEFAULT) -> float:
    """
    Inlet volumetric flow rate (m^3/s) from ideal gas law.

    v0 = F_T0 * R * T0 / P0
    """
    return feed.total_F0 * R * feed.T0 / feed.P0


def superficial_velocity(v: float, A_c: float) -> float:
    """
    Superficial velocity u (m/s) = v / A_c.
    """
    if A_c <= 0:
        return 0.0
    return v / A_c

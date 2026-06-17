"""
Core utility functions for PFR calculations (ideal gas, concentrations, etc.).
"""

from typing import Dict
import math


R_DEFAULT = 8.314  # J / mol / K


def ideal_gas_volumetric_flow(Ft: float, T: float, P: float, Ft0: float, T0: float, P0: float) -> float:
    """
    Volumetric flow rate v (m^3/s) for ideal gas at current conditions.

    v = v0 * (Ft / Ft0) * (T / T0) * (P0 / P)
    """
    if Ft0 <= 0 or P <= 0:
        return 0.0
    # We don't store v0 explicitly; caller should compute ratio
    # To avoid needing v0, we return the factor relative to inlet.
    # But often useful: caller knows v0 or we compute concentration directly.
    # For direct v:
    # We return a scaled value; better to provide concentration helper.
    return (Ft / Ft0) * (T / T0) * (P0 / P)   # this is the *factor*, multiply by v0 outside if known


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


def delta_n(reaction: "Reaction") -> float:
    """
    Change in number of moles for the reaction (sum nu).
    Used for volume change factor.
    """
    return sum(reaction.stoichiometry.values())

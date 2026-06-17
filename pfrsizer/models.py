"""
Data models for PFR reactor design.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional, Literal, Any

from .components import Component, compute_delta_H as _compute_delta_H_from_components


@dataclass
class Reaction:
    """
    Definition of a single gas-phase reaction.

    Attributes:
        stoichiometry: Dict of stoichiometric coefficients. e.g. {'A': -1, 'B': 1, 'C': 1}
                       Negative for reactants. Must include at least one negative.
        k0: Pre-exponential factor. Units depend on reaction order.
            For first order (s^-1), second order (m^3/mol/s), etc.
        E: Activation energy (J/mol).
        delta_H: Heat of reaction at reference T (J/mol). Negative = exothermic.
                 Used for adiabatic energy balance. Sign convention: heat released by rxn is negative delta_H.
        orders: Optional explicit reaction orders per species (power law).
                If None, uses |stoich coeff| for each participating species (common approximation).
                Example: {'A': 1, 'B': 0.5}
        name: Optional label for the reaction.
    """
    stoichiometry: Dict[str, float]
    k0: float
    E: float
    delta_H: float = 0.0
    orders: Optional[Dict[str, float]] = None
    name: str = "rxn"

    def __post_init__(self):
        if not self.stoichiometry:
            raise ValueError("Stoichiometry cannot be empty")
        reactants = [s for s, nu in self.stoichiometry.items() if nu < 0]
        if not reactants:
            raise ValueError("Reaction must have at least one reactant (negative stoich coeff)")

    @property
    def species(self) -> list[str]:
        return list(self.stoichiometry.keys())

    def rate(self, C: Dict[str, float], T: float, R: float = 8.314) -> float:
        """
        Compute reaction rate r (mol / m^3 / s) given concentrations and temperature.

        r = k(T) * product( C_i ** order_i )
        """
        k = self.k0 * pow(2.718281828, -self.E / (R * T))  # exp(-E/RT) without numpy for core
        rate = k
        orders = self.orders if self.orders is not None else {
            sp: abs(nu) for sp, nu in self.stoichiometry.items()
        }
        for sp, ord_ in orders.items():
            conc = max(C.get(sp, 0.0), 0.0)
            rate *= pow(conc, ord_)
        return rate

    def extent_rates(self, r: float) -> Dict[str, float]:
        """Return species generation rates r_j = nu_j * r"""
        return {sp: nu * r for sp, nu in self.stoichiometry.items()}

    def compute_delta_H_from_Hf(self, components: Dict[str, Component]) -> Optional[float]:
        """Compute delta_H using component heats of formation if all are available."""
        return _compute_delta_H_from_components(components, self.stoichiometry)


@dataclass
class Feed:
    """
    Inlet conditions to the reactor.

    All flows are molar flow rates (mol/s).
    T0 in Kelvin.
    P0 in Pascal (Pa). Common: 101325 Pa = 1 atm.
    """
    F0: Dict[str, float]  # species -> molar flow rate (mol/s)
    T0: float = 300.0     # K
    P0: float = 101325.0  # Pa

    @property
    def total_F0(self) -> float:
        return sum(max(f, 0.0) for f in self.F0.values())

    @property
    def species(self) -> list[str]:
        return list(self.F0.keys())

    def mole_fractions(self) -> Dict[str, float]:
        Ft = self.total_F0
        if Ft <= 0:
            return {sp: 0.0 for sp in self.F0}
        return {sp: max(f, 0.0) / Ft for sp, f in self.F0.items()}


@dataclass
class PFRConfig:
    """
    Configuration and operating conditions for the PFR simulation.

    mode:
        - 'isothermal': Temperature is held constant (T = isothermal_T or feed.T0).
                        The solver will compute the heat duty (Q) required to maintain isothermality.
        - 'adiabatic': No heat transfer. Energy balance is integrated → temperature profile.

    isothermal_T: If mode='isothermal' and this is set, the reactor runs at this fixed temperature (K).
                  Otherwise feed.T0 is used.

    pressure_model:
        - 'constant': P = P0 always (default for many homogeneous PFRs)
        - 'simple_drop': Uses dP/dV = -alpha * (volumetric flow ratio) * (P0/P) etc.

    alpha: Pressure drop parameter. See solvers for formulation. 0 = no drop.

    target_X: If provided and >0, integrate until limiting reactant reaches this conversion.
    max_V: Hard upper limit on volume to integrate (safety). m^3

    cross_sectional_area: A_c = pi*D^2/4 (m^2). Used for superficial velocity reporting.
    """
    mode: Literal["isothermal", "adiabatic"] = "isothermal"
    isothermal_T: Optional[float] = None
    pressure_model: Literal["constant", "simple_drop"] = "constant"
    alpha: float = 0.0          # Pressure drop param

    target_X: Optional[float] = None
    max_V: float = 10.0         # m^3 upper bound

    cross_sectional_area: Optional[float] = None  # m^2

    R: float = 8.314            # J/mol/K  (gas constant)


@dataclass
class PFRResult:
    """
    Results from a PFR integration run.
    """
    V: list[float]                  # reactor volume points (m^3)
    X: list[float]                  # conversion of limiting reactant
    T: list[float]                  # temperature (K)
    P: list[float]                  # pressure (Pa)
    F: list[Dict[str, float]]       # molar flows at each point
    r: list[float]                  # reaction rate at each point
    Q: list[float] = field(default_factory=list)   # cumulative heat added to reactor up to V (J/s = W)
    limiting_species: str = ""
    config: Optional[PFRConfig] = None
    feed: Optional[Feed] = None
    reaction: Optional[Reaction] = None
    final_V: float = 0.0
    final_X: float = 0.0
    final_T: float = 0.0
    final_P: float = 0.0
    total_Q: float = 0.0            # total heat duty required for the reactor (W). Negative = heat removed
    success: bool = False
    message: str = ""

    def outlet_flows(self) -> Dict[str, float]:
        if not self.F:
            return {}
        return self.F[-1].copy()

    def outlet_concentrations(self, T: Optional[float] = None, P: Optional[float] = None) -> Dict[str, float]:
        """Approximate outlet concentrations using ideal gas at final (or provided) T, P."""
        from .core import ideal_gas_concentration
        Ft = sum(self.F[-1].values())
        T_use = T if T is not None else self.final_T
        P_use = P if P is not None else self.final_P
        return {
            sp: ideal_gas_concentration(Fi, Ft, T_use, P_use, self.config.R if self.config else 8.314)
            for sp, Fi in self.F[-1].items()
        }

    def heat_duty_summary(self) -> str:
        if not self.Q:
            return "No heat duty data (adiabatic or not computed)"
        sign = "removed" if self.total_Q < 0 else "added"
        return f"Total heat duty: {abs(self.total_Q):.2f} W ({sign}) to maintain conditions"

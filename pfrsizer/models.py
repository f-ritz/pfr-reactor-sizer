"""
Data models for PFR reactor design.

Uses reactor length (z) as the independent variable, following the
design equations from H. Scott Fogler, *Elements of Chemical Reaction Engineering*.

Key relationships:
    V(z) = A_c * z          (volume = cross-sectional area × length)
    dF_i/dz = nu_i * r * A_c  (mole balance per unit length)
    C_i = F_i / v            (concentration = molar flow / volumetric flow)
    v = v0 * (F_T/F_T0) * (T/T0) * (P0/P)  (ideal gas volumetric flow)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional, Literal, Any
import math

from .components import Component, compute_delta_H as _compute_delta_H_from_components


@dataclass
class Reaction:
    """
    Definition of a single gas-phase reaction (power-law kinetics).

    Attributes:
        stoichiometry: Dict of stoichiometric coefficients. e.g. {'A': -1, 'B': 1, 'C': 1}
                       Negative for reactants. Must include at least one negative.
        k0: Pre-exponential factor (Arrhenius) OR the constant k value if E=0 or given-k mode.
            Units depend on reaction order (e.g. s^-1 1st order, m^3/mol/s for 2nd).
        E: Activation energy (J/mol). Set E=0 for temperature-independent (constant) k.
        delta_H: Heat of reaction at reference T (J/mol). Negative = exothermic.
                 Used for adiabatic energy balance. Sign convention: heat released by rxn is negative delta_H.
        orders: Optional explicit reaction orders per species (power law) for FORWARD rate.
                If None, uses |stoich coeff| for reactants only (sensible default).
                Example: {'A': 1, 'B': 0.5}
        name: Optional label for the reaction.
        reversible: If True, rate includes reverse term using Kc.
        Kc: Concentration equilibrium constant (Kc) in (mol/m^3)^(delta_nu) units consistent
            with concentrations in mol/m^3. Used only if reversible=True.
            For delta_nu=0 it is dimensionless. Provide at the relevant T.
    """
    stoichiometry: Dict[str, float]
    k0: float
    E: float
    delta_H: float = 0.0
    orders: Optional[Dict[str, float]] = None
    name: str = "rxn"
    reversible: bool = False
    Kc: Optional[float] = None

    def __post_init__(self):
        if not self.stoichiometry:
            raise ValueError("Stoichiometry cannot be empty")
        reactants = [s for s, nu in self.stoichiometry.items() if nu < 0]
        if not reactants:
            raise ValueError("Reaction must have at least one reactant (negative stoich coeff)")

    @property
    def species(self) -> list[str]:
        return list(self.stoichiometry.keys())

    @property
    def delta_nu(self) -> float:
        """Change in total moles per mole of reaction (sum of stoichiometric coefficients)."""
        return sum(self.stoichiometry.values())

    def rate(self, C: Dict[str, float], T: float, R: float = 8.314) -> float:
        """
        Compute reaction rate r (mol / m^3 / s) given concentrations and temperature.

        For irreversible:   r = k(T) * product (C_i ** order_i)     for reactants
        For reversible:     r = k(T) * ( forward - reverse / Kc )

        k(T) = k0 * exp(-E/(R T))   ; if E=0 then k is constant (=k0)
        """
        k = self.k0 * math.exp(-self.E / (R * T)) if abs(self.E) > 1e-12 else self.k0

        orders = self.orders if self.orders is not None else {
            sp: abs(nu) for sp, nu in self.stoichiometry.items() if nu < 0
        }

        # Forward rate contribution
        fwd = 1.0
        for sp, ord_ in orders.items():
            # Only apply reactant orders for forward (if orders includes products ignore for fwd)
            if self.stoichiometry.get(sp, 0) >= 0:
                continue
            conc = max(C.get(sp, 0.0), 0.0)
            fwd *= pow(conc, ord_)

        if not self.reversible or self.Kc is None or self.Kc <= 0:
            return k * fwd

        # Reverse term: use product stoichiometric magnitudes as reverse orders (elementary-like)
        # If user wants different reverse orders they can be supplied by setting orders for products too,
        # but here we default to |nu| for products.
        rev = 1.0
        for sp, nu in self.stoichiometry.items():
            if nu <= 0:
                continue
            ord_ = abs(nu)
            conc = max(C.get(sp, 0.0), 0.0)
            rev *= pow(conc, ord_)

        return k * (fwd - rev / max(self.Kc, 1e-30))

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

    def volumetric_flow0(self, R: float = 8.314) -> float:
        """Inlet volumetric flow rate (m^3/s) from ideal gas law: v0 = F_T0 * R * T0 / P0"""
        return self.total_F0 * R * self.T0 / self.P0


@dataclass
class PFRConfig:
    """
    Configuration and operating conditions for the PFR simulation.

    The reactor is defined by its geometry (diameter or cross-sectional area)
    and length. The integration variable is reactor length z (m).

    mode:
        - 'isothermal': Temperature is held constant (T = isothermal_T or feed.T0).
                        The solver will compute the heat duty (Q) required to maintain isothermality.
        - 'adiabatic': No heat transfer. Energy balance is integrated → temperature profile.

    isothermal_T: If mode='isothermal' and this is set, the reactor runs at this fixed temperature (K).
                  Otherwise feed.T0 is used.

    pressure_model:
        - 'constant': P = P0 always ("No pressure drop")
        - 'simple_drop': Detailed pressure drop integrated along z ("Calculate detailed pressure drop").
                         Uses dP/dz = -alpha * (Ft/Ft0)*(T/T0)*(P0/P)

    alpha: Pressure drop parameter (units 1/m). See solvers for formulation. 0 = no drop.
           The "detailed" option uses this model.

    target_X: If provided and >0, integrate until limiting reactant reaches this conversion.
    max_L: Hard upper limit on reactor length to integrate (safety). m

    Geometry:
        diameter: Tube inner diameter (m). If provided, cross_sectional_area is computed.
        cross_sectional_area: A_c = pi*D^2/4 (m^2). If both diameter and area are None,
                              defaults to 0.1 m diameter.
    """
    mode: Literal["isothermal", "adiabatic"] = "isothermal"
    isothermal_T: Optional[float] = None
    pressure_model: Literal["constant", "simple_drop"] = "constant"
    alpha: float = 0.0          # Pressure drop param

    target_X: Optional[float] = None
    max_L: float = 100.0        # m (reactor length upper bound)

    # Geometry
    diameter: Optional[float] = None   # m
    cross_sectional_area: Optional[float] = None  # m^2

    R: float = 8.314            # J/mol/K  (gas constant)

    def __post_init__(self):
        # Ensure we have a cross-sectional area
        if self.cross_sectional_area is not None and self.cross_sectional_area > 0:
            pass  # already set
        elif self.diameter is not None and self.diameter > 0:
            self.cross_sectional_area = math.pi * (self.diameter / 2.0) ** 2
        else:
            # Default: 0.1 m diameter
            self.diameter = 0.1
            self.cross_sectional_area = math.pi * (self.diameter / 2.0) ** 2

    @property
    def A_c(self) -> float:
        """Cross-sectional area (m^2)."""
        if self.cross_sectional_area is not None and self.cross_sectional_area > 0:
            return self.cross_sectional_area
        if self.diameter is not None and self.diameter > 0:
            return math.pi * (self.diameter / 2.0) ** 2
        return math.pi * (0.1 / 2.0) ** 2  # default 0.1 m diameter


@dataclass
class PFRResult:
    """
    Results from a PFR integration run.

    The independent variable is reactor length z (m).
    Volume V(z) = A_c * z is also provided.
    """
    z: list[float]                  # reactor length points (m)
    V: list[float]                  # reactor volume points (m^3) = A_c * z
    X: list[float]                  # conversion of limiting reactant
    T: list[float]                  # temperature (K)
    P: list[float]                  # pressure (Pa)
    F: list[Dict[str, float]]       # molar flows at each point
    r: list[float]                  # reaction rate at each point
    Q: list[float] = field(default_factory=list)   # cumulative heat added up to z (J/s = W)
    limiting_species: str = ""
    config: Optional[PFRConfig] = None
    feed: Optional[Feed] = None
    reaction: Optional[Reaction] = None
    final_z: float = 0.0
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

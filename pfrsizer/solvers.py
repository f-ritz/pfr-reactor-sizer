"""
PFR solvers for isothermal and adiabatic gas-phase plug flow reactors.

Integrates the design equations along reactor length z (m), following
H. Scott Fogler, *Elements of Chemical Reaction Engineering*.

Mole balance (per unit length):
    dF_i / dz = nu_i * r * A_c

where A_c is the cross-sectional area (m^2).

Concentrations from ideal gas:
    C_i = F_i / v
    v = v0 * (F_T/F_T0) * (T/T0) * (P0/P)

Energy balance (adiabatic):
    dT/dz = (-r * delta_H * A_c) / sum(F_i * Cp_i)

    This is the standard form from Fogler (Elements of CRE, e.g. eq. 11-3 / 12-1):
    sum Fj Cpj * dT/dV = -r * ΔH_rx
    where ΔH_rx < 0 for exothermic (T rises with X), >0 for endothermic (T falls with X).
    For constant Cp and ΔH, this implies T(X) is linear in conversion X (the "adiabatic operating line").
    (Not S-shaped; the X_e(T) equilibrium curve may be sigmoidal/decreasing.)
    The trajectory vs length z is determined by the coupled mole balance (r varies with T, C).

Pressure drop:
    "No pressure drop" → P constant
    "Detailed pressure drop" → dP/dz = -alpha * (F_T/F_T0) * (T/T0) * (P0/P)
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from typing import Dict, List, Tuple, Optional, Callable
from dataclasses import replace

from .models import Reaction, Feed, PFRConfig, PFRResult
from .core import (
    compute_concentrations_from_v,
    inlet_volumetric_flow,
    R_DEFAULT,
)


def _get_limiting_species(reaction: Reaction, F0: Dict[str, float]) -> str:
    """Identify the limiting reactant based on stoichiometric feed requirement."""
    min_ratio = float("inf")
    limiting = None
    for sp, nu in reaction.stoichiometry.items():
        if nu >= 0:
            continue
        f0 = max(F0.get(sp, 0.0), 0.0)
        ratio = f0 / abs(nu)
        if ratio < min_ratio:
            min_ratio = ratio
            limiting = sp
    if limiting is None:
        for sp, nu in reaction.stoichiometry.items():
            if nu < 0:
                return sp
    return limiting or list(reaction.stoichiometry.keys())[0]


def _build_state_vector(
    F: Dict[str, float],
    T: float,
    P: float,
    species_order: List[str],
    integrate_T: bool,
    integrate_P: bool,
) -> np.ndarray:
    """Pack state into flat array [F0, F1, ..., T?, P?]."""
    state = [max(F.get(sp, 0.0), 0.0) for sp in species_order]
    if integrate_T:
        state.append(T)
    if integrate_P:
        state.append(P)
    return np.array(state, dtype=float)


def _unpack_state(
    y: np.ndarray,
    species_order: List[str],
    integrate_T: bool,
    integrate_P: bool,
    T_default: float,
    P_default: float,
) -> Tuple[Dict[str, float], float, float]:
    """Unpack flat state array back to F_dict, T, P."""
    n = len(species_order)
    F = {sp: float(y[i]) for i, sp in enumerate(species_order)}
    idx = n
    T = float(y[idx]) if integrate_T else T_default
    idx += 1 if integrate_T else 0
    P = float(y[idx]) if integrate_P else P_default
    return F, T, P  # do not clip here; clip for derived quantities in ode to preserve energy balance


def _make_ode(
    reaction: Reaction,
    feed: Feed,
    config: PFRConfig,
    species_order: List[str],
    integrate_T: bool,
    integrate_P: bool,
    Cp: Optional[Dict[str, float]],
) -> Callable:
    """
    Return the ODE function dy/dz = f(z, y) compatible with solve_ivp.

    The independent variable is reactor length z (m).
    All derivatives are per unit length.
    """
    R = config.R
    A_c = config.A_c
    Ft0 = feed.total_F0
    T0 = feed.T0
    P0 = feed.P0
    v0 = inlet_volumetric_flow(feed, R)  # m^3/s

    def ode(z: float, y: np.ndarray) -> np.ndarray:
        F, T, P = _unpack_state(y, species_order, integrate_T, integrate_P, feed.T0, feed.P0)
        Ft = sum(max(f, 0.0) for f in F.values())

        # Use safe values for volumetric flow and concentrations to avoid numerical issues (negative T, v=0)
        # but use the actual T for rate constant k(T) and report the true integrated T for energy balance.
        T_safe = max(T, 1.0)
        P_safe = max(P, 1.0)
        Ft_safe = max(Ft, 1e-12)

        # Volumetric flow rate at current conditions (ideal gas)
        v = v0 * (Ft / max(Ft0, 1e-12)) * (T_safe / max(T0, 1.0)) * (P0 / max(P_safe, 1.0))

        # Concentrations from molar flows and volumetric flow
        C = compute_concentrations_from_v(F, v)

        # Reaction rate (mol/m^3/s) -- use actual T for k(T)
        r = reaction.rate(C, T, R)

        # --- Mole balances: dF_i/dz = nu_i * r * A_c ---
        dF_dict = reaction.extent_rates(r)
        dF = np.array([dF_dict.get(sp, 0.0) * A_c for sp in species_order], dtype=float)

        # --- Energy balance ---
        dT = 0.0
        if integrate_T:
            if Cp:
                sum_FCp = sum(F.get(sp, 0.0) * Cp.get(sp, 0.0) for sp in species_order)
                if sum_FCp > 1e-12:
                    # dT/dz = (-r * delta_H * A_c) / sum(F_i * Cp_i)
                    # This follows directly from adiabatic energy balance (Fogler):
                    # sum(Fj Cpj) dT/dV = -r * delta_H_rx
                    # dT/dz = Ac * dT/dV
                    dT = (-r * reaction.delta_H * A_c) / sum_FCp

        # --- Pressure drop ---
        dP = 0.0
        if integrate_P and config.pressure_model == "simple_drop" and config.alpha != 0.0:
            # Simple pressure drop model (per unit length):
            # dP/dz = -alpha * (Ft / Ft0) * (T / T0) * (P0 / P)
            if P > 0:
                dP = -config.alpha * (Ft / max(Ft0, 1e-12)) * (T / max(T0, 1.0)) * (P0 / P)

        # Assemble derivative vector
        dydz = list(dF)
        if integrate_T:
            dydz.append(dT)
        if integrate_P:
            dydz.append(dP)
        return np.array(dydz, dtype=float)

    return ode


def _make_conversion_event(
    limiting_species: str,
    species_order: List[str],
    F0_lim: float,
    target_X: float,
) -> Callable:
    """Event function that triggers when X_lim = target_X."""
    lim_idx = species_order.index(limiting_species)

    def event(z, y):
        Flim = y[lim_idx]
        X = 1.0 - (Flim / max(F0_lim, 1e-12))
        return X - target_X

    event.terminal = True
    event.direction = 1.0
    return event


def _make_low_rate_event(
    reaction: Reaction,
    species_order: List[str],
    R: float,
    feed: "Feed",
    v0: float,
    T0: float,
    P0: float,
    min_r: float = 1e-8,
) -> Callable:
    """Terminal event when reaction rate becomes very small (tapers off)."""
    def event(z, y):
        n = len(species_order)
        F = {sp: max(float(y[i]), 0.0) for i, sp in enumerate(species_order)}
        T = float(y[n]) if len(y) > n else T0
        P = float(y[n+1]) if len(y) > n+1 else P0

        Ft = sum(F.values())
        if Ft <= 0 or T <= 0 or P <= 0:
            return 1.0
        v = v0 * (Ft / max(feed.total_F0, 1e-12)) * (T / max(T0, 1.0)) * (P0 / max(P, 1.0))
        C = compute_concentrations_from_v(F, v)
        rr = reaction.rate(C, T, R)
        return rr - min_r

    event.terminal = True
    event.direction = -1.0  # crossing from positive to below min_r
    return event


def _reconstruct_profiles(
    sol,
    species_order: List[str],
    integrate_T: bool,
    integrate_P: bool,
    feed: Feed,
    config: PFRConfig,
    reaction: Reaction,
    limiting: str,
    F0_lim: float,
) -> PFRResult:
    """Reconstruct all profiles from the ODE solution."""
    R = config.R
    A_c = config.A_c
    Ft0 = feed.total_F0
    T0 = feed.T0
    P0 = feed.P0
    v0 = inlet_volumetric_flow(feed, R)

    z = sol.t.tolist()

    F_profiles: List[Dict[str, float]] = []
    T_profiles: List[float] = []
    P_profiles: List[float] = []
    r_profiles: List[float] = []
    X_profiles: List[float] = []
    V_profiles: List[float] = []

    for i in range(len(z)):
        y_i = sol.y[:, i]
        F, T, P = _unpack_state(y_i, species_order, integrate_T, integrate_P, feed.T0, feed.P0)
        F_profiles.append(F)
        T_profiles.append(T)
        P_profiles.append(P)

        # Volumetric flow and concentrations (use safe T/P for derived quantities)
        Ft = sum(max(f, 0.0) for f in F.values())
        T_safe = max(T, 1.0)
        P_safe = max(P, 1.0)
        v = v0 * (Ft / max(Ft0, 1e-12)) * (T_safe / max(T0, 1.0)) * (P0 / max(P_safe, 1.0))
        C = compute_concentrations_from_v(F, v)
        r = reaction.rate(C, T, R)  # actual T for k(T)
        r_profiles.append(r)

        # Conversion
        Flim = F.get(limiting, 0.0)
        X = 1.0 - (Flim / max(F0_lim, 1e-12))
        X_profiles.append(max(min(X, 1.0), 0.0))

        # Volume
        V_profiles.append(A_c * z[i])

    success = sol.success
    final_z = z[-1] if z else 0.0
    final_V = V_profiles[-1] if V_profiles else 0.0
    final_X = X_profiles[-1] if X_profiles else 0.0
    final_T = T_profiles[-1] if T_profiles else feed.T0
    final_P = P_profiles[-1] if P_profiles else feed.P0

    msg = sol.message or (
        "Target conversion reached" if (config.target_X and final_X >= config.target_X - 1e-4)
        else "Integration completed"
    )

    return PFRResult(
        z=z,
        V=V_profiles,
        X=X_profiles,
        T=T_profiles,
        P=P_profiles,
        F=F_profiles,
        r=r_profiles,
        limiting_species=limiting,
        config=config,
        feed=feed,
        reaction=reaction,
        final_z=final_z,
        final_V=final_V,
        final_X=final_X,
        final_T=final_T,
        final_P=final_P,
        success=success,
        message=msg,
    )


def solve_pfr_isothermal(
    reaction: Reaction,
    feed: Feed,
    config: Optional[PFRConfig] = None,
) -> PFRResult:
    """
    Solve isothermal PFR (constant T = T0 or config.isothermal_T).

    Integrates along reactor length z (m).
    """
    if config is None:
        config = PFRConfig(mode="isothermal")
    else:
        config = replace(config, mode="isothermal")

    species_order = sorted(set(feed.species) | set(reaction.species))
    limiting = _get_limiting_species(reaction, feed.F0)
    F0_lim = max(feed.F0.get(limiting, 0.0), 0.0)

    integrate_P = (config.pressure_model == "simple_drop")
    integrate_T = False

    y0 = _build_state_vector(feed.F0, feed.T0, feed.P0, species_order, integrate_T, integrate_P)

    ode_fun = _make_ode(reaction, feed, config, species_order, integrate_T, integrate_P, Cp=None)

    events = []
    if config.target_X is not None and 0 < config.target_X <= 1.0:
        events.append(_make_conversion_event(limiting, species_order, F0_lim, config.target_X))

    sol = solve_ivp(
        ode_fun,
        t_span=(0.0, config.max_L),
        y0=y0,
        method="RK45",
        rtol=1e-7,
        atol=1e-9,
        events=events if events else None,
        dense_output=False,
        max_step=0.02 * config.max_L if config.max_L > 0 else 0.05,
    )

    return _reconstruct_profiles(sol, species_order, integrate_T, integrate_P, feed, config, reaction, limiting, F0_lim)


def solve_pfr_adiabatic(
    reaction: Reaction,
    feed: Feed,
    Cp: Dict[str, float],
    config: Optional[PFRConfig] = None,
) -> PFRResult:
    """
    Solve adiabatic PFR. Temperature changes according to energy balance.

    Cp: dict species -> molar heat capacity (J/mol/K). Must cover all relevant species.
    """
    if config is None:
        config = PFRConfig(mode="adiabatic")
    else:
        config = replace(config, mode="adiabatic")

    if not Cp:
        raise ValueError("Cp (heat capacity dict) is required for adiabatic mode")

    species_order = sorted(set(feed.species) | set(reaction.species))
    limiting = _get_limiting_species(reaction, feed.F0)
    F0_lim = max(feed.F0.get(limiting, 0.0), 0.0)

    integrate_P = (config.pressure_model == "simple_drop")
    integrate_T = True

    y0 = _build_state_vector(feed.F0, feed.T0, feed.P0, species_order, integrate_T, integrate_P)

    ode_fun = _make_ode(reaction, feed, config, species_order, integrate_T, integrate_P, Cp=Cp)

    v0_local = inlet_volumetric_flow(feed, config.R)

    events = []
    if config.target_X is not None and 0 < config.target_X <= 1.0:
        events.append(_make_conversion_event(limiting, species_order, F0_lim, config.target_X))

    # Add low-rate event so integration stops when reaction tapers off (prevents running
    # to max_L with negligible r, which can lead to numerical oddities in T).
    events.append(_make_low_rate_event(reaction, species_order, config.R, feed, v0_local, feed.T0, feed.P0))

    sol = solve_ivp(
        ode_fun,
        t_span=(0.0, config.max_L),
        y0=y0,
        method="RK45",
        rtol=1e-7,
        atol=1e-9,
        events=events if events else None,
        max_step=0.02 * config.max_L if config.max_L > 0 else 0.05,
    )

    return _reconstruct_profiles(sol, species_order, integrate_T, integrate_P, feed, config, reaction, limiting, F0_lim)


# =============================================================================
# Heat duty calculation (for isothermal reactors) + unified solver
# =============================================================================

def _compute_heat_duty_profiles(
    V: list[float],
    r: list[float],
    delta_H: float,
) -> tuple[list[float], float]:
    """
    Given reaction rate profile and heat of reaction, compute cumulative heat duty.

    For isothermal operation:
        q_vol = r * delta_H     (W / m³)
    Total Q(V) = ∫ q_vol dV   (W)
    Negative total_Q → heat must be removed to stay isothermal.
    """
    if not V or not r or delta_H is None:
        return [0.0] * len(V or []), 0.0

    V_arr = np.asarray(V, dtype=float)
    r_arr = np.asarray(r, dtype=float)
    q_vol = r_arr * delta_H

    if len(V_arr) < 2:
        Q = np.zeros_like(V_arr)
    else:
        dQ = np.zeros_like(V_arr)
        for i in range(1, len(V_arr)):
            try:
                dQ[i] = np.trapezoid(q_vol[i-1:i+1], V_arr[i-1:i+1])
            except AttributeError:
                dQ[i] = (q_vol[i] + q_vol[i-1]) * (V_arr[i] - V_arr[i-1]) / 2.0
        Q = np.cumsum(dQ)

    total_Q = float(Q[-1]) if len(Q) > 0 else 0.0
    return Q.tolist(), total_Q


def solve_pfr(
    reaction: Reaction,
    feed: Feed,
    Cp: Optional[Dict[str, float]] = None,
    config: Optional[PFRConfig] = None,
) -> PFRResult:
    """
    Recommended unified entry point.

    Supports:
      - Adiabatic: full T profile (requires Cp)
      - Isothermal: fixed T, automatically computes required heat duty Q (W)
        that must be supplied/removed to keep temperature constant.

    The independent variable is reactor length z (m).
    """
    if config is None:
        config = PFRConfig()

    delta_H = getattr(reaction, "delta_H", 0.0) or 0.0

    if config.mode == "adiabatic":
        if not Cp:
            raise ValueError("Cp dictionary is required for adiabatic simulations")
        result = solve_pfr_adiabatic(reaction, feed, Cp, config)
        result.Q = [0.0] * len(result.V)
        result.total_Q = 0.0
        return result

    # Isothermal
    result = solve_pfr_isothermal(reaction, feed, config)

    T_iso = config.isothermal_T if config.isothermal_T is not None else feed.T0
    result.T = [T_iso] * len(result.T) if result.T else [T_iso]
    result.final_T = T_iso

    Q_profile, total_Q = _compute_heat_duty_profiles(result.V, result.r, delta_H)
    result.Q = Q_profile
    result.total_Q = total_Q

    if abs(total_Q) > 1e-6:
        duty_type = "heat removal (cooling required)" if total_Q < 0 else "heat addition (heating required)"
        result.message = (result.message or "") + f" | Heat duty: {abs(total_Q):.2f} W ({duty_type})"
    else:
        result.message = (result.message or "") + " | Isothermal (negligible heat effect)"

    return result

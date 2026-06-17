"""
PFR solvers for isothermal and adiabatic gas-phase plug flow reactors.
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from typing import Dict, List, Tuple, Optional, Callable
from dataclasses import replace

from .models import Reaction, Feed, PFRConfig, PFRResult
from .core import compute_concentrations, R_DEFAULT


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
        # fallback to first reactant
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
    return F, max(T, 1.0), max(P, 1.0)


def _make_ode(
    reaction: Reaction,
    feed: Feed,
    config: PFRConfig,
    species_order: List[str],
    integrate_T: bool,
    integrate_P: bool,
    Cp: Optional[Dict[str, float]],
) -> Callable:
    """Return the ODE function dy/dV = f(V, y) compatible with solve_ivp."""
    R = config.R
    Ft0 = feed.total_F0
    T0 = feed.T0
    P0 = feed.P0

    def ode(V: float, y: np.ndarray) -> np.ndarray:
        F, T, P = _unpack_state(y, species_order, integrate_T, integrate_P, feed.T0, feed.P0)
        Ft = sum(max(f, 0.0) for f in F.values())

        C = compute_concentrations(F, T, P, R)
        r = reaction.rate(C, T, R)

        # Molar flow derivatives
        dF_dict = reaction.extent_rates(r)
        dF = np.array([dF_dict.get(sp, 0.0) for sp in species_order], dtype=float)

        dT = 0.0
        if integrate_T:
            # Adiabatic energy balance (constant Cp assumption)
            if Cp:
                sum_FCp = sum(F.get(sp, 0.0) * Cp.get(sp, 0.0) for sp in species_order)
                if sum_FCp > 1e-12:
                    dT = (-r * reaction.delta_H) / sum_FCp
            # If no Cp provided, dT remains 0 (degenerate case)

        dP = 0.0
        if integrate_P and config.pressure_model == "simple_drop" and config.alpha != 0.0:
            # Simple pressure drop model (volume based, approximate):
            # dP/dV = -alpha * (Ft / Ft0) * (T / T0) * (P0 / P)
            # This is a controllable approximation often used for teaching.
            if P > 0:
                dP = -config.alpha * (Ft / max(Ft0, 1e-12)) * (T / max(T0, 1.0)) * (P0 / P)

        # Assemble derivative
        dydV = list(dF)
        if integrate_T:
            dydV.append(dT)
        if integrate_P:
            dydV.append(dP)
        return np.array(dydV, dtype=float)

    return ode


def _make_conversion_event(
    limiting_species: str,
    species_order: List[str],
    F0_lim: float,
    target_X: float,
) -> Callable:
    """Event function that triggers when X_lim = target_X."""
    lim_idx = species_order.index(limiting_species)

    def event(V, y):
        Flim = y[lim_idx]
        X = 1.0 - (Flim / max(F0_lim, 1e-12))
        return X - target_X

    event.terminal = True
    event.direction = 1.0
    return event


def solve_pfr_isothermal(
    reaction: Reaction,
    feed: Feed,
    config: Optional[PFRConfig] = None,
) -> PFRResult:
    """
    Solve isothermal PFR (constant T = T0).
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
        t_span=(0.0, config.max_V),
        y0=y0,
        method="RK45",
        rtol=1e-6,
        atol=1e-8,
        events=events if events else None,
        dense_output=False,
        max_step=0.05 * config.max_V if config.max_V > 0 else 0.1,
    )

    # Reconstruct profiles
    V = sol.t.tolist()
    n_species = len(species_order)

    F_profiles: List[Dict[str, float]] = []
    T_profiles: List[float] = []
    P_profiles: List[float] = []
    r_profiles: List[float] = []
    X_profiles: List[float] = []

    for i in range(len(V)):
        y_i = sol.y[:, i]
        F, T, P = _unpack_state(y_i, species_order, integrate_T, integrate_P, feed.T0, feed.P0)
        F_profiles.append(F)
        T_profiles.append(T)
        P_profiles.append(P)

        C = compute_concentrations(F, T, P, config.R)
        r = reaction.rate(C, T, config.R)
        r_profiles.append(r)

        Flim = F.get(limiting, 0.0)
        X = 1.0 - (Flim / max(F0_lim, 1e-12))
        X_profiles.append(max(min(X, 1.0), 0.0))

    success = sol.success
    final_V = V[-1] if V else 0.0
    final_X = X_profiles[-1] if X_profiles else 0.0
    final_T = T_profiles[-1] if T_profiles else feed.T0
    final_P = P_profiles[-1] if P_profiles else feed.P0

    msg = sol.message or ("Target conversion reached" if (config.target_X and final_X >= config.target_X - 1e-4) else "Integration completed")

    return PFRResult(
        V=V,
        X=X_profiles,
        T=T_profiles,
        P=P_profiles,
        F=F_profiles,
        r=r_profiles,
        limiting_species=limiting,
        config=config,
        feed=feed,
        reaction=reaction,
        final_V=final_V,
        final_X=final_X,
        final_T=final_T,
        final_P=final_P,
        success=success,
        message=msg,
    )


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

    events = []
    if config.target_X is not None and 0 < config.target_X <= 1.0:
        events.append(_make_conversion_event(limiting, species_order, F0_lim, config.target_X))

    sol = solve_ivp(
        ode_fun,
        t_span=(0.0, config.max_V),
        y0=y0,
        method="RK45",
        rtol=1e-6,
        atol=1e-8,
        events=events if events else None,
        max_step=0.05 * config.max_V if config.max_V > 0 else 0.1,
    )

    V = sol.t.tolist()
    n_species = len(species_order)

    F_profiles: List[Dict[str, float]] = []
    T_profiles: List[float] = []
    P_profiles: List[float] = []
    r_profiles: List[float] = []
    X_profiles: List[float] = []

    for i in range(len(V)):
        y_i = sol.y[:, i]
        F, T, P = _unpack_state(y_i, species_order, integrate_T, integrate_P, feed.T0, feed.P0)
        F_profiles.append(F)
        T_profiles.append(T)
        P_profiles.append(P)

        C = compute_concentrations(F, T, P, config.R)
        r = reaction.rate(C, T, config.R)
        r_profiles.append(r)

        Flim = F.get(limiting, 0.0)
        X = 1.0 - (Flim / max(F0_lim, 1e-12))
        X_profiles.append(max(min(X, 1.0), 0.0))

    success = sol.success
    final_V = V[-1] if V else 0.0
    final_X = X_profiles[-1] if X_profiles else 0.0
    final_T = T_profiles[-1] if T_profiles else feed.T0
    final_P = P_profiles[-1] if P_profiles else feed.P0

    msg = sol.message or ("Target conversion reached" if (config.target_X and final_X >= config.target_X - 1e-4) else "Integration completed")

    return PFRResult(
        V=V,
        X=X_profiles,
        T=T_profiles,
        P=P_profiles,
        F=F_profiles,
        r=r_profiles,
        limiting_species=limiting,
        config=config,
        feed=feed,
        reaction=reaction,
        final_V=final_V,
        final_X=final_X,
        final_T=final_T,
        final_P=final_P,
        success=success,
        message=msg,
    )

"""
Lightweight chemical reaction string parser.

Supports simple forms such as:
    "A -> B"
    "A + B -> C + D"
    "2A + B -> 3C"
    "C2H4 + H2O -> C2H5OH"
    "CH4 + 2 O2 -> CO2 + 2 H2O"

It splits on '->' or '=>' or '<->' (treats reversible as forward for sizing purposes).
Coefficients are optional integers/floats before species names.

Returns stoichiometry dict: species_key -> nu (negative for reactants).
Species keys are normalized (stripped). Case is preserved as typed (user can decide keys).

This is intentionally simple — good enough for PFR design tools.
For very complex reactions, user can still supply the dict manually.
"""

from __future__ import annotations
import re
from typing import Dict, Tuple, List


ARROW_RE = re.compile(r"\s*(?:->|=>|<->|→|⇒)\s*")
TOKEN_RE = re.compile(r"([0-9]*\.?[0-9]+)?\s*([A-Za-z0-9][A-Za-z0-9_\(\)\[\]\-\+\.\s]*)")


def _parse_side(side: str) -> Dict[str, float]:
    """Parse one side of the reaction (reactants or products)."""
    side = side.strip()
    if not side:
        return {}

    # Split on '+' but respect that it's a top-level separator
    tokens = re.split(r"\s*\+\s*", side)
    stoich: Dict[str, float] = {}

    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue

        # Match optional coefficient followed by species name
        m = re.match(r"^([0-9]*\.?[0-9]+)?\s*(.+)$", tok)
        if not m:
            continue

        coeff_str, species = m.groups()
        species = species.strip()

        if not species:
            continue

        # Normalize species key (keep original case, just strip)
        key = species

        coeff = 1.0
        if coeff_str:
            try:
                coeff = float(coeff_str)
            except ValueError:
                coeff = 1.0

        if key in stoich:
            stoich[key] += coeff
        else:
            stoich[key] = coeff

    return stoich


def parse_reaction(reaction_str: str) -> Dict[str, float]:
    """
    Parse a reaction string into a stoichiometry dictionary.

    Example:
        parse_reaction("2A + B -> C + D") -> {'A': -2.0, 'B': -1.0, 'C': 1.0, 'D': 1.0}

    Returns empty dict on failure (caller should validate).
    """
    if not reaction_str or "->" not in reaction_str and "=>" not in reaction_str:
        # Try to be lenient
        pass

    # Normalize arrows
    normalized = ARROW_RE.sub(" -> ", reaction_str)
    parts = normalized.split("->", 1)
    if len(parts) != 2:
        return {}

    left, right = parts[0], parts[1]

    reactants = _parse_side(left)
    products = _parse_side(right)

    stoich: Dict[str, float] = {}

    for sp, coeff in reactants.items():
        stoich[sp] = stoich.get(sp, 0.0) - coeff

    for sp, coeff in products.items():
        stoich[sp] = stoich.get(sp, 0.0) + coeff

    # Remove zero coefficients
    stoich = {k: v for k, v in stoich.items() if abs(v) > 1e-12}

    return stoich


def get_species_list(stoich: Dict[str, float]) -> List[str]:
    return list(stoich.keys())


def pretty_stoich(stoich: Dict[str, float]) -> str:
    """Return a human readable reaction string from the dict."""
    left = []
    right = []
    for sp, nu in stoich.items():
        if nu < 0:
            c = -nu
            left.append(f"{c:g} {sp}" if c != 1 else sp)
        else:
            c = nu
            right.append(f"{c:g} {sp}" if c != 1 else sp)
    return " + ".join(left) + " -> " + " + ".join(right)

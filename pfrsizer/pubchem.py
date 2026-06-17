"""
PubChem integration for automatic lookup of chemical species properties.

Provides molecular weight, formula, and attempts to retrieve basic thermodynamic data
(heat of formation, etc.) when available. Users should always verify/override values
from primary literature for reactor design work.
"""

from __future__ import annotations
from typing import Dict, List, Optional, Any
import warnings

try:
    import pubchempy as pcp
except ImportError:
    pcp = None

# Common fallback / default Cp values (rough, 298K, J/mol/K) for common gases
# These are only used if user does not provide better data.
DEFAULT_CP = {
    "H2": 28.8,
    "O2": 29.4,
    "N2": 29.1,
    "H2O": 33.6,
    "CO": 29.1,
    "CO2": 37.1,
    "CH4": 35.7,
    "C2H4": 43.6,
    "C2H6": 52.5,
    "C3H8": 73.6,
    "NH3": 35.6,
    "HCl": 29.1,
    "Cl2": 33.9,
    "SO2": 39.9,
}

# Very rough average Cp for organics when nothing is known
DEFAULT_CP_GENERIC = 80.0


def _safe_float(val: Any) -> Optional[float]:
    try:
        if val is None:
            return None
        return float(val)
    except (ValueError, TypeError):
        return None


def _extract_heat_of_formation(compound) -> Optional[float]:
    """
    Attempt to extract ΔHf (kJ/mol or J/mol). PubChem data is patchy.
    Returns value in J/mol if found, else None.
    """
    if not hasattr(compound, "experimental_properties"):
        return None

    props = getattr(compound, "experimental_properties", []) or []
    for prop in props:
        # pubchempy returns objects with .name and .value (string often)
        name = str(getattr(prop, "name", "")).lower()
        if "heat of formation" in name or "enthalpy of formation" in name:
            val = _safe_float(getattr(prop, "value", None))
            if val is not None:
                # PubChem often reports in kJ/mol — we convert to J/mol
                # Heuristic: if |val| < 1000 treat as kJ, else J (crude but practical)
                if abs(val) < 2000:
                    val *= 1000.0
                return val
    return None


def _extract_cp(compound) -> Optional[float]:
    """Attempt to pull a Cp value (very rarely present in simple PubChem records)."""
    # PubChem structured Cp data is uncommon via this API. Return None to force user input.
    return None


def lookup_species(identifier: str, namespace: str = "name") -> Dict[str, Any]:
    """
    Lookup a chemical species via PubChem.

    Args:
        identifier: Common name, formula, or synonym (e.g. "water", "ethylene", "CH4", "ethanol")
        namespace: 'name', 'formula', 'cid', etc. (see pubchempy docs)

    Returns:
        dict with keys:
            - query
            - found (bool)
            - name, iupac_name, formula, mw (float g/mol), cid (int or None)
            - delta_hf (J/mol or None)   # attempted
            - cp (J/mol/K or None)       # attempted (often None)
            - error (str or None)
            - raw (optional pubchempy Compound for advanced use)
    """
    result = {
        "query": identifier,
        "found": False,
        "name": identifier,
        "iupac_name": None,
        "formula": None,
        "mw": None,
        "cid": None,
        "delta_hf": None,
        "cp": None,
        "error": None,
        "raw": None,
    }

    if pcp is None:
        result["error"] = "pubchempy not installed"
        return result

    try:
        compounds = pcp.get_compounds(identifier, namespace=namespace, timeout=10)
        if not compounds:
            result["error"] = "No compound found on PubChem"
            return result

        comp = compounds[0]
        result["found"] = True
        result["raw"] = comp

        result["cid"] = getattr(comp, "cid", None)
        result["name"] = getattr(comp, "iupac_name", None) or getattr(comp, "name", identifier) or identifier
        result["iupac_name"] = getattr(comp, "iupac_name", None)
        result["formula"] = getattr(comp, "molecular_formula", None)

        mw_str = getattr(comp, "molecular_weight", None)
        result["mw"] = _safe_float(mw_str)

        # Thermo (best effort)
        result["delta_hf"] = _extract_heat_of_formation(comp)
        result["cp"] = _extract_cp(comp)

        # Provide a sensible default Cp if we have a formula match in our table
        if result["cp"] is None and result["formula"]:
            # Try common keys by formula too
            for key in (result["formula"], result.get("name", "").upper()):
                if key in DEFAULT_CP:
                    result["cp"] = DEFAULT_CP[key]
                    break

        if result["cp"] is None:
            result["cp"] = DEFAULT_CP_GENERIC

    except Exception as exc:
        result["error"] = str(exc)
        warnings.warn(f"PubChem lookup failed for '{identifier}': {exc}")

    return result


def lookup_species_list(identifiers: List[str]) -> Dict[str, Dict[str, Any]]:
    """Batch lookup. Returns mapping of original identifier -> result dict."""
    out = {}
    for ident in identifiers:
        out[ident] = lookup_species(ident)
    return out


def get_mw(name_or_formula: str) -> Optional[float]:
    """Convenience: just return molecular weight (g/mol) or None."""
    res = lookup_species(name_or_formula)
    return res.get("mw") if res.get("found") else None


def suggest_default_cp(species_key: str) -> float:
    """Return a reasonable default constant Cp (J/mol/K)."""
    key = species_key.strip().upper()
    return DEFAULT_CP.get(key, DEFAULT_CP_GENERIC)

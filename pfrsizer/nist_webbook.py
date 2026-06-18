"""
NIST Chemistry WebBook integration for thermodynamic properties.

Provides lookup of heat of formation (Hf) and heat capacity (Cp) for gas-phase species
using the public NIST Webbook (https://webbook.nist.gov/).

Data is approximate / best-effort. Users MUST verify against primary literature
or databases for accurate reactor design (Fogler recommendations).

Returns values in SI: J/mol for Hf, J/mol/K for Cp.
"""

from __future__ import annotations
from typing import Dict, Any, Optional
import re
import warnings

try:
    import requests
except ImportError:
    requests = None


DEFAULT_TIMEOUT = 12


def _safe_float(val: Any) -> Optional[float]:
    if val is None:
        return None
    try:
        # Handle strings like "52.63 kJ/mol" or "-235.1"
        if isinstance(val, str):
            val = val.strip()
            # Remove common units and parentheses
            val = re.sub(r"(?i)(kJ|J|kJ/mol|J/mol|kJ mol-1|kJ/mol-1)", "", val)
            val = val.replace("±", "").replace("~", "").strip()
            # Take first number
            m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", val)
            if m:
                num = float(m.group(0))
                # If original had kJ-ish unit in string, scale
                if "kJ" in (val.lower() + " "):
                    num *= 1000.0
                return num
            return float(val)
        return float(val)
    except (ValueError, TypeError):
        return None


def _get_page(url: str) -> Optional[str]:
    if requests is None:
        return None
    try:
        headers = {"User-Agent": "PFR-Sizer/1.0 (educational reactor design tool)"}
        resp = requests.get(url, headers=headers, timeout=DEFAULT_TIMEOUT)
        if resp.status_code == 200:
            return resp.text
        return None
    except Exception:
        return None


def _extract_nist_hf(html: str) -> Optional[float]:
    """
    Extract ΔfH°gas from NIST gas phase thermochemistry section.
    Looks for patterns like:
      ΔfH°gas   -235.1 kJ/mol   or 52.63 ± ...
    Returns J/mol or None.
    """
    if not html:
        return None

    # Common markers in the thermochemistry table
    patterns = [
        r"ΔfH°gas[^<>\n]*?([-+]?\d+\.?\d*)\s*(?:kJ/mol|kJ mol|J/mol)",
        r"enthalpy of formation[^<>\n]*?gas[^<>\n]*?([-+]?\d+\.?\d*)\s*(kJ|J)",
        r"Δ<sub>f</sub>H°<sub>gas</sub>[^<]*?([-+]?\d+\.?\d*)",
        r"ΔfHgas[^0-9]*?([-+]?\d+\.?\d*)",
    ]
    for p in patterns:
        m = re.search(p, html, re.IGNORECASE | re.DOTALL)
        if m:
            val = _safe_float(m.group(1))
            if val is not None:
                # crude: many NIST are listed in kJ but we try to detect
                # The regex often catches the number; we inspect nearby text for kJ
                context = html[max(0, m.start()-30):m.end()+30].lower()
                if "kj" in context:
                    val *= 1000.0
                return val
    return None


def _extract_nist_cp(html: str) -> Optional[float]:
    """
    Try to find a constant-pressure heat capacity near 298 K or a listed Cp.
    NIST often lists "Cp" in gas phase or Shomate params.
    We take a representative value around room temp if available.
    """
    if not html:
        return None

    # Look for explicit Cp at 298.15K or similar
    # Examples: "Cp  43.56 J/mol*K  " or in tables
    patterns = [
        r"Cp[^<>\n]{0,40}?(\d+\.?\d*)\s*(?:J/mol|J mol|J K-1 mol-1)",
        r"heat capacity[^<>\n]{0,30}gas[^<>\n]*?(\d+\.?\d*)",
        r"298\.15[^<>\n]{0,50}?Cp[^<>\n]*?(\d+\.?\d*)",
        r"Cp°\s*=?\s*(\d+\.?\d*)\s*J",
    ]
    for p in patterns:
        m = re.search(p, html, re.IGNORECASE | re.DOTALL)
        if m:
            val = _safe_float(m.group(1))
            if val is not None and 1 < val < 1000:  # plausible range
                return val
    return None


def lookup_nist(species: str) -> Dict[str, Any]:
    """
    Lookup basic thermo data from NIST Chemistry Webbook.

    Args:
        species: name or formula e.g. "ethylene", "C2H4", "water", "ethanol"

    Returns dict similar to pubchem.lookup_species:
        found, name, formula, mw (sometimes), delta_hf (J/mol), cp (J/mol/K), error, source
    """
    result: Dict[str, Any] = {
        "query": species,
        "found": False,
        "name": species,
        "formula": None,
        "mw": None,
        "delta_hf": None,
        "cp": None,
        "error": None,
        "source": "nist",
    }

    if requests is None:
        result["error"] = "requests not available for NIST lookup"
        return result

    # Build search URL. NIST prefers name or formula.
    # Use the main cbook.cgi endpoint with Name param
    import urllib.parse
    q = urllib.parse.quote(species)
    url = f"https://webbook.nist.gov/cgi/cbook.cgi?Name={q}&Units=SI&cTG=on&cTC=on&cT=on"

    html = _get_page(url)
    if not html:
        # Try a broader search
        url2 = f"https://webbook.nist.gov/cgi/cbook.cgi?Name={q}&Units=SI"
        html = _get_page(url2)

    if not html or "not found" in html.lower() or "no matching" in html.lower():
        result["error"] = "No matching species found on NIST"
        return result

    result["found"] = True

    # Try to grab formula if present
    fm = re.search(r"Formula:\s*([A-Za-z0-9]+)", html, re.IGNORECASE)
    if fm:
        result["formula"] = fm.group(1)

    # Name
    nm = re.search(r"<title>([^<]+)</title>", html)
    if nm:
        result["name"] = nm.group(1).split("- NIST")[0].strip()

    hf = _extract_nist_hf(html)
    if hf is not None:
        result["delta_hf"] = hf

    cp = _extract_nist_cp(html)
    if cp is not None:
        result["cp"] = cp

    # Fallbacks for common species if extraction misses (still better than nothing)
    if result["cp"] is None:
        common_cp = {
            "H2": 28.8, "O2": 29.4, "N2": 29.1, "H2O": 33.6, "CO": 29.1,
            "CO2": 37.1, "CH4": 35.7, "C2H4": 43.6, "C2H5OH": 65.6,
            "C2H6": 52.5, "C3H8": 73.6, "NH3": 35.6,
        }
        key = (result.get("formula") or species or "").upper().replace(" ", "")
        if key in common_cp:
            result["cp"] = common_cp[key]

    if result["delta_hf"] is None and result["cp"] is None:
        # Still mark found if page loaded but we couldn't parse numbers (user can fill)
        pass

    return result


def lookup_nist_list(species_list: list[str]) -> Dict[str, Dict[str, Any]]:
    out = {}
    for s in species_list:
        out[s] = lookup_nist(s)
    return out

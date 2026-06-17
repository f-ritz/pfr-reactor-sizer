"""
Component / species data model with thermodynamic properties.

Each species in the system is represented by a Component that can be populated
from PubChem lookup + user overrides (critical for accurate design).
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Dict, Optional, Any


@dataclass
class Component:
    """
    Thermodynamic and physical properties for one chemical species.

    All energies in J/mol (or J/mol/K for Cp).
    mw in g/mol.
    """
    name: str                          # User-facing key, e.g. "A", "ethylene", "C2H4"
    formula: Optional[str] = None
    mw: Optional[float] = None         # g/mol
    Cp: Optional[float] = None         # constant Cp at average T, J/mol/K
    Hf: Optional[float] = None         # heat of formation at 298K, J/mol
    pubchem_cid: Optional[int] = None
    pubchem_name: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @property
    def has_thermo(self) -> bool:
        return self.Cp is not None or self.Hf is not None


def compute_delta_H(components: Dict[str, Component], stoichiometry: Dict[str, float]) -> Optional[float]:
    """
    Compute heat of reaction from heats of formation.

    delta_H = sum( nu_i * Hf_i )
    Returns None if any required Hf is missing.
    """
    delta = 0.0
    missing = []
    for sp, nu in stoichiometry.items():
        comp = components.get(sp)
        if comp is None or comp.Hf is None:
            missing.append(sp)
            continue
        delta += nu * comp.Hf
    if missing:
        return None
    return delta


def build_components_from_pubchem(lookup_results: Dict[str, Dict[str, Any]]) -> Dict[str, Component]:
    """
    Convert raw pubchem.lookup_species_list() results into Component objects.
    User can (and should) override Cp, Hf etc. afterwards.
    """
    comps: Dict[str, Component] = {}
    for key, data in lookup_results.items():
        if not data.get("found"):
            # Still create a stub so user can fill values manually
            comps[key] = Component(
                name=key,
                formula=data.get("formula"),
                mw=data.get("mw"),
                Cp=data.get("cp"),
            )
            continue

        comps[key] = Component(
            name=data.get("name") or key,
            formula=data.get("formula"),
            mw=data.get("mw"),
            Cp=data.get("cp"),
            Hf=data.get("delta_hf"),
            pubchem_cid=data.get("cid"),
            pubchem_name=data.get("iupac_name") or data.get("name"),
        )
    return comps

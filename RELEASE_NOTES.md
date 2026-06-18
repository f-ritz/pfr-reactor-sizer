# PFR Reactor Sizer — Release Notes

## v1.0.0 — Initial Public Release (2026)

This is the first major public release of the PFR Reactor Sizer.

### Highlights

- Full-featured Windows GUI for sizing non-catalytic gas-phase PFRs (isothermal + adiabatic)
- Reaction parser accepts arbitrary symbolic or real chemical reactions
- Two convenient example loaders:
  - "A → B (unspecified)" — classic textbook-style symbolic reaction
  - "Real chemicals" — ethylene hydration example with realistic species names
- Flexible kinetics input:
  - Arrhenius (k0 + activation energy E)
  - Given constant k (value at the operating temperature)
- Reversible reaction support with user-supplied concentration equilibrium constant Kc
- Dual-source thermodynamic lookup:
  - PubChem (molecular weight, identity)
  - NIST Chemistry Webbook (enthalpy of formation, Cp data)
- Support for additional inert/diluent species in the feed (fully accounted in flow, concentration, and energy balances)
- Conditional user interface:
  - Only relevant fields shown for isothermal vs adiabatic
  - Clear pressure drop choice: "No pressure drop" vs "Calculate detailed pressure drop"
- All calculations performed exclusively in **SI units** with prominent on-screen warnings and documentation
- In-app "Design Equations Help" window explains the source (Fogler) of every parameter and equation
- Modern Tkinter + ttk GUI with embedded interactive matplotlib plots
- Export of full profiles to CSV and plots to PNG
- Packaged as a single standalone Windows .exe (PyInstaller)
- Optional custom icon support (`icon.ico` in root at build time appears in taskbar + window)
- Comprehensive README covering usage (EXE + source), development, units, and equations

### What's New vs Prior Internal Versions

- Complete kinetics choice (given k vs Arrhenius)
- Mode-aware dynamic UI (isothermal/adiabatic)
- Proper two-option pressure model
- NIST + PubChem combined lookup
- Inert feed support polished
- Example loading buttons
- Reversible reactions + Kc
- Icon preparation for EXE
- Extensive parameter descriptions / equation reference
- SI unit enforcement + warnings throughout
- Full README rewrite
- Release notes (this file)

### Known Limitations / Notes

- Cp is treated as constant (no Cp(T) polynomials yet)
- Single reaction only (no parallel/series networks)
- Pressure drop uses the textbook simplified α model (not full Ergun for now)
- NIST and PubChem lookups are best-effort; always override thermo data with trusted values
- Reversible: user is responsible for ensuring target X < Xeq
- For "given k" mode with adiabatic operation, k is held constant (no temperature dependence)

### Installation & Running

**End users (EXE)**: Download the .exe from Releases. Double-click. Done.

**Developers**:
```bash
pip install -e .
python -m pfrsizer.gui
```

**Build EXE** (Windows):
```powershell
python build_exe.py
```

### Credits & References

Built following the design equations and methodology in:
- H. Scott Fogler — *Elements of Chemical Reaction Engineering*
- Octave Levenspiel — *Chemical Reaction Engineering*

---

Thank you for using the PFR Reactor Sizer!

Report issues on GitHub. Feedback greatly appreciated.

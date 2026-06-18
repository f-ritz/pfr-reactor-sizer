# PFR Reactor Sizer — Release Notes

## v1.1.0 — Enhanced Reversible Adiabatic Support + Polish (2026)

This release adds full support for temperature-dependent equilibrium constants (Kc(T)) for reversible reactions in both isothermal and adiabatic modes, along with a dedicated equilibrium plot, calculation fixes, improved validation, and documentation updates. The Windows EXE is ready for distribution.

### Highlights

- **Temperature-dependent Kc for reversible reactions**: Enter Kc0 (at reference temperature Tr) instead of a fixed Kc. Kc(T) is computed via the van't Hoff equation using ΔH_rx:
  ```
  Kc(T) = Kc0 × exp[ (ΔH_rx / R) × (1/T - 1/Tr) ]
  ```
  (Fully integrated into rate law, equilibrium solver X_e(T), and GUI.)
- **New Equilibrium plot**: "Equilibrium Conversion vs T (Kc)" — X_e on Y-axis, T on X-axis (only for reversible reactions with valid Kc/Kc0). Shows how equilibrium limits conversion as temperature changes.
- **Improved T vs X plot**: Now shows both simulated path and the theoretical energy balance line (T vs X). x-axis always extends to X=1.
- **Adiabatic fixes**: Corrected equilibrium curve direction for endo/exothermic cases; added safeguards against unphysical temperatures; low-rate termination for cleaner adiabatic profiles; improved ODE tolerances and step sizing.
- **Better validation**: Clear error messages if required data is missing (e.g., "Missing/invalid Cp for adiabatic...").
- **UI polish**: "Clear All Inputs" button; updated reversible UI fields for Kc0 + Tr; equilibrium plot gracefully shows message when not applicable.
- All prior features retained and refined (SI units, PubChem+NIST, inerts, conditional UI, pressure drop, examples, icon support, etc.).

### Changes & Improvements

- Reversible handling upgraded to support Kc(T) (Fogler Example 11-3 style) while keeping backward compatibility for fixed Kc.
- Equilibrium conversion X_e(T) now correctly uses temperature-dependent Kc in `compute_equilibrium_conversion` and plots.
- Updated plots in both GUI and `plot_profiles()` (standalone).
- Documentation: README and in-app "Design Equations Help" updated for new Kc(T), plots (T on x / X_e on y for equilibrium), and behavior.
- Various robustness fixes from user feedback (equilibrium calc, adiabatic tapering, axis limits, etc.).
- EXE build process remains the same (icon.ico support included).

### Known Limitations

- Still assumes constant Cp (no Cp(T) polynomials).
- Kc(T) assumes constant ΔH (no ΔCp correction).
- Single reaction only (no networks).
- Pressure drop uses simplified α model.
- Equilibrium plot assumes constant P = P0 (pressure drop effects not reflected in X_e).

### Installation & Running the EXE

Download `PFR_Reactor_Sizer.exe` from the Releases page.

Double-click to run. For reversible adiabatic cases:
1. Check "Reversible reaction".
2. Enter Kc0 and Tr (K).
3. Provide ΔH_rx (sign: negative = exothermic).
4. Set Cp values for all species.
5. Choose Adiabatic mode.
6. The X_e vs T plot will appear and the simulation will respect the temperature-dependent equilibrium.

See the in-app "Design Equations Help" and README for full details.

### Credits

Built following Fogler, *Elements of Chemical Reaction Engineering*. Thanks to all testers and feedback providers!

---

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

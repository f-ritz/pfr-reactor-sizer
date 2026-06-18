# pfr-reactor-sizer

**PFR Reactor Sizer** — A practical, self-contained tool for **sizing and analyzing non-catalytic gas-phase plug-flow reactors (PFRs)**.

All calculations are performed in **strict SI units**. See the prominent warning inside the application and the Units section below.

## Key Features

- Type **any reaction string** — symbolic (`A + B -> C`) or real chemicals (`C2H4 + H2O -> C2H5OH`, `2 NO + O2 -> 2 NO2`, etc.)
- **Two example loaders**:
  - Symbolic A → B (unspecified chemistry)
  - Real chemicals (ethylene + water → ethanol)
- **Kinetics options**:
  - Arrhenius: enter k0 + E
  - Given constant k (at operating temperature) — E forced to 0
- **Reversible reactions**: toggle + enter concentration equilibrium constant Kc
- **Thermodynamic data lookup**:
  - Automatic search of **PubChem** (MW + basic info) **and NIST Chemistry Webbook** (Hf, Cp when available)
  - Always verify/override with literature values
- **Additional inert / diluent feeds**: freely add inerts; they participate in total flow, dilution, and adiabatic energy balance
- **Isothermal vs Adiabatic** with only the relevant input fields shown for the selected mode
- **Pressure drop**:
  - Assume no pressure drop (P = constant = P0)
  - Calculate detailed pressure drop (Fogler-style α model)
- **Isothermal**: automatic heat duty Q calculation (how much heat to add/remove)
- **Adiabatic**: full temperature profile from energy balance
- Full profiles: X(z), T(z), P(z), F_i(z), r(z), V(z)
- Modern native-looking **Tkinter GUI** (ttk) + embedded matplotlib plots + CSV export
- **Single-file Windows EXE** (PyInstaller) — runs on any Windows PC without Python
- Parameter help: "Design Equations Help" button explains the origin of every variable (Fogler-based design equations)
- All units **SI** with clear warnings in the app

## Quick Start (GUI — Recommended)

```bash
pip install -r requirements.txt
python -m pfrsizer.gui
```

Or after `pip install -e .`:
```bash
pfrsizer-gui
```

**IMPORTANT (SI Units)**: Every numerical input and all internal calculations use SI units:
- Flows: mol/s
- Lengths / diameter: m
- Volume: m³
- T: K
- P: Pa (strict SI only; enter e.g. 101325 for 1 atm)
- Energies: J/mol and J/(mol·K)
- E: J/mol

The GUI displays a bright warning banner and labels reinforce the units.

### Basic GUI Workflow
1. (Recommended) Click one of the **Load Example** buttons (symbolic A→B or real chemicals).
2. Or type your reaction (e.g. `A -> B` or `C2H4 + H2O -> C2H5OH`) and click **Parse**.
3. Click **Lookup All (PubChem + NIST)** — then review/override Cp and Hf with good literature values.
4. Optionally add extra inert species with the "Add extra (inert)" control and set their feed flows.
5. Choose **Kinetics model**: Arrhenius or Given constant k.
6. For reversible reactions, check the box and supply Kc0 (at Tr) and Tr (K). ΔH_rx is also required for Kc(T). The GUI will show the X_e vs T equilibrium curve (T on x-axis, X_e on y-axis).
7. Select **Isothermal** or **Adiabatic** — only the fields needed for that mode are shown.
8. Choose pressure drop option (no drop vs calculate detailed).
9. Fill feed conditions (T0, P0), target X, diameter, safety max length.
10. Click **▶ Calculate / Size Reactor** (or Ctrl+Enter).
11. Review summary + plots (including the new T vs X graph, especially useful for adiabatic cases). Export CSV or PNG.

The GUI remains responsive during network lookups and ODE integration.

## Building a Standalone Windows EXE

On a Windows machine:

```powershell
pip install pyinstaller
python build_exe.py
```

- Produces `dist/PFR_Reactor_Sizer.exe` (single-file, windowed app).
- Copy to any Windows PC — no Python or dependencies required on the target machine.
- The EXE is large (~150–250 MB) because it bundles the complete scientific Python stack.
- **Custom icon**: Place an `icon.ico` file in the project root directory **before** running `build_exe.py`. The icon will automatically be embedded and will appear in the Windows taskbar and in the top-left corner of the application window.
- If no `icon.ico` is present, the EXE is built without a custom icon.

## Goals (original + expanded)

- Size tubular PFRs (compute required volume for a target conversion)
- Explore composition, temperature, and pressure profiles along the reactor
- Support gas-phase behavior (mole change on reaction, volume change with T and P)
- Handle both **isothermal** and **adiabatic** operation
- Simple pressure drop model (optional)
- Easy to use from Python or via command line

## Full Feature List

- **Reaction input**: free-form parser supporting integer and decimal coefficients and any species names
- **Two one-click example sets** (symbolic A→B and a realistic multi-species reaction)
- **Kinetics flexibility**:
  - Arrhenius form (k0 + activation energy E)
  - Constant-k mode (user provides the value of k at the operating temperature)
- **Reversible reactions**: enter Kc0 (at reference Tr in K) for temperature-dependent Kc via van't Hoff (Fogler ex. 11-3 style):
  ```
  Kc(T) = Kc0 × exp[(ΔH_rx / R) × (1/T - 1/Tr)]
  ```
  Includes dedicated "Equilibrium Conversion vs T" plot (T on x-axis, X_e on y-axis). Works for both isothermal and adiabatic.
- **Dual database lookup**: PubChem (MW, identity) + NIST Chemistry Webbook (enthalpies of formation, heat capacities)
- **Inert / diluent support**: add any number of extra non-reacting species to the feed; automatically included in mole balances, volume change, and energy balance
- **Operating mode**:
  - Isothermal (constant T, computed heat duty Q)
  - Adiabatic (full coupled energy balance, T profile)
  - UI dynamically shows/hides fields appropriate to the selected mode
- **Pressure options**:
  - No pressure drop (P = P0 throughout)
  - Detailed pressure drop calculation using the lumped α model (Fogler formulation)
- **Geometry**: specify diameter (m); length and volume are computed from the design equation integration
- **Design equation basis**: all equations follow standard mole and energy balances in Fogler (dF_i/dz = ν_i r A_c, energy balance, ideal-gas variable-density flow)
- **Parameter help**: in-app "Design Equations Help" window explains the origin and units of every field
- Full profile plots in GUI: X(z), T(z), T vs X (energy balance), molar flows, rate + for reversible: X_e vs T (equilibrium conversion vs temperature, T on x, X_e on y) + CSV export
- "Clear All Inputs" button to reset the form
- Validation errors for missing required data (e.g., Cp for adiabatic)
- Professional Windows GUI with live status, threading for long operations
- Completely self-contained EXE distributable

See the in-app help for the exact equations used for each variable.

## Installation (for development / running from source)

```bash
# From the project root
pip install -r requirements.txt

# Or install in editable mode (recommended for development)
pip install -e .
```

Core dependencies: `numpy`, `scipy`, `matplotlib`. `rich` is used for nicer tables in the CLI.

## Quick Start (CLI)

```bash
# Run built-in examples
pfrsizer example A_to_B_isothermal
pfrsizer example A_to_2B_adiabatic
pfrsizer example pressure_drop
```

## Quick Start (Python)

```python
from pfrsizer import Reaction, Feed, PFRConfig
from pfrsizer.solvers import solve_pfr
from pfrsizer.plot import plot_profiles

# 1. Define the reaction: A -> B   (first order)
rxn = Reaction(
    stoichiometry={"A": -1, "B": 1},
    k0=0.05,          # s^-1
    E=0.0,            # J/mol (irrelevant for isothermal)
    delta_H=0.0,      # J/mol
    orders={"A": 1},
    name="A -> B"
)

# 2. Inlet conditions
feed = Feed(
    F0={"A": 1.0, "B": 0.0},
    T0=350.0,                 # K
    P0=5 * 101325,            # Pa (strict SI; use full Pa values in GUI)
)

# 3. Configuration
cfg = PFRConfig(
    target_X=0.80,            # stop at 80% conversion
    max_L=50.0,               # safety upper bound on length (m)
    diameter=0.1,             # m
    pressure_model="constant",
)

# 4. Solve (unified entry point)
result = solve_pfr(rxn, feed, config=cfg)

print("Required volume:", round(result.final_V, 5), "m³")
print("Final conversion:", round(result.final_X, 4))

# 5. Visualize
plot_profiles(result)
```

### Adiabatic Example

```python
rxn = Reaction(
    stoichiometry={"A": -1, "B": 2},
    k0=8.5e8,
    E=65000,           # J/mol
    delta_H=-72000,    # exothermic
    orders={"A": 1},
)

feed = Feed(F0={"A": 1.0, "Inert": 0.4}, T0=330.0, P0=3*101325)
Cp = {"A": 115.0, "B": 82.0, "Inert": 29.5}   # J / mol / K

cfg = PFRConfig(mode="adiabatic", target_X=0.8, max_L=4.0, diameter=0.1)
result = solve_pfr(rxn, feed, Cp=Cp, config=cfg)

print("Outlet T:", round(result.final_T, 1), "K")
plot_profiles(result)
```

## How It Works

The program integrates the **design equation** (mole balance) along reactor **length z** (m), where volume V = A_c × z and A_c = πD²/4:

```
dFⱼ / dz = νⱼ · r · A_c
```

For gas phase, concentrations are obtained via the ideal gas law:

```
Cⱼ = yⱼ · P / (R T)
```

Volumetric flow (and thus concentrations) changes with total moles, temperature, and pressure.

**Isothermal**: T = T₀ (constant). Heat duty Q is computed by integrating the energy term.

**Adiabatic** (no heat exchange):

```
Σ(Fⱼ Cpⱼ) dT/dz = −r · ΔH_rx · A_c
```

For constant Cp and ΔH, this produces a **linear T(X)** relationship (the adiabatic operating line). Temperature vs. Conversion (T vs X) is plotted in the GUI.

Pressure drop:

- **No pressure drop** → P = P0 constant (recommended for many lab-scale / low ΔP cases)
- **Detailed pressure drop** → integrates dP/dz = −α × (F_T/F_T0) × (T/T0) × (P0/P)

α is a user-supplied lumped parameter (units 1/m). For packed beds you can derive α from the Ergun equation. The integration uses length z as the independent variable, so concentration changes from pressure are captured. The GUI includes a "T vs Conversion" plot that shows the energy balance line.

## Units — CRITICAL (Strict SI Only)

**The entire application (GUI + core) performs every calculation using SI units exclusively.**

The GUI shows a permanent warning banner and labels on every field.

| Quantity                    | Unit                  |
|-----------------------------|-----------------------|
| Molar flow F                | mol/s                 |
| Length z, L, diameter D     | m                     |
| Volume V                    | m³                    |
| Temperature T / T0          | K                     |
| Pressure P / P0             | Pa (strict SI only)   |
| Concentration C             | mol/m³                |
| Rate r                      | mol/(m³·s)            |
| ΔH_rx , Hf                  | J/mol                 |
| Cp                          | J/(mol·K)             |
| E (activation energy)       | J/mol                 |
| k0 / given k                | depends on order      |
| α (pressure drop parameter) | 1/m                   |
| Kc (reversible)             | (mol/m³)^(Δν) via Kc0 at Tr | temperature-dependent using van't Hoff with ΔH_rx |


R = 8.314 J mol⁻¹ K⁻¹ (fixed).

**P0 must be in Pa (SI)**. Example: 101325 = 1 atm, 506625 = 5 atm. The GUI input field now labels it clearly as P0 (Pa) with no automatic conversion.

Mis-entering units (e.g. typing atm values) will give physically meaningless results. All calculations are strict SI.

## Common Reaction Types

| Reaction type     | Volume change? | Adiabatic interesting? |
|-------------------|----------------|------------------------|
| A → B             | No             | Only if heat effects   |
| A → 2B            | Yes (expansion)| Strong T rise if exo   |
| 2A → B            | Contraction    | Cooling if endo        |

The code automatically accounts for mole change via total molar flow in concentrations and (if enabled) pressure drop.

## Roadmap / Future Ideas

- Multiple reactions (parallel/series networks)
- More complete pressure drop (Ergun equation with particle diameter, porosity, etc.)
- Better heat capacity handling (Cp(T) polynomials)
- Coolant temperature profiles (non-adiabatic with Ua)
- Diameter + length calculations + velocity constraints
- Streamlit web UI for interactive reactor design
- Export of results (CSV, Excel)
- More validation examples from Fogler / Levenspiel

## Project Structure

```
pfrsizer/
├── __init__.py
├── models.py       # Reaction, Feed, PFRConfig, PFRResult
├── core.py         # Ideal gas helpers
├── solvers.py      # The numerical heart (solve_ivp)
├── plot.py         # matplotlib visualization
└── cli.py          # Command line interface
examples/
tests/ (planned)
```

## Development

```bash
pip install -e ".[dev]"
# or manually: pip install pytest black ruff
```

## How to Use the EXE (No Python Required)

1. Download the latest `PFR_Reactor_Sizer.exe` from the GitHub Releases.
2. Double-click — it may take a few seconds to start the first time (large bundled Python runtime).
3. The icon (if you supplied `icon.ico` at build time) will appear in the taskbar and window.
4. Use exactly as the GUI instructions above. No installation or Python needed on the target machine.
5. All the same warnings about SI units apply.

## For Developers — Code Setup & Compilation

```bash
# Clone
git clone ...
cd pfr-reactor-sizer

# Install in editable mode (includes all runtime deps)
pip install -e .

# Or minimal
pip install -r requirements.txt
```

Run GUI from source:
```bash
python -m pfrsizer.gui
# or
python run_gui.py
```

Run CLI examples:
```bash
pfrsizer example A_to_B_isothermal
```

Build the EXE (Windows only):
```powershell
pip install pyinstaller
python build_exe.py
# Output: dist/PFR_Reactor_Sizer.exe
```

To include a custom taskbar/window icon:
- Put `icon.ico` (recommended 256×256 or multi-resolution) in the repo root.
- Run `build_exe.py` again. The build script auto-detects it.

Lint / style (optional):
```bash
pip install -e ".[dev]"
ruff check .
black .
```

## Design Equations

The program is based on the standard plug-flow reactor design equations (mole and energy balances) from Fogler, *Elements of Chemical Reaction Engineering* (primarily Chapters 1–4, 8, 11–12). The independent variable is reactor **length z (m)** (not volume directly). Volume is computed as:

```
V = A_c × z
A_c = π D² / 4
```

### 1. Mole Balance (Design Equation)
```
dF_i / dz = ν_i · r · A_c
```
- `F_i`: molar flow rate of species i (mol/s)
- `ν_i`: stoichiometric coefficient (negative for reactants)
- `r`: reaction rate (mol/m³/s)
- `A_c`: cross-sectional area (m²)

This is equivalent to the volume form `dF_i / dV = ν_i · r` via `dV = A_c dz`.

### 2. Concentrations and Volumetric Flow (Ideal Gas, Variable Density)
```
v = v0 · (F_T / F_T0) · (T / T0) · (P0 / P)
C_i = F_i / v
```
or equivalently:
```
C_i = y_i · P / (R T)     where y_i = F_i / F_T
```
- `v`: volumetric flow rate (m³/s)
- `F_T`: total molar flow rate
- `y_i`: mole fraction
- `R = 8.314` J/mol·K

### 3. Reaction Rate Law
**Arrhenius (or constant k)**:
```
k(T) = k0 · exp(−E / (R T))     (or constant k if E = 0 or "given k" mode)
```

**Irreversible power-law**:
```
r = k(T) · ∏ C_j^{order_j}     (over reactants)
```

**Reversible** (when enabled):
```
r = k(T) · (fwd − rev / Kc(T))
```
where `fwd` uses reactant terms and `rev` uses product terms (default orders = |ν|).

### 4. Equilibrium Constant Kc(T) (Reversible Reactions)
```
Kc(T) = Kc0 · exp[(ΔH_rx / R) · (1/T − 1/Tr)]
```
- Enter `Kc0` (at reference temperature `Tr`) + `ΔH_rx` (van't Hoff integration, constant ΔH assumption).
- Used to compute equilibrium conversion `X_e(T)` and net rate.
- The GUI plots **X_e vs T** (T on x-axis, X_e on y-axis) for reversible cases.

### 5. Energy Balance
**Isothermal** (T = constant):
```
Q = ∫_0^V r · ΔH_rx dV
```
(Heat duty required to hold temperature constant; positive = heat added.)

**Adiabatic** (Q = 0):
```
Σ (F_j Cp_j) dT/dz = − r · ΔH_rx · A_c
```
- For constant Cp and ΔH, this integrates to a **linear T(X)** relationship (the adiabatic operating/energy balance line):
  ```
  T(X) ≈ T0 + X · (−ΔH_rx · F_{A0} / Σ F_{j0} Cp_j)
  ```
- The GUI plots **T vs X** (shows both the simulated path and the theoretical energy balance line). x-axis always runs to X = 1.
- Inerts contribute to the heat capacity flow `Σ F_j Cp_j`.

### 6. Pressure Drop (Detailed Mode)
```
dP/dz = −α · (F_T / F_T0) · (T / T0) · (P0 / P)
```
- α is a user-supplied lumped parameter (1/m).
- "No pressure drop" keeps P = P0 constant.

### 7. Conversion
```
X = 1 − F_lim / F0_lim     (for the limiting reactant)
```

### Additional Notes
- All profiles are generated by simultaneously integrating the mole and energy (adiabatic) balances using `scipy.integrate.solve_ivp`.
- For reversible adiabatic operation, the maximum achievable conversion is limited by the intersection of the energy balance line T(X) and the equilibrium curve X_e(T).
- Assumptions: ideal gas, constant Cp (unless you override per species), single reaction, power-law kinetics.
- See the in-app **"Design Equations Help"** button for more details, variable definitions, and Fogler references.

All derivations follow standard treatments in Fogler, *Elements of Chemical Reaction Engineering* (5th ed.).

## References & Theory

- H. Scott Fogler, *Elements of Chemical Reaction Engineering*, 5th ed.
- Octave Levenspiel, *Chemical Reaction Engineering*, 3rd ed.
- Standard mole/energy balances for ideal-gas PFRs

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

Contributions and feedback welcome. This tool is intended for learning, rapid scoping, and educational use in chemical reaction engineering. Always validate numbers against primary sources and pilot data before using for real equipment design.

**Remember: SI units only — the software tells you this at every opportunity.**

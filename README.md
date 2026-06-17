# pfr-reactor-sizer

A complete tool for **designing non-catalytic gas-phase plug-flow reactors (PFRs)**.

**Key features (as requested):**
- Type **any chemical reaction** (e.g. `C2H4 + H2O -> C2H5OH` or `A + B -> C`)
- **Automatic lookup** of chemicals in the **PubChem** database (MW, formula, rough Cp)
- Supply your own kinetics from papers (k0, E, reaction orders, ΔH or Hf)
- **Isothermal** mode with **automatic heat duty calculation** (how much heat you must add or remove to stay isothermal)
- **Adiabatic** mode with full temperature profile
- Full concentration, flow, T, P, X, and heat duty **profiles**
- Modern **Windows-looking GUI** (Tkinter + ttk)
- Packageable as a single **.exe** for any Windows PC (no Python required)

## Quick Start (GUI — Recommended)

```bash
pip install -r requirements.txt
python -m pfrsizer.gui
```

Or after `pip install -e .`:
```bash
pfrsizer-gui
```

In the GUI:
1. Type a reaction string and click **Parse**.
2. Click **Lookup All in PubChem** (or per-species buttons).
3. Fill/edit **Cp** and **Hf** (J/mol/K and J/mol) with literature values — PubChem data is only a starting point.
4. Enter feed molar flows (mol/s), T0, P0.
5. Enter kinetics (k0, E, ΔH_rx).
6. Choose **Isothermal** (program will calculate required heat duty) or **Adiabatic**.
7. Set target conversion and click the big **Calculate** button.
8. View summary numbers + interactive embedded plots. Export CSV of all profiles.

The GUI stays responsive during the (sometimes slow) PubChem calls and ODE solves. Export CSV of all profiles.

## Building a Standalone Windows EXE

On a Windows machine:

```powershell
pip install pyinstaller
python build_exe.py
```

- Produces `dist/PFR_Reactor_Sizer.exe` (single file, windowed).
- Can be copied to any other Windows PC — no Python or libraries needed.
- Size is large (~150-250 MB) because it bundles the full scientific Python stack + matplotlib.
- Edit `build_exe.py` (or the generated `.spec`) to add an icon.

## Goals (original + expanded)

- Size tubular PFRs (compute required volume for a target conversion)
- Explore composition, temperature, and pressure profiles along the reactor
- Support gas-phase behavior (mole change on reaction, volume change with T and P)
- Handle both **isothermal** and **adiabatic** operation
- Simple pressure drop model (optional)
- Easy to use from Python or via command line

## Features

- Type **any reaction string** (real chemicals or symbolic)
- **PubChem lookup** for species (MW, formula, basic properties)
- User-supplied kinetics + thermo from papers
- **Isothermal** with full **heat duty calculation** (Q in/out required to hold T constant)
- **Adiabatic** with temperature profile
- Gas-phase volume change, ideal gas concentrations, optional pressure drop
- Full profiles + summary
- Windows Tkinter GUI (looks native) + embedded live plots
- Export CSV profiles
- PyInstaller-ready for standalone EXE

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
from pfrsizer.solvers import solve_pfr_isothermal, solve_pfr_adiabatic
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
    P0=5 * 101325,            # Pa (5 atm)
)

# 3. Configuration
cfg = PFRConfig(
    target_X=0.80,            # stop at 80% conversion
    max_V=50.0,               # safety upper bound (m^3)
    pressure_model="constant",
)

# 4. Solve
result = solve_pfr_isothermal(rxn, feed, cfg)

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

cfg = PFRConfig(mode="adiabatic", target_X=0.8, max_V=4.0)
result = solve_pfr_adiabatic(rxn, feed, Cp, cfg)

print("Outlet T:", round(result.final_T, 1), "K")
plot_profiles(result)
```

## How It Works

The program integrates the **design equation** (mole balance) along reactor volume V:

```
dFⱼ / dV = νⱼ · r
```

For gas phase, concentrations are obtained via the ideal gas law:

```
Cⱼ = yⱼ · P / (R T)
```

Volumetric flow (and thus concentrations) changes with total moles, temperature, and pressure.

**Isothermal**: T = T₀ (constant)

**Adiabatic** (no heat exchange):

```
Σ(Fⱼ Cpⱼ) dT/dV = -r · ΔH_rx
```

Pressure drop (when enabled):

A simple volume-based model is provided:

```
dP/dV = -α · (F_T / F_T0) · (T / T0) · (P0 / P)
```

You control the magnitude with the `alpha` parameter in `PFRConfig`.

## Units (Important)

All calculations use **consistent SI-derived units**:

| Quantity          | Unit          |
|-------------------|---------------|
| Molar flow F      | mol/s         |
| Volume V          | m³            |
| Temperature T     | K             |
| Pressure P        | Pa            |
| Concentration C   | mol/m³        |
| Rate r            | mol/(m³·s)    |
| ΔH_rx             | J/mol         |
| Cp (heat capacity)| J/(mol·K)     |
| Activation energy E | J/mol       |
| k0 (rate constant pre-factor) | depends on order |

**R = 8.314 J/mol/K** (configurable on `PFRConfig`).

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

## References & Theory

- H. Scott Fogler, *Elements of Chemical Reaction Engineering*
- Octave Levenspiel, *Chemical Reaction Engineering*
- Standard mole balance + energy balance derivations for PFRs

## License

MIT

---

Contributions and feedback welcome. This is a work-in-progress tool focused on learning and practical sizing of homogeneous gas-phase PFRs.

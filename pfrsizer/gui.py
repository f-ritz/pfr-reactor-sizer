"""
PFR Reactor Sizer - Tkinter GUI (modern Windows look via ttk)

Run with:
    python -m pfrsizer.gui
    python run_gui.py
    python pfrsizer/gui.py     # also works (direct script execution)
"""

from __future__ import annotations

import sys
from pathlib import Path

# --- Robust support for running the GUI file directly ----------------------
# This makes `python gui.py`, `python pfrsizer/gui.py`, running from inside
# the pfrsizer folder, or even from an IDE "Run" button all work, as long as
# the pfrsizer package directory is somewhere above this file.
#
# We locate the directory that contains `pfrsizer/__init__.py` and put it on
# sys.path so that `from pfrsizer.xxx import ...` succeeds.
def _ensure_pfrsizer_on_path():
    here = Path(__file__).resolve()
    for ancestor in [here.parent] + list(here.parents):
        if (ancestor / "pfrsizer" / "__init__.py").is_file():
            if str(ancestor) not in sys.path:
                sys.path.insert(0, str(ancestor))
            return
    # Fallback (should rarely be needed)
    fallback = here.parent.parent
    if str(fallback) not in sys.path:
        sys.path.insert(0, str(fallback))

_ensure_pfrsizer_on_path()
# ---------------------------------------------------------------------------

import math
import threading
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from typing import Dict, Optional, Any

import numpy as np

# Matplotlib TkAgg embedding
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Core pfrsizer (use absolute imports so direct execution always works)
from pfrsizer.models import Reaction, Feed, PFRConfig, PFRResult
from pfrsizer.solvers import solve_pfr
from pfrsizer.reaction_parser import parse_reaction, pretty_stoich
from pfrsizer.pubchem import lookup_species, lookup_species_list
from pfrsizer.nist_webbook import lookup_nist, lookup_nist_list
from pfrsizer.components import build_components_from_pubchem, Component
from pfrsizer.plot import plot_profiles


class PFRSizerApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("PFR Reactor Sizer")
        self.root.geometry("1200x820")
        self.root.minsize(1000, 700)

        # Set custom window/taskbar icon if icon.ico is available (works in dev and frozen EXE)
        self._set_window_icon()

        # State
        self.current_result: Optional[PFRResult] = None
        self._species_vars: Dict[str, Dict[str, tk.StringVar]] = {}
        self._feed_vars: Dict[str, tk.StringVar] = {}

        # Input variables
        self.reaction_str = tk.StringVar(value="A -> B")
        self.k0_var = tk.StringVar(value="0.05")
        self.E_var = tk.StringVar(value="0.0")
        self.deltaH_var = tk.StringVar(value="0.0")
        self.orders_var = tk.StringVar(value="")  # e.g. "A:1"

        self.kinetics_model_var = tk.StringVar(value="arrhenius")  # "arrhenius" or "given_k"

        self.T0_var = tk.StringVar(value="350.0")
        self.P0_var = tk.StringVar(value="101325")  # SI units: Pa (1 atm)

        self.mode_var = tk.StringVar(value="isothermal")
        self.targetX_var = tk.StringVar(value="0.80")
        self.diameter_var = tk.StringVar(value="0.10")
        self.maxL_var = tk.StringVar(value="50.0")
        self.iso_T_var = tk.StringVar(value="")  # optional override

        self.pressure_model_var = tk.StringVar(value="none")  # "none" or "detailed"
        self.alpha_var = tk.StringVar(value="0.05")

        self.reversible_var = tk.BooleanVar(value=False)
        self.Kc_var = tk.StringVar(value="")  # Kc for reversible reactions

        # UI
        self._build_ui()

        # Seed an initial reaction
        self.root.after(150, self._initial_seed)

    def _initial_seed(self):
        """Parse the default reaction so the UI is ready."""
        try:
            self.on_parse_reaction()
        except Exception:
            pass

    def _set_window_icon(self):
        """Set the application icon for the window title bar and taskbar.
        This is needed in addition to the PyInstaller EXE icon= setting.
        Bundled via datas in the build.
        """
        import os
        try:
            if getattr(sys, 'frozen', False):
                # Running inside PyInstaller bundle (onefile extracts to _MEIPASS)
                icon_path = os.path.join(getattr(sys, '_MEIPASS', ''), 'icon.ico')
            else:
                # Development: locate the project root (dir containing pfrsizer/ package)
                # using the same logic as _ensure_pfrsizer_on_path, then look for icon.ico there.
                icon_path = None
                here = Path(__file__).resolve()
                for ancestor in [here.parent] + list(here.parents):
                    if (ancestor / "pfrsizer" / "__init__.py").is_file():
                        candidate = ancestor / "icon.ico"
                        if candidate.is_file():
                            icon_path = str(candidate)
                        break
                if icon_path is None:
                    icon_path = "icon.ico"

            if icon_path and os.path.isfile(icon_path):
                self.root.iconbitmap(icon_path)
                # Improve taskbar icon on Windows (prevents generic icon)
                try:
                    import ctypes
                    appid = "pfr.reactor.sizer.v1"
                    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(appid)
                except Exception:
                    pass
        except Exception:
            # Icon is optional - don't crash the app
            pass

    def _build_ui(self):
        # Top bar
        top = ttk.Frame(self.root, padding=8)
        top.pack(fill="x")

        title = ttk.Label(
            top,
            text="PFR Reactor Sizer  •  Non-catalytic Gas-Phase Plug Flow Reactor",
            font=("Segoe UI", 14, "bold"),
        )
        title.pack(side="left")

        self.btn_run = ttk.Button(
            top,
            text="▶ Calculate / Size Reactor",
            command=self.on_calculate,
            style="Accent.TButton",
        )
        self.btn_run.pack(side="right", padx=6)

        ttk.Button(
            top,
            text="Design Equations Help",
            command=self.show_design_equations_help,
        ).pack(side="right", padx=6)

        ttk.Button(
            top,
            text="Clear All Inputs",
            command=self.clear_all,
        ).pack(side="right", padx=6)

        # Try to style the accent button (safe if theme doesn't have it)
        try:
            style = ttk.Style(self.root)
            style.configure("Accent.TButton", font=("Segoe UI", 11, "bold"))
        except Exception:
            pass

        # Prominent SI units warning banner (requirement #11)
        si_banner = ttk.Label(
            self.root,
            text="⚠ ALL CALCULATIONS USE STRICT SI UNITS: mol/s  |  m  |  m³  |  K  |  Pa (e.g. 101325 for 1 atm)  |  J/mol  |  J/(mol·K).  "
                 "P0 must be entered in Pa. All other inputs must be SI-consistent. Verify your data!",
            font=("Segoe UI", 9, "bold"),
            foreground="#b33",
            background="#fff3cd",
            padding=4,
        )
        si_banner.pack(fill="x", padx=8, pady=(0, 4))

        # Main horizontal split
        main = ttk.PanedWindow(self.root, orient="horizontal")
        main.pack(fill="both", expand=True, padx=8, pady=4)

        # ===== LEFT: Inputs =====
        left_frame = ttk.Frame(main, padding=6)
        main.add(left_frame, weight=1)

        # Reaction frame
        rxn_frame = ttk.LabelFrame(left_frame, text="Reaction", padding=8)
        rxn_frame.pack(fill="x", pady=4)

        ttk.Label(rxn_frame, text="Reaction string (e.g.  C2H4 + H2O -> C2H5OH   or   A + B -> C )").grid(
            row=0, column=0, columnspan=6, sticky="w", pady=(0, 2)
        )

        ttk.Entry(rxn_frame, textvariable=self.reaction_str, width=55).grid(
            row=1, column=0, columnspan=4, sticky="ew", pady=2
        )
        ttk.Button(rxn_frame, text="Parse", command=self.on_parse_reaction).grid(row=1, column=4, padx=2)

        # Example loaders (requirement #6)
        ttk.Button(rxn_frame, text="Load Ex: A → B (unspecified)", command=self.load_example_symbolic).grid(
            row=2, column=0, padx=2, pady=2, sticky="w"
        )
        ttk.Button(rxn_frame, text="Load Ex: Real chemicals", command=self.load_example_real).grid(
            row=2, column=1, padx=2, pady=2, sticky="w"
        )

        ttk.Button(rxn_frame, text="Lookup All (PubChem + NIST)", command=self.on_lookup_all).grid(
            row=2, column=3, columnspan=2, pady=2, sticky="w"
        )

        # Dynamic species properties area
        self.species_frame = ttk.Frame(rxn_frame)
        self.species_frame.grid(row=3, column=0, columnspan=4, sticky="ew", pady=4)

        ttk.Label(
            rxn_frame,
            text="Edit Cp and Hf from literature. PubChem gives good MW + rough defaults.",
            font=("Segoe UI", 8),
        ).grid(row=4, column=0, columnspan=4, sticky="w")

        # Feed frame
        feed_frame = ttk.LabelFrame(left_frame, text="Feed Conditions", padding=8)
        feed_frame.pack(fill="x", pady=4)

        ttk.Label(feed_frame, text="T0 (K)").grid(row=0, column=0, sticky="w")
        ttk.Entry(feed_frame, textvariable=self.T0_var, width=10).grid(row=0, column=1, padx=4)

        ttk.Label(feed_frame, text="P0 (Pa)").grid(row=0, column=2, sticky="w")
        ttk.Entry(feed_frame, textvariable=self.P0_var, width=12).grid(row=0, column=3, padx=4)

        # Dynamic feed flows row label
        self.feed_flows_label = ttk.Label(feed_frame, text="Molar feed flows F0 (mol/s) — include inerts/diluents here too (they have nu=0):")
        self.feed_flows_label.grid(row=1, column=0, columnspan=3, sticky="w", pady=(6, 2))

        # Quick add extra species (inerts, diluents, etc.) #5
        self.extra_species_var = tk.StringVar(value="")
        ttk.Label(feed_frame, text="Add extra (inert):").grid(row=1, column=3, sticky="e")
        ttk.Entry(feed_frame, textvariable=self.extra_species_var, width=12).grid(row=1, column=4, padx=2)
        ttk.Button(feed_frame, text="Add", command=self._add_extra_species).grid(row=1, column=5, padx=2)

        self.feed_frame = ttk.Frame(feed_frame)
        self.feed_frame.grid(row=2, column=0, columnspan=6, sticky="ew")

        # Kinetics frame
        kin_frame = ttk.LabelFrame(left_frame, text="Kinetics (power-law)  —  SI units: k0 or k depends on order (e.g. s⁻¹)", padding=8)
        kin_frame.pack(fill="x", pady=4)

        # Kinetics model selection (#1)
        ttk.Label(kin_frame, text="Kinetics model:").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(
            kin_frame, text="Arrhenius (k0, E)", variable=self.kinetics_model_var, value="arrhenius",
            command=self._on_kinetics_model_change
        ).grid(row=0, column=1, sticky="w")
        ttk.Radiobutton(
            kin_frame, text="Given constant k (at T)", variable=self.kinetics_model_var, value="given_k",
            command=self._on_kinetics_model_change
        ).grid(row=0, column=2, sticky="w")

        self.k_label = ttk.Label(kin_frame, text="k0 (or k)")
        self.k_label.grid(row=1, column=0, sticky="w")
        self.k0_entry = ttk.Entry(kin_frame, textvariable=self.k0_var, width=12)
        self.k0_entry.grid(row=1, column=1)

        self.e_label = ttk.Label(kin_frame, text="E (J/mol)")
        self.e_label.grid(row=1, column=2, sticky="w")
        self.e_entry = ttk.Entry(kin_frame, textvariable=self.E_var, width=12)
        self.e_entry.grid(row=1, column=3)

        # Reversible (#7)
        ttk.Checkbutton(
            kin_frame, text="Reversible reaction", variable=self.reversible_var,
            command=self._on_reversible_change
        ).grid(row=2, column=0, sticky="w", pady=(4,0))

        ttk.Label(kin_frame, text="Kc (conc. equil. const, mol/m³ units)").grid(row=2, column=1, sticky="w")
        self.kc_entry = ttk.Entry(kin_frame, textvariable=self.Kc_var, width=14)
        self.kc_entry.grid(row=2, column=2, columnspan=2, sticky="w")

        ttk.Label(kin_frame, text="ΔH_rx (J/mol)").grid(row=3, column=0, sticky="w")
        ttk.Entry(kin_frame, textvariable=self.deltaH_var, width=12).grid(row=3, column=1)

        ttk.Label(kin_frame, text="Orders (optional, A:1 B:0.5)").grid(row=3, column=2, sticky="w")
        ttk.Entry(kin_frame, textvariable=self.orders_var, width=18).grid(row=3, column=3)

        ttk.Button(kin_frame, text="Use ΔH from Hf values", command=self._update_deltaH_from_Hf_if_possible).grid(
            row=4, column=0, columnspan=2, pady=4, sticky="w"
        )

        # Initial state of fields
        self.root.after(10, self._on_kinetics_model_change)
        self.root.after(10, self._on_reversible_change)

        # Operation / Config
        op_frame = ttk.LabelFrame(left_frame, text="Operation & Reactor  (SI: K, m, Pa, mol/s)", padding=8)
        op_frame.pack(fill="x", pady=4)

        # Mode (#2) - conditional display handled via trace
        ttk.Label(op_frame, text="Mode:").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(op_frame, text="Isothermal", variable=self.mode_var, value="isothermal",
                        command=self._on_mode_change).grid(row=0, column=1, sticky="w")
        ttk.Radiobutton(op_frame, text="Adiabatic", variable=self.mode_var, value="adiabatic",
                        command=self._on_mode_change).grid(row=0, column=2, sticky="w")

        ttk.Label(op_frame, text="Target X (0-1)").grid(row=1, column=0, sticky="w")
        ttk.Entry(op_frame, textvariable=self.targetX_var, width=8).grid(row=1, column=1)

        ttk.Label(op_frame, text="Diameter (m)  [A_c=πD²/4]").grid(row=1, column=2, sticky="w")
        ttk.Entry(op_frame, textvariable=self.diameter_var, width=8).grid(row=1, column=3)

        ttk.Label(op_frame, text="Max length safety (m)").grid(row=2, column=0, sticky="w")
        ttk.Entry(op_frame, textvariable=self.maxL_var, width=8).grid(row=2, column=1)

        # These will be conditionally gridded in _on_mode_change
        self.iso_label = ttk.Label(op_frame, text="Isothermal T override (K, optional)")
        self.iso_entry = ttk.Entry(op_frame, textvariable=self.iso_T_var, width=10)

        # Pressure drop options (#3)
        ttk.Label(op_frame, text="Pressure drop:").grid(row=3, column=0, sticky="w", pady=(6,0))
        ttk.Radiobutton(op_frame, text="No pressure drop (P=P0 constant)", variable=self.pressure_model_var,
                        value="none", command=self._on_pressure_model_change).grid(row=3, column=1, columnspan=2, sticky="w", pady=(6,0))
        ttk.Radiobutton(op_frame, text="Calculate detailed pressure drop", variable=self.pressure_model_var,
                        value="detailed", command=self._on_pressure_model_change).grid(row=4, column=1, columnspan=2, sticky="w")

        self.alpha_label = ttk.Label(op_frame, text="alpha (1/m) [Fogler]")
        self.alpha_label.grid(row=5, column=0, sticky="w")
        self.alpha_entry = ttk.Entry(op_frame, textvariable=self.alpha_var, width=8)
        self.alpha_entry.grid(row=5, column=1, sticky="w")

        # Initial conditional UI setup
        self.root.after(20, self._on_mode_change)
        self.root.after(20, self._on_pressure_model_change)

        # ===== RIGHT: Results + Plots =====
        right_frame = ttk.Frame(main, padding=6)
        main.add(right_frame, weight=2)

        # Results summary
        res_frame = ttk.LabelFrame(right_frame, text="Results Summary", padding=6)
        res_frame.pack(fill="both", expand=False)

        self.summary_text = tk.Text(res_frame, height=9, wrap="word", font=("Consolas", 10))
        self.summary_text.pack(fill="both", expand=True)

        # Plots
        plot_frame = ttk.LabelFrame(right_frame, text="Profiles", padding=4)
        plot_frame.pack(fill="both", expand=True, pady=6)

        self.fig = Figure(figsize=(11, 5.5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        # Matplotlib toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()

        # Bottom action buttons
        btns = ttk.Frame(right_frame)
        btns.pack(fill="x", pady=4)

        ttk.Button(btns, text="Export CSV Profiles", command=self.on_export_csv).pack(side="left", padx=4)
        ttk.Button(btns, text="Save Plot PNG", command=self.on_save_plot).pack(side="left", padx=4)
        ttk.Button(btns, text="Clear Results", command=self.clear_results).pack(side="left", padx=4)

        # Status bar
        self.status_var = tk.StringVar(value="Ready (SI UNITS ONLY). P0 in Pa (e.g. 101325). Parse reaction or use example buttons, lookup (PubChem+NIST), set kinetics/conditions, Calculate.")
        ttk.Label(self.root, textvariable=self.status_var, relief="sunken", padding=4).pack(fill="x")

    # -----------------------
    # Reaction parsing + species UI
    # -----------------------
    def on_parse_reaction(self):
        rxn_str = self.reaction_str.get().strip()
        if not rxn_str:
            return

        stoich = parse_reaction(rxn_str)
        if not stoich:
            messagebox.showerror(
                "Parse Error",
                "Could not parse the reaction string.\nExample: 'C2H4 + H2O -> C2H5OH' or 'A + B -> C'",
            )
            return

        species = list(stoich.keys())
        self._rebuild_species_ui(species)
        self._rebuild_feed_ui(species)
        self._update_deltaH_from_Hf_if_possible()

        # Prefill a reasonable default orders string using reactants only
        if not self.orders_var.get().strip():
            reactants = [sp for sp, nu in stoich.items() if nu < 0]
            if reactants:
                self.orders_var.set(" ".join(f"{sp}:1" for sp in reactants))

    def _rebuild_species_ui(self, species: list[str]):
        # Clear previous
        for w in self.species_frame.winfo_children():
            w.destroy()
        self._species_vars.clear()

        if not species:
            return

        # Header row
        headers = ["Species", "MW (g/mol)", "Cp (J/mol·K)", "Hf (J/mol)", "Lookup"]
        for col, h in enumerate(headers):
            ttk.Label(self.species_frame, text=h, font=("Segoe UI", 9, "bold")).grid(
                row=0, column=col, padx=3, pady=1
            )

        for i, sp in enumerate(species):
            mw_var = tk.StringVar(value="")
            cp_var = tk.StringVar(value="")
            hf_var = tk.StringVar(value="")

            self._species_vars[sp] = {"mw": mw_var, "cp": cp_var, "hf": hf_var}

            ttk.Label(self.species_frame, text=sp).grid(row=i + 1, column=0, padx=3, sticky="w")
            ttk.Entry(self.species_frame, textvariable=mw_var, width=12).grid(row=i + 1, column=1, padx=2)
            ttk.Entry(self.species_frame, textvariable=cp_var, width=12).grid(row=i + 1, column=2, padx=2)
            ttk.Entry(self.species_frame, textvariable=hf_var, width=12).grid(row=i + 1, column=3, padx=2)

            ttk.Button(
                self.species_frame,
                text="Lookup",
                width=8,
                command=lambda s=sp: self._lookup_one_species(s),
            ).grid(row=i + 1, column=4, padx=4)

    def _rebuild_feed_ui(self, species: list[str]):
        # Preserve previous values only for overlapping species from previous edit of same reaction.
        # Do NOT auto-carry species/extras from completely different previous runs/reactions
        # (that was causing "species from other runs" in the results summary).
        old_values = {sp: var.get() for sp, var in self._feed_vars.items()}

        for w in self.feed_frame.winfo_children():
            w.destroy()
        self._feed_vars.clear()

        # Only the species from the current reaction/parse (extras must be re-added if desired)
        all_species = species

        cols = 0
        for i, sp in enumerate(all_species):
            col = i % 4
            row = (i // 4) * 2

            default = old_values.get(sp, "1.0" if i == 0 else "0.0")
            var = tk.StringVar(value=default)
            self._feed_vars[sp] = var

            ttk.Label(self.feed_frame, text=f"{sp}").grid(row=row, column=col * 2, sticky="e", padx=2)
            ttk.Entry(self.feed_frame, textvariable=var, width=9).grid(row=row, column=col * 2 + 1, padx=2)

    def _add_extra_species(self):
        text = self.extra_species_var.get().strip()
        if not text:
            return
        extras = [s.strip() for s in text.replace(",", " ").split() if s.strip()]
        if not extras:
            return

        # Rebuild species props UI with extras added (without clearing user edits)
        current = list(self._species_vars.keys())
        for ex in extras:
            if ex not in current:
                current.append(ex)

        self._rebuild_species_ui(current)
        self._rebuild_feed_ui(current)  # will merge
        self.extra_species_var.set("")

    def on_lookup_all(self):
        species = list(self._species_vars.keys())
        if not species:
            messagebox.showinfo("No species", "Parse a reaction first.")
            return

        self.root.config(cursor="watch")
        self.root.update()

        try:
            # PubChem first
            pc_lookups = lookup_species_list(species)
            # NIST for better thermo (Hf, Cp) when real chemicals
            nist_lookups = lookup_nist_list(species)

            for key in species:
                if key not in self._species_vars:
                    continue
                vars_ = self._species_vars[key]

                pc = pc_lookups.get(key, {})
                ni = nist_lookups.get(key, {})

                # MW prefers PubChem (more reliable)
                mw = pc.get("mw") if pc.get("mw") else ni.get("mw")
                if mw is not None:
                    vars_["mw"].set(f"{mw:.4g}")

                # Prefer NIST thermo if present (often more complete for Hf/Cp), else PubChem
                cp = ni.get("cp") if ni.get("cp") is not None else pc.get("cp")
                hf = ni.get("delta_hf") if ni.get("delta_hf") is not None else pc.get("delta_hf")

                if cp is not None:
                    vars_["cp"].set(f"{cp:.4g}")
                if hf is not None:
                    vars_["hf"].set(f"{hf:.4g}")
        finally:
            self.root.config(cursor="")
            self._update_deltaH_from_Hf_if_possible()

    def _lookup_one_species(self, species_key: str):
        # Try PubChem first for MW, then NIST for thermo (or combine)
        pc = lookup_species(species_key)
        ni = lookup_nist(species_key)

        vars_ = self._species_vars.get(species_key)
        if not vars_:
            return

        mw = pc.get("mw") if pc.get("mw") else ni.get("mw")
        if mw is not None:
            vars_["mw"].set(f"{mw:.4g}")

        cp = ni.get("cp") if ni.get("cp") is not None else pc.get("cp")
        hf = ni.get("delta_hf") if ni.get("delta_hf") is not None else pc.get("delta_hf")

        if cp is not None:
            vars_["cp"].set(f"{cp:.4g}")
        if hf is not None:
            vars_["hf"].set(f"{hf:.4g}")

        if not (pc.get("found") or ni.get("found")):
            messagebox.showwarning("Not found", f"Neither PubChem nor NIST found '{species_key}'")

        self._update_deltaH_from_Hf_if_possible()

    def _collect_components(self) -> Dict[str, Component]:
        comps: Dict[str, Component] = {}
        for sp, vars_ in self._species_vars.items():
            def _f(key: str, default: Optional[float] = None) -> Optional[float]:
                raw = vars_.get(key, tk.StringVar(value="")).get().strip()
                if not raw:
                    return default
                try:
                    return float(raw)
                except ValueError:
                    return default

            mw = _f("mw")
            cp = _f("cp", 80.0)
            hf = _f("hf")

            comps[sp] = Component(name=sp, mw=mw, Cp=cp, Hf=hf)
        return comps

    def _update_deltaH_from_Hf_if_possible(self):
        try:
            comps = self._collect_components()
            stoich = parse_reaction(self.reaction_str.get())
            if not stoich:
                return
            delta = None
            # Use the helper
            from pfrsizer.components import compute_delta_H
            delta = compute_delta_H(comps, stoich)
            if delta is not None:
                self.deltaH_var.set(f"{delta:.2f}")
        except Exception:
            pass

    # -----------------------
    # Conditional UI handlers + examples
    # -----------------------
    def _on_kinetics_model_change(self):
        """Toggle E field visibility and labels for given-k vs Arrhenius (#1)."""
        model = self.kinetics_model_var.get()
        if model == "given_k":
            self.e_label.grid_remove()
            try:
                self.e_entry.grid_remove()
            except Exception:
                pass
            self.k_label.config(text="k (constant at T)")
            # Force E=0 for constant k semantics
            self.E_var.set("0.0")
        else:
            try:
                self.e_label.grid()
                self.e_entry.grid()
            except Exception:
                self.e_label.grid(row=1, column=2, sticky="w")
                self.e_entry.grid(row=1, column=3)
            self.k_label.config(text="k0 (or k)")

    def _on_reversible_change(self):
        """Show/hide Kc entry based on reversible flag (#7)."""
        if self.reversible_var.get():
            self.kc_entry.grid()
        else:
            try:
                self.kc_entry.grid_remove()
            except Exception:
                pass
            self.Kc_var.set("")

    def _on_mode_change(self):
        """Show only relevant fields for isothermal vs adiabatic (#2)."""
        mode = self.mode_var.get()
        if mode == "isothermal":
            self.iso_label.grid(row=2, column=2, sticky="w")
            self.iso_entry.grid(row=2, column=3)
        else:
            try:
                self.iso_label.grid_remove()
                self.iso_entry.grid_remove()
            except Exception:
                pass
            self.iso_T_var.set("")

    def _on_pressure_model_change(self):
        """Show alpha only for detailed pressure drop (#3)."""
        model = self.pressure_model_var.get()
        if model == "detailed":
            try:
                self.alpha_label.grid()
                self.alpha_entry.grid()
            except Exception:
                self.alpha_label.grid(row=5, column=0, sticky="w")
                self.alpha_entry.grid(row=5, column=1, sticky="w")
        else:
            for w in (getattr(self, 'alpha_label', None), getattr(self, 'alpha_entry', None)):
                try:
                    if w: w.grid_remove()
                except Exception:
                    pass
            self.alpha_var.set("0.0")

    def load_example_symbolic(self):
        """Load simple A -> B unspecified chemicals example (#6)."""
        self.reaction_str.set("A -> B")
        self.on_parse_reaction()
        self.k0_var.set("0.05")
        self.E_var.set("0.0")
        self.kinetics_model_var.set("arrhenius")
        self._on_kinetics_model_change()
        self.deltaH_var.set("0.0")
        self.orders_var.set("A:1")
        self.reversible_var.set(False)
        self.Kc_var.set("")
        self._on_reversible_change()

        self.T0_var.set("350.0")
        self.P0_var.set("506625")  # 5 atm in Pa (SI)
        self.mode_var.set("isothermal")
        self._on_mode_change()
        self.targetX_var.set("0.80")
        self.diameter_var.set("0.10")
        self.maxL_var.set("50.0")
        self.iso_T_var.set("")
        self.pressure_model_var.set("none")
        self._on_pressure_model_change()
        self.alpha_var.set("0.0")

        # Thermodynamic properties for A and B (letter reactants DO need these for adiabatic runs
        # and heat duty calculations; they are placeholders - override with real values or Lookup):
        # Cp (J/mol·K), Hf (J/mol at 298K), MW (g/mol)
        # Example values used here:
        # A: Cp=100, Hf=0, MW=28
        # B: Cp=80, Hf=-20000, MW=28
        # deltaH=0 here (isothermal demo); for endothermic adiabatic test set deltaH>0 e.g. 40000
        for sp in list(self._species_vars.keys()):
            if sp == "A":
                self._species_vars[sp]["cp"].set("100.0")
                self._species_vars[sp]["hf"].set("0.0")
                self._species_vars[sp]["mw"].set("28.0")
            elif sp == "B":
                self._species_vars[sp]["cp"].set("80.0")
                self._species_vars[sp]["hf"].set("-20000.0")
                self._species_vars[sp]["mw"].set("28.0")
        self.status_var.set("Loaded symbolic A -> B example. All values SI. Adjust k/Cp/Hf as needed.")

    def load_example_real(self):
        """Load real chemical example: ethylene hydration to ethanol (#6)."""
        self.reaction_str.set("C2H4 + H2O -> C2H5OH")
        self.on_parse_reaction()
        # Typical values (illustrative, verify with literature before use!)
        self.k0_var.set("1.2e8")   # example pre-factor
        self.E_var.set("65000")    # J/mol ~65 kJ/mol
        self.kinetics_model_var.set("arrhenius")
        self._on_kinetics_model_change()
        self.deltaH_var.set("-46000")  # approx exothermic J/mol
        self.orders_var.set("C2H4:1 H2O:0")
        self.reversible_var.set(False)
        self.Kc_var.set("")
        self._on_reversible_change()

        self.T0_var.set("500.0")   # K ~227 C
        self.P0_var.set("3039750")  # 30 atm in Pa (SI)
        self.mode_var.set("isothermal")
        self._on_mode_change()
        self.targetX_var.set("0.6")
        self.diameter_var.set("0.08")
        self.maxL_var.set("30.0")
        self.iso_T_var.set("500")
        self.pressure_model_var.set("none")
        self._on_pressure_model_change()

        # Seed realistic Cp/Hf (will be improved by lookup)
        defaults = {
            "C2H4": {"mw": "28.05", "cp": "43.6", "hf": "52500"},
            "H2O": {"mw": "18.015", "cp": "33.6", "hf": "-241800"},
            "C2H5OH": {"mw": "46.07", "cp": "65.6", "hf": "-235100"},
        }
        for sp, vals in defaults.items():
            if sp in self._species_vars:
                self._species_vars[sp]["mw"].set(vals["mw"])
                self._species_vars[sp]["cp"].set(vals["cp"])
                self._species_vars[sp]["hf"].set(vals["hf"])

        self.status_var.set("Loaded real example C2H4 + H2O -> C2H5OH (illustrative numbers). Use Lookup + verify literature values!")

    # -----------------------
    # Input collection + solve
    # -----------------------
    def _collect_inputs(self) -> tuple[Reaction, Feed, Dict[str, float], PFRConfig]:
        stoich = parse_reaction(self.reaction_str.get().strip())
        if not stoich:
            raise ValueError("Invalid reaction string")

        k0 = float(self.k0_var.get())
        E = float(self.E_var.get())
        delta_H = float(self.deltaH_var.get())

        # Kinetics model handling (#1)
        if self.kinetics_model_var.get() == "given_k":
            # Constant k at operating temperature: force E=0 so k(T)=k0
            E = 0.0

        orders = None
        ord_text = self.orders_var.get().strip()
        if ord_text:
            orders = {}
            for pair in ord_text.replace(",", " ").split():
                if ":" not in pair:
                    continue
                k, v = pair.split(":", 1)
                orders[k.strip()] = float(v)

        # Reversible handling (#7)
        rev = self.reversible_var.get()
        kc = None
        if rev:
            kc_text = self.Kc_var.get().strip()
            if kc_text:
                try:
                    kc = float(kc_text)
                except ValueError:
                    raise ValueError("Invalid Kc value for reversible reaction")

        rxn = Reaction(
            stoichiometry=stoich,
            k0=k0,
            E=E,
            delta_H=delta_H,
            orders=orders,
            name=self.reaction_str.get().strip(),
            reversible=rev,
            Kc=kc,
        )

        # Feed
        F0: Dict[str, float] = {}
        for sp, var in self._feed_vars.items():
            F0[sp] = float(var.get())

        T0 = float(self.T0_var.get())

        # P0 must be in Pa (SI units). No automatic atm conversion.
        try:
            P0 = float(self.P0_var.get())
        except ValueError:
            P0 = 101325.0

        feed = Feed(F0=F0, T0=T0, P0=P0)

        # Components -> Cp for adiabatic + deltaH already handled
        comps = self._collect_components()
        Cp = {sp: (c.Cp if c.Cp is not None else 80.0) for sp, c in comps.items()}

        # Config
        mode = self.mode_var.get()
        target_x = float(self.targetX_var.get() or 0.0)
        diameter = float(self.diameter_var.get() or 0.1)
        max_L = float(self.maxL_var.get() or 50.0)

        iso_T = None
        if self.iso_T_var.get().strip():
            try:
                iso_T = float(self.iso_T_var.get())
            except ValueError:
                iso_T = None

        # Pressure model (#3)
        pmodel = self.pressure_model_var.get()
        pressure_model = "constant"
        alpha = 0.0
        if pmodel == "detailed":
            pressure_model = "simple_drop"  # uses the existing detailed-ish alpha model (Fogler form)
            try:
                alpha = float(self.alpha_var.get())
            except ValueError:
                alpha = 0.05

        cfg = PFRConfig(
            mode=mode,
            isothermal_T=iso_T,
            pressure_model=pressure_model,
            alpha=alpha,
            target_X=target_x if 0 < target_x <= 1.0 else None,
            max_L=max_L,
            diameter=diameter,
        )

        return rxn, feed, Cp, cfg

    def on_calculate(self):
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("1.0", "Solving PFR... (GUI stays responsive)")

        def worker():
            try:
                rxn, feed, Cp, cfg = self._collect_inputs()

                # Validation: required thermo data (especially for adiabatic / heat calcs)
                comps = self._collect_components()
                if cfg.mode == "adiabatic":
                    missing = [sp for sp, c in comps.items() if c.Cp is None or c.Cp <= 0]
                    if missing:
                        raise ValueError(
                            f"Missing/invalid Cp (J/mol·K) for adiabatic: {missing}. "
                            "Cp is REQUIRED for the energy balance dT/dz. "
                            "Enter values in the species table or click Lookup (PubChem/NIST)."
                        )
                    if abs(rxn.delta_H) < 1e-6:
                        # Not fatal (no temp change), but note it
                        pass
                # For any run using heat (isothermal Q or adiabatic), deltaH should be set (0 is allowed)
                # but warn in status if using defaults blindly? here we raise only for critical

                # Choose solve path
                if cfg.mode == "adiabatic":
                    result = solve_pfr(rxn, feed, Cp=Cp, config=cfg)
                else:
                    result = solve_pfr(rxn, feed, config=cfg)

                self.root.after(0, lambda: self._display_result(result))
            except Exception as exc:
                self.root.after(0, lambda: messagebox.showerror("Solve Error", str(exc)))

        threading.Thread(target=worker, daemon=True).start()

    def _display_result(self, result: PFRResult):
        self.current_result = result

        # Summary text
        self.summary_text.delete("1.0", "end")
        lines = []
        rxn_name = pretty_stoich(result.reaction.stoichiometry) if result.reaction else ""
        lines.append(f"Reaction: {rxn_name}")
        lines.append(f"Mode: {result.config.mode if result.config else 'unknown'}")
        lines.append(f"Required Length: {result.final_z:.6g} m")
        lines.append(f"Required Volume: {result.final_V:.6g} m³   (D={result.config.diameter if result.config else '?'} m)")
        lines.append(f"Final Conversion X: {result.final_X:.5f}")
        lines.append(f"Final T: {result.final_T:.2f} K")
        if result.config and result.config.mode == 'adiabatic' and result.final_X > 0.01 and result.feed:
            dT_dX = (result.final_T - result.feed.T0) / result.final_X
            lines.append(f"  (adiabatic dT/dX ~ {dT_dX:.1f} K ; linear in X per energy balance)")
        lines.append(f"Final P: {result.final_P:.1f} Pa  ({result.final_P/101325:.4f} atm)")
        if result.total_Q:
            duty = "removed" if result.total_Q < 0 else "added"
            lines.append(f"Heat duty (isothermal): {abs(result.total_Q):.2f} W ({duty})")
        lines.append(f"Message: {result.message}")
        lines.append("")
        lines.append("Outlet molar flows (mol/s):")
        for sp, f in result.outlet_flows().items():
            lines.append(f"  {sp}: {f:.6g}")

        self.summary_text.insert("1.0", "\n".join(lines))

        # Plots
        self._update_plots(result)

        self.status_var.set(f"Done. V = {result.final_V:.5g} m³, X = {result.final_X:.4f}")

    def _update_plots(self, result: PFRResult):
        self.fig.clear()

        z = np.array(result.z)
        if len(z) == 0:
            self.canvas.draw()
            return

        # Replicate a useful subset of plot_profiles but inline for Tk control
        # Added T vs X (key for adiabatic: should be linear per energy balance)
        nrows = 2
        ncols = 3
        axes = self.fig.subplots(nrows, ncols)
        axes = axes.flatten()
        ax_idx = 0

        # Conversion
        ax = axes[ax_idx]
        ax_idx += 1
        ax.plot(z, result.X, "b-", linewidth=2)
        ax.set_xlabel("z (m)")
        ax.set_ylabel("X")
        ax.set_title("Conversion")
        ax.grid(True, alpha=0.3)

        # Temperature (always show, useful for both modes)
        ax = axes[ax_idx]
        ax_idx += 1
        ax.plot(z, result.T, "r-", linewidth=2)
        if result.feed:
            ax.axhline(result.feed.T0, color="gray", ls="--", label=f"T0={result.feed.T0:.1f}")
            ax.legend(fontsize=8)
        ax.set_xlabel("z (m)")
        ax.set_ylabel("T (K)")
        ax.set_title("Temperature Profile")
        ax.grid(True, alpha=0.3)

        # NEW: Temperature as a function of conversion (T vs X)
        # For adiabatic: should be (approx) straight line per energy balance
        # For isothermal: horizontal line
        ax = axes[ax_idx]
        ax_idx += 1
        ax.plot(result.X, result.T, "g-", linewidth=2)
        if result.feed:
            ax.axhline(result.feed.T0, color="gray", ls="--", label=f"T0={result.feed.T0:.1f}")
            ax.legend(fontsize=8)
        ax.set_xlabel("X")
        ax.set_ylabel("T (K)")
        ax.set_title("T vs Conversion")
        ax.grid(True, alpha=0.3)

        # Molar flows
        ax = axes[ax_idx]
        ax_idx += 1
        if result.F:
            species = sorted(result.F[0].keys())
            for sp in species:
                Fi = [f.get(sp, 0.0) for f in result.F]
                ax.plot(z, Fi, linewidth=1.6, label=sp)
            ax.legend(loc="best", fontsize=8, ncol=2)
        ax.set_xlabel("z (m)")
        ax.set_ylabel("Fi (mol/s)")
        ax.set_title("Molar Flows")
        ax.grid(True, alpha=0.3)

        # Rate
        ax = axes[ax_idx]
        ax_idx += 1
        ax.plot(z, result.r, "m-", linewidth=2)
        ax.set_xlabel("z (m)")
        ax.set_ylabel("r (mol/m³/s)")
        ax.set_title("Reaction Rate")
        ax.grid(True, alpha=0.3)

        # Hide the 6th subplot (2x3 grid)
        if ax_idx < len(axes):
            axes[ax_idx].axis("off")

        self.fig.suptitle(
            f"{result.reaction.name if result.reaction else ''}  |  final V={result.final_V:.4g} m³, X={result.final_X:.4f}",
            fontsize=11,
        )
        self.fig.tight_layout(rect=[0, 0.02, 1, 0.96])
        self.canvas.draw()

    # -----------------------
    # Export / utility
    # -----------------------
    def on_export_csv(self):
        if not self.current_result:
            messagebox.showinfo("No data", "Run a calculation first.")
            return

        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv")],
        )
        if not path:
            return

        res = self.current_result
        import csv

        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            # Header
            header = ["z_m", "V_m3", "X", "T_K", "P_Pa"]
            if res.F:
                for s in sorted(res.F[0].keys()):
                    header.append(f"F_{s}")
            writer.writerow(header)

            for i in range(len(res.z)):
                row = [
                    res.z[i],
                    res.V[i],
                    res.X[i],
                    res.T[i],
                    res.P[i],
                ]
                if res.F:
                    for s in sorted(res.F[0].keys()):
                        row.append(res.F[i].get(s, 0.0))
                writer.writerow(row)

        self.status_var.set(f"Exported: {path}")

    def on_save_plot(self):
        if not self.current_result:
            messagebox.showinfo("No data", "Run a calculation first.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG", "*.png")],
        )
        if not path:
            return
        self.fig.savefig(path, dpi=150, bbox_inches="tight")
        self.status_var.set(f"Plot saved: {path}")

    def clear_results(self):
        self.current_result = None
        self.summary_text.delete("1.0", "end")
        self.fig.clear()
        self.canvas.draw()
        self.status_var.set("Cleared.")

    def clear_all(self):
        """Clear all inputs, species, feed, kinetics, options, and results."""
        # Reset core input vars to sensible defaults (SI)
        self.reaction_str.set("A -> B")
        self.k0_var.set("0.05")
        self.E_var.set("0.0")
        self.deltaH_var.set("0.0")
        self.orders_var.set("")
        self.kinetics_model_var.set("arrhenius")
        self.T0_var.set("350.0")
        self.P0_var.set("101325")
        self.mode_var.set("isothermal")
        self.targetX_var.set("0.80")
        self.diameter_var.set("0.10")
        self.maxL_var.set("50.0")
        self.iso_T_var.set("")
        self.pressure_model_var.set("none")
        self.alpha_var.set("0.05")
        self.reversible_var.set(False)
        self.Kc_var.set("")
        self.extra_species_var.set("")

        # Force no stale old feed values so rebuild uses clean defaults
        self._feed_vars.clear()

        # Rebuild UIs from the (now default) reaction - only current species
        self.on_parse_reaction()

        # Apply conditional UI states
        self._on_kinetics_model_change()
        self._on_reversible_change()
        self._on_mode_change()
        self._on_pressure_model_change()

        # Clear any previous results
        self.clear_results()
        self.status_var.set("All inputs cleared.")

    def show_design_equations_help(self):
        """Show explanations of where each variable/parameter comes from (#10)."""
        win = tk.Toplevel(self.root)
        win.title("PFR Design Equations & Parameter Reference (SI units)")
        win.geometry("820x620")

        txt = tk.Text(win, wrap="word", font=("Consolas", 9))
        txt.pack(fill="both", expand=True, padx=8, pady=8)

        help_text = """PFR REACTOR SIZER — DESIGN EQUATIONS REFERENCE
(Based primarily on H. Scott Fogler, Elements of Chemical Reaction Engineering, 5th ed. Chapters 1–4, 8, 11–12)

ALL QUANTITIES IN THE SOFTWARE ARE IN STRICT SI UNITS:
- Molar flow F: mol/s
- Reactor length z, diameter D, length L: m
- Volume V = (π D²/4) * z : m³
- Temperature T, T0: K
- Pressure P, P0: Pa   (SI only; e.g. 101325 for 1 atm, 506625 for 5 atm)
- Concentration C: mol/m³   (ideal gas: C_i = (F_i / v) or y_i * P / (R T) )
- Rate r: mol/(m³ s)
- ΔH_rx, Hf : J/mol
- Cp : J/(mol·K)
- E (activation energy): J/mol
- k0 or given k : consistent with order (s⁻¹, m³ mol⁻¹ s⁻¹, ...)

CORE MOLE BALANCE (design equation) — independent var = length z (m):
  dF_i / dz = ν_i * r * A_c     where A_c = π D² / 4   (Fogler 1-3 / 4-1 per length)
  Equivalent in volume: dF_i / dV = ν_i * r

CONCENTRATIONS (gas phase, variable volume):
  v  = v0 * (F_T / F_T0) * (T / T0) * (P0 / P)     (Fogler eq 4-23, ideal gas expansion)
  C_i = F_i / v     or     C_i = y_i * P / (R T)   with R = 8.314 J/mol/K

REACTION RATE (power law):
  k(T) = k0 * exp(−E / (R T))     [Arrhenius]
  If "given constant k" selected: k = user value, E forced = 0 (no T dependence)
  r = k(T) * ∏ C_j ^ order_j     (for irreversible)
  For reversible (when checked):
     r = k(T) * [ ∏_react (C^o)  −  (1/Kc) * ∏_prod (C^o) ]
  Kc is user-supplied concentration equilibrium constant at relevant T (units (mol/m³)^Δν )

ENERGY BALANCE:
  Isothermal mode: T = T0 (or override). Heat duty Q required is computed as:
     Q = ∫ r * ΔH_rx * A_c  d z    (integrated over volume)   (Fogler 11-3, 12)
  Adiabatic mode (no Q):
     Σ (F_j Cp_j) dT/dz = − r * ΔH_rx * A_c     (Fogler 11-1 / 12-1 simplified, constant Cp)
  For constant Cp/ΔH, this implies the adiabatic operating line T = T0 + X * (-ΔH_rx / ΣCp) is LINEAR in conversion X.
  (Not S-shaped; X_e(T) equilibrium curve may appear sigmoidal.) The profile vs z depends on kinetics.

PRESSURE DROP:
  "No pressure drop": P(z) = P0 constant
  "Detailed pressure drop": uses simplified model (Fogler 4-33 type)
     dP/dz = − α * (F_T / F_T0) * (T / T0) * (P0 / P)
  α (alpha) is a lumped parameter you fit or estimate from Ergun / pipe friction for your system.
  Resulting length and volume account for changing P affecting concentrations.

INERTS:
  Inerts/diluents are species present in feed but not in stoichiometry (ν=0).
  They contribute to total Ft, affect dilution, volume change factor, and (for adiabatic) heat capacity flow.
  Add them with the "Add extra (inert)" control or include in feed values. They are carried through all balances.

TARGET CONVERSION:
  Integration stops when X_limiting = target_X (event driven solve_ivp).
  X = 1 − (F_lim / F_lim0)

GEOMETRY:
  Diameter D (m) → A_c = π(D/2)²   → V = A_c * z
  You size for required L (or V) to reach target X.

OTHER:
  - When real species (not A/B/C) are used, PubChem + NIST lookup populate MW, rough Cp, Hf.
  - Always override Cp and Hf with literature values for your T range.
  - For reversible reactions be sure X_target < X_eq at the conditions (software does not auto-check).
  - The ODEs are integrated vs z (length) for physical meaning (velocity, residence).

See also core.py, solvers.py, models.py docstrings for equation references.
"""

        txt.insert("1.0", help_text)
        txt.config(state="disabled")

        ttk.Button(win, text="Close", command=win.destroy).pack(pady=4)

    def on_exit(self):
        self.root.destroy()


def main():
    root = tk.Tk()
    # Set icon as early as possible
    try:
        # quick early set (will be reinforced inside app)
        import os
        if getattr(sys, 'frozen', False):
            ip = os.path.join(getattr(sys, '_MEIPASS', ''), 'icon.ico')
        else:
            ip = None
            here = Path(__file__).resolve()
            for ancestor in [here.parent] + list(here.parents):
                if (ancestor / "pfrsizer" / "__init__.py").is_file():
                    cand = ancestor / "icon.ico"
                    if cand.is_file():
                        ip = str(cand)
                    break
            if not ip:
                ip = "icon.ico"
        if ip and os.path.isfile(ip):
            root.iconbitmap(ip)
    except Exception:
        pass

    app = PFRSizerApp(root)
    # Keyboard shortcut hint
    root.bind("<Control-Return>", lambda e: app.on_calculate())
    root.protocol("WM_DELETE_WINDOW", app.on_exit)
    root.mainloop()


if __name__ == "__main__":
    main()

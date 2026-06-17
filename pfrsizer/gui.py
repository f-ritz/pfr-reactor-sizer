"""
PFR Reactor Sizer - Tkinter GUI (Windows-friendly)

Main application window for entering any chemical reaction, looking up
species in PubChem, providing kinetics, choosing adiabatic or isothermal
(with automatic heat duty calculation), running the PFR simulation,
and viewing all profiles.

Designed to be packaged into a standalone Windows .exe via PyInstaller.
"""

from __future__ import annotations
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import threading
from typing import Dict, Optional

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import numpy as np

# Import core package (works whether run as module or script)
try:
    from . import (
        Reaction, Feed, PFRConfig, PFRResult,
        solve_pfr,
        lookup_species, lookup_species_list,
        parse_reaction, pretty_stoich,
        Component, build_components_from_pubchem,
    )
    from .plot import plot_profiles
except ImportError:
    # Fallback when running gui.py directly
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from pfrsizer import (
        Reaction, Feed, PFRConfig, PFRResult,
        solve_pfr,
        lookup_species, lookup_species_list,
        parse_reaction, pretty_stoich,
        Component, build_components_from_pubchem,
    )
    from pfrsizer.plot import plot_profiles


class PFRSizerApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("PFR Reactor Sizer")
        root.geometry("1150x780")
        root.minsize(980, 650)

        # Theming - makes it look more like a modern Windows app
        style = ttk.Style()
        try:
            style.theme_use("vista")   # or "winnative" on older Windows
        except tk.TclError:
            pass

        self._species_vars: Dict[str, Dict[str, tk.StringVar]] = {}  # per-species editable fields
        self._feed_vars: Dict[str, tk.StringVar] = {}
        self.current_result: Optional[PFRResult] = None
        self.current_components: Dict[str, Component] = {}

        self._build_ui()

    # ------------------------------------------------------------------
    # UI Construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        # Top control bar
        top = ttk.Frame(self.root, padding=8)
        top.pack(fill="x")

        ttk.Label(top, text="PFR Reactor Sizer  •  Non-catalytic Gas-Phase Plug Flow Reactor", font=("Segoe UI", 14, "bold")).pack(side="left")

        btn_run = ttk.Button(top, text="▶ Calculate / Size Reactor", command=self.on_calculate, style="Accent.TButton")
        btn_run.pack(side="right", padx=6)

        # Main container with left input panel and right results
        main = ttk.PanedWindow(self.root, orient="horizontal")
        main.pack(fill="both", expand=True, padx=8, pady=4)

        # LEFT: Inputs
        left_frame = ttk.Frame(main, padding=6)
        main.add(left_frame, weight=1)

        # --- Reaction section ---
        rxn_frame = ttk.LabelFrame(left_frame, text="Reaction", padding=8)
        rxn_frame.pack(fill="x", pady=4)

        ttk.Label(rxn_frame, text="Reaction string (e.g.  C2H4 + H2O -> C2H5OH   or   A + B -> C )").grid(row=0, column=0, columnspan=4, sticky="w")

        self.reaction_str = tk.StringVar(value="A -> B")
        ttk.Entry(rxn_frame, textvariable=self.reaction_str, width=55).grid(row=1, column=0, columnspan=3, sticky="ew", pady=2)

        ttk.Button(rxn_frame, text="Parse", command=self.on_parse_reaction).grid(row=1, column=3, padx=4)
        ttk.Button(rxn_frame, text="Lookup All in PubChem", command=self.on_lookup_all).grid(row=2, column=0, pady=4)

        # Species property table (dynamic rows)
        self.species_frame = ttk.Frame(rxn_frame)
        self.species_frame.grid(row=3, column=0, columnspan=4, sticky="ew", pady=4)

        ttk.Label(rxn_frame, text="Edit Cp and Hf from literature. PubChem gives good MW + rough defaults.", font=("Segoe UI", 8)).grid(row=4, column=0, columnspan=4, sticky="w")

        # --- Feed section ---
        feed_frame = ttk.LabelFrame(left_frame, text="Feed Conditions (molar flow rates)", padding=8)
        feed_frame.pack(fill="x", pady=4)

        self.feed_frame = ttk.Frame(feed_frame)
        self.feed_frame.pack(fill="x")

        ttk.Label(feed_frame, text="T0 (K)").grid(row=1, column=0, sticky="w", padx=2)
        self.T0_var = tk.StringVar(value="400.0")
        ttk.Entry(feed_frame, textvariable=self.T0_var, width=10).grid(row=1, column=1, padx=2)

        ttk.Label(feed_frame, text="P0 (Pa or atm)").grid(row=1, column=2, sticky="w", padx=2)
        self.P0_var = tk.StringVar(value="101325")
        ttk.Entry(feed_frame, textvariable=self.P0_var, width=12).grid(row=1, column=3, padx=2)

        # --- Kinetics ---
        kin_frame = ttk.LabelFrame(left_frame, text="Kinetics (power-law, Arrhenius)", padding=8)
        kin_frame.pack(fill="x", pady=4)

        ttk.Label(kin_frame, text="k0").grid(row=0, column=0)
        self.k0_var = tk.StringVar(value="0.05")
        ttk.Entry(kin_frame, textvariable=self.k0_var, width=12).grid(row=0, column=1)

        ttk.Label(kin_frame, text="E (J/mol)").grid(row=0, column=2)
        self.E_var = tk.StringVar(value="0")
        ttk.Entry(kin_frame, textvariable=self.E_var, width=12).grid(row=0, column=3)

        ttk.Label(kin_frame, text="ΔH_rx (J/mol)   [auto if Hf available]").grid(row=1, column=0, columnspan=2, sticky="w")
        self.deltaH_var = tk.StringVar(value="-45000")
        ttk.Entry(kin_frame, textvariable=self.deltaH_var, width=14).grid(row=1, column=2, columnspan=2, sticky="w")

        ttk.Label(kin_frame, text="Orders (e.g. A:1  or leave blank for |stoich|)").grid(row=2, column=0, columnspan=4, sticky="w")
        self.orders_var = tk.StringVar(value="")
        ttk.Entry(kin_frame, textvariable=self.orders_var, width=50).grid(row=3, column=0, columnspan=4, sticky="ew")

        # --- Operating Conditions ---
        op_frame = ttk.LabelFrame(left_frame, text="Operating Conditions", padding=8)
        op_frame.pack(fill="x", pady=4)

        self.mode_var = tk.StringVar(value="isothermal")
        ttk.Radiobutton(op_frame, text="Isothermal (T constant — program calculates required heat duty)", variable=self.mode_var, value="isothermal").pack(anchor="w")
        ttk.Radiobutton(op_frame, text="Adiabatic (no heat transfer — temperature will change)", variable=self.mode_var, value="adiabatic").pack(anchor="w")

        ttk.Label(op_frame, text="Isothermal target T (K, blank = use T0)").grid(row=2, column=0, sticky="w")
        self.iso_T_var = tk.StringVar(value="")
        ttk.Entry(op_frame, textvariable=self.iso_T_var, width=12).grid(row=2, column=1, sticky="w")

        ttk.Label(op_frame, text="Target conversion X (0-1)").grid(row=3, column=0, sticky="w")
        self.targetX_var = tk.StringVar(value="0.80")
        ttk.Entry(op_frame, textvariable=self.targetX_var, width=10).grid(row=3, column=1)

        ttk.Label(op_frame, text="Max volume safety limit (m³)").grid(row=3, column=2, sticky="w")
        self.maxV_var = tk.StringVar(value="50")
        ttk.Entry(op_frame, textvariable=self.maxV_var, width=10).grid(row=3, column=3)

        # Pressure drop (simple)
        self.pressure_drop_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(op_frame, text="Include simple pressure drop (alpha)", variable=self.pressure_drop_var).grid(row=4, column=0, columnspan=2, sticky="w")
        ttk.Label(op_frame, text="alpha").grid(row=4, column=2)
        self.alpha_var = tk.StringVar(value="0.0")
        ttk.Entry(op_frame, textvariable=self.alpha_var, width=10).grid(row=4, column=3)

        # RIGHT: Results + Plots
        right_frame = ttk.Frame(main, padding=6)
        main.add(right_frame, weight=2)

        # Summary box
        sum_frame = ttk.LabelFrame(right_frame, text="Results Summary", padding=8)
        sum_frame.pack(fill="x", pady=4)

        self.summary_text = tk.Text(sum_frame, height=7, width=70, font=("Consolas", 9))
        self.summary_text.pack(fill="x")
        self.summary_text.insert("1.0", "Enter a reaction, feed, kinetics and conditions, then click Calculate.")

        # Plot area
        plot_frame = ttk.LabelFrame(right_frame, text="Profiles (embedded Matplotlib)", padding=4)
        plot_frame.pack(fill="both", expand=True, pady=4)

        self.plot_notebook = ttk.Notebook(plot_frame)
        self.plot_notebook.pack(fill="both", expand=True)

        # We will create dynamic tabs for plots
        self.plot_frames: Dict[str, tk.Frame] = {}
        self.canvases: Dict[str, FigureCanvasTkAgg] = {}

        # Bottom buttons
        bottom = ttk.Frame(right_frame)
        bottom.pack(fill="x", pady=4)

        ttk.Button(bottom, text="Export Profiles (CSV)", command=self.on_export_csv).pack(side="left", padx=4)
        ttk.Button(bottom, text="Save Current Plot", command=self.on_save_plot).pack(side="left", padx=4)
        ttk.Button(bottom, text="Clear", command=self.clear_results).pack(side="left", padx=4)

        # Initial species from default reaction
        self.root.after(150, self.on_parse_reaction)

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------
    def on_parse_reaction(self):
        rxn_str = self.reaction_str.get().strip()
        if not rxn_str:
            return

        stoich = parse_reaction(rxn_str)
        if not stoich:
            messagebox.showerror("Parse Error", "Could not parse the reaction string.\nExample: 'C2H4 + H2O -> C2H5OH' or 'A + B -> C'")
            return

        self._rebuild_species_ui(list(stoich.keys()))
        self._rebuild_feed_ui(list(stoich.keys()))

        # Try to auto-fill deltaH if user has Hf values later
        self._update_deltaH_from_Hf_if_possible()

    def on_lookup_all(self):
        """Lookup every species via PubChem and populate MW, Cp, Hf where possible."""
        species = list(self._species_vars.keys())
        if not species:
            messagebox.showinfo("No species", "Parse a reaction first.")
            return

        self.root.config(cursor="watch")
        self.root.update()

        try:
            lookups = lookup_species_list(species)
            comps = build_components_from_pubchem(lookups)

            for key, comp in comps.items():
                if key in self._species_vars:
                    vars_ = self._species_vars[key]
                    if comp.mw is not None:
                        vars_["mw"].set(f"{comp.mw:.4g}")
                    if comp.Cp is not None:
                        vars_["cp"].set(f"{comp.Cp:.4g}")
                    if comp.Hf is not None:
                        vars_["hf"].set(f"{comp.Hf:.4g}")
                    if comp.formula:
                        vars_["formula"].set(comp.formula)

            self.current_components = comps
            self._update_deltaH_from_Hf_if_possible()
            messagebox.showinfo("PubChem", "Lookup complete. Review and override Cp / Hf with literature values!")
        except Exception as e:
            messagebox.showerror("PubChem Error", str(e))
        finally:
            self.root.config(cursor="")

    def _rebuild_species_ui(self, species: list[str]):
        for w in self.species_frame.winfo_children():
            w.destroy()
        self._species_vars.clear()

        headers = ["Species", "Formula", "MW (g/mol)", "Cp (J/mol·K)", "Hf (J/mol)", "Lookup"]
        for col, h in enumerate(headers):
            ttk.Label(self.species_frame, text=h, font=("Segoe UI", 9, "bold")).grid(row=0, column=col, padx=3, pady=1)

        for row, sp in enumerate(species, start=1):
            mw_var = tk.StringVar(value="")
            cp_var = tk.StringVar(value="80")
            hf_var = tk.StringVar(value="")
            formula_var = tk.StringVar(value="")

            self._species_vars[sp] = {
                "mw": mw_var, "cp": cp_var, "hf": hf_var, "formula": formula_var
            }

            ttk.Label(self.species_frame, text=sp).grid(row=row, column=0, padx=3)
            ttk.Entry(self.species_frame, textvariable=formula_var, width=10).grid(row=row, column=1)
            ttk.Entry(self.species_frame, textvariable=mw_var, width=10).grid(row=row, column=2)
            ttk.Entry(self.species_frame, textvariable=cp_var, width=10).grid(row=row, column=3)
            ttk.Entry(self.species_frame, textvariable=hf_var, width=12).grid(row=row, column=4)

            btn = ttk.Button(self.species_frame, text="PubChem", width=8,
                             command=lambda s=sp: self._lookup_one_species(s))
            btn.grid(row=row, column=5, padx=2)

        self._update_deltaH_from_Hf_if_possible()

    def _lookup_one_species(self, species_key: str):
        res = lookup_species(species_key)
        if not res.get("found"):
            messagebox.showwarning("Not found", f"PubChem did not find '{species_key}'")
            return

        vars_ = self._species_vars.get(species_key)
        if vars_:
            if res.get("mw"):
                vars_["mw"].set(f"{res['mw']:.4g}")
            if res.get("cp"):
                vars_["cp"].set(f"{res['cp']:.4g}")
            if res.get("delta_hf"):
                vars_["hf"].set(f"{res['delta_hf']:.4g}")
            if res.get("formula"):
                vars_["formula"].set(res["formula"])

        self._update_deltaH_from_Hf_if_possible()

    def _rebuild_feed_ui(self, species: list[str]):
        for w in self.feed_frame.winfo_children():
            w.destroy()
        self._feed_vars.clear()

        ttk.Label(self.feed_frame, text="Species", font=("Segoe UI", 9, "bold")).grid(row=0, column=0)
        ttk.Label(self.feed_frame, text="F0 (mol/s)", font=("Segoe UI", 9, "bold")).grid(row=0, column=1)

        for i, sp in enumerate(species, start=1):
            ttk.Label(self.feed_frame, text=sp).grid(row=i, column=0, sticky="w", padx=4)
            var = tk.StringVar(value="1.0" if i == 1 else "0.0")
            self._feed_vars[sp] = var
            ttk.Entry(self.feed_frame, textvariable=var, width=12).grid(row=i, column=1, padx=4)

    def _collect_components(self) -> Dict[str, Component]:
        comps: Dict[str, Component] = {}
        for sp, vars_ in self._species_vars.items():
            try:
                mw = float(vars_["mw"].get()) if vars_["mw"].get().strip() else None
            except ValueError:
                mw = None
            try:
                cp = float(vars_["cp"].get()) if vars_["cp"].get().strip() else 80.0
            except ValueError:
                cp = 80.0
            try:
                hf = float(vars_["hf"].get()) if vars_["hf"].get().strip() else None
            except ValueError:
                hf = None

            comps[sp] = Component(name=sp, mw=mw, Cp=cp, Hf=hf)
        return comps

    def _update_deltaH_from_Hf_if_possible(self):
        comps = self._collect_components()
        stoich = parse_reaction(self.reaction_str.get())
        if not stoich:
            return
        rxn = Reaction(stoichiometry=stoich, k0=1, E=0, delta_H=0)  # dummy
        delta = rxn.compute_delta_H_from_Hf(comps)
        if delta is not None:
            self.deltaH_var.set(f"{delta:.4g}")

    def _collect_inputs(self) -> tuple[Reaction, Feed, Dict[str, float], PFRConfig]:
        # Reaction
        stoich = parse_reaction(self.reaction_str.get())
        if not stoich:
            raise ValueError("Invalid reaction string")

        try:
            k0 = float(self.k0_var.get())
            E = float(self.E_var.get())
            delta_H = float(self.deltaH_var.get())
        except ValueError:
            raise ValueError("k0, E and ΔH must be numbers")

        orders = None
        ord_text = self.orders_var.get().strip()
        if ord_text:
            orders = {}
            for pair in ord_text.split():
                if ":" in pair:
                    k, v = pair.split(":", 1)
                    orders[k.strip()] = float(v)

        rxn = Reaction(
            stoichiometry=stoich,
            k0=k0,
            E=E,
            delta_H=delta_H,
            orders=orders,
            name=self.reaction_str.get()
        )

        # Feed
        F0 = {}
        for sp, var in self._feed_vars.items():
            F0[sp] = float(var.get())

        try:
            T0 = float(self.T0_var.get())
            P0_str = self.P0_var.get().strip()
            P0 = float(P0_str) if P0_str.replace(".", "", 1).isdigit() else float(P0_str) * 101325
        except Exception:
            T0, P0 = 400.0, 101325.0

        feed = Feed(F0=F0, T0=T0, P0=P0)

        # Cp from current species table
        comps = self._collect_components()
        Cp = {sp: (c.Cp or 80.0) for sp, c in comps.items()}

        # Config
        mode = self.mode_var.get()
        try:
            target_x = float(self.targetX_var.get())
        except:
            target_x = None
        try:
            max_v = float(self.maxV_var.get())
        except:
            max_v = 50.0

        iso_T = None
        if self.iso_T_var.get().strip():
            try:
                iso_T = float(self.iso_T_var.get())
            except:
                pass

        alpha = 0.0
        if self.pressure_drop_var.get():
            try:
                alpha = float(self.alpha_var.get())
            except:
                alpha = 0.05

        cfg = PFRConfig(
            mode=mode,
            isothermal_T=iso_T,
            target_X=target_x if 0 < (target_x or 0) < 1 else None,
            max_V=max_v,
            pressure_model="simple_drop" if alpha > 0 else "constant",
            alpha=alpha,
        )

        return rxn, feed, Cp, cfg

    def on_calculate(self):
        try:
            rxn, feed, Cp, cfg = self._collect_inputs()
        except Exception as e:
            messagebox.showerror("Input Error", str(e))
            return

        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("1.0", "Solving PFR... (GUI will stay responsive)")

        def worker():
            try:
                result = solve_pfr(rxn, feed, Cp=Cp, config=cfg)
                self.root.after(0, lambda: self._display_result(result))
            except Exception as exc:
                self.root.after(0, lambda: messagebox.showerror("Solve Error", str(exc)))

        threading.Thread(target=worker, daemon=True).start()

    def _display_result(self, result: PFRResult):
        self.current_result = result

        # Summary text
        self.summary_text.delete("1.0", "end")
        lines = []
        lines.append(f"Reaction: {pretty_stoich(result.reaction.stoichiometry) if result.reaction else ''}")
        lines.append(f"Mode: {result.config.mode if result.config else 'unknown'}")
        lines.append(f"Required Volume: {result.final_V:.6g} m³")
        lines.append(f"Achieved Conversion (limiting {result.limiting_species}): {result.final_X:.5f}")
        lines.append(f"Outlet Temperature: {result.final_T:.2f} K")
        lines.append(f"Outlet Pressure: {result.final_P:.1f} Pa ({result.final_P/101325:.4f} atm)")
        if result.Q:
            duty = result.total_Q
            lines.append(f"Total Heat Duty: {duty:.2f} W  ({'cooling required' if duty < 0 else 'heating required' if duty > 0 else 'neutral'})")
        lines.append(f"Success: {result.success}  — {result.message}")
        self.summary_text.insert("1.0", "\n".join(lines))

        # Plots
        self._update_plots(result)

    def _update_plots(self, result: PFRResult):
        # Clear old tabs
        for tab in self.plot_notebook.tabs():
            self.plot_notebook.forget(tab)
        self.canvases.clear()

        V = np.array(result.V)

        def add_plot_tab(title: str, plot_func):
            frame = ttk.Frame(self.plot_notebook)
            self.plot_notebook.add(frame, text=title)

            fig = Figure(figsize=(6.5, 3.8), dpi=100)
            ax = fig.add_subplot(111)
            plot_func(ax, V, result)
            ax.grid(True, alpha=0.3)
            ax.set_xlabel("V (m³)")

            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill="both", expand=True)
            NavigationToolbar2Tk(canvas, frame).pack()
            self.canvases[title] = canvas

        def plot_conv(ax, V, res):
            ax.plot(V, res.X, "b-", lw=2)
            ax.set_ylabel("Conversion X")
            ax.set_title("Conversion Profile")

        def plot_T(ax, V, res):
            ax.plot(V, res.T, "r-", lw=2)
            ax.axhline(res.feed.T0 if res.feed else res.T[0], ls="--", color="gray", label="T0")
            ax.set_ylabel("T (K)")
            ax.set_title("Temperature Profile")
            ax.legend(fontsize=8)

        def plot_flows(ax, V, res):
            if not res.F:
                return
            species = list(res.F[0].keys())
            for sp in species:
                Fi = [f.get(sp, 0.0) for f in res.F]
                ax.plot(V, Fi, lw=1.6, label=sp)
            ax.set_ylabel("Fi (mol/s)")
            ax.set_title("Molar Flow Profiles")
            ax.legend(fontsize=7, ncol=2)

        def plot_heat(ax, V, res):
            if res.Q:
                ax.plot(V, np.array(res.Q)/1000.0, "m-", lw=2)
                ax.set_ylabel("Cumulative Q (kW)")
                ax.set_title("Heat Duty Profile (isothermal)")
            else:
                ax.text(0.5, 0.5, "Adiabatic (Q=0)", ha="center", transform=ax.transAxes)

        def plot_rate(ax, V, res):
            ax.plot(V, res.r, "g-", lw=2)
            ax.set_ylabel("r (mol/m³/s)")
            ax.set_title("Reaction Rate Profile")

        add_plot_tab("Conversion", plot_conv)
        add_plot_tab("Temperature", plot_T)
        add_plot_tab("Flows", plot_flows)
        add_plot_tab("Heat Duty", plot_heat)
        add_plot_tab("Rate", plot_rate)

    def on_export_csv(self):
        if not self.current_result:
            messagebox.showinfo("No data", "Run a calculation first.")
            return

        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not path:
            return

        res = self.current_result
        try:
            import csv
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                header = ["V_m3", "X", "T_K", "P_Pa"] + [f"F_{s}" for s in (res.F[0].keys() if res.F else [])] + ["r", "Q_W"]
                writer.writerow(header)
                for i in range(len(res.V)):
                    row = [res.V[i], res.X[i], res.T[i], res.P[i]]
                    if res.F:
                        row += [res.F[i].get(s, 0) for s in res.F[0].keys()]
                    row += [res.r[i], res.Q[i] if res.Q and i < len(res.Q) else 0.0]
                    writer.writerow(row)
            messagebox.showinfo("Exported", f"Saved profiles to {path}")
        except Exception as e:
            messagebox.showerror("Export failed", str(e))

    def on_save_plot(self):
        if not self.canvases:
            messagebox.showinfo("No plot", "No plots available.")
            return
        # Save the currently visible tab's figure
        current = self.plot_notebook.tab(self.plot_notebook.select(), "text")
        canvas = self.canvases.get(current)
        if canvas:
            path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG", "*.png")])
            if path:
                canvas.figure.savefig(path, dpi=150, bbox_inches="tight")
                messagebox.showinfo("Saved", f"Plot saved to {path}")

    def clear_results(self):
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("1.0", "Results cleared.")
        for tab in self.plot_notebook.tabs():
            self.plot_notebook.forget(tab)
        self.current_result = None


def main():
    root = tk.Tk()
    app = PFRSizerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()

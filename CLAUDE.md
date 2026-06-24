# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

A two-script ANUGA hydrodynamic experiment that models stormwater inlets (pits/grates/lintels) on a sloping surface-flow plane, plus a Tkinter GUI for inspecting the resulting hydrographs.

- `stormwater_inlet_simulation.py` — runs the 2D shallow-water simulation and writes one `hydrograph_<Asset_ID>.csv` per inlet, plus `sloped_inlet_experiment.sww`.
- `stormwater_inlet_viewer.py` — standalone Tkinter dashboard that scans a folder for those CSVs and renders four diagnostic subplots per inlet.

The two scripts are coupled only by the CSV schema (see below), not by imports.

## Commands

```bash
# Run the simulation (requires a working ANUGA install; regenerates the hydrograph_*.csv and .sww files)
python stormwater_inlet_simulation.py

# Launch the hydrograph viewer GUI (needs a display / X server for Tkinter)
python stormwater_inlet_viewer.py

# Run the tests (pytest; fast unit tests + the slower domain-based integration tests)
python -m pytest
```

Tests use **pytest** (plain functions, `assert`, `pytest.approx`, `@pytest.mark.parametrize`). The simulation only runs under `if __name__ == "__main__"` (via `run_experiment()`), so the module can be imported for testing without evolving a domain. `test_inlet_hydraulics.py` covers the pure asset/hydraulics logic and the CSV-contract columns (no domain built). `test_inlet_operator_integration.py` builds a real ANUGA pond and (a) asserts the live `update_Q` discharge equals the weir/orifice equation evaluated with each catalogued grate's operational area & perimeter (including blockage derating), and (b) runs a flat, reflective-walled pond drained by one inlet, checking the domain's water loss equals the operator's captured volume (to machine precision) while the draining depth carries the inlet through both regimes. `test_parallel_inlet_mass_balance.py` is an MPI test that follows the anuga_core convention: the pytest function shells out via `anuga.mpicmd` to re-run the same file under `mpiexec -np 2` (skipped if `mpi4py`/`mpiexec` are absent), and that worker (`__main__`) distributes the pond, drains it with a `Depth_driven_parallel_inlet_operator`, and asserts the global mass balance and regime crossing. There is no build step, linter config, or package manifest. Dependencies (`anuga`, `numpy`, `pandas`, `matplotlib`, `tkinter`) are assumed present in the environment; ANUGA in particular is not pip-trivial and is the main setup hurdle.

## Architecture

### Simulation script
The hydraulics live in an asset-library + operator design:

- `Inlet_specification` holds geometry (`clear_area`, `effective_perimeter`) and a `blockage_factor`; `operational_area`/`operational_perimeter` derate the geometry by blockage. `INLET_LIBRARY` is the catalog of named standard inlet types keyed by spec string.
- `Depth_driven_inlet_operator` is the core. It must subclass the Inlet_operator **class** from `anuga.structures.inlet_operator` — note `anuga.Inlet_operator` (the top-level name) is a factory *function* and is **not** subclassable. It overrides `update_Q(t)` to return the capture discharge as a **negative** value (water leaving the domain), computed from the inlet's region-averaged depth (`self.inlet.get_average_depth()`) via a **dual-regime capture law**: weir flow (`C_w * P * depth^1.5`) below a precomputed transition depth `d_trans`, orifice flow (`C_o * A * sqrt(2g*depth)`) above it. The capture law + hydrograph logging are factored into a shared `_Depth_driven_capture_mixin`: pure `@staticmethod`s `transition_depth(A, P, C_w, C_o, g)` and `capture_discharge(depth, A, P, C_w, C_o, g, d_trans)` (unit-testable without an ANUGA domain), plus `_capture_Q`/`_log_capture` helpers. `d_trans = (C_o*A*sqrt(2g)) / (C_w*P)` — the exponent-1 weir/orifice crossover (an earlier `**(2/3)` was wrong: it made the law discontinuous and non-monotonic). The parent `Inlet_operator` distributes the discharge over the inlet region and enforces mass balance — the operator no longer touches the stage array directly. `__call__` calls `super().__call__()` then logs a hydrograph record (using `applied_Q` for realised capture, and region-averaged momentum × `sqrt(area)` to estimate approach flow and hence bypass).
- `Depth_driven_parallel_inlet_operator` is the MPI-safe sibling (subclasses `anuga.parallel.parallel_inlet_operator.Parallel_Inlet_operator`, shares the same mixin). On a distributed domain `get_average_depth()` is per-subdomain (local), so it uses the collective `get_global_average_depth()`/`get_global_*` reductions. Because `Parallel_Inlet_operator.__call__` runs `update_Q` on the **master rank only** (then broadcasts the volume), the global depth is sampled collectively in `__call__` (all ranks) and stashed for the master-only `update_Q` — calling a collective inside `update_Q` would deadlock. `Stormwater_inlet_network.add_inlet` auto-selects serial vs parallel from `domain.parallel` and forwards `**operator_kwargs` (e.g. `master_proc`/`procs`). The parallel import is guarded (`_HAVE_PARALLEL`) so the module still imports without ANUGA's parallel stack. See `docs/HYDRAULICS.md` → "Running in parallel (MPI)".
- `Stormwater_inlet_network` registers inlets, building a small circular `anuga.Region(center=[x,y], radius=...)` footprint per pit (default radius 1.5 m), owns the per-asset capture logs, and exposes `to_dataframe()` for export. To inspect ANUGA's `Inlet_operator`/`Inlet` API, see `/home/steve/anuga_core/anuga/structures/`.

Runtime flow (under `if __name__ == "__main__"`, via `run_experiment()`): `build_domain()` makes a 100×20 m domain via `create_domain_from_regions`, sets a 1% slope elevation, dry initial stage, Dirichlet inflow on `left`, Transmissive `right`, Reflective walls; `build_network()` places the 6 `PIT_PLACEMENTS` down the channel centerline; `domain.evolve(yieldstep=10, finaltime=120)`; then `print_summary()` prints a table and dumps per-inlet CSVs.

### Key physical constants / knobs
- Weir/orifice coefficients `C_w=1.66`, `C_o=0.67` (operator defaults; HEC-22 metric grate-inlet values).
- `constant_inflow_stage = 1.35` and the `elevation_function` slope set the steady approach depth.
- `maximum_triangle_area=2.0` controls mesh resolution; pit capture resolves to a single triangle, so mesh density at the pit location materially affects results.

### CSV contract (the integration point)
Both scripts depend on this exact column set, written by `Depth_driven_inlet_operator.capture_log` and validated by the viewer's `required_headers`:

```
Time_s, Depth_m, Approach_Q_cms, Captured_Q_cms, Bypass_Q_cms, Cum_Captured_m3, Cum_Bypassed_m3
```

Flows are stored in m³/s (cms) and converted to L/s only at display time. **If you change a column name in the simulation, update `required_headers` in `stormwater_inlet_viewer.py` or the GUI will reject the file.**

### Viewer
`Stormwater_inlet_viewer_app` (in `stormwater_inlet_viewer.py`) is a single-class Tkinter app: a sidebar lists CSVs in the chosen directory, and selecting one calls `generate_hydraulic_plots` to draw four stacked matplotlib axes (approach-vs-capture, accumulated volumes, flows+depth twin-axis, and a time-colored depth-vs-discharge hysteresis plot) embedded via `FigureCanvasTkAgg`. A **View** menu switches between this per-inlet "Pit Hydrograph" and a "Combined Hydrograph" (`show_combined_hydrograph`) that sums `Captured_Q_cms`/`Bypass_Q_cms` across every CSV in the folder onto a common `Time_s` axis and plots instantaneous flows (L/s, left axis) with cumulative volumes (m³, right twin axis). `self.current_view` ("pit"/"combined") tracks which is showing so the plot-font slider re-renders the right one. The menu bar is built from `tk.Label`+`tk_popup` (not native menus, which ignore fonts on the user's Linux/Wayland theme).

# Simple SW Inlets

A small ANUGA hydrodynamic experiment that models stormwater inlets
(pits / grates / lintels) capturing water from a 2-D shallow-water surface
flow, plus a Tkinter GUI for inspecting the resulting hydrographs.

The simulation places a row of inlets down a sloping channel, drives a steady
inflow over them, and records how much water each inlet captures versus how much
bypasses it. Each inlet captures water according to a **dual-regime weir/orifice
capture law** driven by the local ponded depth. See
[`docs/HYDRAULICS.md`](docs/HYDRAULICS.md) for the equations and their basis.

---

## Contents

| File | Purpose |
|------|---------|
| `py_ANUGA_Simple_SW_WORKING.py` | Runs the 2-D shallow-water simulation; writes one `hydrograph_<Asset_ID>.csv` per inlet plus `sloped_inlet_experiment.sww`. |
| `py_Plot_Pit_Depth_Q2.py` | Standalone Tkinter dashboard that scans a folder for the hydrograph CSVs and renders four diagnostic plots per inlet. |
| `test_inlet_hydraulics.py` | Fast unit tests for the asset/hydraulics logic and the CSV schema (no ANUGA domain built). |
| `test_inlet_operator_integration.py` | Integration tests that evolve/query a real ANUGA domain (capture law ↔ live `update_Q`, and mass balance). |
| `docs/HYDRAULICS.md` | The capture-law theory, equations, transition depth, and references. |

The two main scripts are coupled **only** by the CSV schema (below), not by
imports.

---

## Requirements

- Python 3
- [`anuga`](https://github.com/anuga-community/anuga_core) (the main setup
  hurdle — not pip-trivial)
- `numpy`, `pandas`, `matplotlib`
- `tkinter` (for the viewer GUI; needs a display / X server)
- `pytest` (to run the tests)

ANUGA is assumed already installed in the environment; there is no package
manifest or build step for this project.

---

## Usage

### Run the simulation

```bash
python py_ANUGA_Simple_SW_WORKING.py
```

This builds the domain, registers the six inlets, evolves for 120 s, prints a
results summary, and writes `hydrograph_*.csv` (one per inlet) and
`sloped_inlet_experiment.sww` into the working directory.

The runtime is guarded by `if __name__ == "__main__"`, so importing the module
(e.g. from a notebook or the tests) does **not** start a simulation. You can
also drive it programmatically:

```python
import py_ANUGA_Simple_SW_WORKING as sim

network = sim.run_experiment(yieldstep=10, finaltime=120, write_csv=True)
df = network.to_dataframe("Pit_01_SmallGrate")
```

### Launch the hydrograph viewer

```bash
python py_Plot_Pit_Depth_Q2.py
```

A sidebar lists the hydrograph CSVs found in the chosen directory; selecting one
draws four stacked plots:

1. Approach vs captured discharge over time
2. Accumulated captured / bypassed volumes
3. Flows with ponded depth on a twin axis
4. A time-coloured depth-vs-discharge hysteresis loop

The viewer is HiDPI-aware (it reads the configured `Xft.dpi` rather than the
DPI Tk reports, which is unreliable on Wayland/XWayland) and has live
**UI Font Size** and **Plot Font Size** sliders. Per-machine slider values and
window geometry are saved to `~/.stormwater_inlet_viewer.json`.

### Run the tests

```bash
python -m pytest
```

Tests use **pytest**. `test_inlet_hydraulics.py` is fast (no domain);
`test_inlet_operator_integration.py` builds a small ANUGA pond and is a little
slower but still runs in well under a second.

---

## The experiment

`build_domain()` creates a 100 m × 20 m channel meshed with
`maximum_triangle_area = 2.0`, with:

- a **1 % slope** bed: `elevation = 1.0 − 0.01·x` (falls left → right),
- an initially dry stage (`stage = elevation`), Manning friction `0.015`,
- a **Dirichlet inflow** on the left at `stage = 1.35` (≈ 0.35 m deep entering
  at `x = 0`), a **Transmissive** right boundary, and **Reflective** side walls.

`build_network()` then places six inlets along the centreline (`y = 10 m`), each
with a circular sampling footprint (default radius 1.5 m):

| Asset ID | x (m) | Spec |
|----------|-------|------|
| Pit_01_SmallGrate | 15 | Grate_600x600 |
| Pit_02_LargeGrate | 28 | Grate_900x900 |
| Pit_03_ShortLintel | 42 | Lintel_1.2m |
| Pit_04_LongLintel | 56 | Lintel_2.4m |
| Pit_05_ComboSmall | 70 | Combo_1.2m_G600 |
| Pit_06_ComboLarge | 85 | Combo_2.4m_G900 |

`domain.evolve(yieldstep=10, finaltime=120)` runs the simulation, then
`print_summary()` prints the steady-state table and dumps the CSVs.

> **Note:** capture resolves over the inlet footprint, which can be a small
> number of triangles, so mesh density (`maximum_triangle_area`) at the pit
> location materially affects results.

---

## Code overview

The hydraulics use an asset-library + operator design (see
[`docs/HYDRAULICS.md`](docs/HYDRAULICS.md) for the physics):

- **`Inlet_specification`** — geometry of a standard inlet (`clear_area`,
  `effective_perimeter`) plus a `blockage_factor`. `operational_area` and
  `operational_perimeter` derate the geometry by blockage. `INLET_LIBRARY` is
  the catalogue of named standard inlets.
- **`Depth_driven_inlet_operator(anuga.Inlet_operator)`** — the core. Its
  `update_Q(t)` returns the capture discharge as a **negative** value (water
  leaving the domain), computed from the inlet's region-averaged depth via the
  dual-regime law. The pure hydraulics are factored into the static methods
  `transition_depth(...)` and `capture_discharge(...)` so they can be tested
  without a domain. The parent `Inlet_operator` distributes that discharge over
  the inlet region and enforces mass balance.
- **`Stormwater_inlet_network`** — registers inlets onto the domain (one small
  circular `anuga.Region` per pit), owns the per-asset capture logs, and exposes
  `to_dataframe()` for export.

The viewer is a single class, **`Stormwater_inlet_viewer_app`**, in
`py_Plot_Pit_Depth_Q2.py`.

---

## CSV contract (the integration point)

Both scripts depend on this exact column set, written by the operator's capture
log and validated by the viewer's `required_headers`:

```
Time_s, Depth_m, Approach_Q_cms, Captured_Q_cms, Bypass_Q_cms, Cum_Captured_m3, Cum_Bypassed_m3
```

Flows are stored in m³/s (cms) and converted to L/s only at display time.
The exported per-inlet DataFrame additionally prepends an `Asset_ID` column.

> If you rename a column in the simulation, update `required_headers` in
> `py_Plot_Pit_Depth_Q2.py` or the GUI will reject the file.

---

## Generated / artifact files

Running the simulation or tests produces files in the working directory that are
**outputs, not source**: `hydrograph_*.csv`, `sloped_inlet_experiment.sww`, and
the test ponds (`inlet_integration_test.sww`, etc.). They can be deleted and
regenerated at any time.

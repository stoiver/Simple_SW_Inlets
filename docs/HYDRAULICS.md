# Inlet capture hydraulics

This document describes how `Depth_driven_inlet_operator` decides how much water
an inlet captures, and the basis for the equations and coefficients.

## Overview

Each inlet captures water as a function of the **ponded depth** `d` over the
inlet. The operator's `use_max_depth` argument (default True, from the module
`USE_MAX_DEPTH` constant) selects whether `d` is the maximum depth over the inlet
footprint or ANUGA's region average (`inlet.get_average_depth()`). A grate / pit
in a sag captures water in two regimes:

- **Weir flow** at shallow depths — water spills over the *perimeter* of the
  opening.
- **Orifice flow** at greater depths — the opening is submerged and flow is
  limited by the *clear area*.

The capture discharge is the **smaller** of the two at any given depth, which is
the standard HEC-22 procedure: *"evaluate both equations and use the lower
value."*

## Equations (SI / metric)

Let

- `d` = ponded depth over the inlet (m),
- `A` = clear opening area (m²) = `Inlet_specification.operational_area`,
- `P` = effective perimeter (m) = `Inlet_specification.operational_perimeter`,
- `g` = 9.81 m/s²,
- `C_w` = weir coefficient = **1.66** (metric; 3.0 in US customary units),
- `C_o` = orifice coefficient = **0.67**.

**Weir flow:**

```
Q_weir = C_w · P · d^1.5
```

**Orifice flow:**

```
Q_orifice = C_o · A · √(2 · g · d)
```

`A` and `P` are the *operational* values, i.e. the clear geometry derated by the
inlet's `blockage_factor`:

```
operational_area      = clear_area          · (1 − blockage_factor)
operational_perimeter = effective_perimeter · (1 − blockage_factor)
```

So a 50 % blocked grate captures half as much (linear in `A` for orifice flow,
linear in `P` for weir flow).

## Transition depth

Rather than evaluating `min(Q_weir, Q_orifice)` every step, the operator
precomputes the depth at which the two laws are equal and switches between them.
Setting `Q_weir = Q_orifice` and solving:

```
C_w · P · d^1.5 = C_o · A · √(2g) · d^0.5
        (divide both sides by d^0.5)
C_w · P · d      = C_o · A · √(2g)

  d_trans = (C_o · A · √(2g)) / (C_w · P)
```

The `d^0.5` factors cancel, so the exponent is **1** — `d_trans` is simply that
ratio.

This switch is **exactly equivalent** to taking the minimum of the two laws.
Their ratio is

```
Q_weir / Q_orifice = [C_w·P / (C_o·A·√(2g))] · d = d / d_trans
```

so `Q_weir < Q_orifice` precisely when `d < d_trans`. Therefore:

```
              ⎧ 0                        if d ≤ 1e-4   (negligible)
Q_capture =   ⎨ C_w · P · d^1.5          if d < d_trans   (weir)
              ⎩ C_o · A · √(2g·d)        if d ≥ d_trans   (orifice)
```

`Q_capture` is returned from `update_Q(t)` as a **negative** discharge (water
leaving the domain). A degenerate inlet with `P ≤ 0` has `d_trans = 0`, i.e. it
is always in the orifice regime.

> **History:** an earlier version used `d_trans = (ratio)^(2/3)`, which is
> wrong: it placed the switch away from the true crossover, making the capture
> law **discontinuous and non-monotonic** (it would capture *more* at 0.2 m than
> at 0.4 m). The unit tests `test_law_is_continuous_at_transition` and
> `test_capture_increases_with_depth` guard against this.

## Where this lives in the code

In `stormwater_inlets.py`:

```python
class Depth_driven_inlet_operator(Inlet_operator):

    @staticmethod
    def transition_depth(A, P, C_w, C_o, g=9.81):
        if P <= 0:
            return 0.0
        return (C_o * A * np.sqrt(2 * g)) / (C_w * P)

    @staticmethod
    def capture_discharge(depth, A, P, C_w, C_o, g=9.81, d_trans=None):
        if depth <= 1e-4:
            return 0.0
        if d_trans is None:
            d_trans = Depth_driven_inlet_operator.transition_depth(A, P, C_w, C_o, g)
        if depth < d_trans:
            return C_w * P * (depth ** 1.5)      # Weir flow
        return C_o * A * np.sqrt(2 * g * depth)  # Orifice flow
```

`transition_depth` and `capture_discharge` are pure functions of their
arguments (no ANUGA domain needed), which is what the unit tests in
`test_inlet_hydraulics.py` exercise. `update_Q` simply calls them with the live
inlet depth and the spec's operational geometry, and the integration tests in
`test_inlet_operator_integration.py` confirm the live `update_Q` output matches
these equations for every catalogued grate.

## Worked example

For `Grate_600x600` (`clear_area = 0.21 m²`, `effective_perimeter = 2.40 m`,
no blockage), with `C_w = 1.66`, `C_o = 0.67`, `g = 9.81`:

```
d_trans = (0.67 · 0.21 · √(2·9.81)) / (1.66 · 2.40) ≈ 0.156 m
```

- At `d = 0.10 m` (< `d_trans`, weir):
  `Q = 1.66 · 2.40 · 0.10^1.5 ≈ 0.126 m³/s`
- At `d = 0.30 m` (> `d_trans`, orifice):
  `Q = 0.67 · 0.21 · √(2·9.81·0.30) ≈ 0.341 m³/s`

## Running in parallel (MPI)

On a single (serial) domain, `inlet.get_average_depth()` covers the whole inlet
and the capture law is correct as written. On a **distributed** domain
(`anuga.distribute(domain)`, run under `mpirun`/`mpiexec`) the inlet footprint
can be split across ranks, and `get_average_depth()` is then **per-subdomain
(local)** — feeding the capture law only this rank's slice of the pond. The
*global* average is needed instead.

`Depth_driven_parallel_inlet_operator` handles this. It subclasses ANUGA's
`Parallel_Inlet_operator` and shares the exact same capture hydraulics via the
`_Depth_driven_capture_mixin`, differing only in *how the inlet state is
sampled*:

- It uses the collective `get_global_average_depth()` (and the `get_global_*`
  momentum/area reductions for the hydrograph log) instead of the local ones.
- **`update_Q` runs on the master rank only.** `Parallel_Inlet_operator.__call__`
  computes the discharge on the master and broadcasts the resulting volume to the
  other ranks (which wait in a receive). The `get_global_*` calls are collective
  and must be entered by *every* rank, so they are sampled in `__call__` (all
  ranks) and the master-only `update_Q` reads the stashed global depth. Calling a
  collective reduction *inside* `update_Q` would deadlock the waiting ranks.

`Stormwater_inlet_network.add_inlet` selects the operator automatically from
`domain.parallel` (serial vs distributed). ANUGA does **not** auto-discover which
ranks hold an inlet, so for a footprint spanning multiple ranks pass the
participating ranks through:

```python
network.add_inlet("Pit_01", x, y, "Grate_600x600",
                  master_proc=0, procs=list(range(numprocs)))
```

This is exercised by `test_parallel_inlet_mass_balance.py`, which shells out to
`mpiexec -np 2` (anuga_core convention): the global water loss matches the
operator's captured volume to machine precision (~3e-15), confirming the global
reduction, the master-only `update_Q`, and the collective sampling are correct.

## Inlet catalogue (`INLET_LIBRARY`)

| Key | clear_area (m²) | effective_perimeter (m) |
|-----|-----------------|--------------------------|
| Grate_600x600 | 0.21 | 2.40 |
| Grate_900x900 | 0.48 | 3.60 |
| Lintel_1.2m | 0.18 | 1.20 |
| Lintel_2.4m | 0.36 | 2.40 |
| Combo_1.2m_G600 | 0.39 | 3.00 |
| Combo_2.4m_G900 | 0.84 | 5.10 |

### Where these numbers come from

These are **representative standard-inlet values**, hard-coded in
`INLET_LIBRARY` — they are *given* geometry, not quantities the model derives.
For reference, this is how such figures are arrived at:

- **Clear area** is the gross footprint times the grate's open-area ratio (the
  fraction not blocked by bars). A "600×600" grate has a gross plan area of
  `0.600 m × 0.600 m = 0.36 m²`; the catalogue's `0.21 m²` therefore implies an
  open ratio of `0.21 / 0.36 ≈ 58 %`. Real bar grates run roughly 50–80 % open
  depending on bar spacing, so the precise value comes from the manufacturer's
  grate geometry (sum of the open slot areas), not from the nominal size alone.
- **Effective perimeter** here is the full nominal outer perimeter,
  `4 × 0.600 m = 2.40 m`. Note HEC-22's convention is to *disregard the side
  against the curb* (3 sides) for a kerb-side grate; using all four sides
  corresponds to a grate **not** against a curb — e.g. the mid-channel /
  sag placement used in this experiment.

The `blockage_factor` is applied on top of these to give
`operational_area` / `operational_perimeter`; it models *additional* clogging
(debris, leaves) beyond the grate's built-in open ratio.

## References

The weir/orifice formulation, the metric coefficients `C_w = 1.66` and
`C_o = 0.67`, and the "use the lower discharge" transition rule follow FHWA
HEC-22 (Urban Drainage Design Manual) and its derivatives:

- TxDOT Hydraulic Design Manual — *Gutter and Inlet Equations*:
  <https://www.txdot.gov/manuals/des/hyd/chapter-10--storm-drains/section-6--gutter-and-inlet-equations.html>
- FHWA HEC-22, 3rd edition — *Urban Drainage Design Manual*:
  <https://www.fhwa.dot.gov/engineering/hydraulics/pubs/10009/10009.pdf>

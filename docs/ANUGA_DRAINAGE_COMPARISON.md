# Simple_SW_Inlets vs. anuga_drainage

An analysis of how this project relates to
[`anuga_drainage`](https://github.com/anuga-community/anuga_drainage) (the
ANUGA↔SWMM/pipedream coupling package), whether the two should be combined, and
a concrete diff of their inlet-exchange flux laws.

## TL;DR

- They share the **same physics family** (weir/orifice exchange through an ANUGA
  `Inlet_operator`) but use **different formulations, coefficients, and
  regime-selection logic** — they are *not* drop-in compatible.
- This project models **one-way capture** to an implicit sink; `anuga_drainage`
  models **two-way coupling** to a real 1D pipe network (capture *and* surcharge).
- Recommendation: **don't merge wholesale.** Keep this repo as the lightweight,
  dependency-free, teaching/standalone tool; cross-pollinate selectively (see
  [Recommendation](#recommendation)).

## What each project is

| | Simple_SW_Inlets (this repo) | anuga_drainage |
|---|---|---|
| Purpose | Self-contained inlet-capture experiment / teaching tool | Production 2D↔1D coupled urban drainage |
| 1D network | None — captured water leaves the domain (implicit infinite sink) | SWMM (`pyswmm`) or pipedream (`pipedream_solver`) |
| Exchange | One-way capture (always extraction) | Two-way: capture **and** surcharge back onto the surface |
| Flux law | `capture_discharge` (HEC-22 metric weir/orifice, depth-driven) | `calculate_Q` (Leandro & Martins 2016, head-driven) |
| Distinctive extras | Named inlet **catalogue** (`clear_area`/`effective_perimeter`), **blockage** derating, TOML config, **Tkinter hydrograph viewer**, serial+MPI operators | `.inp` parsing, `Coupler` + SWMM/pipedream backends, volume/mass-balance tracking, packaging (src/pyproject), RTD docs |
| Home | `stoiver` (personal) | `anuga-community` (org) |

Conceptually this project is the **pipe-free special case** of `anuga_drainage`:
capture with the 1D head held at/below the inlet bed, so flow is always inward
and never surcharges.

## The flux laws, side by side

**Simple_SW_Inlets** — `capture_discharge(depth, A, P, C_w, C_o, g)` with
`C_w=1.66`, `C_o=0.67` (see [`HYDRAULICS.md`](HYDRAULICS.md)):

```
d < d_trans :  Q = C_w · P · d^1.5                 (weir)
d ≥ d_trans :  Q = C_o · A · √(2g·d)               (orifice)
d_trans     =  (C_o·A·√2g)/(C_w·P)                 # weir/orifice crossover
```

- Inputs: surface ponded depth `d`, orifice area `A`, weir perimeter `P`.
- Output: a non-negative magnitude; `update_Q` returns `−Q` (water always leaves).
- The two branches are two fits to **one** capture-vs-depth curve; the code takes
  the lower (switch at `d_trans`), so capture **plateaus** at depth.

**anuga_drainage** — `calculate_Q(head1D, depth2D, bed2D, length_weir L,
area_manhole A, cw=0.67, co=0.67)`:

```
head1D < bed2D                   :  Q =  cw · L · depth2D · √(2g·depth2D)         (free weir in)
bed2D ≤ head1D < bed2D+depth2D   :  Q =  co · A · √(2g·(depth2D+bed2D−head1D))    (orifice in)
head1D > bed2D+depth2D           :  Q = −co · A · √(2g·(head1D−bed2D−depth2D))    (surcharge out)
```

- Inputs: 1D pipe hydraulic head `head1D`, surface depth `depth2D`, inlet bed
  elevation `bed2D`, weir length `L`, manhole area `A`.
- Output: **signed** — positive = surface→pipe, negative = pipe→surface.
- The three branches are three **physical states** of the surface↔pipe
  connection, selected by the **pipe head** (not a depth crossover).

Geometry mapping: `P ↔ length_weir`, `A ↔ area_manhole`, `d ↔ depth2D`.

## Three substantive differences

### 1. Weir coefficient — differ by ~1.8×
Rewriting drainage's weir as a `P·d^1.5` law gives an effective coefficient
`cw·√(2g) = 0.67·4.429 = 2.968`, versus this project's `1.66`. Same shape,
different calibration (HEC-22 lumped coefficient vs Leandro & Martins'
explicit-√(2g) form). **Ratio 1.788.**

> Note: drainage's effective weir coefficient (2.97) is also higher than the
> textbook sharp-crested value (~1.77, the `⅔·Cd·√2g` form), so it's worth
> checking before adopting wholesale.

### 2. Orifice — identical, in the empty-pipe limit
Both orifice branches are `0.67·A·√(2g·Δh)` with the same `Co`. This project's
driving head is the surface depth `d`; drainage's is the head difference
`(depth2D + bed2D − head1D)`. When the pipe water level sits at the inlet bed
(`head1D = bed2D`), that difference *is* `d`, and the two coincide **exactly**.
So this project's capture law is precisely drainage's orifice branch with an
empty pipe.

### 3. Regime selection means different things
- **This project**: `min(weir, orifice)` vs **depth** → capture curve plateaus.
- **drainage**: weir / orifice / surcharge selected by **pipe head**; in "weir
  mode" (drawn-down pipe) capture keeps growing as `d^1.5`, uncapped.

### 4. Direction
This project is one-way (always extraction). drainage is signed and can
**surcharge** water back onto the surface — impossible here without a pipe head.

## Worked numbers

For `Grate_600x600` (`A = 0.21 m²`, `P = L = 2.40 m`), `g = 9.81`,
this project's `d_trans = 0.156 m`:

```
 d (m) | Simple Q  (branch) | drain weir | drain orifice (empty pipe)
---------------------------------------------------------------------
 0.050 |   0.0445    weir   |   0.0796   |        0.1394
 0.100 |   0.1260    weir   |   0.2252   |        0.1971
 0.156 |   0.2465    orif   |   0.4407   |        0.2465   ← crossover
 0.300 |   0.3414    orif   |   1.1704   |        0.3414
 0.600 |   0.4827    orif   |   3.3103   |        0.4827
```

- The **orifice** column matches this project's `Q` exactly once past `d_trans`
  (0.2465, 0.3414, 0.4827) — confirming difference #2.
- The **weir** column is `1.788×` this project's weir values (0.0445→0.0796,
  0.1260→0.2252) — confirming difference #1.
- This project's `Q` plateaus (≈0.48 at 0.6 m); drainage's weir branch does not
  (3.31 at 0.6 m) — confirming difference #3.

(Reproduce with the small script in the commit/PR that added this doc, or by
evaluating the two formulas above.)

## Recommendation

A full merge would dissolve this project's value (simple, dependency-light, no
SWMM/pipedream, followable). Prefer **selective cross-pollination**:

1. **Contribute the additive pieces to `anuga_drainage`** (it lacks them):
   the inlet **asset catalogue + blockage** model, and the **hydrograph viewer**.
2. **Align the physics**: adopt a single weir convention (or document the
   divergence), and treat this project's `capture_discharge` as the
   `head1D ≤ bed2D` degenerate of `calculate_Q`.
3. If one repo is desired, fold this in as an **"ANUGA-only inlet capture"
   example/module inside `anuga_drainage`** (it already has an `examples/` dir
   and an ANUGA-only Boyd reference) — not the reverse.

A shared flux kernel, if pursued, should be the **general** `calculate_Q`, with
this project's law recovered by fixing the pipe head — because the two laws do
not decompose the same way (`min(weir,orifice)` vs head-state selection).

> Combining is also a **governance/scope** decision (personal vs `anuga-community`,
> teaching vs production), not only a code one.

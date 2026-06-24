"""Integration tests that evolve / query a real ANUGA domain.

Unlike test_inlet_hydraulics (pure logic, no domain), these build an actual
ANUGA domain, so they are slower. Two things are checked:

  * The discharge the operator actually feeds to ANUGA (``update_Q``) equals the
    weir / orifice equation evaluated with the grate's own operational area and
    perimeter, for every catalogued grate and for a blocked grate. This ties the
    live operator output to the hydraulics and to the asset geometry.

  * Mass balance -- in a closed (reflective-walled) pond the inlet is the only
    sink, so the domain's water loss must equal the volume the operator reports
    capturing, and the draining depth must carry the inlet through both the
    orifice and weir regimes.

Run with::

    python -m pytest test_inlet_operator_integration.py
"""

import math

import pytest

import anuga

import stormwater_inlets as sim

OP = sim.Depth_driven_inlet_operator
# HEC-22 metric coefficients (the operator defaults); asserting against these
# directly means a coefficient change would be caught.
C_W, C_O, G = 1.66, 0.67, 9.81

SIDE = 5.0          # square pond side length (m)
CENTER = [SIDE / 2, SIDE / 2]
RADIUS = 1.0        # inlet footprint radius (m)


def build_pond(initial_depth, side=SIDE, n=10):
    """Flat, frictionless square pond of uniform depth with closed walls."""
    domain = anuga.rectangular_cross_domain(n, n, len1=side, len2=side)
    domain.set_name("inlet_integration_test")
    domain.set_store(False)                  # don't leave a .sww artifact behind
    domain.set_quantity("elevation", 0.0)
    domain.set_quantity("friction", 0.0)
    domain.set_quantity("stage", initial_depth)   # flat bed => uniform depth
    reflective = anuga.Reflective_boundary(domain)
    domain.set_boundary({"left": reflective, "right": reflective,
                         "top": reflective, "bottom": reflective})
    return domain


def operator_at_uniform_depth(domain, spec, depth):
    """Set a uniform pond depth and return (operator, sampled avg depth).

    The operators never interact (no evolve is run), so a single domain can be
    reused across checks by just resetting the stage.
    """
    domain.set_quantity("stage", depth)
    region = anuga.Region(domain, center=CENTER, radius=RADIUS)
    op = OP(domain, region, spec)
    return op, op.inlet.get_average_depth()


@pytest.fixture(scope="module")
def pond():
    """A reusable pond for the (evolve-free) update_Q equation checks."""
    return build_pond(0.30)


# --------------------------------------------------------------------------- #
# use_max_depth is a per-operator argument (not a hard-wired global)           #
# --------------------------------------------------------------------------- #

def test_use_max_depth_argument_flows_to_operator(pond):
    pond.set_quantity("stage", 0.30)
    net = sim.Stormwater_inlet_network(pond)
    op_default = net.add_inlet("d", CENTER[0], CENTER[1], "Grate_600x600", radius=RADIUS)
    op_avg = net.add_inlet("a", CENTER[0], CENTER[1], "Grate_600x600",
                           radius=RADIUS, use_max_depth=False)
    assert op_default.use_max_depth is True          # default from USE_MAX_DEPTH
    assert op_avg.use_max_depth is False             # per-inlet override honored
    # On a uniform pond max == average, so both sample the same depth.
    assert op_default._sample_depth() == pytest.approx(op_avg._sample_depth())


# --------------------------------------------------------------------------- #
# update_Q matches the weir/orifice law with the grate's area & perimeter      #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("key", list(sim.INLET_LIBRARY))
def test_orifice_Q_uses_spec_area_for_every_grate(pond, key):
    spec = sim.INLET_LIBRARY[key]
    A, P = spec.operational_area, spec.operational_perimeter
    d_trans = OP.transition_depth(A, P, C_W, C_O, G)
    op, depth = operator_at_uniform_depth(pond, spec, 2.0 * d_trans)
    assert depth > op.d_trans                       # confirm orifice regime
    assert op.update_Q(0.0) == pytest.approx(-(C_O * A * math.sqrt(2 * G * depth)),
                                             abs=1e-9)


@pytest.mark.parametrize("key", list(sim.INLET_LIBRARY))
def test_weir_Q_uses_spec_perimeter_for_every_grate(pond, key):
    spec = sim.INLET_LIBRARY[key]
    A, P = spec.operational_area, spec.operational_perimeter
    d_trans = OP.transition_depth(A, P, C_W, C_O, G)
    op, depth = operator_at_uniform_depth(pond, spec, 0.5 * d_trans)
    assert depth < op.d_trans                       # confirm weir regime
    assert op.update_Q(0.0) == pytest.approx(-(C_W * P * depth ** 1.5), abs=1e-9)


def test_update_Q_is_negative_for_extraction(pond):
    spec = sim.INLET_LIBRARY["Grate_600x600"]
    op, _ = operator_at_uniform_depth(pond, spec, 0.30)
    assert op.update_Q(0.0) < 0.0


@pytest.mark.parametrize("regime,factor", [("weir", 0.5), ("orifice", 2.0)])
def test_blockage_derates_area_and_perimeter_in_Q(pond, regime, factor):
    base = sim.INLET_LIBRARY["Grate_600x600"]
    # 50% blockage halves both clear_area and effective_perimeter; the transition
    # depth depends only on A/P (unchanged here), so the same pond depth stays in
    # the same regime for blocked and unblocked operators.
    blocked = sim.Inlet_specification(
        "Blocked", base.clear_area, base.effective_perimeter, blockage_factor=0.5)
    d_trans = OP.transition_depth(base.operational_area,
                                  base.operational_perimeter, C_W, C_O, G)
    depth = factor * d_trans
    op_full, _ = operator_at_uniform_depth(pond, base, depth)
    op_blocked, _ = operator_at_uniform_depth(pond, blocked, depth)
    # Q scales linearly with A (orifice) and P (weir), so halving both halves Q.
    assert op_blocked.update_Q(0.0) == pytest.approx(0.5 * op_full.update_Q(0.0),
                                                     abs=1e-9)


# --------------------------------------------------------------------------- #
# Mass balance: domain water loss == operator's reported capture              #
# --------------------------------------------------------------------------- #

def test_capture_matches_domain_loss_and_spans_both_regimes():
    initial_depth = 0.30          # > d_trans (~0.156 m) => starts in orifice regime
    domain = build_pond(initial_depth)

    spec = sim.INLET_LIBRARY["Grate_600x600"]
    capture_log = []
    region = anuga.Region(domain, center=CENTER, radius=RADIUS)
    op = OP(domain, region, spec, capture_log=capture_log, label="Pit")

    volume_before = domain.get_water_volume()
    for _t in domain.evolve(yieldstep=2.0, finaltime=40.0):
        pass
    volume_after = domain.get_water_volume()

    removed = volume_before - volume_after

    # 1. Something drained, but the pond is not fully empty.
    assert removed > 0.0
    assert volume_after > 0.0, "pond drained completely; pick a shorter finaltime"

    # 2. Mass balance: the inlet is the only sink, so loss == captured.
    tol = max(1e-6, 1e-3 * removed)
    assert op.total_volume_captured == pytest.approx(removed, abs=tol)
    assert -op.total_applied_volume == pytest.approx(removed, abs=tol)

    # 3. Requested ~ applied: nothing was clamped (we stop while still wet).
    assert -op.total_requested_volume == pytest.approx(removed, abs=tol)

    # 4. Both regimes exercised as the pond drained.
    depths = [r["Depth_m"] for r in capture_log]
    assert max(depths) > op.d_trans, "never sampled the orifice regime"
    assert min(depths) < op.d_trans, "never dropped into the weir regime"

    # 5. Capture is non-negative and actually happened.
    assert all(r["Captured_Q_cms"] >= 0.0 for r in capture_log)
    assert max(r["Captured_Q_cms"] for r in capture_log) > 0.0

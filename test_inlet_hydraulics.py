"""Unit tests for the importable parts of stormwater_inlets.

These exercise the pure asset/hydraulics logic and the CSV-contract columns
without building a live ANUGA domain, so they run fast. Importing the module
under test still imports ANUGA (a dependency of the operator class), but the
simulation itself only runs under ``if __name__ == "__main__"`` and is not
triggered here.

Run with::

    python -m pytest test_inlet_hydraulics.py
"""

import math
import os

import pytest

import stormwater_inlets as sim
import stormwater_inlet_simulation as sim_run   # the experiment script (PIT_PLACEMENTS)

_CONFIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config")

OP = sim.Depth_driven_inlet_operator
C_W, C_O, G = 1.66, 0.67, 9.81
# Representative inlet geometry (Grate_600x600 clear values).
A, P = 0.21, 2.40


# --------------------------------------------------------------------------- #
# Inlet_specification                                                         #
# --------------------------------------------------------------------------- #

def test_no_blockage_uses_full_geometry():
    spec = sim.Inlet_specification("S", clear_area=0.5, effective_perimeter=3.0)
    assert spec.operational_area == pytest.approx(0.5)
    assert spec.operational_perimeter == pytest.approx(3.0)


def test_blockage_derates_area_and_perimeter():
    spec = sim.Inlet_specification("S", 0.5, 3.0, blockage_factor=0.25)
    assert spec.operational_area == pytest.approx(0.5 * 0.75)
    assert spec.operational_perimeter == pytest.approx(3.0 * 0.75)


def test_full_blockage_zeroes_geometry():
    spec = sim.Inlet_specification("S", 0.5, 3.0, blockage_factor=1.0)
    assert spec.operational_area == pytest.approx(0.0)
    assert spec.operational_perimeter == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# transition_depth                                                            #
# --------------------------------------------------------------------------- #

def test_zero_perimeter_returns_zero():
    assert OP.transition_depth(0.5, 0.0, C_W, C_O) == 0.0


def test_weir_and_orifice_agree_at_transition_depth():
    d = OP.transition_depth(A, P, C_W, C_O, G)
    weir = C_W * P * d ** 1.5
    orifice = C_O * A * math.sqrt(2 * G * d)
    assert weir == pytest.approx(orifice)


# --------------------------------------------------------------------------- #
# capture_discharge                                                           #
# --------------------------------------------------------------------------- #

def cap(depth):
    return OP.capture_discharge(depth, A, P, C_W, C_O, G)


def test_negligible_depth_captures_nothing():
    assert cap(0.0) == 0.0
    assert cap(1e-5) == 0.0


def test_weir_regime_below_transition():
    depth = OP.transition_depth(A, P, C_W, C_O, G) * 0.5
    assert cap(depth) == pytest.approx(C_W * P * depth ** 1.5)


def test_orifice_regime_above_transition():
    depth = OP.transition_depth(A, P, C_W, C_O, G) * 2.0
    assert cap(depth) == pytest.approx(C_O * A * math.sqrt(2 * G * depth))


def test_law_is_continuous_at_transition():
    d_trans = OP.transition_depth(A, P, C_W, C_O, G)
    assert cap(d_trans - 1e-6) == pytest.approx(cap(d_trans + 1e-6), abs=1e-4)


def test_capture_increases_with_depth():
    caps = [cap(d) for d in (0.05, 0.1, 0.2, 0.4, 0.8)]
    assert caps == sorted(caps)
    assert all(b > a for a, b in zip(caps, caps[1:]))


def test_precomputed_d_trans_matches_internal():
    d_trans = OP.transition_depth(A, P, C_W, C_O, G)
    depth = 0.15
    with_arg = OP.capture_discharge(depth, A, P, C_W, C_O, G, d_trans)
    without_arg = OP.capture_discharge(depth, A, P, C_W, C_O, G)
    assert with_arg == pytest.approx(without_arg, abs=1e-12)


# --------------------------------------------------------------------------- #
# INLET_LIBRARY                                                               #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("key", ["Grate_600x600", "Grate_900x900", "Lintel_1.2m",
                                 "Lintel_2.4m", "Combo_1.2m_G600", "Combo_2.4m_G900"])
def test_expected_spec_present(key):
    assert key in sim.INLET_LIBRARY


@pytest.mark.parametrize("key", list(sim.INLET_LIBRARY))
def test_specs_have_positive_geometry(key):
    spec = sim.INLET_LIBRARY[key]
    assert spec.clear_area > 0.0
    assert spec.effective_perimeter > 0.0


# --------------------------------------------------------------------------- #
# Network data reporting (domain-free)                                        #
# --------------------------------------------------------------------------- #

CSV_COLUMNS = ["Asset_ID", "Time_s", "Depth_m", "Approach_Q_cms",
               "Captured_Q_cms", "Bypass_Q_cms",
               "Cum_Inflow_m3", "Cum_Captured_m3", "Cum_Bypassed_m3"]


def test_unknown_spec_raises_keyerror():
    # add_inlet validates the spec key before touching the domain, so a
    # placeholder domain is enough to reach the check.
    net = sim.Stormwater_inlet_network(domain=None)
    with pytest.raises(KeyError):
        net.add_inlet("Pit_X", 1.0, 2.0, "NoSuchSpec")


def test_to_dataframe_empty_when_no_log():
    net = sim.Stormwater_inlet_network(domain=None)
    assert net.to_dataframe("missing").empty


def test_to_dataframe_matches_csv_contract():
    net = sim.Stormwater_inlet_network(domain=None)
    net.logs["Pit_A"] = [{
        "Time_s": 10.0, "Depth_m": 0.3, "Approach_Q_cms": 0.5,
        "Captured_Q_cms": 0.2, "Bypass_Q_cms": 0.3,
        "Cum_Inflow_m3": 5.0, "Cum_Captured_m3": 2.0, "Cum_Bypassed_m3": 3.0,
    }]
    df = net.to_dataframe("Pit_A")
    assert list(df.columns) == CSV_COLUMNS
    assert df.iloc[0]["Asset_ID"] == "Pit_A"


# --------------------------------------------------------------------------- #
# TOML config loaders                                                         #
# --------------------------------------------------------------------------- #

def test_load_inlet_library_roundtrip(tmp_path):
    p = tmp_path / "lib.toml"
    p.write_text(
        "[inlets.MyGrate]\n"
        "clear_area = 0.5\n"
        "effective_perimeter = 3.0\n"
    )
    lib = sim.load_inlet_library(str(p))
    assert set(lib) == {"MyGrate"}
    spec = lib["MyGrate"]
    assert isinstance(spec, sim.Inlet_specification)
    assert spec.clear_area == pytest.approx(0.5)
    assert spec.effective_perimeter == pytest.approx(3.0)


def test_load_inlet_library_missing_key_raises(tmp_path):
    p = tmp_path / "bad.toml"
    p.write_text("[inlets.Broken]\nclear_area = 0.5\n")   # no effective_perimeter
    with pytest.raises(ValueError):
        sim.load_inlet_library(str(p))


def test_load_pit_placements_roundtrip(tmp_path):
    p = tmp_path / "pits.toml"
    p.write_text(
        '[[pits]]\nid = "P1"\nx = 1.0\ny = 2.0\nspec = "MyGrate"\n'
        'radius = 2.0\nblockage = 0.3\n'
    )
    pits = sim.load_pit_placements(str(p))
    assert len(pits) == 1
    assert pits[0] == {"id": "P1", "x": 1.0, "y": 2.0, "spec": "MyGrate",
                       "radius": 2.0, "blockage": 0.3}


def test_load_pit_placements_missing_required_raises(tmp_path):
    p = tmp_path / "pits.toml"
    p.write_text('[[pits]]\nid = "P1"\nx = 1.0\n')   # missing y, spec
    with pytest.raises(ValueError):
        sim.load_pit_placements(str(p))


def test_shipped_config_files_match_builtins():
    """The example TOMLs in config/ should reproduce the built-in library/keys."""
    lib = sim.load_inlet_library(os.path.join(_CONFIG_DIR, "inlet_library.toml"))
    assert set(lib) == set(sim.INLET_LIBRARY)
    for name, spec in lib.items():
        assert spec.clear_area == pytest.approx(sim.INLET_LIBRARY[name].clear_area)
        assert spec.effective_perimeter == pytest.approx(
            sim.INLET_LIBRARY[name].effective_perimeter)

    pits = sim.load_pit_placements(os.path.join(_CONFIG_DIR, "pit_placements.toml"))
    assert [p["id"] for p in pits] == [p["id"] for p in sim_run.PIT_PLACEMENTS]


def test_network_uses_provided_library():
    """add_inlet resolves spec keys against the network's library, not the global."""
    custom = {"OnlyThis": sim.Inlet_specification("OnlyThis", 0.1, 1.0)}
    net = sim.Stormwater_inlet_network(domain=None, library=custom)
    with pytest.raises(KeyError):                 # built-in key absent from custom lib
        net.add_inlet("p", 0.0, 0.0, "Grate_600x600")

"""Sloped-plane stormwater inlet experiment.

A runnable, follow-along example built on the reusable toolkit in
``stormwater_inlets.py``: it puts a row of inlets down a 1%-sloping channel,
drives a steady inflow over them, evolves the 2-D shallow-water flow, and writes
one ``hydrograph_<Asset_ID>.csv`` per inlet plus ``sloped_inlet_experiment.sww``.

Run it::

    python stormwater_inlet_simulation.py
    python stormwater_inlet_simulation.py --library config/inlet_library.toml \
                                          --placements config/pit_placements.toml
    python stormwater_inlet_simulation.py --finaltime 60 --yieldstep 5

The inlet specs and placements below are the built-in defaults; pass TOML files
to override them (see config/*.toml).
"""

import anuga

import stormwater_inlets as si


# --------------------------------------------------------------------------- #
# Experiment definition                                                       #
# --------------------------------------------------------------------------- #

# Six pits down the channel centreline (y = 10 m). "radius" is the circular
# footprint (m) over which depth is sampled; "blockage" derates the inlet's clear
# area & perimeter (0.0 = clear .. 1.0 = fully blocked, optional, default 0.0).
PIT_PLACEMENTS = [
    {"id": "Pit_01_SmallGrate", "x": 15.0, "y": 10.0, "spec": "Grate_600x600",   "radius": 1.5, "blockage": 0.0},
    {"id": "Pit_02_LargeGrate", "x": 28.0, "y": 10.0, "spec": "Grate_900x900",   "radius": 1.5, "blockage": 0.2},
    {"id": "Pit_03_ShortLintel", "x": 42.0, "y": 10.0, "spec": "Lintel_1.2m",    "radius": 1.5, "blockage": 0.4},
    {"id": "Pit_04_LongLintel",  "x": 56.0, "y": 10.0, "spec": "Lintel_2.4m",    "radius": 1.5, "blockage": 0.5},
    {"id": "Pit_05_ComboSmall",  "x": 70.0, "y": 10.0, "spec": "Combo_1.2m_G600", "radius": 1.5, "blockage": 0.6},
    {"id": "Pit_06_ComboLarge",  "x": 85.0, "y": 10.0, "spec": "Combo_2.4m_G900", "radius": 1.5, "blockage": 0.8},
]


def build_domain():
    """Create the 100 x 20 m sloped-plane domain with inflow/outflow boundaries."""
    print("Setting up ANUGA Sloped Plane Domain...")

    # 100 m long, 20 m wide rectangular corridor (clockwise boundary polygon).
    bounding_polygon = [[0.0, 0.0], [100.0, 0.0], [100.0, 20.0], [0.0, 20.0]]
    boundary_tags = {'bottom': [0], 'right': [1], 'top': [2], 'left': [3]}

    domain = anuga.create_domain_from_regions(
        bounding_polygon, boundary_tags=boundary_tags,
        maximum_triangle_area=2.0)
    domain.set_name("sloped_inlet_experiment")

    # 1% bed slope falling left -> right; dry start; concrete-rough friction.
    domain.set_quantity('elevation', function=lambda x, y: 1.0 - (x * 0.01))
    domain.set_quantity('stage', expression='elevation')
    domain.set_quantity('friction', 0.015)

    # Steady inflow on the left (stage 1.35 => ~0.35 m deep at x=0), free outflow
    # on the right, solid walls top and bottom.
    domain.set_boundary({
        'left':   anuga.Dirichlet_boundary([1.35, 0.0, 0.0]),
        'right':  anuga.Transmissive_boundary(domain),
        'top':    anuga.Reflective_boundary(domain),
        'bottom': anuga.Reflective_boundary(domain),
    })
    return domain


def build_network(domain, pit_placements, library=None):
    """Register the configured inlets onto the domain and return the network."""
    network = si.Stormwater_inlet_network(domain, library=library)
    for pit in pit_placements:
        network.add_inlet(pit["id"], pit["x"], pit["y"], pit["spec"],
                          blockage_factor=pit.get("blockage", 0.0),
                          radius=pit.get("radius", 1.5),
                          use_max_depth=pit.get("use_max_depth", si.USE_MAX_DEPTH))
    print(f"Registered {len(network.inlets)} inlets.")
    return network


def print_summary(network, pit_placements, write_csv=True):
    """Print the steady-state results table and optionally dump per-inlet CSVs."""
    header = (f"{'Asset ID':<18} | {'Type':<15} | {'Block':<5} | {'Depth (m)':<9} | "
              f"{'Q_In (L/s)':<10} | {'Bypass (L/s)':<12} | "
              f"{'Cum_Captured (m3)':<17} | {'Cum_Inflow (m3)'}")
    width = len(header)

    print("\n" + "=" * width)
    print(f"{'STORMWATER INLET EXPERIMENT RESULTS SUMMARY':^{width}}")
    print("=" * width)
    print(header)
    print("-" * width)

    for pit in pit_placements:
        df = network.to_dataframe(pit["id"])
        if df.empty:
            continue
        final = df.iloc[-1]   # steady-state final row
        print(f"{pit['id']:<18} | {pit['spec']:<15} | {pit.get('blockage', 0.0):5.2f} | "
              f"{final['Depth_m']:9.3f} | {final['Captured_Q_cms'] * 1000.0:10.1f} | "
              f"{final['Bypass_Q_cms'] * 1000.0:12.1f} | "
              f"{final['Cum_Captured_m3']:17.2f} | {final['Cum_Inflow_m3']:.2f}")
        if write_csv:
            df.to_csv(f"hydrograph_{pit['id']}.csv", index=False)

    print("=" * width)
    if write_csv:
        print("Individual hydrograph log CSVs have been saved to your workspace.")


def run_experiment(yieldstep=10, finaltime=120, write_csv=True,
                   library_path=None, placements_path=None):
    """Build the domain, register inlets, evolve, and report results.

    library_path / placements_path: optional TOML files overriding the built-in
    INLET_LIBRARY / PIT_PLACEMENTS.
    """
    library = si.load_inlet_library(library_path) if library_path else None
    placements = si.load_pit_placements(placements_path) if placements_path else PIT_PLACEMENTS
    if library_path:
        print(f"Loaded {len(library)} inlet specs from {library_path}")
    if placements_path:
        print(f"Loaded {len(placements)} pit placements from {placements_path}")

    domain = build_domain()
    network = build_network(domain, placements, library=library)

    print("\nStarting simulation loop...")
    for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        print(f"Simulation Time: {t:.1f}s")
        domain.report_water_volume_statistics()

    print_summary(network, placements, write_csv=write_csv)
    return network


def _parse_args(argv=None):
    import argparse
    p = argparse.ArgumentParser(
        description="Run the ANUGA stormwater-inlet experiment.")
    p.add_argument("--library", metavar="TOML",
                   help="inlet asset catalogue TOML (default: built-in INLET_LIBRARY)")
    p.add_argument("--placements", metavar="TOML",
                   help="pit placements TOML (default: built-in PIT_PLACEMENTS)")
    p.add_argument("--finaltime", type=float, default=120.0,
                   help="simulation end time in seconds (default: 120)")
    p.add_argument("--yieldstep", type=float, default=10.0,
                   help="reporting interval in seconds (default: 10)")
    return p.parse_args(argv)


if __name__ == "__main__":
    args = _parse_args()
    run_experiment(yieldstep=args.yieldstep, finaltime=args.finaltime,
                   library_path=args.library, placements_path=args.placements)

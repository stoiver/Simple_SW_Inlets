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
import argparse
import stormwater_inlets as si


# --------------------------------------------------------------------------- #
# Experiment definition                                                       #
# --------------------------------------------------------------------------- #

# Six pits down the channel centreline (y = 10 m). "radius" is the circular
# footprint (m) over which depth is sampled; "blockage" derates the inlet's clear
# area & perimeter (0.0 = clear .. 1.0 = fully blocked, optional, default 0.0);
# "use_max_depth" drives the capture law from the max footprint depth (True) or
# the region-averaged depth (False).
blockage = 0.0          # applied to every pit below
use_max_depth = True   # applied to every pit below

PIT_PLACEMENTS = [
    {"id": "Pit_01_SmallGrate", "x": 15.0, "y": 10.0, "spec": "Grate_600x600",   "radius": 1.5, "blockage": blockage, "use_max_depth": use_max_depth},
    {"id": "Pit_02_LargeGrate", "x": 28.0, "y": 10.0, "spec": "Grate_900x900",   "radius": 1.5, "blockage": blockage, "use_max_depth": use_max_depth},
    {"id": "Pit_03_ShortLintel", "x": 42.0, "y": 10.0, "spec": "Lintel_1.2m",    "radius": 1.5, "blockage": blockage, "use_max_depth": use_max_depth},
    {"id": "Pit_04_LongLintel",  "x": 56.0, "y": 10.0, "spec": "Lintel_2.4m",    "radius": 1.5, "blockage": blockage, "use_max_depth": use_max_depth},
    {"id": "Pit_05_ComboSmall",  "x": 70.0, "y": 10.0, "spec": "Combo_1.2m_G600", "radius": 1.5, "blockage": blockage, "use_max_depth": use_max_depth},
    {"id": "Pit_06_ComboLarge",  "x": 85.0, "y": 10.0, "spec": "Combo_2.4m_G900", "radius": 1.5, "blockage": blockage, "use_max_depth": use_max_depth},
]



# Everything below runs only when the script is executed directly, so importing
# this module (e.g. to reuse PIT_PLACEMENTS / INLET_LIBRARY) does not start a run.
if __name__ == "__main__":

    # --- Command-line options ----------------------------------------------
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
    args = p.parse_args()

    yieldstep = args.yieldstep
    finaltime = args.finaltime
    library_path = args.library
    placements_path = args.placements
    write_csv = True   # set False to skip writing hydrograph CSVs

    # --- Build the 100 x 20 m sloped-plane domain --------------------------
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

    # --- Register inlets (TOML files override the built-in defaults) --------
    library = si.load_inlet_library(library_path) if library_path else None
    placements = si.load_pit_placements(placements_path) if placements_path else PIT_PLACEMENTS
    if library_path:
        print(f"Loaded {len(library)} inlet specs from {library_path}")
    if placements_path:
        print(f"Loaded {len(placements)} pit placements from {placements_path}")

    network = si.build_network(domain, placements, library=library)

    # --- Evolve and report -------------------------------------------------
    print("\nStarting simulation loop...")
    for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        domain.print_timestepping_statistics()
        domain.report_water_volume_statistics()

    si.print_summary(network, placements, write_csv=write_csv)


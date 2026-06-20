import os
import anuga
import numpy as np
import pandas as pd
from anuga.structures.inlet_operator import Inlet_operator

# ==========================================
# 1. CORE STORMWATER ASSET LIBRARY CLASSES
# ==========================================

class Inlet_specification:
    """Defines the geometric parameters of a standard inlet asset."""
    def __init__(self, name, clear_area, effective_perimeter, blockage_factor=0.0):
        self.name = name
        self.clear_area = clear_area  
        self.effective_perimeter = effective_perimeter  
        self.blockage_factor = blockage_factor  

    @property
    def operational_area(self):
        return self.clear_area * (1.0 - self.blockage_factor)

    @property
    def operational_perimeter(self):
        return self.effective_perimeter * (1.0 - self.blockage_factor)


# Structural Catalog Definitions
INLET_LIBRARY = {
    "Grate_600x600": Inlet_specification("Grate_600x600", 0.21, 2.40),
    "Grate_900x900": Inlet_specification("Grate_900x900", 0.48, 3.60),
    "Lintel_1.2m":   Inlet_specification("Lintel_1.2m",   0.18, 1.20),
    "Lintel_2.4m":   Inlet_specification("Lintel_2.4m",   0.36, 2.40),
    "Combo_1.2m_G600": Inlet_specification("Combo_1.2m_G600", 0.39, 3.00),
    "Combo_2.4m_G900": Inlet_specification("Combo_2.4m_G900", 0.84, 5.10)
}


class Depth_driven_inlet_operator(Inlet_operator):
    """Stormwater inlet implemented as an ANUGA Inlet_operator.

    The capture discharge is a function of the inlet's (region-averaged) ponded
    depth, using dual-regime weir/orifice hydraulics, and is returned from
    update_Q() as a negative discharge (water leaving the domain). The parent
    Inlet_operator distributes that discharge over the inlet region and enforces
    mass balance, so we never touch the stage array directly.
    """
    def __init__(self, domain, region, spec, capture_log=None,
                 C_w=1.66, C_o=0.67, label=None, **kwargs):
        Inlet_operator.__init__(self, domain, region, Q=0.0, label=label, **kwargs)

        self.spec = spec
        self.C_w = C_w
        self.C_o = C_o
        self.g = 9.81

        # Pre-calculate the weir->orifice transition depth (where both laws agree)
        A = self.spec.operational_area
        P = self.spec.operational_perimeter
        self.d_trans = self.transition_depth(A, P, self.C_w, self.C_o, self.g)

        self.capture_log = capture_log
        self.total_volume_captured = 0.0
        self.total_volume_bypassed = 0.0

    @staticmethod
    def transition_depth(A, P, C_w, C_o, g=9.81):
        """Ponded depth at which the weir and orifice laws give equal discharge.

        Returns 0.0 for a degenerate inlet (no perimeter), so capture is always
        in the orifice regime there.
        """
        if P <= 0:
            return 0.0
        # Equating C_w*P*d^1.5 == C_o*A*sqrt(2g*d) and solving for d: the
        # d^0.5 factors cancel, leaving d = (C_o*A*sqrt(2g)) / (C_w*P).
        return (C_o * A * np.sqrt(2 * g)) / (C_w * P)

    @staticmethod
    def capture_discharge(depth, A, P, C_w, C_o, g=9.81, d_trans=None):
        """Positive capture discharge (m3/s) for a ponded depth (dual-regime law).

        Weir flow (C_w * P * depth^1.5) below the transition depth, orifice flow
        (C_o * A * sqrt(2g*depth)) above it. Negligible depths capture nothing.
        Pure function of the arguments — no ANUGA domain required — so it carries
        the core hydraulics and is what the unit tests exercise.
        """
        if depth <= 1e-4:
            return 0.0
        if d_trans is None:
            d_trans = Depth_driven_inlet_operator.transition_depth(A, P, C_w, C_o, g)
        if depth < d_trans:
            return C_w * P * (depth ** 1.5)              # Weir flow
        return C_o * A * np.sqrt(2 * g * depth)          # Orifice flow

    def update_Q(self, t):
        """Capture discharge as a function of inlet state (negative => extraction).

        t is unused: Q depends on the current ponded depth, not on time.
        """
        depth = self.inlet.get_average_depth()
        A = self.spec.operational_area
        P = self.spec.operational_perimeter
        return -self.capture_discharge(depth, A, P, self.C_w, self.C_o,
                                       self.g, self.d_trans)

    def __call__(self):
        # Parent applies the discharge (with mass-balance clamping) and sets self.applied_Q.
        Inlet_operator.__call__(self)

        if self.capture_log is None:
            return

        dt = self.domain.get_timestep()
        depth = self.inlet.get_average_depth()

        # applied_Q is negative for extraction; report the realised capture as positive.
        captured_Q = max(0.0, -self.applied_Q)

        # Approach discharge across the inlet footprint, from the region-averaged
        # specific discharge and a representative width (sqrt of the inlet area).
        uh = self.inlet.get_average_xmom()
        vh = self.inlet.get_average_ymom()
        specific_discharge = np.sqrt(uh ** 2 + vh ** 2)
        width = np.sqrt(self.inlet.get_area())
        Q_approach = specific_discharge * width

        bypass_Q = max(0.0, Q_approach - captured_Q)

        self.total_volume_captured += captured_Q * dt
        self.total_volume_bypassed += bypass_Q * dt

        self.capture_log.append({
            "Time_s": self.domain.get_time(),
            "Depth_m": depth,
            "Approach_Q_cms": Q_approach,
            "Captured_Q_cms": captured_Q,
            "Bypass_Q_cms": bypass_Q,
            "Cum_Captured_m3": self.total_volume_captured,
            "Cum_Bypassed_m3": self.total_volume_bypassed
        })


class Stormwater_inlet_network:
    """Manages collections of point-based inlets and handles data reporting."""
    def __init__(self, domain):
        self.domain = domain
        self.inlets = {}
        self.logs = {}

    def add_inlet(self, asset_id, x, y, spec_key, blockage_factor=0.0, radius=1.5):
        if spec_key not in INLET_LIBRARY:
            raise KeyError(f"Asset spec '{spec_key}' not found.")

        base_spec = INLET_LIBRARY[spec_key]
        spec = Inlet_specification(base_spec.name, base_spec.clear_area, base_spec.effective_perimeter, blockage_factor)

        # Define the inlet footprint as a small circular region around the pit.
        region = anuga.Region(self.domain, center=[x, y], radius=radius)

        self.logs[asset_id] = []
        operator = Depth_driven_inlet_operator(self.domain, region, spec,
                                            capture_log=self.logs[asset_id], label=asset_id)
        self.inlets[asset_id] = operator
        return operator

    def to_dataframe(self, asset_id):
        if asset_id not in self.logs or not self.logs[asset_id]:
            return pd.DataFrame()
        df = pd.DataFrame(self.logs[asset_id])
        df.insert(0, "Asset_ID", asset_id)
        return df


# ==========================================
# 2. RUNTIME SIMULATION EXPERIMENT SETUP
# ==========================================

# Place 6 pits spaced out down the sloping channel at y = 10.0m (mid-channel).
# "radius" sets the circular inlet-footprint region (m) over which depth is sampled.
PIT_PLACEMENTS = [
    {"id": "Pit_01_SmallGrate", "x": 15.0, "y": 10.0, "spec": "Grate_600x600",   "radius": 1.5},
    {"id": "Pit_02_LargeGrate", "x": 28.0, "y": 10.0, "spec": "Grate_900x900",   "radius": 1.5},
    {"id": "Pit_03_ShortLintel", "x": 42.0, "y": 10.0, "spec": "Lintel_1.2m",    "radius": 1.5},
    {"id": "Pit_04_LongLintel",  "x": 56.0, "y": 10.0, "spec": "Lintel_2.4m",    "radius": 1.5},
    {"id": "Pit_05_ComboSmall",  "x": 70.0, "y": 10.0, "spec": "Combo_1.2m_G600", "radius": 1.5},
    {"id": "Pit_06_ComboLarge",  "x": 85.0, "y": 10.0, "spec": "Combo_2.4m_G900", "radius": 1.5}
]


def build_domain():
    """Create the 100x20 m sloped-plane domain with inflow/outflow boundaries."""
    print("Setting up ANUGA Sloped Plane Domain...")

    # Create a 100m long, 20m wide domain corridor clockwise
    bounding_polygon = [[0.0, 0.0], [100.0, 0.0], [100.0, 20.0], [0.0, 20.0]]

    # Explicit sequential list indices matching the clockwise boundary polygon edges
    boundary_tags = {
        'bottom': [0],
        'right': [1],
        'top': [2],
        'left': [3]
    }

    # Generate the domain mesh structure correctly
    domain = anuga.create_domain_from_regions(
        bounding_polygon,
        boundary_tags=boundary_tags,
        maximum_triangle_area=2.0
    )

    domain.set_name("sloped_inlet_experiment")

    # Define a sloping bed terrain (1% gradient sloping down from Left to Right)
    def elevation_function(x, y):
        return 1.0 - (x * 0.01)

    domain.set_quantity('elevation', function=elevation_function)
    domain.set_quantity('stage', expression='elevation')  # Initially dry channel
    domain.set_quantity('friction', 0.015)                # Concrete rough lining profile

    # Establish boundary conditions
    constant_inflow_stage = 1.35  # Results in roughly ~0.35m deep water entering at x=0
    bd_inflow = anuga.Dirichlet_boundary([constant_inflow_stage, 0.0, 0.0])
    bd_outflow = anuga.Transmissive_boundary(domain)
    bd_wall = anuga.Reflective_boundary(domain)

    domain.set_boundary({
        'left': bd_inflow,
        'right': bd_outflow,
        'top': bd_wall,
        'bottom': bd_wall
    })

    return domain


def build_network(domain, pit_placements=PIT_PLACEMENTS):
    """Register the configured inlets onto the domain and return the network."""
    network = Stormwater_inlet_network(domain)
    for pit in pit_placements:
        network.add_inlet(pit["id"], pit["x"], pit["y"], pit["spec"],
                          blockage_factor=0.0, radius=pit.get("radius", 1.5))
    print(f"Registered {len(network.inlets)} distinct configurations in path line.")
    return network


def print_summary(network, pit_placements=PIT_PLACEMENTS, write_csv=True):
    """Print the steady-state results table and optionally dump per-inlet CSVs."""
    print("\n" + "="*70)
    print(f"{'STORMWATER INLET EXPERIMENT RESULTS SUMMARY':^70}")
    print("="*70)
    print(f"{'Asset ID':<18} | {'Type':<15} | {'Depth (m)':<9} | {'Q_In (L/s)':<10} | {'Bypass (L/s)'}")
    print("-"*70)

    for pit in pit_placements:
        asset_id = pit["id"]
        df = network.to_dataframe(asset_id)

        if not df.empty:
            # Extract the steady-state final row values
            final_row = df.iloc[-1]
            depth = final_row["Depth_m"]
            q_cap_lps = final_row["Captured_Q_cms"] * 1000.0  # Convert to Litres/sec
            q_byp_lps = final_row["Bypass_Q_cms"] * 1000.0   # Convert to Litres/sec

            print(f"{asset_id:<18} | {pit['spec']:<15} | {depth:9.3f} | {q_cap_lps:10.1f} | {q_byp_lps:11.1f}")

            if write_csv:
                # Save out to CSV cleanly
                filename = f"hydrograph_{asset_id}.csv"
                df.to_csv(filename, index=False)

    print("="*70)
    if write_csv:
        print("Individual hydrograph log CSVs have been saved to your workspace.")


def run_experiment(yieldstep=10, finaltime=120, write_csv=True):
    """Build the domain, register inlets, evolve, and report results."""
    domain = build_domain()
    network = build_network(domain)

    print("\nStarting simulation loop...")
    # Run for `finaltime` seconds to allow conditions to stabilise across the slope
    for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
        print(f"Simulation Time: {t:.1f}s")

    print_summary(network, write_csv=write_csv)
    return network


if __name__ == "__main__":
    run_experiment()

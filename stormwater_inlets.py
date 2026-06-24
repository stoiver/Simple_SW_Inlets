"""Reusable toolkit for depth-driven stormwater inlet modelling in ANUGA.

This module holds the building blocks — the inlet asset specification + catalogue,
the TOML loaders, the depth-driven capture operators (serial and MPI), and the
inlet network manager. It contains no experiment-specific setup; see
``stormwater_inlet_simulation.py`` for a runnable example that uses these pieces.

The capture hydraulics (dual-regime weir/orifice) are documented in
``docs/HYDRAULICS.md``.
"""

import anuga
import numpy as np
import pandas as pd

try:
    import tomllib            # Python 3.11+
except ModuleNotFoundError:   # pragma: no cover - fallback for older Pythons
    import tomli as tomllib
# NB: import the Inlet_operator *class* (subclassable). The top-level
# anuga.Inlet_operator is a factory *function* and cannot be subclassed.
from anuga.structures.inlet_operator import Inlet_operator

# The parallel operator base is only needed for MPI (distributed-domain) runs and
# pulls in ANUGA's parallel stack; degrade gracefully if it is unavailable so the
# serial path still imports.
try:
    from anuga.parallel.parallel_inlet_operator import Parallel_Inlet_operator
    _HAVE_PARALLEL = True
except Exception:
    Parallel_Inlet_operator = object
    _HAVE_PARALLEL = False

# Default for the operators' `use_max_depth` argument: drive the capture law from
# the max ponded depth over the inlet footprint (True) or the region-averaged
# depth (False). Max is more representative when the inlet sits in a local
# depression; both agree for a uniform pond. Override per inlet via the operator
# / add_inlet argument; this constant only sets the default.
USE_MAX_DEPTH = True

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


def load_inlet_library(path):
    """Load an inlet asset catalogue from a TOML file.

    Expected layout (one table per named inlet)::

        [inlets.Grate_600x600]
        clear_area = 0.21
        effective_perimeter = 2.40

    Returns a {name: Inlet_specification} dict suitable for
    ``Stormwater_inlet_network(domain, library=...)``.
    """
    with open(path, "rb") as f:
        data = tomllib.load(f)

    inlets = data.get("inlets", {})
    if not inlets:
        raise ValueError(f"No [inlets.*] tables found in {path}")

    library = {}
    for name, props in inlets.items():
        try:
            library[name] = Inlet_specification(
                name, props["clear_area"], props["effective_perimeter"])
        except KeyError as e:
            raise ValueError(
                f"Inlet '{name}' in {path} is missing required key {e}") from e
    return library


def load_pit_placements(path):
    """Load pit placements from a TOML file.

    Expected layout (an array of tables)::

        [[pits]]
        id = "Pit_01_SmallGrate"
        x = 15.0
        y = 10.0
        spec = "Grate_600x600"
        radius = 1.5      # optional
        blockage = 0.0    # optional

    Returns a list of placement dicts (the same shape as a hard-coded placement
    list).
    """
    with open(path, "rb") as f:
        data = tomllib.load(f)

    pits = data.get("pits", [])
    if not pits:
        raise ValueError(f"No [[pits]] entries found in {path}")

    required = ("id", "x", "y", "spec")
    for i, pit in enumerate(pits):
        missing = [k for k in required if k not in pit]
        if missing:
            raise ValueError(f"Pit #{i} in {path} is missing keys: {missing}")
    return pits


class _Depth_driven_capture_mixin:
    """Shared weir/orifice capture hydraulics and hydrograph logging.

    Mixed into both the serial and parallel operators so the physics
    (transition_depth / capture_discharge) and the capture-log record are defined
    once. The host operator supplies the ponded depth and momentum samples —
    per-subdomain (local) in serial, globally reduced across ranks in parallel.
    """

    def _init_capture(self, spec, C_w, C_o, capture_log, use_max_depth=USE_MAX_DEPTH):
        self.spec = spec
        self.C_w = C_w
        self.C_o = C_o
        self.g = 9.81
        self.use_max_depth = use_max_depth

        # Pre-calculate the weir->orifice transition depth (where both laws agree)
        A = self.spec.operational_area
        P = self.spec.operational_perimeter
        self.d_trans = self.transition_depth(A, P, self.C_w, self.C_o, self.g)

        self.capture_log = capture_log
        self.total_volume_inflow = 0.0
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
            d_trans = _Depth_driven_capture_mixin.transition_depth(A, P, C_w, C_o, g)
        if depth < d_trans:
            return C_w * P * (depth ** 1.5)              # Weir flow
        return C_o * A * np.sqrt(2 * g * depth)          # Orifice flow

    def _capture_Q(self, depth):
        """Signed capture discharge (negative => extraction) for a ponded depth."""
        A = self.spec.operational_area
        P = self.spec.operational_perimeter
        return -self.capture_discharge(depth, A, P, self.C_w, self.C_o,
                                       self.g, self.d_trans)

    def _log_capture(self, depth, uh, vh, area, dt):
        """Append one hydrograph record from the sampled inlet state.

        applied_Q (set by the parent operator) is negative for extraction; the
        realised capture is reported as positive. Approach discharge is estimated
        from the region-averaged specific discharge and a representative width
        (sqrt of the inlet area), and bypass is whatever the approach flow exceeds
        the captured flow by.
        """
        captured_Q = max(0.0, -self.applied_Q)
        specific_discharge = np.sqrt(uh ** 2 + vh ** 2)
        width = np.sqrt(area)
        Q_approach = specific_discharge * width
        bypass_Q = max(0.0, Q_approach - captured_Q)

        self.total_volume_inflow += Q_approach * dt
        self.total_volume_captured += captured_Q * dt
        self.total_volume_bypassed += bypass_Q * dt

        self.capture_log.append({
            "Time_s": self.domain.get_time(),
            "Depth_m": depth,
            "Approach_Q_cms": Q_approach,
            "Captured_Q_cms": captured_Q,
            "Bypass_Q_cms": bypass_Q,
            "Cum_Inflow_m3": self.total_volume_inflow,
            "Cum_Captured_m3": self.total_volume_captured,
            "Cum_Bypassed_m3": self.total_volume_bypassed
        })


class Depth_driven_inlet_operator(_Depth_driven_capture_mixin, Inlet_operator):
    """Stormwater inlet implemented as a (serial) ANUGA Inlet_operator.

    The capture discharge is a function of the inlet's region-averaged ponded
    depth, using dual-regime weir/orifice hydraulics, and is returned from
    update_Q() as a negative discharge (water leaving the domain). The parent
    Inlet_operator distributes that discharge over the inlet region and enforces
    mass balance, so we never touch the stage array directly.

    The depth/momentum samples here are per-subdomain (local). For an MPI run on
    a distributed domain use Depth_driven_parallel_inlet_operator instead.
    """
    def __init__(self, domain, region, spec, capture_log=None,
                 C_w=1.66, C_o=0.67, label=None, use_max_depth=USE_MAX_DEPTH,
                 **kwargs):
        Inlet_operator.__init__(self, domain, region, Q=0.0, label=label, **kwargs)
        self._init_capture(spec, C_w, C_o, capture_log, use_max_depth)

    def _sample_depth(self):
        """Ponded depth driving the capture law: max over the footprint if
        self.use_max_depth, else the region average."""
        if self.use_max_depth:
            depths = self.inlet.get_depths()
            return float(np.max(depths)) if len(depths) > 0 else 0.0
        return self.inlet.get_average_depth()

    def update_Q(self, t):
        """Capture discharge from the current ponded depth (t unused)."""
        return self._capture_Q(self._sample_depth())

    def __call__(self):
        # Parent applies the discharge (with mass-balance clamping) and sets applied_Q.
        Inlet_operator.__call__(self)
        if self.capture_log is None:
            return
        self._log_capture(self._sample_depth(),
                          self.inlet.get_average_xmom(),
                          self.inlet.get_average_ymom(),
                          self.inlet.get_area(),
                          self.domain.get_timestep())


class Depth_driven_parallel_inlet_operator(_Depth_driven_capture_mixin,
                                           Parallel_Inlet_operator):
    """MPI-safe counterpart of Depth_driven_inlet_operator.

    Two parallel subtleties drive the differences from the serial class:

    * ``Inlet.get_average_depth()`` is per-subdomain (local). When the inlet
      footprint straddles ranks the capture law needs the *global* average, so
      this class uses ``get_global_average_depth()`` (and the ``get_global_*``
      momentum/area reductions for the hydrograph log).
    * ``Parallel_Inlet_operator.__call__`` invokes ``update_Q`` on the master
      rank only (it then broadcasts the resulting volume). The ``get_global_*``
      reductions are collective and must be entered by every rank, so they are
      sampled in ``__call__`` (all ranks) and the master-only ``update_Q`` reads
      the stashed global depth — calling a collective inside ``update_Q`` would
      deadlock the non-master ranks waiting in the broadcast.

    Note: ANUGA does not auto-discover which ranks hold the inlet; for a footprint
    spanning multiple ranks pass ``procs=[...]`` (and ``master_proc``) through, per
    ANUGA's parallel-inlet conventions.
    """
    def __init__(self, domain, region, spec, capture_log=None,
                 C_w=1.66, C_o=0.67, label=None, use_max_depth=USE_MAX_DEPTH,
                 **kwargs):
        Parallel_Inlet_operator.__init__(self, domain, region, Q=0.0,
                                         label=label, **kwargs)
        self._init_capture(spec, C_w, C_o, capture_log, use_max_depth)
        self._global_depth = 0.0

    def _sample_global_depth(self):
        """Global ponded depth across all ranks holding the inlet (collective).

        Honors self.use_max_depth via an MPI max-reduction of each rank's local
        max, else uses the global average. Must be entered by every rank.
        """
        if self.use_max_depth:
            from mpi4py import MPI
            depths = self.inlet.get_depths()
            local_max = float(np.max(depths)) if len(depths) > 0 else 0.0
            return float(MPI.COMM_WORLD.allreduce(local_max, op=MPI.MAX))
        return self.inlet.get_global_average_depth()

    def update_Q(self, t):
        # Master-only; uses the global depth gathered collectively in __call__.
        return self._capture_Q(self._global_depth)

    def __call__(self):
        # Collective: every rank must enter this so the reduction completes.
        self._global_depth = self._sample_global_depth()
        # Master computes Q via update_Q and broadcasts the volume to the others.
        Parallel_Inlet_operator.__call__(self)

        if self.capture_log is None:
            return
        # get_global_* are collective -> sample on every rank before guarding.
        uh = self.inlet.get_global_average_xmom()
        vh = self.inlet.get_global_average_ymom()
        area = self.inlet.get_global_area()
        if self.myid != self.master_proc:
            return
        self._log_capture(self._global_depth, uh, vh, area,
                          self.domain.get_timestep())


class Stormwater_inlet_network:
    """Manages collections of point-based inlets and handles data reporting."""
    def __init__(self, domain, library=None):
        self.domain = domain
        # Catalogue of named inlet specs to resolve spec keys against; defaults to
        # the built-in INLET_LIBRARY but can be a catalogue loaded from a file.
        self.library = library if library is not None else INLET_LIBRARY
        self.inlets = {}
        self.logs = {}

    def add_inlet(self, asset_id, x, y, spec_key, blockage_factor=0.0, radius=1.5,
                  use_max_depth=USE_MAX_DEPTH, **operator_kwargs):
        if spec_key not in self.library:
            raise KeyError(f"Asset spec '{spec_key}' not found.")

        base_spec = self.library[spec_key]
        spec = Inlet_specification(base_spec.name, base_spec.clear_area, base_spec.effective_perimeter, blockage_factor)

        # Define the inlet footprint as a small circular region around the pit.
        region = anuga.Region(self.domain, center=[x, y], radius=radius)

        # Pick the serial or MPI-safe operator based on whether the domain has
        # been distributed. Extra kwargs (e.g. master_proc/procs for parallel)
        # are forwarded to the operator.
        if getattr(self.domain, "parallel", False):
            if not _HAVE_PARALLEL:
                raise RuntimeError(
                    "Domain is distributed but ANUGA's parallel inlet operator "
                    "could not be imported.")
            operator_cls = Depth_driven_parallel_inlet_operator
        else:
            operator_cls = Depth_driven_inlet_operator

        self.logs[asset_id] = []
        operator = operator_cls(self.domain, region, spec,
                                capture_log=self.logs[asset_id], label=asset_id,
                                use_max_depth=use_max_depth, **operator_kwargs)
        self.inlets[asset_id] = operator
        return operator

    def to_dataframe(self, asset_id):
        if asset_id not in self.logs or not self.logs[asset_id]:
            return pd.DataFrame()
        df = pd.DataFrame(self.logs[asset_id])
        df.insert(0, "Asset_ID", asset_id)
        return df

"""MPI test for Depth_driven_parallel_inlet_operator.

Following the anuga_core parallel-test convention, this single file is both:

  * the pytest entry point -- ``test_parallel_inlet_mass_balance`` shells out via
    ``anuga.mpicmd`` to re-run this file under ``mpiexec``; and
  * the MPI worker -- when executed under mpiexec (``__main__``), ``_run_parallel``
    builds a distributed pond, drains it with one parallel inlet, and asserts the
    global water loss equals the operator's captured volume while the inlet
    crosses both flow regimes.

The worker exits non-zero if any assertion fails, so the shelled-out command's
return code is what the pytest assertion checks.

Run with::

    python -m pytest test_parallel_inlet_mass_balance.py
    # or directly:  mpiexec -np 2 python -m mpi4py test_parallel_inlet_mass_balance.py
"""

import os
import shutil
import sys

import pytest

try:
    import mpi4py  # noqa: F401
    _HAVE_MPI4PY = True
except ImportError:
    _HAVE_MPI4PY = False

NPROCS = 2
SIDE = 5.0
INITIAL_DEPTH = 0.30        # > d_trans (~0.156 m) => starts in the orifice regime


# --------------------------------------------------------------------------- #
# pytest entry point: shell out to mpiexec                                     #
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not _HAVE_MPI4PY, reason="requires the mpi4py module")
@pytest.mark.skipif(shutil.which("mpiexec") is None, reason="requires mpiexec")
def test_parallel_inlet_mass_balance():
    import anuga
    cmd = anuga.mpicmd(os.path.abspath(__file__), numprocs=NPROCS)
    assert os.system(cmd) == 0, f"parallel run failed: {cmd}"


# --------------------------------------------------------------------------- #
# MPI worker: runs once per rank under mpiexec                                 #
# --------------------------------------------------------------------------- #

def _run_parallel():
    import warnings
    warnings.simplefilter("ignore")

    # Ensure the module under test is importable regardless of how mpiexec
    # set up sys.path.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    import anuga
    from anuga import distribute, myid, numprocs, finalize
    import stormwater_inlet_simulation as sim

    # Build the (flat, reflective-walled) pond on rank 0, then distribute it.
    if myid == 0:
        domain = anuga.rectangular_cross_domain(10, 10, len1=SIDE, len2=SIDE)
        domain.set_name("parallel_inlet_mass_balance")
        domain.set_store(False)                 # no .sww artifact
        domain.set_quantity("elevation", 0.0)
        domain.set_quantity("friction", 0.0)
        domain.set_quantity("stage", INITIAL_DEPTH)
    else:
        domain = None

    domain = distribute(domain)
    reflective = anuga.Reflective_boundary(domain)
    domain.set_boundary({"left": reflective, "right": reflective,
                         "top": reflective, "bottom": reflective})

    # Footprint may straddle ranks -> let every rank participate in the inlet.
    network = sim.Stormwater_inlet_network(domain)
    op = network.add_inlet("Pit", SIDE / 2, SIDE / 2, "Grate_600x600",
                           radius=1.0, master_proc=0,
                           procs=list(range(numprocs)))

    # get_water_volume() is a global (MPI_Allreduce) quantity.
    volume_before = domain.get_water_volume()
    for _t in domain.evolve(yieldstep=2.0, finaltime=20.0):
        pass
    volume_after = domain.get_water_volume()

    # The master rank holds the authoritative applied_Q and the hydrograph log.
    if myid == op.master_proc:
        removed = volume_before - volume_after

        assert type(op).__name__ == "Depth_driven_parallel_inlet_operator", \
            "network did not select the parallel operator on a distributed domain"
        assert removed > 0.0, "nothing drained"
        assert volume_after > 0.0, "pond drained completely; shorten finaltime"

        # Global mass balance: the inlet is the only sink.
        assert abs(removed - op.total_volume_captured) < 1e-6, \
            f"capture {op.total_volume_captured} != domain loss {removed}"
        assert abs(removed + op.total_applied_volume) < 1e-6

        # Both regimes exercised as the (global) depth drained past d_trans.
        depths = [r["Depth_m"] for r in network.logs["Pit"]]
        assert max(depths) > op.d_trans, "never sampled the orifice regime"
        assert min(depths) < op.d_trans, "never dropped into the weir regime"

        print(f"parallel inlet mass-balance OK on {numprocs} ranks: "
              f"removed={removed:.4f} captured={op.total_volume_captured:.4f} "
              f"mass_err={abs(removed - op.total_volume_captured):.2e}")

    finalize()


if __name__ == "__main__":
    _run_parallel()

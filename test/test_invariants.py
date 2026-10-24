"""Invariants of the architecture: the times since events never go back; without resources to limit it, the root
system does not depend on the time step."""
import pytest

from openalea.archisimple.drawing import totals
from openalea.archisimple.parameter_sets import WHEAT
from openalea.archisimple.run import run


def test_clocks_never_go_back():
    result = run(WHEAT, n_steps=96, record_every=1)
    for before, after in zip(result.steps[:-1], result.steps[1:]):
        common = [v for v in before.vertices if v in after.vertices]
        for name in ("thermal_time_since_primordium_formation", "actual_time_since_primordium_formation"):
            assert all(after.props[name][v] >= before.props[name][v] for v in common), name


def test_the_time_step_does_not_matter():
    """Hourly and 15-min steps over 5 days: the same root system's length within 2 %."""
    lengths = {}
    for time_step in (3600, 900):
        n_steps = int(5 * 86400 / time_step)
        lengths[time_step] = totals(run(WHEAT, n_steps=n_steps, time_step=time_step,
                                        record_every=n_steps).steps[-1])["length"]
    assert lengths[3600] == pytest.approx(lengths[900], rel=0.02)

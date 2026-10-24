"""The root systems' geometry, drawn by rhizodep's turtle (openalea.archisimple.drawing): the laterals leave their
bearer at their drawn angles, and no plant leans to one side (the seminal and nodal roots leave the base all around)."""
import numpy as np
import pytest

from openalea.archisimple.drawing import positions
from openalea.archisimple.parameter_sets import WHEAT
from openalea.archisimple.run import run

SEEDS = (1, 2, 3, 4)


@pytest.fixture(scope="module")
def records():
    return [run(dict(WHEAT, random_choice=seed), n_steps=240, record_every=240).steps[-1] for seed in SEEDS]


@pytest.mark.parametrize("gravitropism", [0., 0.06])
def test_no_lean(records, gravitropism):
    """The elements' mean horizontal position relative to the plant's extent, averaged over seeds, is near 0 (the roots
    leaving the base all in one plane gave 0.23 to 0.33 in x)."""
    leans = []
    for record in records:
        ends = np.array([end for _, end in positions(record, gravitropism=gravitropism).values()])
        leans.append(ends[:, :2].mean(axis=0) / np.abs(ends[:, :2]).max())
    assert (np.abs(np.mean(leans, axis=0)) < 0.15).all(), np.mean(leans, axis=0)


def test_laterals_leave_their_bearer_at_their_angles(records):
    """Without gravitropism, a lateral's first element makes its drawn angle_down with its bearer (45° on a
    first-order root, 70° beyond, ± 10°)."""
    for record in records[:2]:
        xyz = positions(record, gravitropism=0.)
        p = record.props
        heading = {v: (e - s) / np.linalg.norm(e - s) for v, (s, e) in xyz.items() if p["length"][v] > 0}
        laterals = [v for v in heading if record.edge_type[v] == "+" and record.parent[v] in heading
                    and p["root_order"][v] == 2]
        assert laterals
        angles = np.degrees(np.arccos(np.clip([heading[v] @ heading[record.parent[v]] for v in laterals], -1., 1.)))
        drawn = np.array([p["angle_down"][v] for v in laterals])
        np.testing.assert_allclose(angles, drawn, atol=1e-6)

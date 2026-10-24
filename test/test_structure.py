"""The seedling's structure in the MPG (an Axis and an Organ per first-order root, the seminal and nodal roots
branched on the base), the elements' lineage keys, the directions the roots leave the base in."""
import numpy as np
import pytest

from openalea.archisimple.architecture import TYPE_BASE
from openalea.archisimple.drawing import positions
from openalea.archisimple.parameter_sets import WHEAT
from openalea.archisimple.run import run, scene


@pytest.fixture(scope="module")
def seedling():
    parameters = dict(WHEAT, n_seminal_roots=3, n_adventitious_roots=4)
    return scene(parameters).populations[0].data_structure, parameters


def test_grouping(seedling):
    """One Axis (and Organ) per first-order root: the primary, each seminal, each nodal."""
    ds, p = seedling
    g = ds.mtg
    s = g.scales
    axes = [v for v in g.vertices(scale=s.Axis) if not g.property("isanchor").get(v)]
    organs = [v for v in g.vertices(scale=s.Organ) if not g.property("isanchor").get(v)]
    n_roots = 1 + (p["n_seminal_roots"] - 1) + p["n_adventitious_roots"]
    assert len(axes) == len(organs) == n_roots
    for axis in axes:
        assert len([o for o in organs if g.complex_at_scale(o, s.Axis) == axis]) == 1


def test_branching_on_the_base(seedling):
    """Each seminal or nodal root's apex is a RootApex branched (+) on the base segment, in an Organ of its own; the
    primary's Organ is its `<` chain: the base, the apex."""
    ds, p = seedling
    g = ds.mtg
    s, labels = g.scales, g.labels
    ids = np.asarray(ds.entity_ids("node")).tolist()
    types = dict(zip(ids, np.asarray(ds.get("type")).tolist()))
    base, = [v for v in ids if types[v] == TYPE_BASE]
    primary = g.complex_at_scale(base, s.Organ)
    branched = [c for c in g.children(base) if g.edge_type(c) == "+"]
    assert len(branched) == p["n_seminal_roots"] - 1 + p["n_adventitious_roots"]
    for apex in branched:
        assert g.label(apex) == labels.SubOrgan.RootApex
        assert g.complex_at_scale(apex, s.Organ) != primary
    chain = sorted(v for v in ids if g.complex_at_scale(v, s.Organ) == primary)
    assert [g.label(v) for v in chain] == [labels.SubOrgan.RootSegment, labels.SubOrgan.RootApex]


def test_keys_are_distinct(seedling):
    """Each element its own lineage key, and so its own draws."""
    ds, _ = seedling
    keys = np.asarray(ds.get("key"))
    assert (keys >= 0).all() and np.unique(keys).size == keys.size


def test_the_first_order_roots_leave_the_base_all_around():
    """The seminal and nodal roots leave the base in their drawn azimuth, not all in one plane: after 3 days, their
    first elements point to both sides in x and in y."""
    record = run(dict(WHEAT, emergence_delay=0., ER=1e-4), n_steps=72, record_every=72).steps[-1]
    xyz = positions(record)
    base = [v for v in record.vertices if record.parent[v] is None][0]
    first = [v for v in record.vertices if record.parent[v] == base and record.edge_type[v] == "+"
             and record.props["length"][v] > 0]
    directions = np.array([xyz[v][1] - xyz[v][0] for v in first])
    assert len(first) >= 4
    assert (directions[:, 0] > 0).any() and (directions[:, 0] < 0).any()
    assert (directions[:, 1] > 0).any() and (directions[:, 1] < 0).any()

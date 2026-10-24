"""A population in one MPG: plants of different scenarios grown together, each identical to the same plant grown
alone, bit for bit: a plant's development and draws (keyed by lineage) do not depend on its neighbours. The plants of
one population differ by their numeric parameters (stored per plant); a switch (a boolean option) is the
population's."""
import pytest

from openalea.archisimple.parameter_sets import WHEAT
from openalea.archisimple.run import Record, run

PLANTS = [dict(WHEAT),
          dict(WHEAT, D_ini=5e-4, random_choice=3),
          dict(WHEAT, n_adventitious_roots=10, EL=4e-4, IPD=6e-3),
          dict(WHEAT, D_ini=9e-4, n_seminal_roots=3, CVDD=0.3, random_choice=2)]


def by_key(record):
    """The record keyed by the elements' lineage keys (their ids differ between a population and a plant alone)."""
    key = {v: int(record.props["key"][v]) for v in record.vertices}
    props = {name: {key[v]: value for v, value in values.items()} for name, values in record.props.items()}
    for name in ("axis_apex_id", "lateral_primordium_index"):
        props[name] = {key[v]: key.get(int(a), -1) for v, a in record.props[name].items()}
    return Record(vertices=sorted(key.values()), parent={key[v]: key.get(p) for v, p in record.parent.items()},
                  edge_type={key[v]: e for v, e in record.edge_type.items()}, props=props, plant=record.plant)


@pytest.fixture(scope="module")
def population():
    return run(PLANTS, n_steps=120, record_every=120)


@pytest.mark.parametrize("plant", range(len(PLANTS)))
def test_each_plant_is_the_plant_alone(population, plant):
    alone = by_key(run(PLANTS[plant], n_steps=120, record_every=120).steps[-1])
    together = by_key(population.steps[-1][plant])
    assert together.vertices == alone.vertices and len(together.vertices) > 50
    assert together.parent == alone.parent and together.edge_type == alone.edge_type
    assert together.props == alone.props
    assert together.plant == alone.plant

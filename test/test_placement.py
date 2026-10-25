"""The roots' placement as state (RootArchitecture._placement, on metafspm's follow_parents): its positions are
drawing.positions' of the same records, which equal the former drawing's PlantGL turtle (kept here as the reference)
to 1e-10 m, with and without gravitropism (without, the former drawing used a numpy turtle mirrored in x from PlantGL's;
the placement is PlantGL's in both cases); plants placed at their position, depth and rotation."""
import numpy as np
import pytest

from openalea.archisimple.drawing import positions
from openalea.archisimple.parameter_sets import PEA, WHEAT
from openalea.archisimple.run import plant_model, record, run
from openalea.metafspm.scene.population import Planting, planting
from openalea.metafspm.scene.scene import Scene

NAMES = ("x1", "y1", "z1", "x2", "y2", "z2")


# --- The former drawing (rhizodep's plot_mtg by PlantGL's turtle), the reference ---

def pglturtle_positions(snapshot, origin=(0., 0., 0.), gravitropism=0.06):
    """{vid: (start, end)}: the elements' end points, by the turtle above, from *origin* (the plant's position)."""
    from openalea.plantgl.all import PglTurtle, Vector3
    props = snapshot.props
    children = {vid: [] for vid in snapshot.vertices}
    roots = []
    for vid in snapshot.vertices:
        parent = snapshot.parent[vid]
        (children[parent].append(vid) if parent in children else roots.append(vid))
    turtle = PglTurtle()
    turtle.down(180)                                   # rhizodep's plot_mtg: the graph upside down
    turtle.tropism = (0, 0, -1)
    start_state = (Vector3(*origin), turtle.getHeading(), turtle.getUp())
    radius = props.get("original_radius", props["radius"])
    out = {}
    for root in roots:
        base_radius = radius.get(root) or 1.
        stack = [(root, start_state)]
        while stack:
            vid, (position, heading, up) = stack.pop()
            turtle.move(position)
            turtle.setHead(heading, up)
            azimuth = props.get("azimuth", {}).get(vid, 0.) or 0.
            if azimuth:
                turtle.rollL(azimuth)
            turtle.down(props["angle_down"].get(vid, 0.) or 0.)
            turtle.rollL(props["angle_roll"].get(vid, 0.) or 0.)
            turtle.elasticity = gravitropism * ((radius.get(vid) or 0.) / base_radius)
            turtle.F(props["length"].get(vid, 0.) or 0.)
            end = turtle.getPosition()
            out[vid] = (np.array([position.x, position.y, position.z]), np.array([end.x, end.y, end.z]))
            state = (Vector3(end), turtle.getHeading(), turtle.getUp())
            for child in children[vid]:
                stack.append((child, state))
    return out


@pytest.fixture(scope="module")
def records():
    return {name: run(dict(parameters, random_choice=2), n_steps=240, record_every=240).steps[-1]
            for name, parameters in (("wheat", WHEAT), ("pea", PEA))}


@pytest.mark.parametrize("species", ["wheat", "pea"])
@pytest.mark.parametrize("gravitropism", [0., 0.06])
def test_as_the_former_drawing(records, species, gravitropism):
    pytest.importorskip("openalea.plantgl.all")
    record = records[species]
    new, old = positions(record, gravitropism=gravitropism), pglturtle_positions(record, gravitropism=gravitropism)
    assert len(new) > 100
    for vid, (start, end) in old.items():
        np.testing.assert_allclose(new[vid][0], start, atol=1e-10)
        np.testing.assert_allclose(new[vid][1], end, atol=1e-10)


def test_the_state_is_the_drawing():
    """Three plants, each at its position, depth and rotation: every element's start and end (the step's) are
    positions' of its plant's record."""
    plants = [Planting(0., 0., -0.02, 0.), Planting(0.1, 0.05, -0.03, 90.), Planting(0.2, -0.05, 0., 215.)]
    scenarios = [dict(WHEAT, random_choice=k + 1, gravitropism=0.06 * k) for k in range(3)]
    scene = Scene(planting(plants, plant_model(), plant_parameters=scenarios), time_step=3600)
    for _ in range(120):
        scene.run()
    for k, plant in enumerate(plants):
        state = record(scene, plant=k)
        drawn = positions(state, origin=(plant.x, plant.y, plant.z), gravitropism=0.06 * k, rotation=plant.rotation)
        for vid in state.vertices:
            start, end = (np.array([state.props[name][vid] for name in names]) for names in (NAMES[:3], NAMES[3:]))
            np.testing.assert_allclose(start, drawn[vid][0], atol=1e-12)
            np.testing.assert_allclose(end, drawn[vid][1], atol=1e-12)
        base = [v for v in state.vertices if state.parent[v] is None][0]
        assert [state.props[name][base] for name in NAMES[:3]] == pytest.approx([plant.x, plant.y, plant.z])


def test_triangles_and_no_placement():
    scene = Scene(planting([Planting()], plant_model(), plant_parameters=[dict(WHEAT)]), time_step=3600)
    for _ in range(48):
        scene.run()
    growth = scene.populations[0].instance.growth
    triangles = growth.triangles()
    assert len(triangles) > 0 and triangles.vertices[..., 2].max() <= 1e-12          # below the soil surface
    skipped = Scene(planting([Planting()], plant_model(), plant_parameters=[dict(WHEAT, placement=0)]),
                    time_step=3600)
    for _ in range(48):
        skipped.run()
    assert not np.any(skipped.populations[0].data_structure.get("z2"))

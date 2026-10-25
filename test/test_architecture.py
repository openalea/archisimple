"""RootArchitecture, ArchiSimple's root architecture, by its properties: the satisfaction coefficient (unlimited growth
is the potential one; half the demand offered, half the increments), the root system's development."""
import numpy as np
import pytest

from openalea.archisimple.architecture import RootArchitecture, TYPE_EMERGED
from openalea.archisimple.parameter_sets import WHEAT
from openalea.archisimple.run import run, run_stage, scene


def stepped(built, n_steps, struct_mass_input=None):
    """*n_steps* steps of a scene, stage by stage; the arrays after the potential and after the actual growth of the
    last step."""
    ds = built.populations[0].data_structure
    for step in range(n_steps):
        if struct_mass_input is not None:
            ds.set("struct_mass_input", np.full(len(ds.entity_ids("Plant")), struct_mass_input))
        for stage in ("potential", "actual", "segmentation", "postsegmentation"):
            run_stage(built, stage)
            if step == n_steps - 1 and stage == "potential":
                potential = {name: np.array(ds.get(name)) for name in (
                    "length", "radius", "initial_length", "initial_radius", "potential_length", "potential_radius",
                    "struct_mass", "volume", "root_tissue_density", "type")}
            if step == n_steps - 1 and stage == "actual":
                actual = {name: np.array(ds.get(name)) for name in ("length", "radius", "struct_mass_produced")}
                return potential, actual


def test_unlimited_growth_is_the_potential_one():
    potential, actual = stepped(scene(WHEAT), 72)
    grows = potential["potential_length"] > potential["initial_length"]
    assert grows.any()
    np.testing.assert_allclose(actual["length"][grows], potential["potential_length"][grows], rtol=1e-14)


def test_half_satisfied_growth_takes_half_the_increments():
    built = scene(WHEAT)
    stepped(built, 71)
    ds = built.populations[0].data_structure
    run_stage(built, "potential")
    growth = built.populations[0].instance.growth
    arrays = {name: np.array(ds.get(name)) for name in ("type", "initial_radius", "initial_length", "potential_radius",
                                                         "potential_length", "root_tissue_density")}
    demand = np.sum((growth._volume(arrays["type"], arrays["potential_radius"], arrays["potential_length"])
                     - growth._volume(arrays["type"], arrays["initial_radius"], arrays["initial_length"]))
                    * arrays["root_tissue_density"])
    ds.set("struct_mass_input", np.full(1, demand / 2.))
    assert growth.satisfaction_coefficient(
        arrays["type"], arrays["initial_radius"], arrays["initial_length"], arrays["potential_radius"],
        arrays["potential_length"], arrays["root_tissue_density"])[0] == pytest.approx(0.5, rel=1e-12)
    run_stage(built, "actual")
    increments = np.array(ds.get("length")) - arrays["initial_length"]
    expected = (arrays["potential_length"] - arrays["initial_length"]) * 0.5
    grows = expected > 0
    np.testing.assert_allclose(increments[grows], expected[grows], rtol=1e-9)
    assert np.sum(np.array(ds.get("struct_mass_produced"))) == pytest.approx(demand / 2., rel=1e-9)


def test_the_root_system_develops():
    result = run(WHEAT, n_steps=240, record_every=240)
    first, last = result.steps[0], result.steps[-1]
    assert len(last.vertices) > 3 * len(first.vertices)                    # segmentation
    laterals = [v for v in last.vertices if last.props["root_order"][v] > 1]
    assert laterals                                                         # primordia
    assert any(last.props["type"][v] == TYPE_EMERGED for v in laterals)     # emergences


def test_a_scenario_of_another_component_is_refused():
    """A parameter RootArchitecture does not have (a RootGrowth one) is an error, not ignored."""
    with pytest.raises(ValueError, match="growth_duration_rule"):
        scene(dict(WHEAT, growth_duration_rule="probabilistic"))

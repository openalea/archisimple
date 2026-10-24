"""The reviewed fixes of RhizoDep's root growth that are the architecture's, each
by its property."""
import numpy as np
import pytest

from openalea.archisimple.architecture import TYPE_DEAD, TYPE_JUST_DEAD, TYPE_SEMINAL_BEFORE_EMERGENCE
from openalea.archisimple.parameter_sets import WHEAT
from openalea.archisimple.run import run, scene

# Elements dying within the runs: their life durations shortened, the main roots' durations not extended
SHORT_LIVED = dict(WHEAT, LDs=WHEAT["LDs"] / 200, main_roots_growth_extender=1)


def test_each_seminal_root_its_own_draws():
    seedling = run(WHEAT, n_steps=0).steps[0]
    radii = {seedling.props["radius"][v] for v in seedling.vertices
             if seedling.props["type"][v] == TYPE_SEMINAL_BEFORE_EMERGENCE}
    assert len(radii) == WHEAT["n_seminal_roots"] - 1


def test_actual_time_at_zero_temperature():
    growth = scene(WHEAT).populations[0].instance.growth
    assert growth._actual_time(np.array([3600.]), np.array([0.]))[0] == 0.
    assert growth._actual_time(np.array([3600.]), np.array([0.5]))[0] == 7200.


def test_determinate_laterals():
    """The laterals' durations are GDs·d², without the main roots' extender: determinate laterals (Pagès et al. 2014,
    Eq. 2)."""
    last = run(WHEAT, n_steps=120, record_every=120).steps[-1]
    laterals = [v for v in last.vertices if last.props["root_order"][v] > 1 and last.props["label"][v] == 2]
    assert laterals
    for v in laterals[:20]:
        diameter = 2. * last.props["original_radius"][v]
        assert last.props["growth_duration"][v] == pytest.approx(WHEAT["GDs"] * diameter ** 2, rel=1e-9)


def test_laterals_by_their_diameter():
    """Laterals form when their diameter (not radius) reaches Dmin."""
    parameters = dict(WHEAT, D_ini=2.4e-4, Dmin=1.6e-4, CVDD=0.3)
    last = run(parameters, n_steps=168, record_every=168).steps[-1]
    radii = [last.props["original_radius"][v] for v in last.vertices if last.props["root_order"][v] > 1
             and last.props["label"][v] == 2]
    assert radii
    assert all(2 * r >= parameters["Dmin"] for r in radii)


def test_life_duration_exponent():
    """ArchiSimple's life duration, LDs·d²·density (the parameter sets' exponent, 2)."""
    ds = scene(WHEAT).populations[0].data_structure
    radius, life = np.asarray(ds.get("radius")), np.asarray(ds.get("life_duration"))
    assert life[0] == pytest.approx(WHEAT["LDs"] * (2 * radius[0]) ** 2 * WHEAT["new_root_tissue_density"], rel=1e-12)


def test_the_dead_have_no_living_mass():
    last = run(SHORT_LIVED, n_steps=240, record_every=240).steps[-1]
    dead = [v for v in last.vertices if last.props["type"][v] in (TYPE_DEAD, TYPE_JUST_DEAD)
            and last.props["struct_mass"][v] > 0]
    assert dead
    assert all(last.props["living_struct_mass"][v] == 0 for v in dead)

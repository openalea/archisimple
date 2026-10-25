"""
The component run alone, and its records: to test it, draw it and compare it.

    run = run(parameters, n_steps=240)             RootArchitecture alone, one plant (a list of scenarios: a population)
    run.scene                                      the metafspm Scene
    run.steps[i]                                   Records after run.recorded_steps[i] steps (0: the seedling)
    record(scene, plant=0)                         one plant's state: its elements by MPG id, their parents and edges,
                                                   every declared node variable, the MPG label (2: an apex, 1: a segment)
    run_stage(scene, "segmentation")               only that stage's steps of the component

A subclass (rhizodep's RootGrowth) is run the same way, with its inputs given (*given*: {name: a number, an array or a
function of the step and the node ids}).
"""
from dataclasses import dataclass, field

import numpy as np

from openalea.metafspm.coupling.composite_wrapper import CompositeModel
from openalea.metafspm.coupling.translator import Translator
from openalea.metafspm.scene.population import Planting, planting
from openalea.metafspm.scene.scene import Scene
from openalea.archisimple.architecture import RootArchitecture

LABEL_SEGMENT, LABEL_APEX = 1, 2


@dataclass
class Record:
    vertices: list
    parent: dict
    edge_type: dict
    props: dict
    plant: dict
    meta: dict = field(default_factory=dict)


@dataclass
class Run:
    scene: Scene
    steps: list = field(default_factory=list)
    recorded_steps: list = field(default_factory=list)


_MODELS = {}


def plant_model(component=RootArchitecture):
    """A plant model (a CompositeModel) of *component* (RootArchitecture or a subclass) alone, nothing to couple; its
    other inputs given by run()."""
    if component not in _MODELS:
        def __init__(self, data_structure, time_step):
            self.growth = component(data_structure=data_structure, time_step=time_step)
            self.declare_data_and_couple_components(translator_path=Translator(), components=(self.growth,))

        def run(self):
            self.growth()

        _MODELS[component] = type(f"{component.__name__}Only", (CompositeModel,), dict(
            initiators=(component,), __init__=__init__, run=run, __doc__=plant_model.__doc__))
    return _MODELS[component]


def scene(parameters=None, time_step=3600, component=RootArchitecture, spacing=0.1):
    """A Scene of *component*: one plant, or one per element of *parameters* when it is a list (a population in one
    MPG, each plant its own scenario), *spacing* m apart along x."""
    plants = parameters if isinstance(parameters, list) else [parameters]
    table = planting([Planting(x=spacing * i) for i in range(len(plants))], plant_model(component),
                     plant_parameters=[dict(p or {}) for p in plants], xrange=max(1., spacing * len(plants)), yrange=1.)
    return Scene(table, time_step=time_step)


def record(scene, plant=0, decode=None) -> Record:
    """The *plant*-th plant's state (in the planting table's order). *decode*: {variable: {code: value}}, values
    written instead of the codes (e.g. a legacy's strings)."""
    population = scene.populations[0]
    ds, growth = population.data_structure, population.instance.components[0]
    g = ds.mtg
    owner = np.asarray(ds.owner("Plant"))
    ids = np.asarray(ds.entity_ids("node"))[owner == plant].tolist()
    decode = decode or {}
    props = {}
    for name, spec in growth._variable_specs.items():
        if spec.variable_type in ("state_variable", "input") and ds.location(name) == "node":
            values = np.asarray(ds.get(name))[owner == plant]
            codes = decode.get(name)
            props[name] = ({v: codes[int(x)] for v, x in zip(ids, values)} if codes is not None
                           else {v: x.item() for v, x in zip(ids, values)})
    apex = g.labels.SubOrgan.RootApex
    props["label"] = {v: LABEL_APEX if g.label(v) == apex else LABEL_SEGMENT for v in ids}
    for name in type(growth).ELEMENT_PARAMETERS:
        values = np.asarray(growth.parameter_values(name))[owner == plant]
        props[name] = {v: x.item() for v, x in zip(ids, values)}
    known = set(ids)
    parent, edge_type = {}, {}
    for v in ids:
        mpg_parent = g.parent(v)
        parent[v] = mpg_parent if mpg_parent in known else None
        edge_type[v] = g.edge_type(v) if parent[v] is not None else ""
    totals = {"total_living_struct_mass": float(np.asarray(ds.get("total_living_struct_mass"))[plant])}
    return Record(vertices=sorted(ids), parent=parent, edge_type=edge_type, props=props, plant=totals)


def records(scene, parameters, decode=None):
    """The scene's state: a Record, or one per plant when *parameters* gives a population (a list)."""
    if isinstance(parameters, list):
        return [record(scene, plant, decode) for plant in range(len(parameters))]
    return record(scene, decode=decode)


def write_given(scene, step, name, value):
    """A given input on every element: a number, an array, or a function of the step and the node ids."""
    ds = scene.populations[0].data_structure
    if value is None or not ds.has(name):              # not given, or not an input of this component
        return
    if callable(value):
        value = value(step, np.asarray(ds.entity_ids("node")))
    ds.set(name, np.broadcast_to(np.asarray(value, dtype=float), (ds.n_nodes(),)))


def run(parameters=None, n_steps=24, time_step=3600, record_every=1, record_states=True, component=RootArchitecture,
        given=None, setup=None, decode=None) -> Run:
    """*n_steps* steps of *component*; the records every *record_every* steps (and the last). *given*: inputs written
    before each step ({name: value}, write_given); *setup*: a function of the scene, called once it is built."""
    built = scene(parameters, time_step, component)
    if setup is not None:
        setup(built)
    out = Run(scene=built)
    if record_states:
        out.steps.append(records(built, parameters, decode))
        out.recorded_steps.append(0)
    for step in range(n_steps):
        for name, value in (given or {}).items():
            write_given(built, step, name, value)
        built.run()
        if record_states and ((step + 1) % record_every == 0 or step == n_steps - 1):
            out.steps.append(records(built, parameters, decode))
            out.recorded_steps.append(step + 1)
    return out


def run_stage(scene, stage):
    """Run only the component's steps of one stage (metafspm's run_stage: "potential", "actual", "segmentation",
    "postsegmentation")."""
    scene.populations[0].instance.growth.run_stage(stage)

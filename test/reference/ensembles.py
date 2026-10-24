"""
Ensembles over seeds: the statistical reference of the architecture, against which its changes are validated.

    python test/reference/ensembles.py        writes test/reference/ensembles.json, from the code checked out

The reference was computed on rhizodep's tag port-exact, where RootArchitecture reproduced RhizoDep's original root
growth bit for bit. Each configuration: a parameter set, 10 days, 8 seeds; the totals after the run, their mean and
standard deviation. Its element counts include the zero-length supports RhizoDep's seedling then had (SUPPORTS).
"""
import json
import os

import numpy as np

from openalea.archisimple.architecture import RootArchitecture
from openalea.archisimple.drawing import totals
from openalea.archisimple.parameter_sets import MAIZE, PEA
from openalea.archisimple.run import run

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = range(1, 9)
DAYS = 10
CONFIGURATIONS = {"pea": PEA, "maize": MAIZE}
METRICS = ("length", "struct_mass", "vertices", "apices", "lateral_apices", "mean_radius", "first_order_radius_spread")


def supports(name):
    """The zero-length supports the reference's seedlings had (one per seminal root but the primary, one per nodal
    root), which the elements counted since no longer include."""
    parameters = CONFIGURATIONS[name]
    return (parameters.get("n_seminal_roots", RootArchitecture.n_seminal_roots) - 1
            + parameters.get("n_adventitious_roots", RootArchitecture.n_adventitious_roots))


def metrics(record):
    t = totals(record)
    radii = [record.props["radius"][v] for v in record.vertices if record.props["length"][v] > 0]
    first = [record.props["original_radius"][v] for v in record.vertices
             if record.props["root_order"][v] == 1 and record.props["label"][v] == 2]
    return dict(length=t["length"], struct_mass=t["struct_mass"], vertices=t["vertices"], apices=t["apices"],
                lateral_apices=t["lateral_apices"], mean_radius=float(np.mean(radii)),
                first_order_radius_spread=float(np.std(first)))


def ensemble(name, seeds=SEEDS, days=DAYS):
    """{metric: [value per seed]} of a configuration."""
    values = {metric: [] for metric in METRICS}
    for seed in seeds:
        result = run(dict(CONFIGURATIONS[name], random_choice=seed), n_steps=24 * days, record_every=24 * days)
        for metric, value in metrics(result.steps[-1]).items():
            values[metric].append(value)
    return values


if __name__ == "__main__":
    reference = {}
    for name in CONFIGURATIONS:
        values = ensemble(name)
        reference[name] = {metric: dict(mean=float(np.mean(v)), std=float(np.std(v, ddof=1)), values=v)
                           for metric, v in values.items()}
        print(name, {m: round(reference[name][m]["mean"], 4) for m in METRICS})
    with open(os.path.join(HERE, "ensembles.json"), "w") as file:
        json.dump(dict(source="archisimple", days=DAYS, seeds=list(SEEDS), configurations=reference), file, indent=1)

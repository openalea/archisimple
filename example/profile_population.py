"""
A population in one MPG against its plants grown one by one: the population's plants are in the same arrays, each
step one call for all; grown one by one, each plant is its own Scene. Each plant is the same either way
(test/test_population.py).

    left     the time of the run against the number of plants, one by one and as one population
    right    the time per plant

The plants differ by their scenario (random_choice: their own draws, so their own root systems). Run from
archisimple/example:

    python profile_population.py [--plants 1 2 4 8 16 32] [--days D] [--species NAME]    figure in ./outputs/
"""
import argparse
import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from openalea.archisimple.parameter_sets import SPECIES
from openalea.archisimple.run import run

OUTPUTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


def timed(function, *args, **keywords):
    start = time.perf_counter()
    function(*args, **keywords)
    return time.perf_counter() - start


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plants", type=int, nargs="*", default=[1, 2, 4, 8, 16, 32])
    parser.add_argument("--days", type=float, default=7.)
    parser.add_argument("--species", default="wheat (Soissons)", choices=list(SPECIES))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    if not args.show:
        matplotlib.use("Agg")
    os.makedirs(OUTPUTS, exist_ok=True)
    n_steps = int(round(args.days * 24))
    one_by_one, together = [], []
    for n in args.plants:
        plants = [dict(SPECIES[args.species], random_choice=1 + k) for k in range(n)]
        one_by_one.append(sum(timed(run, plant, n_steps=n_steps, record_states=False) for plant in plants))
        together.append(timed(run, plants, n_steps=n_steps, record_states=False))
        print(f"{n:3d} plants: one by one {one_by_one[-1]:6.1f} s, one population {together[-1]:6.1f} s "
              f"({one_by_one[-1] / together[-1]:.1f}× faster)")
    n = np.array(args.plants)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for label, seconds, style in (("one by one", one_by_one, "ko-"), ("one population", together, "bs-")):
        axes[0].plot(n, seconds, style, label=label)
        axes[1].plot(n, np.array(seconds) / n, style, label=label)
    for ax, ylabel in zip(axes, ("time of the run (s)", "time per plant (s)")):
        ax.set_xscale("log", base=2)
        ax.set_xlabel("plants")
        ax.set_ylabel(ylabel)
        ax.legend()
    axes[0].set_yscale("log")
    fig.suptitle(f"{args.species}, {args.days:g} days: {n[-1]} plants, one population "
                 f"{one_by_one[-1] / together[-1]:.1f}× faster")
    fig.tight_layout()
    fig.savefig(os.path.join(OUTPUTS, f"population_{args.days:g}days.png"), dpi=120)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

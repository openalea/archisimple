"""
ArchiSimple's calibrated species side by side (openalea.archisimple.parameter_sets: wheat, pea, maize, peach), grown
alone: above, their elements coloured by root order; below, by radius (widths and colours on one logarithmic scale for
all), where the dicotyledons' secondary growth shows (SGC 0.9: pea, peach), the monocotyledons having none (SGC 0).

    python species.py [--days D] [--maize-days D] [--show]          figure in ./outputs/species.png

Maize's nodal roots make it the largest (about 37 000 elements at 20 days): it grows 20 days by default, the others
30.
"""
import argparse
import os
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from openalea.archisimple.drawing import draw_root_system, totals
from openalea.archisimple.parameter_sets import SPECIES
from openalea.archisimple.run import run

OUTPUTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--days", type=float, default=30.)
    parser.add_argument("--maize-days", type=float, default=20.)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    if not args.show:
        matplotlib.use("Agg")
    os.makedirs(OUTPUTS, exist_ok=True)
    records, days = {}, {}
    for species, parameters in SPECIES.items():
        days[species] = args.maize_days if species == "maize" else args.days
        n_steps = int(round(days[species] * 24))
        start = time.perf_counter()
        records[species] = run(parameters, n_steps=n_steps, record_every=n_steps).steps[-1]
        print(f"{species}: {days[species]:g} days, {len(records[species].vertices)} elements, "
              f"{time.perf_counter() - start:.1f} s")
    radii_of = {species: np.array([r.props["radius"][v] for v in r.vertices]) for species, r in records.items()}
    smallest = min(radii[radii > 0].min() for radii in radii_of.values())
    largest = max(radii.max() for radii in radii_of.values())
    fig, axes = plt.subplots(2, len(SPECIES), figsize=(17, 11))
    for column, (species, record) in enumerate(records.items()):
        orders = draw_root_system(axes[0, column], record, color_by="root_order", linewidth_by_radius=False,
                                  vmin=1, vmax=3)
        radii = draw_root_system(axes[1, column], record, color_by="radius")
        radii.set_norm(LogNorm(smallest, largest))
        radii.set_linewidths(0.3 + 4. * np.sqrt(radii.get_array() / largest))
        t = totals(record)
        axes[0, column].set_title(f"{species}, {days[species]:g} days\n{t['length']:.1f} m, "
                                  f"{t['struct_mass'] * 1e3:.0f} mg, {t['apices']} apices", fontsize=9)
        positive = radii_of[species][radii_of[species] > 0]
        axes[1, column].set_title(f"radii from {positive.min() * 1e3:.2f} to {positive.max() * 1e3:.2f} mm",
                                  fontsize=9)
        axes[1, column].set_xlabel("x (m)")
    for row in axes:
        row[0].set_ylabel("z (m)")
    fig.colorbar(orders, ax=axes[0].tolist(), label="root order", shrink=0.8)
    fig.colorbar(radii, ax=axes[1].tolist(), label="radius (m)", shrink=0.8)
    fig.savefig(os.path.join(OUTPUTS, "species.png"), dpi=110)
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()

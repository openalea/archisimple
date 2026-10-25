"""
Root systems drawn from a record (openalea.archisimple.run.record, or any record of that form: rhizodep's records of its
legacy), to look at them and compare them.

    xyz = positions(snapshot)                {vid: (start, end)}, in m
    draw_root_system(ax, snapshot, color_by="type")
    totals(snapshot)                         {"length": ..., "struct_mass": ..., "apices": ..., "lateral_apices": ...}

The geometry is rhizodep's (tools.plot_mtg and get_root_visitor), PlantGL's turtle reproduced without PlantGL
(metafspm's follow_parents, as RootArchitecture's _placement step): from the plant's position, pointing down, each
element turns down by angle_down then rolls by angle_roll (a root leaving the base first rolls by its azimuth), bends
towards the vertical by its gravitropism (gravitropism · its original radius over the base's, Wheat-BRIDGES' 0.06 by
default), and goes `length` forward.
"""
import numpy as np

from openalea.metafspm.coupling.geometry import follow_parents, roll_right

LABEL_APEX = 2
TYPE_NAMES = {1: "base", 2: "seminal support (RhizoDep's legacy)", 3: "seminal before emergence", 4: "adventitious support (legacy)",
              5: "adventitious before emergence", 6: "before emergence", 7: "emerged", 8: "stopped",
              9: "just stopped", 10: "dead", 11: "just dead", 12: "nodule"}


def positions(snapshot, origin=(0., 0., 0.), gravitropism=0.06, rotation=0.):
    """{vid: (start, end)}: the elements' end points, as RootArchitecture's _placement places them (metafspm's
    follow_parents: PlantGL's turtle, without PlantGL), from *origin* (the plant's position) rolled by *rotation* (°).
    For records of any form (a run's, or rhizodep's legacy's)."""
    props = snapshot.props
    vids = list(snapshot.vertices)
    index = {vid: k for k, vid in enumerate(vids)}
    parent = np.array([index.get(snapshot.parent[vid], -1) for vid in vids], dtype=np.int64)
    children = {k: [] for k in range(len(vids))}
    for k, p in enumerate(parent.tolist()):
        if p >= 0:
            children[p].append(k)
    order, stack = [], [k for k in range(len(vids)) if parent[k] < 0][::-1]
    while stack:
        k = stack.pop()
        order.append(k)
        stack.extend(children[k][::-1])

    def values(name, default=None):
        given = props.get(name, {}) if default is None else props.get(name, default)
        return np.array([float(given.get(vid, 0.) or 0.) for vid in vids])

    radius = values("original_radius") if "original_radius" in props else values("radius")
    root_of = np.arange(len(vids))
    for k in order:
        if parent[k] >= 0:
            root_of[k] = root_of[parent[k]]
    reference = np.where(radius[root_of] > 0., radius[root_of], 1.)
    heading, up = np.array([0., 0., -1.]), np.array([-1., 0., 0.])
    roll_right(heading, up, float(rotation))
    n = len(vids)
    start, _, _, end = follow_parents(
        np.asarray(order, dtype=np.int64), parent, np.tile(np.asarray(origin, dtype=float), (n, 1)),
        np.tile(heading, (n, 1)), np.tile(up, (n, 1)), values("azimuth"), values("angle_down"), values("angle_roll"),
        gravitropism * radius / reference, values("length"), np.array([0., 0., -1.]))
    return {vid: (start[k], end[k]) for k, vid in enumerate(vids)}


def draw_root_system(ax, snapshot, color_by="type", cmap="viridis", plane=(0, 2), linewidth_by_radius=True,
                     xyz=None, vmin=None, vmax=None, origin=(0., 0., 0.)):
    """The elements as lines projected on *plane* (axes indices, x-z by default), coloured by *color_by* (a
    per-vertex property); returns the LineCollection (for a colorbar)."""
    from matplotlib.collections import LineCollection

    xyz = positions(snapshot, origin) if xyz is None else xyz
    vids = [vid for vid in snapshot.vertices if (snapshot.props["length"].get(vid, 0.) or 0.) > 0.]
    i, j = plane
    lines = [((xyz[v][0][i], xyz[v][0][j]), (xyz[v][1][i], xyz[v][1][j])) for v in vids]
    values = np.array([float(snapshot.props[color_by].get(v, np.nan)) for v in vids])
    radii = np.array([snapshot.props["radius"].get(v, 0.) for v in vids])
    widths = 0.4 + 3. * radii / radii.max() if linewidth_by_radius and len(radii) and radii.max() > 0 else 0.8
    collection = LineCollection(lines, array=values, cmap=cmap, linewidths=widths)
    if vmin is not None or vmax is not None:
        collection.set_clim(vmin, vmax)
    ax.add_collection(collection)
    ax.set_aspect("equal")
    ax.autoscale_view()
    return collection


def totals(snapshot):
    p = snapshot.props
    apices = [v for v in snapshot.vertices if p["label"].get(v) == LABEL_APEX]
    return {
        "length": sum(p["length"].get(v, 0.) for v in snapshot.vertices),
        "struct_mass": sum(p["struct_mass"].get(v, 0.) for v in snapshot.vertices),
        "living_struct_mass": snapshot.plant.get("total_living_struct_mass"),
        "vertices": len(snapshot.vertices),
        "apices": len(apices),
        "lateral_apices": sum(1 for v in apices if p["root_order"].get(v, 1) > 1),
    }

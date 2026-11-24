"""Closed convex solids as triangle meshes, outward CCW, with exact analytic references.

Every solid is a union of disjoint axis-aligned boxes, so its reference volume,
centroid, inertia tensor and second-moment tensor are known in closed form.
"""
import numpy as np

# ---------- box mesh generator (outward-facing, CCW when seen from outside) ----------

def box_mesh(lo, hi):
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    v = np.array([
        [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],  # z0 face
        [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],  # z1 face
    ])
    f = np.array([
        [0, 2, 1], [0, 3, 2],   # -z
        [4, 5, 6], [4, 6, 7],   # +z
        [0, 1, 5], [0, 5, 4],   # -y
        [2, 3, 7], [2, 7, 6],   # +y
        [0, 4, 7], [0, 7, 3],   # -x
        [1, 2, 6], [1, 6, 5],   # +x
    ])
    return v, f


def merge(*parts):
    vs, fs, off = [], [], 0
    for v, f in parts:
        vs.append(v)
        fs.append(f + off)
        off += len(v)
    return np.vstack(vs), np.vstack(fs)


# Three models, three different mass distributions.  Boxes in each model have
# disjoint interiors, so the union is a valid closed solid.
#
#   A  heavy slab + small corner pedestal   -> centroid well off origin,
#                                             3D-tilted principal axes
#   B  L-plate with a tall chimney           -> Ixx << Iyy < Izz, axes tilted
#                                             in all three planes
#   C  flat L-bracket in the xy plane        -> plan-view principal axes tilted
#                                             ~45 deg away from world axes
MODEL_BOXES = {
    'A': [([0, 0, 0], [4, 3, 1]),
          ([3, 2, 1], [4, 3, 2])],
    'B': [([0, 0, 0], [6, 2, 1]),
          ([0, 2, 0], [2, 6, 1]),
          ([6, 0, 0], [7, 2, 3])],
    'C': [([0, 0, 0], [6, 1, 1]),
          ([0, 1, 0], [1, 5, 1])],
}


def gen_three():
    return {name: merge(*[box_mesh(lo, hi) for lo, hi in parts])
            for name, parts in MODEL_BOXES.items()}


def write_obj(path, v, f):
    with open(path, "w") as fp:
        for p in v:
            fp.write(f"v {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}\n")
        for t in f:
            fp.write(f"f {t[0] + 1} {t[1] + 1} {t[2] + 1}\n")


# ---------- analytic reference ----------

def box_inertia_about_own_center(dims, rho):
    dx, dy, dz = dims
    m = rho * dx * dy * dz
    return np.diag([m * (dy * dy + dz * dz) / 12.0,
                    m * (dx * dx + dz * dz) / 12.0,
                    m * (dx * dx + dy * dy) / 12.0])


def reference(boxes, rho=1.0):
    """Volume, centroid, inertia about the centroid, second moment about centroid."""
    vols = [rho * np.prod(np.subtract(hi, lo)) for lo, hi in boxes]
    cents = [0.5 * (np.asarray(lo, float) + np.asarray(hi, float)) for lo, hi in boxes]
    V = sum(vols)
    c = sum(v * p for v, p in zip(vols, cents)) / V
    I = np.zeros((3, 3))
    S = np.zeros((3, 3))
    for v, p, (lo, hi) in zip(vols, cents, boxes):
        d = p - c
        dims = np.subtract(hi, lo)
        # second moment of this box about its own centre
        S0 = np.diag(v * dims * dims / 12.0)
        I += box_inertia_about_own_center(dims, rho)
        I += v * (np.dot(d, d) * np.eye(3) - np.outer(d, d))
        S += S0 + v * np.outer(d, d)
    return V, c, I, S

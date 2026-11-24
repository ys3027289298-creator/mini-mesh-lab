#!/usr/bin/env python3
"""Generate three test meshes with different vertex distributions.

  m1_grid.obj    uniform 6x6x6 grid + 10 exact duplicate vertices
  m2_cluster.obj 12 tight clusters (sigma << step) -> heavy collisions
  m3_soup.obj    random triangle soup, every triangle emitted twice
                 (exact duplicate corners, like an STL round-trip)

Usage: gen_models.py <outdir>
"""
import random
import sys
import os

def write_obj(path, verts, faces):
    with open(path, "w") as f:
        for v in verts:
            f.write("v %.17g %.17g %.17g\n" % v)
        for a, b, c in faces:
            f.write("f %d %d %d\n" % (a + 1, b + 1, c + 1))

def m1_grid(out):
    verts, faces = [], []
    n = 6
    idx = {}
    for i in range(n):
        for j in range(n):
            for k in range(n):
                idx[(i, j, k)] = len(verts)
                verts.append((0.5 * i, 0.5 * j, 0.5 * k))
    # surface quads -> triangles
    def quad(p, q, r, s):
        faces.append((idx[p], idx[q], idx[r]))
        faces.append((idx[p], idx[r], idx[s]))
    for a in range(n - 1):
        for b in range(n - 1):
            quad((0, a, b), (0, a + 1, b), (0, a + 1, b + 1), (0, a, b + 1))
            quad((n - 1, a, b), (n - 1, a, b + 1), (n - 1, a + 1, b + 1), (n - 1, a + 1, b))
            quad((a, 0, b), (a + 1, 0, b), (a + 1, 0, b + 1), (a, 0, b + 1))
            quad((a, n - 1, b), (a, n - 1, b + 1), (a + 1, n - 1, b + 1), (a + 1, n - 1, b))
            quad((a, b, 0), (a + 1, b, 0), (a + 1, b + 1, 0), (a, b + 1, 0))
            quad((a, b, n - 1), (a, b + 1, n - 1), (a + 1, b + 1, n - 1), (a + 1, b, n - 1))
    # 10 exact duplicates, referenced by extra triangles
    rng = random.Random(1)
    for _ in range(10):
        src = rng.randrange(len(verts))
        dup = len(verts)
        verts.append(verts[src])
        a, b = rng.randrange(n ** 3), rng.randrange(n ** 3)
        if a != src and b != src and a != b:
            faces.append((src, a, b))
            faces.append((dup, a, b))
    write_obj(os.path.join(out, "m1_grid.obj"), verts, faces)

def m2_cluster(out):
    rng = random.Random(2)
    verts, faces = [], []
    centers = [(2.0 * (i % 4), 2.0 * ((i // 4) % 2), 2.0 * (i // 8)) for i in range(12)]
    for c in centers:
        base = len(verts)
        for _ in range(8):
            verts.append(tuple(c[d] + rng.gauss(0.0, 1e-4) for d in range(3)))
        # fan triangles inside the cluster (degenerate after weld)
        for t in range(1, 7):
            faces.append((base, base + t, base + t + 1))
    # a few cross-cluster triangles (survive the weld)
    for i in range(0, 12, 3):
        faces.append((i * 8, ((i + 1) % 12) * 8, ((i + 2) % 12) * 8))
    write_obj(os.path.join(out, "m2_cluster.obj"), verts, faces)

def m3_soup(out):
    rng = random.Random(3)
    verts, faces = [], []
    tris = []
    for _ in range(60):
        t = []
        for _ in range(3):
            t.append(len(verts))
            verts.append((rng.random(), rng.random(), rng.random()))
        tris.append(tuple(t))
    # duplicate every triangle's corners exactly (STL-style soup)
    for t in tris:
        dup = []
        for vi in t:
            dup.append(len(verts))
            verts.append(verts[vi])
        faces.append(t)
        faces.append(tuple(dup))
    write_obj(os.path.join(out, "m3_soup.obj"), verts, faces)

if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    m1_grid(out)
    m2_cluster(out)
    m3_soup(out)
    print("wrote m1_grid.obj m2_cluster.obj m3_soup.obj to", out)

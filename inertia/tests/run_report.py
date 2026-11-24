"""Run the C solver on the three models and check both criteria together.

For each model, verify against the closed-form reference (inertia/tests/boxes.py):
  1. volume and centroid
  2. the three principal moments of inertia about the centroid
  3. the angle between each principal axis and the matching shape principal
     direction (eigenvectors of the second-moment tensor, ordered by size)

The same table is printed with the un-shifted ("origin") tensor, so a model
whose volume matches but whose axis angle is large is visibly the case where
the parallel-axis term was added but the axes were not recomputed.
"""
import json
import subprocess
import sys
import numpy as np
import boxes

MODEL_NAMES = ['A', 'B', 'C']


def shape_axes(S):
    """Second-moment eigenvectors, ordered by descending moment (long axis first)."""
    w, V = np.linalg.eigh(S)
    order = np.argsort(w)[::-1]
    return w[order], V[:, order].T


def axis_angles(shape, axes):
    """Angle from each shape axis to the matching inertia axis, degrees."""
    return np.array([np.degrees(np.arccos(min(1.0, abs(np.dot(s, a))))) for s, a in zip(shape, axes)])


def main():
    solver = sys.argv[1] if len(sys.argv) > 1 else '../inertia'
    models = boxes.gen_three()
    rows = []
    ok = True

    for name in MODEL_NAMES:
        obj = f'/tmp/model_{name}.obj'
        boxes.write_obj(obj, *models[name])
        out = subprocess.run([solver, obj], capture_output=True, text=True)
        if out.returncode != 0:
            print(f"model {name}: solver exited {out.returncode}: {out.stderr.strip()}")
            return 1
        r = json.loads(out.stdout)

        V_ref, c_ref, I_ref, S_ref = boxes.reference(boxes.MODEL_BOXES[name])
        lam_ref, ax_ref = shape_axes(S_ref)
        ev_ref = np.linalg.eigvalsh(I_ref)

        # 1. volume + centroid
        dV = abs(r['volume'] - V_ref)
        dc = np.linalg.norm(np.array(r['centroid']) - c_ref)
        # 2. principal moments
        lam_c = np.array(r['principal_moments'])
        dl = np.max(np.abs(lam_c - ev_ref))
        # 3. axis angles.  Compare the solver's inertia principal axes (about
        #    the centroid) with the reference inertia axes, and separately the
        #    solver's shape axes with the reference shape axes (sign- and
        #    order-insensitive alignment deficit, because arccos near 1 turns
        #    float noise into microdegrees).
        ax_in = shape_axes(np.array(r['inertia_com']))[1]
        _, ax_in_ref = shape_axes(I_ref)
        inertia_def = max(1.0 - abs(np.dot(a, b)) for a, b in zip(ax_in, ax_in_ref))
        ax_s = shape_axes(np.array(r['second_moment']))[1]
        shape_def = max(1.0 - abs(np.dot(a, b)) for a, b in zip(ax_s, ax_ref))

        a = r['axis_angles_deg']
        good = dV < 1e-9 * V_ref and dc < 1e-9 and dl < 1e-9 * max(1.0, ev_ref[-1]) \
            and max(a) < 1e-4

        rows.append(dict(
            model=name, volume=r['volume'], V_ref=V_ref, dV=dV,
            centroid=r['centroid'], c_ref=c_ref.tolist(), dc=dc,
            moments=r['principal_moments'], moments_ref=ev_ref.tolist(), dlam=dl,
            angles=r['axis_angles_deg'], angle_max=r['axis_angle_max_deg'],
            shape_def_vs_ref=shape_def, inertia_def_vs_ref=inertia_def,
            origin_angles=r['origin_axis_angles_deg'],
            origin_angle_max=r['origin_axis_angle_max_deg'],
            centroid_offset=r['centroid_offset'],
            closed=r['closed'], good=bool(good)))

        ok = ok and good
        ok = ok and shape_def < 1e-12 and inertia_def < 1e-12

    hdr = (f"{'model':>5} {'V':>9} {'c':>28} {'I1':>10} {'I2':>10} {'I3':>10}"
           f" {'ax1':>7} {'ax2':>7} {'ax3':>7} {'axmax':>7}")
    print(hdr)
    for r in rows:
        c = r['centroid']
        print(f"{r['model']:>5} {r['volume']:9.4f}"
              f" ({c[0]:7.4f},{c[1]:7.4f},{c[2]:7.4f}) "
              f"{r['moments'][0]:10.5f} {r['moments'][1]:10.5f} {r['moments'][2]:10.5f}"
              f" {r['angles'][0]:7.3f} {r['angles'][1]:7.3f} {r['angles'][2]:7.3f}"
              f" {r['angle_max']:7.4f}")

    print("\nreference check")
    for r in rows:
        print(f"{r['model']}: volume {r['volume']:.17g} (ref {r['V_ref']:.17g}, d={r['dV']:.2e})"
              f"  centroid_err={r['dc']:.2e}  moment_err={r['dlam']:.2e}"
              f"  inertia_axes_deficit={r['inertia_def_vs_ref']:.2e}  shape_axes_deficit={r['shape_def_vs_ref']:.2e}  closed={r['closed']}")

    print("\nwhat the un-shifted (origin) tensor would have reported")
    for r in rows:
        print(f"{r['model']}: centroid offset |c|={r['centroid_offset']:.4f}"
              f"  origin axis angles={np.round(r['origin_angles'], 3).tolist()}"
              f"  max={r['origin_angle_max']:.3f} deg")

    with open('results.json', 'w') as fp:
        json.dump(rows, fp, indent=2)

    print("models_ok=%d" % (1 if ok else 0))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

"""Validate the model generators and the analytic reference before trusting them."""
import numpy as np
import boxes


def signed_volume_moments(v, f):
    """Volume, first moment, and integral of x x^T over a closed triangle mesh."""
    V6 = 0.0
    m1 = np.zeros(3)
    P = np.zeros((3, 3))   # integral of x x^T dV
    for (ia, ib, ic) in f:
        a, b, c = v[ia], v[ib], v[ic]
        vol6 = np.dot(a, np.cross(b, c))
        V6 += vol6
        m1 += vol6 * (a + b + c)
        u = a + b + c
        S = np.outer(a, a) + np.outer(b, b) + np.outer(c, c)
        P += vol6 / 120.0 * (S + np.outer(u, u))
    vol = V6 / 6.0
    return vol, m1 / 24.0, P


def main():
    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"{name}={'PASS' if cond else 'FAIL'}")
        ok = ok and bool(cond)

    for name, (v, f) in boxes.gen_three().items():
        V_mesh, m1, P = signed_volume_moments(v, f)
        c_mesh = m1 / V_mesh
        S_mesh = P - V_mesh * np.outer(c_mesh, c_mesh)
        I_mesh = np.trace(S_mesh) * np.eye(3) - S_mesh

        V_ref, c_ref, I_ref, S_ref = boxes.reference(boxes.MODEL_BOXES[name])
        check(f"{name}_volume_positive", V_mesh > 0)
        check(f"{name}_volume_vs_reference", abs(V_mesh - V_ref) < 1e-9 * V_ref)
        check(f"{name}_centroid_vs_reference", np.allclose(c_mesh, c_ref, atol=1e-9))
        check(f"{name}_inertia_vs_reference", np.allclose(I_mesh, I_ref, atol=1e-9))
        check(f"{name}_second_moment_vs_reference", np.allclose(S_mesh, S_ref, atol=1e-9))

        ev = np.linalg.eigvalsh(I_ref)
        check(f"{name}_inertia_spd", np.all(ev > 0))
        check(f"{name}_triangle_ineq", ev[0] + ev[1] >= ev[2] - 1e-9)
        # principal axes must be well separated, otherwise "axis angle" is
        # meaningless: the eigenvector is arbitrary inside a degenerate subspace
        ratios = ev[1:] / ev[:-1]
        check(f"{name}_principal_axes_well_separated", np.all(ratios > 1.05))

        # the corner-origin bug (no parallel-axis shift) must change the tensor
        I_corner = np.trace(P) * np.eye(3) - P
        check(f"{name}_corner_wrongness_is_visible",
              np.linalg.norm(I_corner - I_ref) > 0.1 * np.linalg.norm(I_ref))
        print(f"  V={V_ref:.6f} c={np.round(c_ref, 6)}")
        print(f"  I eigenvalues={np.round(ev, 6)}")

    print("models_ok=%d" % (1 if ok else 0))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

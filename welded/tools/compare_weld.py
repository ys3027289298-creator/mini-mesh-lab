#!/usr/bin/env python3
"""Compare two weld_det outputs element-wise (cross-machine check).

Usage: compare_weld.py <out_A> <out_B> <step>

Prints: vertex sequence match rate, welded vertex count agreement,
face sequence agreement. On vertex mismatch, lists the tied vertex
pairs (same quantized key) around the first divergence, with their
coordinates and the quantization step.
"""
import sys

def load(path):
    vs, fs = [], []
    for line in open(path):
        if line.startswith("v "):
            vs.append(tuple(float(t) for t in line.split()[1:4]))
        elif line.startswith("f "):
            fs.append(tuple(int(t) for t in line.split()[1:4]))
    return vs, fs

def main():
    pa, pb, step = sys.argv[1], sys.argv[2], float(sys.argv[3])
    va, fa = load(pa)
    vb, fb = load(pb)
    n = max(len(va), len(vb))
    # element-wise vertex comparison (exact tuple equality)
    match = 0
    first_bad = None
    for i in range(min(len(va), len(vb))):
        if va[i] == vb[i]:
            match += 1
        elif first_bad is None:
            first_bad = i
    rate = match / n if n else 1.0
    print("vertex_seq_match_rate=%.6f (%d/%d)" % (rate, match, n))
    print("vertex_count A=%d B=%d equal=%d" % (len(va), len(vb), len(va) == len(vb)))
    print("face_seq_equal=%d (A=%d B=%d)" % (fa == fb, len(fa), len(fb)))
    if first_bad is not None:
        print("first divergence at welded index %d:" % first_bad)
        for tag, vs in (("A", va), ("B", vb)):
            if first_bad < len(vs):
                v = vs[first_bad]
                q = tuple(round(c / step) for c in v)
                print("  %s[%d] coord=%r qkey=%r step=%g" % (tag, first_bad, v, q, step))
        # list all quantized-key ties inside each output near the divergence
        for tag, vs in (("A", va), ("B", vb)):
            seen = {}
            for i, v in enumerate(vs):
                q = tuple(round(c / step) for c in v)
                if q in seen:
                    print("  tie in %s: welded %d and %d share qkey=%r step=%g"
                          % (tag, seen[q], i, q, step))
                else:
                    seen[q] = i
    return 0 if (rate == 1.0 and fa == fb and len(va) == len(vb)) else 1

if __name__ == "__main__":
    sys.exit(main())

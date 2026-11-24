import numpy as np
from quant_analysis import ciede2000
from de2000_ref import de2000_ref

PATH = "../output/ciede2000_test_data.txt"
rows = []
with open(PATH) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        # 1-indexed cols: 3-5 L1a1b1, 8 hp1, 16 dE, 17-19 L2a2b2, 22 hp2
        L1, a1, b1 = map(float, p[2:5])
        hp1_ref = float(p[7])
        dE_ref = float(p[15].split()[0])  # dE field is "2.0425 2" (merged set tag)
        L2, a2, b2 = map(float, p[16:19])
        hp2_ref = float(p[21])
        rows.append(((L1, a1, b1), (L2, a2, b2), dE_ref, hp1_ref, hp2_ref))

worst = 0.0
for i, (l1, l2, ref, _, _) in enumerate(rows, 1):
    g1 = float(ciede2000(np.array(l1), np.array(l2)))
    g2 = float(de2000_ref(np.array(l1), np.array(l2)))
    worst = max(worst, abs(g1 - ref))
    flag = "ok " if abs(g1 - ref) < 1e-3 else "FAIL"
    print(f"{i:2d} ref={ref:8.4f} mine={g1:8.4f} refimpl={g2:8.4f} {flag}")
print(f"\n{len(rows)} rows, worst |mine-ref| = {worst:.5f}")
print("ALL PASS" if worst < 1e-3 else "MISMATCH REMAINS")

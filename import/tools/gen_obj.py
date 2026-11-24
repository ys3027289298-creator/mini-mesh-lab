#!/usr/bin/env python3
"""gen_obj.py OUT N — 生成 N×N 网格的 OBJ，顶点按面角点重复输出（不共享），
模拟真实导出器产物，让焊接阶段有实际工作量。用法: python3 -I gen_obj.py out.obj 200"""
import sys

def main():
    out, n = sys.argv[1], int(sys.argv[2])
    with open(out, "w") as f:
        w = f.write
        idx = 1
        for i in range(n):
            for j in range(n):
                x0, y0 = i * 0.1, j * 0.1
                x1, y1 = x0 + 0.1, y0 + 0.1
                # 每个 quad 两个三角形，6 个角点全部重复写出
                for (x, y) in ((x0, y0), (x1, y0), (x1, y1),
                               (x0, y0), (x1, y1), (x0, y1)):
                    w(f"v {x:.6f} {y:.6f} 0.000000\n")
                w(f"f {idx} {idx+1} {idx+2}\n")
                w(f"f {idx+3} {idx+4} {idx+5}\n")
                idx += 6

if __name__ == "__main__":
    main()

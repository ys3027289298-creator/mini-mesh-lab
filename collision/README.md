# collision

碰撞检测代理流水线，C 实现，每帧四个阶段：

1. **运动** `scene_advance`：随机游走平移 + 小角度四元数旋转，wrap 回场景盒
2. **代理生成** `proxy_build`：每个物体把 `BODY_VERTS`(=256) 个局部网格顶点做
   四元数旋转+平移后拟合 AABB proxy —— O(n)，常数大、带真实内存流量
3. **粗检测** `broadphase`：所有 proxy 两两 AABB 重叠测试，O(n²)，输出候选对
4. **精检测** `narrowphase`：候选对上做精确包围球距离判定，O(候选对数)

## 构建

```
make test     # 单元自检（10 个用例）
make bench    # 压测驱动：10 → 2000，分阶段计时
```

压测建议绑核运行，避免调度漂移：

```
taskset -c 0 ./bench > ../artifacts/bench_mesh.csv 2> ../artifacts/bench_mesh.log
gcc -O2 -std=c11 -DCOLLISION_NO_MAIN src/bench_ext.c src/collision.c -o bench_ext -lm
taskset -c 0 ./bench_ext > ../artifacts/bench_ext.csv   # 2500–4000 外推点
```

## 实验设计

- 场景边长 L=∛n，密度固定：增长只反映算法阶数，不反映拥挤度
- 每物体 256 顶点、局部点云 ±0.03，球半径 ~0.052，场景稀疏（候选对 ≈0–5/帧），
  精检测不构成噪声
- 每点 2 个预热 block + 9 个计时 block，帧数自适应（小 n 多帧、大 n 少帧，
  单 block 目标 ~30ms），block 内取每帧耗时，再取 block 中位数；
  整体重复 3 次独立跑取中位
- 计时用 `CLOCK_MONOTONIC`，checksum 汇聚防止循环被优化掉

## 分析与报告（`tools/`）

```
PYTHONPATH=<user-site> python3 -P tools/analyze.py ../artifacts/bench_mesh*.csv
PYTHONPATH=<user-site> python3 -P tools/make_report.py ../artifacts
```

`-P` 把 cwd 排除出 sys.path（防同名模块注入），`PYTHONPATH` 显式指定依赖位置。
产出 `curve_loglog.png`、`exponents.png`、`analysis.json`、离线自包含 `report.html`
（内嵌数据、十字线悬浮、明暗双主题、数据表）。

## 关键结论（Ryzen 9700X 单核，gcc -O2）

| 区间 | 总耗时阶数 | 瓶颈 | 证据 |
|---|---|---|---|
| n ≲ 460 | ≈ n^1.02 | 代理生成（256 顶点/物体的旋转拟合） | 占帧时间 85–95%；顶点池在 n≈16 溢 L1d、n≈341 溢 L2 |
| 460 → ~1600 | 阶数从 1 漂到 1.5 | 过渡带 | 分段拟合二次项首次可测量（n≈460，恰在 L2 边界） |
| n ≳ 1700 | → n^≈2.1 | 粗检测暴力两两 AABB | n≈1619 粗检测追平代理生成；n=2000 时 200 万次测试/帧，占 55% |

- **拐点：n≈1600（瓶颈易主点；过渡带 460–1700）**
- 精检测全程 < 0.001ms，可忽略
- **60 FPS（16.7ms 帧预算）下最多约 3800 个物体**（实测插值 3830，
  二次拟合 3913；3500 个 14.2ms、4000 个 17.9ms）
- 扩容方向：把粗检测换成空间划分（均匀网格 / SAP 排序扫描），优化代理生成
  对 n>1600 的场景已无意义

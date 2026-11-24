# proxy

`make test` 构建并跑完整三档基准与 21 项自检（约 2–3 秒），样例输出见
[`sample_run.txt`](sample_run.txt)。

几百个可破坏物件不再逐个做碰撞检测：先按中心距离用并查集把邻近物件合成一个
碰撞代理（每档代理数稳定在 ~16–17，而物件数 71 / 253 / 387），每个代理在
**AABB** 与 **凸包**（monotone-chain）之间二选一。

## 模型（确定性、可复现）

- 2D 俯视碰撞世界；每个物件是一个带朝向的盒子（OBB），三档用同一组 9 个
  细长斜置"物件堆"，只改内部网格间距（3.2 / 1.75 / 1.35）得到低/中/高密度。
- 每帧 512 条从场地四边射入的碰撞探针（模拟弹道/查询），种子固定。
- 真值（oracle）= 探针与全部 OBB 逐个求交。
- 聚类半径 `CLUSTER_R = 4.2`，并查集单链接合并。

## 两个判据（同时满足才用凸包）

设一个代理含 `n` 个物件：

- **凹度判据**：实际凹度
  `C = 1 − 凸包内被某个 OBB 覆盖的采样点 / 凸包内采样点总数`
  （28×28 规则网格，落在凸包外的点拒绝）。越空越凹，`C → 1`；越实心 `C → 0`。
  判据：`C <= 0.60`（足够接近凸形才值得用凸包）。
- **紧密程度判据**：
  `T = 物件平均最近邻中心距 / 物件平均全宽`。排列越密 `T` 越小。
  判据：`T <= 1.80`。

```c
use_hull = concavity_gate(C, 0.60, +1)   /* C <= 0.60 */
        && tightness_gate(T, 1.80);      /* T <= 1.80 */
```

`n < 3` 的小簇凸包与盒子等价，一律走 AABB。

## "漏判"的口径（重要）

AABB 和凸包都由全部成员角点构造，二者都是成员几何的保守包围体，所以在窄相
预算无限时**几何漏判恒为 0**（自检 `conservative_proxies_no_geometric_FN`
对三档 × 两种判据连跑 6 帧验证）。包围体只可能多判、不可能少判。

因此这里的"碰撞漏判"建模真实引擎的另一机制：**每帧窄相 OBB 测试预算
`NARROW_BUDGET = 8192` 次**。松散的 AABB 让探针在空白角区命中代理、对整簇
做无效窄相扫描；预算耗尽后，当前及后续探针被丢弃，若 oracle 证明其中有真实
命中，则记一次漏判。计时路径不设预算，测的是完整碰撞负载。

## 三组九个数（正确凹度方向，本机一次样例运行）

| 密度 | 物件数 | 代理数量 | 单帧碰撞耗时 (µs) | 漏判次数 (/120 帧) |
|---|---:|---:|---:|---:|
| low  | 71  | **17** | **71.6**  | **0** |
| mid  | 253 | **16** | **81.7**  | **0** |
| high | 387 | **17** | **99.0**  | **0** |

机器可读行：

```
NINE_proxies=17,16,17
NINE_collision_us_per_frame=71.6,81.7,99.0
NINE_misses=0,0,0
```

同次运行的逐个物件基线耗时为 162.3 / 439.8 / 757.6 µs，代理方案加速
2.27× / 5.39× / 7.66×。代理数与漏判数是确定性的；µs 用 best-of-5×120 帧
计时，随机器/负载有 ±15% 量级波动。

## 实际凹度 × 所选代理形状（题目要求的并列对照）

| 密度 | 实际 C | 实际 T | C 门限 | T 门限 | 正确规则选择 (hull/box 数) | 反接凹度后选择 |
|---|---:|---:|---|---|---|---|
| low  | **0.762** | 2.57 | fail | fail | **bbox** (0/17) | bbox (0/17) |
| mid  | **0.630** | 1.38 | fail | pass | **bbox** (0/16，AND 拦截) | hull (9/7) |
| high | **0.472** | 1.07 | pass | pass | **hull** (9 hull / 8 小簇 box) | **bbox** (0/17) |

凹度随密度严格下降 `0.762 > 0.630 > 0.472`，门限 0.60 正好卡在中/高之间。

## 诊断：高密度漏判最多 ⇒ 凹度判据反了

程序内置反接对照，只把凹度门限从 `C <= 0.60` 改成 `C >= 0.60`
（紧密门限与 AND 结构不动）：

| 密度 | hull/box | 漏判 (/120 帧) | 单帧最差 |
|---|---|---:|---:|
| low  | 0/17  | 0 | 0 |
| mid  | 9/7   | 0 | 0 |
| high | 0/17  | **1301** | **48** |

反接后高密度档把最低凹度（0.472，本该用 hull）错判成 bbox，漏判 1301 次、
显著为三档最多，复现题述特征：

```
diagnostic_inverted_high_density_is_worst=1
diagnostic_correct_high_density_is_worst=0
```

因果链有窄相需求量佐证（基准第 0 帧，正确 / 反接）：高密度
`6739 / 8846` 次——正确规则低于 8192 预算，反接成 bbox 后超过预算；
中档 `6056 / 4820`、低档 `2357 / 2357`。即"凹度算反 → 形状选错 →
包围体变松 → 无效窄相吃光帧预算 → 漏判"。把不等号改回 `C <=` 后三档漏判
全部归零。

## 自检项

| 测试 | 含义 |
|---|---|
| `hull_square_4verts` / `hull_square_area_one` | monotone-chain 凸包基本正确性 |
| `ray_hits_box` / `ray_misses_box` | OBB 射线求交正负例 |
| `criteria_both_pass_hull` | 两判据同时满足才出 hull |
| `criteria_high_concavity_box` / `criteria_loose_arrangement_box` | 任一判据失败退回 bbox |
| `criteria_are_AND` | 两个都失败仍为 bbox |
| `inverted_flips_concavity_gate` | 反接对照只翻转凹度方向 |
| `concavity_sparse_gt_filled` | 受控稀疏环的凹度严格大于实心块 |
| `concavity_gate_direction_sane` | 凹度门限方向本身正确（实心过、稀疏不过） |
| `tier_concavity_monotone` / `tier_tightness_monotone` | 三档 C、T 随密度单调 |
| `proxy_count_compresses_objects` | 每档代理数少于物件数一半 |
| `high_density_correct_picks_hulls` | 高密度 9 个物件堆全部选中 hull |
| `conservative_proxies_no_geometric_FN` | 无限预算下代理几何漏判为 0 |
| `correct_budget_meets_zero_misses` | 正确判据在预算内三档零漏判 |
| `diag_inverted_high_density_worst` | 反接后高密度漏判三档最多 |
| `diag_correct_rule_not_high_worst` | 正确规则无此特征 |
| `inverted_box_wastes_narrow_phase` | 反接高密度窄相需求高出 10% 以上 |
| `timing_measured` | 计时路径确实执行且 sink 非零 |

## 文件

- `src/proxy.c` — 全部实现、基准与自检（无外部依赖，仅 libm）
- `Makefile` — `make` / `make test` / `make clean`
- `sample_run.txt` — 一次完整 `make test` 输出样例

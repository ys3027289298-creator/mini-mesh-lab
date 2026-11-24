# normal

`make test` 跑三套自检：

1. `normal` —— 原有的端到端自检（现经由拆分后的 U4 入口）。
2. `unit_tests` —— U1..U4 四个单元的分支级单测。
3. `bit_identity` —— 与拆分前实现（`tests/normal_ref.c`，git HEAD 原版的改符号副本）逐位比对。

## 拆分结构与依赖方向

```
U1 face_normal       (只用 sub/cross 等几何原语)
        │
        ▼
U2 accumulate_face   (调用 U1 的结果做散射累加)
        │
        ▼
U4 vertex_normals    (编排：空面守卫、清零、U1+U2 循环、U3 循环)
        ▲
        │
U3 normalize_vertex  (只用 len；与 U1/U2 无横向调用)
```

边只有 `U1→U2→U4` 和 `U3→U4`，单向无环。

## 逐位一致的保证方式

- 所有浮点表达式（含 `cross`、`dot` 内的乘加次序、`/= l` 的分量次序）逐字保留。
- 不引入任何每面临时缓冲：U2 仍对调用方那块 `out` 原地散射累加，面顺序、
  顶点顺序、分量顺序与原版完全相同，因此 float 加法的结合次序不变。
- `tests/bit_identity.c` 用两个阶段共 40 万组随机网格（含次正规数、带符号零、Inf、
  NaN、1e-12 阈值以下的极小坐标）对两个独立编译单元做比对。有限数分量必须逐位相同；
  NaN 分量只比对 NaN 属性——NaN 载荷经 SSE 运算的传播 C 标准不保证，同一份单体源码
  在不同内联决策下载荷就已经不同，源码层重构无法对此承诺。验证矩阵：gcc/clang ×
  -O0/-O1/-O2/-O3/-Os，加 `-ffloat-store -fno-inline -ffp-contract=off
  -fexcess-precision=standard`，以及 ASan+UBSan，全部通过。

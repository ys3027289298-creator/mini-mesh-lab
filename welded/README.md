# welded

`make test` 跑自检。

## weld_det：确定性顶点合并

`src/welded.c` 里的 `mesh_weld` 用 eps 距离贪心合并，"距离 ≤ eps" 不是等价关系（不可传递），分组依赖扫描顺序和浮点收缩（FMA），跨机器结果可能不同。

`weld_det` 按确定的全序合并：量化键 `q = llround(coord/step)`（int64），按 `(qx, qy, qz, 原始索引)` 字典序排序，量化键相同才合并，原始索引作唯一决胜键。前置条件：坐标有限、step > 0、|coord/step| < 2^62。

```
./weld_det mesh.obj 1e-3 > out.txt   # stdout: 规范化的 v/f 序列, stderr: 统计
```

跨机器验证：`tools/gen_models.py <dir>` 生成三个顶点分布不同的模型，两台机器（或两种工具链）各跑一次，`tools/compare_weld.py outA outB 1e-3` 逐项比较顶点序列与面索引。

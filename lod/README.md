# lod

LOD 链落位硬约束实验：每一级顶点必须落在上一级的原始顶点或边上，
另以相邻级最大表面偏差（采样对称 Hausdorff）为第二判据，并测量切级时
UV 贴图位移（256 纹素）。

```sh
make test     # 打印三模型 × 三级 × lock/free 对照与自检
make csv      # 输出 18 行逐指标 CSV
```

详见 `../artifacts/lod_report.md`，原始输出 `../artifacts/lod_run.txt`、
`../artifacts/lod_metrics.csv`、`../artifacts/lod_metrics.json`。

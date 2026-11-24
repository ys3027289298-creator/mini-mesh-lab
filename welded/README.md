# welded

`make test` 跑自检（含「导出偶发丢面」四类分支的最小复现 `br_*`）。
`make sanitize` 用 ASan/UBSan 跑同一套自检。

「导出后偶发少一块面」的逐条排查、最小输入、判定依据与修复说明见
[`REPORT.md`](REPORT.md)。结论：绕向冲突与索引互异的几何退化面在本模块没有
剔除分支；字面重复/越界索引在插入处被正确拒绝。真正会误丢合法面的是旧贪心
焊接对顶点顺序敏感、会让两个相距超过 ε 的角共用一个号，已改成 complete-linkage
（保证同号即距离 ≤ε，最小复现面在任意顶点顺序下都保留）；另修复了 `mesh_weld`
重映射时未校验面索引的越界读，并增加 `dropped_invalid` / `dropped_degenerate`
计数。

# mini-mesh-lab

网格处理基础模块，C 实现，每目录独立构建可自测。

三个模块共用 [`contract/mml_contract.h`](contract/mml_contract.h) 中定义的非持有 mesh 视图/缓冲契约；所有权、局部索引、无效值和演进规则见 [`contract/CONTRACT.md`](contract/CONTRACT.md)。

| 目录 | 内容 |
|---|---|
| `welded/` | 顶点焊接、面片重映射、退化面剔除 |
| `normal/` | 面积加权顶点法线、零长度累加器统计 |
| `bbox/` | 轴对齐包围盒、薄片检测 |
| `pipeline/` | 焊接 → 法线 → 包围盒的最小串联与演进兼容性检查 |

一键运行模块自测和串联场景：

```sh
make test
```

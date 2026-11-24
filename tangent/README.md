# tangent

从 UV 梯度重建每顶点切线框架：T、B、N 正交且右手系（`cross(T,B)=+N`），两条同时验。

- `B` 由 `cross(N,T)` 生成，正交性与手性按构造成立；原始 `dP/dv` 只用于检测镜像 UV。
- 无切向 UV 梯度的顶点（UV 塌陷、极点、退化 UV 三角形）退化到 UV 轴近似并打 `TAN_FALLBACK`。
- 镜像 UV 布局无法在保持 T 沿 `dP/du` 的同时满足右手系，翻转 B 满足手性并打 `TAN_UV_MIRRORED`。

`make test` 跑自检，打印三个 UV 布局模型的 正交性残差 / 右手占比 / 与 UV 主轴夹角，
并给出未正交化对照与 UV 各向异性比。

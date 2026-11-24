# mini-mesh-lab 共用网格契约

公共 C 契约在 [`mml_contract.h`](mml_contract.h)。三个模块对外都使用同一组类型：

| 模块 | 契约入口 |
|---|---|
| `welded` | `mml_weld_mesh(const mml_mesh_view*, const mml_weld_options*, mml_mesh_buffer*, mml_weld_report*)` |
| `normal` | `mml_compute_vertex_normals(const mml_mesh_view*, const mml_normal_options*, mml_vec3*, uint32_t, mml_normal_report*)` |
| `bbox` | `mml_bbox_compute(const mml_mesh_view*, const mml_bbox_options*, mml_box*)` |

模块内部仍可有自己的数学辅助函数，但跨模块边界不再传递各模块私有的 `V3`、`Tri` 或 `Mesh`。

## 1. 内存所有权

契约是**同步、非持有、调用方拥有内存**的描述符契约：

1. `mml_mesh_view` 只保存借用指针。调用期间，位置、索引、法线数组必须有效。
2. `mml_mesh_buffer` 也不拥有数组；`positions`、`triangles` 等都由调用方分配，并通过 `vertex_capacity`、`triangle_capacity` 说明容量（元素数，不是字节数）。
3. 被调用模块：
   - 不 `free`、`realloc` 或长期保存调用方指针；
   - 不在函数返回后保留引用；
   - 只写输出数组中声明容量以内的元素；
   - 内部临时内存必须在返回前释放。
4. `mml_weld_options::old_to_new`、法线输出数组、`mml_box`、report 同样归调用方所有。
5. 一次调用的输入和输出存储不应重叠；例如不要把同一个三角形数组同时作为 `in->triangles` 和 `out.triangles`。
6. 如果将来需要异步或堆对象所有权，应新增显式 owning 类型和 `release` 函数，不能把当前 view 偷偷改造成持有指针的对象。

典型生命周期就是 `pipeline/minimal_chain.c`：数组在 `main` 栈上，weld 写入调用方提供的数组，随后用这些数组构造一个只读 `mml_mesh_view` 传给 normal 和 bbox。

## 2. 索引是局部索引

- 索引类型固定为 `mml_index`（`uint32_t`）。
- 所有索引都是**零基、局部于当前 `mml_mesh_view` / `mml_mesh_buffer`**：
  - `mml_triangle::a/b/c` 索引同一个 mesh 的 `positions`；
  - `mml_weld_options::old_to_new[i]` 把 weld 输入的局部顶点 `i` 映射到 weld 输出的局部顶点；
  - 不存在跨 mesh、场景或文件的“全局顶点索引”。
- 全局 ID、DCC 名称、场景节点 ID 等属于应用层。需要时应另传 side-car 数组或扩展结构，不能复用局部索引。
- 合法索引范围是 `0 .. MML_INVALID_INDEX - 1`。

## 3. 无效值和错误表示

| 情况 | 契约表示 |
|---|---|
| 无效/缺失索引 | `MML_INVALID_INDEX`，即 `UINT32_MAX` |
| 函数调用失败 | 返回负值 `mml_status`，例如 `MML_E_NULL`、`MML_E_CAPACITY`、`MML_E_INDEX` |
| 成功但有非致命结果 | 返回 `MML_OK`，具体计数放 report，例如退化面数量、零累加器法线数量 |
| 可选指针/派生通道不存在 | 指针为 `NULL`，且对应 flag 不置位；当前有 `MML_MESH_HAS_NORMALS` |
| 无效顶点法线 | 全零向量 `{0,0,0}`，同时 report 给出数量；有效法线为单位向量 |
| 空/无效包围盒 | `mml_box::valid == 0`；此时 `min/max` 被清零且无意义，不依赖 `±1e30` 哨兵 |
| 非法输入浮点 | `NaN/Inf` 位置或法线会被拒绝为 `MML_E_NONFINITE` |
| 空 mesh 算法无结果 | normal/bbox 返回 `MML_E_EMPTY`；weld 空 mesh 可以成功 |
| 退化三角形 | 索引仍必须越界检查通过；是否删除/计为无效由对应模块和 report 决定 |

法线采用 `cross(b-a, c-a)` 的绕序；即三角形 `(a,b,c)` 按右手定则定义正面。坐标单位由调用场景统一，契约不做单位换算。

## 4. 演进时如何加字段而不破坏调用

所有跨模块结构都以 `mml_header` 开头：

```c
typedef struct mml_header {
    uint32_t size;        /* 调用方结构的 sizeof，单位字节 */
    uint32_t abi_version;
} mml_header;
```

规则：

1. 调用方用 `MML_STRUCT_INIT(&obj, Type)` 清零，并填入当前头文件中的 `sizeof(Type)` 和 ABI 版本。
2. 被调用方只要求 `header.size >=` 它认识的**最小前缀长度**，绝不要求 `size == sizeof(当前实现里的结构)`。因此更大的新调用方结构不会被旧实现仅因 size 不同而拒绝。
3. 新字段只能追加在结构尾部，不能插入中间、重排、改类型或改变旧字段语义。
4. 新字段必须是可选的：`0`、`NULL`、默认枚举值应表示“旧行为”。
5. 读取新字段前用 `MML_HAS_FIELD(&header, Type, field)` 判断调用方的 `size` 是否覆盖该字段；较短的旧调用方结构不能直接读取尾部字段。
6. 不兼容变更（改旧字段类型/语义、必填新数据等）必须提高 `abi_version` 或新增 `mml_*_v2` 函数；旧函数继续保留。新主版本实现可以显式兼容旧版本前缀。
7. 更大的实验性/非连续扩展适合做成带类型和大小的扩展链，而不是无限膨胀核心结构。

`normal` 的 `zero_epsilon` 就是按这个规则放置的尾部可调字段：零表示默认阈值。`pipeline/evolution_check.c` 验证了三种情况：

- 旧调用方只声明 `mml_header` 大小：调用成功，尾部阈值被忽略，使用默认值；
- 当前调用方提供完整字段：阈值生效；
- 调用方 size 更大（模拟未来追加字段）：当前实现接受，不访问自己不认识的尾部。

这保证的是 C 源码/二进制边界的平滑演进；如果所有模块始终在同一仓库里一起编译，这些检查主要是在强制边界清晰，而不是日常运行所必需。

## 5. 最小串联场景

运行：

```sh
make test
```

只跑串联和演进检查：

```sh
make -C pipeline test
```

场景使用 5 个输入顶点、3 个三角形：

- 顶点 3 是顶点 0 的精确重复；
- 顶点 4 是顶点 1 在 `1e-3` 容差内的近重复；
- `(3,4,2)` weld 后与 `(0,1,2)` 等价；
- `(0,3,1)` weld 后变成 `(0,0,1)`，作为退化面剔除。

实际串联结果：

```text
weld: 5 -> 3 vertices, 3 -> 2 triangles
map(local old -> local new): 0->0 1->1 2->2 3->0 4->1
normal[0]=(0.000000 0.000000 1.000000)
bbox min=(0.000 0.000 0.000) max=(1.000 1.000 0.000) sliver=1
chain_ok=1
evolution_ok=1
```

平面三角形的 z 范围为零，因此 bbox 按既有 `smallest_extent < ratio * diagonal` 规则报告 `sliver=1`；这正是 bbox 模块参与串联并给出判断的证据。

## 6. 什么时候契约会成为负担

这个契约不是免费的。以下场景中，它可能比直接传一个固定内部结构更慢或更繁琐：

- **单一模块、短期实验或脚本式程序**：没有第二个调用方，header、容量检查和 init 宏只是样板。
- **所有代码始终原子提交、一起编译**：不存在独立二进制/版本漂移时，size/version 协商的实际收益很小。
- **极致性能的逐顶点调用**：如果把契约函数放进热点内层、每个顶点调一次，校验和描述符开销可能不可接受；应保持当前的批量 mesh API，或另做特化快速路径。
- **数据布局不匹配**：契约固定为 `float xyz`、紧凑 `uint32_t` 三角形。若上游是 double、SoA/strided 缓冲、16 位索引或 GPU buffer，仍需转换；契约只统一“标准 CPU 边界格式”，不能消除所有布局差异。
- **需要全局身份、场景图或跨 mesh 关联**：局部索引刻意不解决这些问题；硬把全局 ID 塞进 `mml_index` 会制造错误。
- **异步、GPU 或长生命周期共享内存**：当前模型是同步借用，无法表达引用计数、等待 fence、显存所有权或延迟释放。
- **需求仍在快速变化的早期原型**：严格的追加式 ABI 规则可能让简单修改显得僵硬；这时可以把契约限制在真正跨模块的边界，模块内部继续快速演进。
- **大量可选扩展和版本分支**：如果每个调用都带很多扩展块，验证、兼容矩阵和测试成本会上升；应只在出现真实第二个版本/调用方时引入复杂度。

简言之：当“多个模块/调用方拥有不同内部结构或构建节奏”是真实问题时，这个契约有价值；当系统其实是一个紧耦合、同生共死、布局恰好一致的单体时，它本身就是额外抽象。

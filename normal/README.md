# normal

`make test` 跑自检。

## 面片尺寸相关的法线平滑阈值实验

复现：

```sh
python3 normal/threshold_lab.py
```

输出写入 `artifacts/normal-threshold/`：三张渲染图、`nine_numbers.csv`、
`metrics.json` 和区间图。

### 1. 两个不等式与可用区间

只统计两个面共享的内部边。设相邻单位面法线夹角为

\[
\theta_e=\arccos(n_{e,1}\cdot n_{e,2}),
\]

并约定 \(\theta_e\ge \tau\) 时为硬边。令

\[
L(h)=\max_{e\in\text{同曲率平滑边}}\theta_e(h),\qquad
U(h)=\min_{e\in\text{真折边}}\theta_e(h).
\]

要做到“同曲率边全部归入平滑组”且“真折边不被误平滑”，必须同时满足

\[
\boxed{\tau>L(h)}
\tag{1：下限}
\]

\[
\boxed{\tau\le U(h)}
\tag{2：上限}
\]

所以零误差可用区间是

\[
\boxed{L(h)<\tau\le U(h)}.
\]

交集非空的条件就是 \(L(h)<U(h)\)。若改用“大于阈值才是硬边”的约定，
两个等号位置互换，但交集判据不变。

无噪声时，曲率 \(k\)、弧长/面片尺寸 \(h\) 与相邻法线转角的关系为

\[
L_0(h)=2\arcsin\!\left(\frac{kh}{2}\right)
      \simeq kh+O((kh)^3),
\qquad
U_0(h)=\theta_f,
\]

其中 \(\theta_f\) 是固定折边角。因此理想区间为

\[
kh<\tau\le\theta_f.
\]

取区间中点附近作为阈值，

\[
\tau_0^*(h)=\frac{kh+\theta_f}{2},
\qquad
\boxed{\frac{d\tau_0^*}{dh}=\frac{k}{2}}.
\]

也就是说，固定物理曲率时，阈值的曲率项应与 \(h\) 同比例变化；\(h\) 缩小到
\(1/3\)，这项角度余量也应缩小到约 \(1/3\)。若阈值直接贴着下限取
\(\tau\simeq kh\)，则 \(d\tau/dh\simeq k\)。单位是
\(\mathrm{rad}/\text{长度单位}\)，换成角度要乘 \(180/\pi\)。

注意：如果只是把同一个已剖分模型整体均匀缩放 \(s\) 倍，则 \(h=sh_0\)、
\(k=k_0/s\)，所以 \(kh=k_0h_0\) 为常数；沿这种缩放路径
\(\boxed{d\tau/dh=0}\)。同一网格实例化为大石块和小碎石时，纯几何角度阈值
本来不应随世界尺寸改变。小碎石“全是硬边”更可能来自网格密度变化、固定坐标
误差/量化、远端浮点抵消或屏幕空间走样，而不是均匀缩放本身。

### 2. 固定坐标误差使小面片区间闭合

为复现题目所说的最小组区间为空，实验给每个顶点的 \(z\) 坐标加入固定的
独立高斯误差 \(\sigma_z=0.006\)（可解释为扫描精度或导出量化；它不随 \(h\)
缩小）。面法线的斜率误差约为 \(O(\sigma_z/h)\)，因此采用界

\[
L(h)\simeq kh+B\frac{\sigma_z}{h},\qquad
U(h)\simeq\theta_f-C\frac{\sigma_z}{h}.
\]

相应导数为

\[
\boxed{\frac{dL}{dh}=k-\frac{B\sigma_z}{h^2}},
\qquad
\boxed{\frac{dU}{dh}=\frac{C\sigma_z}{h^2}}.
\]

在小面片、误差主导的区段，\(dL/dh<0\)：\(h\) 越缩小，下限反而以约
\(B\sigma_z/h^2\) 的速率上升；同时上限随 \(h\) 缩小而下降。区间宽度约为

\[
W(h)=\theta_f-kh-(B+C)\frac{\sigma_z}{h},
\]

故最小组为空的条件是

\[
\boxed{kh+(B+C)\frac{\sigma_z}{h}\ge\theta_f}.
\]

这时不存在任何全局标量角度阈值能同时满足两个要求；单纯调高阈值只会开始
抹平折边，单纯调低则会把噪声平滑边判成硬边。

### 3. 三组九个数

三个“同类折叠岩面”使用同一曲面

\[
z=0.5kx^2+0.344|x|+0.03\sin(3y)x^2+0.05y^2+\epsilon_z,
\quad k=0.5,
\]

只是网格分辨率按 3 倍加密：\(6\times4\)、\(18\times12\)、\(54\times36\)。
全局阈值固定为 \(\tau_g=20^\circ\)。

“法线方差”主表采用与边判定直接对应的共边法线对方差：

\[
V_N=\operatorname{mean}_{\text{smooth edges}}
\frac{\|n_1-n_2\|^2}{4}
=\operatorname{mean}\sin^2\!\left(\frac{\theta_e}{2}\right).
\]

小角度下 \(V_N\simeq\overline{\theta_e^2}/4\)。

| 面片组 | 平均边长 \(h\) | 硬边数 | 被误平滑的折边数 | 法线方差 \(V_N\) |
|---|---:|---:|---:|---:|
| large  | 0.4438 | **4** | **0** | **0.001431** |
| medium | 0.1463 | **12** | **0** | **0.002587** |
| small  | 0.0491 | **1087** | **5** | **0.017839** |

这就是要求的三组九个数（后三列）。若改用全局面法线协方差迹
\(\frac1F\sum_i\|n_i-\bar n\|^2\)，三组分别是
\(0.265028,\ 0.267278,\ 0.293366\)；结构化数据同时保存了两种定义。

![large render](../artifacts/normal-threshold/render_large.png)

![medium render](../artifacts/normal-threshold/render_medium.png)

![small render](../artifacts/normal-threshold/render_small.png)

### 4. 两个不等式的实际取值并列

| 组 | \(L=\max\theta_{\text{smooth}}\) | \(U=\min\theta_{\text{crease}}\) | 实际可用区间 | \(\tau_g=20^\circ\) |
|---|---:|---:|---|---|
| large  | \(11.289^\circ\) | \(45.178^\circ\) | \((11.289^\circ,45.178^\circ]\) | 两条件都成立 |
| medium | \(18.793^\circ\) | \(27.808^\circ\) | \((18.793^\circ,27.808^\circ]\) | 两条件都成立，但下限仅剩 \(1.207^\circ\) 余量 |
| small  | \(60.745^\circ\) | \(1.438^\circ\) | **空集** | 两条件都不成立 |

最小组的两个实际要求并列是

\[
\boxed{\tau>60.745^\circ}
\qquad\text{对比}\qquad
\boxed{\tau\le1.438^\circ}.
\]

左侧要求阈值至少 \(60.745^\circ\)，右侧要求阈值至多 \(1.438^\circ\)，
上下限倒挂 \(59.307^\circ\)，交集为空。代入全局 \(20^\circ\)：

- 归组条件：\(20^\circ>60.745^\circ\)，不成立，阈值比下限低 \(40.745^\circ\)；
- 护边条件：\(20^\circ\le1.438^\circ\)，不成立，阈值比上限高 \(18.562^\circ\)。

**先破的是“阈值要足够大以归组”这一侧。** medium 组中它距全局阈值只剩
\(20-18.793=1.207^\circ\)，而护边侧还有 \(27.808-20=7.808^\circ\)。

空区间的两个极端操作点把两种视觉故障都量化了：

- 若强行满足下限，取 \(\tau>L=60.745^\circ\)，所有平滑边被接受，但
  **36 条真折边中有 34 条被误平滑**；
- 若强行满足上限，取 \(\tau\le U=1.438^\circ\)，36 条折边全部保住，但
  **5706 条平滑边中有 5278 条被判硬边**（总硬边数 5314），这正是“小碎石全是硬边”
  的极端情形。主表的 \(\tau_g=20^\circ\) 位于两个矛盾界之间，因此同时出现
  1056 条误硬平滑边和 5 条误平滑折边。

用三点 \(L\simeq ah+b/h\)、\(U\simeq\phi-c/h\) 拟合实测极值，得到

\[
L(h)\simeq0.1554h+\frac{0.05106}{h},\qquad
U(h)\simeq0.8288-\frac{0.04034}{h}\quad(\mathrm{rad}),
\]

其 \(L=U\) 交叉点约为 \(h_*\approx0.113\)，位于 medium \(h=0.146\) 与
small \(h=0.049\) 之间。固定 \(\tau_g=20^\circ\) 时，拟合的下限交叉点约
\(h_L\approx0.157\)，上限交叉点约 \(h_U\approx0.084\)，同样说明归组侧先破。
在 \(h=0.0491\) 处，实测拟合导数约为

\[
\frac{dL}{dh}\approx-1203^\circ/\text{长度单位},\qquad
\frac{dU}{dh}\approx+957^\circ/\text{长度单位}.
\]

即沿 \(h\) 减小方向，下限上升得比上限下降更快。极值和样本数有关，因此这两个
大导数应用于诊断，不应直接当成生产环境的稳定斜率；稳定估计应使用分位数或
已知的 \(\sigma_z\)。

### 5. 理想对照

关闭位置误差后，三组的实际极值区间分别为

| \(h\) | \(L_0=\max\theta_{\text{smooth}}\) | \(U_0=\min\theta_{\text{crease}}\) |
|---:|---:|---:|
| 0.4438 | \(8.703^\circ\) | \(45.396^\circ\) |
| 0.1460 | \(3.073^\circ\) | \(40.402^\circ\) |
| 0.0485 | \(1.049^\circ\) | \(38.725^\circ\) |

最小组区间并不空，反而最宽。对这组无噪声极值做 \(L_0\simeq a_0h\) 拟合得到
\(a_0\simeq0.345\ \mathrm{rad}/\text{长度单位}\approx19.8^\circ/\text{长度单位}\)，
所以贴下限的尺寸斜率约为 \(19.8^\circ\) 每模型长度单位，区间中点阈值的斜率约为
\(9.9^\circ\) 每模型长度单位。

这给出最终判断：

1. **理想曲率离散化**：小面片的相邻法线差更小，低角度全局阈值应更容易平滑，
   不会导致“小碎石全是硬边”；此时阈值可按 \(d\tau/dh\sim k\) 随 \(h\) 线性调整。
2. **固定绝对误差/量化或数值抵消**：小面片法线误差按 \(1/h\) 放大，下限上升、
   上限下降，可能出现题目所说的空区间。
3. **空区间下不存在可补救的单一全局阈值**。应先做尺度感知的降噪/焊接、扩大
   法线估计邻域，或直接保留折边/smoothing-group 标签；不能只调一个全局角度值。

![interval plot](../artifacts/normal-threshold/threshold_intervals.png)

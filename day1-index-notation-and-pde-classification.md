# Day 1：指标记号与 PDE 分类

> 日期：2026-11-03（周二）
> 主线：流体与数值、光学物理、AI4Science/HPC
> 今日目标：把向量/张量的指标记号和 PDE 分类这两块地基打牢，并理解它们如何决定后续数值方法的选择。

## 1. 向量与张量的指标记号

### 1.1 Einstein 求和约定

指标重复即自动求和：

$$
a\cdot b = a_i b_i,\qquad
(Av)_i = A_{ij}v_j,\qquad
(AB)_{ij} = A_{ik}B_{kj}.
$$

### 1.2 Kronecker delta 与 Levi-Civita 符号

- Kronecker delta：$\delta_{ij}$，等价于单位矩阵，作用是"替换指标"，$\delta_{ij}a_j = a_i$。
- Levi-Civita 符号：$\varepsilon_{ijk}$，完全反对称三阶张量。

$$
\varepsilon_{ijk}=
\begin{cases}
+1, & (i,j,k)\text{ 是 }(1,2,3)\text{ 的偶排列},\\
-1, & (i,j,k)\text{ 是 }(1,2,3)\text{ 的奇排列},\\
0, & \text{有任意两个指标相同}.
\end{cases}
$$

非零分量只有 6 个：

| $(i,j,k)$ | $\varepsilon_{ijk}$ |
| --- | --- |
| $(1,2,3)$、$(2,3,1)$、$(3,1,2)$ | $+1$ |
| $(3,2,1)$、$(1,3,2)$、$(2,1,3)$ | $-1$ |

### 1.3 常见算子的指标形式

| 算子 | 指标形式 |
| --- | --- |
| 梯度 | $(\nabla\phi)_i=\partial_i\phi$ |
| 散度 | $\nabla\cdot v=\partial_i v_i$ |
| 旋度 | $(\nabla\times v)_i=\varepsilon_{ijk}\partial_j v_k$ |
| 拉普拉斯 | $\nabla^2\phi=\partial_i\partial_i\phi$ |
| 对流项 | $(\mathbf v\cdot\nabla)\mathbf v$ 的第 $i$ 分量 $=v_j\partial_j v_i$ |

### 1.4 N-S 方程的指标形式

$$
\rho\Big(\partial_t v_i + v_j\partial_j v_i\Big)=-\partial_i p+\partial_j\tau_{ij}+\rho f_i,
$$

牛顿流体的粘性应力张量：

$$
\tau_{ij}=\mu(\partial_j v_i+\partial_i v_j)-\tfrac23\mu(\partial_k v_k)\delta_{ij}.
$$

注意对流项的两种写法：守恒形式 $\partial_j(v_j v_i)$ 与非守恒形式 $v_j\partial_j v_i$ 相差 $v_i(\partial_j v_j)$。不可压时两者相等，可压时不等；有限体积法（FVM）的守恒性正来自守恒形式。

### 1.5 Levi-Civita 符号的几何含义

$\varepsilon_{ijk}$ 的本质是三维空间的**定向（orientation）**编码器，即"有向体积的单位元"。

- 有向体积（混合积）：$a\cdot(b\times c)=\varepsilon_{ijk}a_i b_j c_k$。
- 叉积：$(a\times b)_i=\varepsilon_{ijk}a_j b_k$。
- 旋度：$(\nabla\times v)_i=\varepsilon_{ijk}\partial_j v_k$。
- 行列式：$\det A=\varepsilon_{ijk}A_{1i}A_{2j}A_{3k}$。

二维直觉：$\varepsilon_{ij}$ 对应矩阵 $\begin{pmatrix}0&1\\-1&0\end{pmatrix}$，即"逆时针旋转 90 度"；$\varepsilon_{ij}a_i b_j=a_1b_2-a_2b_1$ 是平行四边形的有向面积。符号即方向。

$\varepsilon_{ijk}$ 是**赝张量**（轴矢量来源）：在坐标反射 $x\to -x$ 下多乘一个 $\det=-1$。涡量、角速度、叉积、磁场在镜像下行为反常，光学里的手性、旋光性都源于此。

### 1.6 缩并恒等式

$$
\varepsilon_{ijk}\varepsilon_{ilm}=\delta_{jl}\delta_{km}-\delta_{jm}\delta_{kl},
$$

进一步缩并：

$$
\varepsilon_{ijk}\varepsilon_{ijm}=2\delta_{km},\qquad \varepsilon_{ijk}\varepsilon_{ijk}=6.
$$

用途：推导所有向量恒等式，例如 $\nabla\times(\nabla\times v)=\nabla(\nabla\cdot v)-\nabla^2 v$。

## 2. PDE 分类

### 2.1 二阶线性 PDE 与判别式

对二阶线性 PDE $A u_{xx}+2B u_{xy}+C u_{yy}+\text{低阶项}=f$，看判别式 $B^2-AC$：

| 判别式 | 类型 | 典型方程 |
| --- | --- | --- |
| $B^2-AC<0$ | 椭圆 | Laplace/Poisson、Helmholtz |
| $B^2-AC=0$ | 抛物 | 热/扩散方程 |
| $B^2-AC>0$ | 双曲 | 波动方程、对流方程 |

### 2.2 三类 PDE 的物理特征

| 类型 | 信息传播 | 边界条件 | 数值方法倾向 |
| --- | --- | --- | --- |
| 椭圆 | 边界全局决定内部，无时间方向 | 全边界 | 直接解大型稀疏线性系统、多重网格 |
| 抛物 | 单向时间方向，空间瞬时全局但衰减 | 初值 + 边界 | 常需隐式/ADI |
| 双曲 | 沿特征线有限速度传播，可有间断 | 初值 + 边界 | 迎风、FVM、Godunov/Riemann |

### 2.3 多维推广：系数矩阵特征值

判别式 $B^2-AC$ 只适用于两个自变量。多维时把方程写成

$$
\sum_{i,j}a_{ij}\frac{\partial^2 u}{\partial x_i\partial x_j}+\text{低阶项}=0,
$$

看对称系数矩阵 $A=(a_{ij})$ 的特征值符号（惯性指数）：

| 类型 | 特征值符号 |
| --- | --- |
| 椭圆 | 全部同号（正定/负定） |
| 抛物 | 有零特征值，其余同号（半定、退化） |
| 双曲 | 有一个与其余符号相反（惯性指数 $(1,d-1)$） |

### 2.4 热方程为什么是抛物

一维：$u_t=\nu u_{xx}$ 写成 $u_{xx}-\tfrac{1}{\nu}u_t=0$，二阶项只有 $u_{xx}$，故 $B=0,\ C=0$，判别式 $B^2-AC=0$，抛物。

$d$ 维：$u_t=\nu\nabla^2 u$，自变量 $(x_1,\dots,x_d,t)$，二阶系数矩阵为 $\mathrm{diag}(\nu,\dots,\nu,0)$：$d$ 个正特征值 + 1 个零特征值，半正定、退化，仍是抛物。**结论与空间维数无关**。

抛物型的本质：一个单向时间方向 + 空间瞬时全局传播。

对比：非定常热方程 $\partial_t u=\nu\nabla^2 u$ 是抛物；稳态热方程 $\nabla^2 u=0$ 是椭圆。抛物 = 椭圆 + 一个时间方向。

### 2.5 混合型：N-S 方程拆解

不可压 N-S 逐项标注类型：

$$
\underbrace{\partial_t u}_{\text{抛物}} + \underbrace{(u\cdot\nabla)u}_{\text{双曲}} = \underbrace{-\nabla p}_{\text{椭圆}} + \underbrace{\nu\nabla^2 u}_{\text{抛物}},\qquad \nabla\cdot u=0.
$$

- 抛物部分：$\partial_t u=\nu\nabla^2 u$（热方程，非定常粘性扩散）。
- 双曲部分：$\partial_t u+(u\cdot\nabla)u=0$（无粘 Burgers/Euler，有限速度、可产生激波）。
- 椭圆部分：对 N-S 取散度得到压力 Poisson 方程
  $$
  \nabla^2 p=-\nabla\cdot[(u\cdot\nabla)u],
  $$
  压力瞬时全局约束，无时间方向。

雷诺数 $Re=UL/\nu$ 是"切换旋钮"：$Re\to0$ 粘性主导（近椭圆 Stokes 流），$Re$ 中等是抛物+椭圆，$Re\to\infty$ 对流主导（趋近双曲 Euler）。

对求解器的直接推论：

1. 完整求解器 = 抛物时间推进 + 椭圆 Poisson 压力修正（投影法、SIMPLE/PISO 的本质）。
2. 高 Re 需迎风、限制器、Godunov/Riemann 处理双曲行为。
3. 低 Re 可用中心差分，主要矛盾是扩散稳定性 $\Delta t\le\Delta x^2/(2\nu)$。
4. 压力 Poisson 是性能瓶颈，是多重网格、Krylov、并行域分解的主战场。

### 2.6 Poisson 方程与 Laplace 方程

Poisson 方程：$\nabla^2 u=f$，$f$ 为源项。当 $f=0$ 时为 Laplace 方程 $\nabla^2 u=0$。

拉普拉斯算子 $\nabla^2 u=\nabla\cdot(\nabla u)=\sum_i\partial_i^2 u$，衡量某点值与邻域平均值的偏差。

调和函数（Laplace 解）的均值性质：任一点的值等于其球面邻域上的平均值，故内部无局部极值，完全由边界决定。

典型物理来源：

| 物理背景 | Poisson 方程 | 源项含义 |
| --- | --- | --- |
| 稳态热传导 | $\nabla^2 T=f$ | 热源 |
| 静电势 | $\nabla^2\phi=-\rho/\varepsilon_0$ | 电荷密度 |
| 不可压流压力 | $\nabla^2 p=-\rho\nabla\cdot[(u\cdot\nabla)u]$ | 速度场散度驱动 |
| 势流 | $\nabla^2\phi=0$（速度势） | 无源 |
| 二维流函数 | $\nabla^2\psi=-\omega$ | 涡量 |

数值上，离散 Poisson 方程即解大型稀疏线性系统 $Au=b$，这是第 2 天数值线性代数的主题。

## 3. 动手练习：三类 PDE 的典型解

用 matplotlib 画三类 PDE 的典型解（标题用英文以避免中文字体乱码）：

```python
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

# Elliptic: Laplace steady solution
x = y = np.linspace(0, 1, 200)
X, Y = np.meshgrid(x, y)
U = np.sin(np.pi*X)*np.sinh(np.pi*Y)/np.sinh(np.pi)
ax = axes[0]
cf = ax.contourf(X, Y, U, levels=20, cmap='coolwarm')
ax.contour(X, Y, U, levels=10, colors='k', linewidths=.5)
ax.set_title('Elliptic: Laplace steady solution')
ax.set_aspect('equal')
fig.colorbar(cf, ax=ax)

# Parabolic: heat equation diffusion
xd = np.linspace(-6, 6, 400)
nu = 1.0
ax = axes[1]
for t, c in [(0.1, 'r'), (0.4, 'b'), (1.0, 'g')]:
    u = np.exp(-xd**2/(4*nu*t))/np.sqrt(4*np.pi*nu*t)
    ax.plot(xd, u, color=c, label=f't={t}')
ax.set_title('Parabolic: heat equation diffusion')
ax.legend()

# Hyperbolic: traveling wave
ax = axes[2]
c0 = 1.0
for t, c in [(0, 'r'), (1, 'b'), (2, 'g')]:
    u = np.exp(-((xd-c0*t)**2)/0.8**2)
    ax.plot(xd, u, color=c, label=f't={t}')
ax.set_title('Hyperbolic: traveling wave')
ax.legend()

plt.tight_layout()
plt.show()
```

若需显示中文，在 import 后设置：

```python
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
```

## 4. 今日自测清单（闭卷检验）

1. 推导 $\varepsilon_{ijk}\varepsilon_{ilm}=\delta_{jl}\delta_{km}-\delta_{jm}\delta_{kl}$。
2. 写出 N-S 对流项守恒形式 $\partial_j(v_j v_i)$，并证明它与 $v_j\partial_j v_i$ 相差 $v_i\,\partial_j v_j$。
3. 判断类型：Poisson、热方程、波动方程、稳态不可压 N-S 的压力方程、Helmholtz、无粘/有粘 Burgers。
4. 解释为什么显式解 1D 热方程受 $\Delta x^2$ 限制，而 1D 对流只受 $\Delta x$ 限制。
5. 用自己的话说明"椭圆问题中，内部一点的值依赖全部边界条件"。

## 5. 今日关键结论速记

- 凡出现"旋转、叉积、定向、体积、手性"的公式，背后几乎都有 $\varepsilon_{ijk}$；凡出现"投影、替换、内积"的，是 $\delta_{ij}$。
- 抛物 = 椭圆 + 一个时间方向；稳态问题去掉时间项就退化为椭圆。
- N-S 是混合型：抛物（时间+粘性）+ 双曲（对流）+ 椭圆（压力 Poisson），雷诺数决定主导类型。

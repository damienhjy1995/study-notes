# -*- coding: utf-8 -*-
"""Day 2 动手：Jacobi / Gauss-Seidel / SOR 解 1D Poisson 方程，对比收敛速度。

方程：-u''(x) = f(x),  u(0) = u(1) = 0
三点中心差分离散后，内点线性系统为 A u = b，其中
    A 主对角为 2，次对角为 -1，b = h^2 f(x)。

只用 numpy + Pillow，不依赖 matplotlib / scipy。
"""
import numpy as np
from PIL import Image, ImageDraw, ImageFont


def A_mul(u):
    """A 与向量 u 相乘（A 主对角 2、次对角 -1，边界为 0）。"""
    Au = 2.0 * u.copy()
    Au[:-1] -= u[1:]
    Au[1:] -= u[:-1]
    return Au


def jacobi(b, tol=1e-8, maxiter=500000):
    m = b.size
    u = np.zeros(m)
    for k in range(maxiter):
        u_new = np.empty(m)
        u_new[0] = 0.5 * (b[0] + u[1])
        u_new[-1] = 0.5 * (b[-1] + u[-2])
        u_new[1:-1] = 0.5 * (b[1:-1] + u[:-2] + u[2:])
        u = u_new
        if np.max(np.abs(b - A_mul(u))) < tol:
            return k + 1, u
    return maxiter, u


def gauss_seidel(b, tol=1e-8, maxiter=500000):
    m = b.size
    u = np.zeros(m)
    for k in range(maxiter):
        u[0] = 0.5 * (b[0] + u[1])
        for i in range(1, m - 1):
            u[i] = 0.5 * (b[i] + u[i - 1] + u[i + 1])
        u[-1] = 0.5 * (b[-1] + u[-2])
        if np.max(np.abs(b - A_mul(u))) < tol:
            return k + 1, u
    return maxiter, u


def sor(b, omega, tol=1e-8, maxiter=500000):
    m = b.size
    u = np.zeros(m)
    for k in range(maxiter):
        u[0] += omega * (0.5 * (b[0] + u[1]) - u[0])
        for i in range(1, m - 1):
            u[i] += omega * (0.5 * (b[i] + u[i - 1] + u[i + 1]) - u[i])
        u[-1] += omega * (0.5 * (b[-1] + u[-2]) - u[-1])
        if np.max(np.abs(b - A_mul(u))) < tol:
            return k + 1, u
    return maxiter, u


def optimal_omega(n):
    h = 1.0 / n
    rho_j = np.cos(np.pi * h)  # Jacobi 谱半径
    return 2.0 / (1.0 + np.sqrt(1.0 - rho_j ** 2))


# ---------------- 主程序 ----------------
ns = [8, 16, 32, 64, 128]
results = {"jacobi": [], "gs": [], "sor": []}

for n in ns:
    h = 1.0 / n
    x = np.linspace(0, 1, n + 1)[1:-1]
    f = np.pi ** 2 * np.sin(np.pi * x)  # 真解 u = sin(pi x)
    b = h ** 2 * f

    kj, uj = jacobi(b)
    kg, ug = gauss_seidel(b)
    ks, us = sor(b, optimal_omega(n))

    u_true = np.sin(np.pi * x)
    results["jacobi"].append(kj)
    results["gs"].append(kg)
    results["sor"].append(ks)
    print(f"n={n:3d}  Jacobi={kj:6d}  GS={kg:6d}  SOR={ks:5d}  "
          f"| SOR err={np.max(np.abs(us-u_true)):.2e}")


# ---------------- 用 Pillow 画 log-log 收敛对比图 ----------------
def _font(path, size):
    try:
        return ImageFont.truetype(path, size, index=0)
    except Exception:
        return ImageFont.load_default()


YAHEI = "C:/Windows/Fonts/msyh.ttc"
F_LABEL = _font(YAHEI, 20)
F_TICK = _font(YAHEI, 16)

W, H = 900, 620
img = Image.new("RGB", (W, H), "white")
d = ImageDraw.Draw(img)

lx, ly = 90, 520  # 左下角
rx, ry = 850, 70  # 右上角
d.rectangle([lx, ry, rx, ly], outline=(40, 40, 40), width=2)

lxv = np.log10(ns)
all_iters = [v for k in results.values() for v in k]
ymin, ymax = np.floor(np.log10(min(all_iters))), np.ceil(np.log10(max(all_iters)))


def px(v):
    return lx + (np.log10(v) - lxv[0]) / (lxv[-1] - lxv[0]) * (rx - lx)


def py(v):
    return ly - (np.log10(v) - ymin) / (ymax - ymin) * (ly - ry)


colors = {"jacobi": (192, 57, 43), "gs": (31, 111, 178), "sor": (30, 132, 73)}
labels = {"jacobi": "Jacobi", "gs": "Gauss-Seidel", "sor": "SOR"}
for key in results:
    pts = [(px(v), py(v)) for v in results[key]]
    d.line(pts, fill=colors[key], width=3)
    for xp, yp in pts:
        d.ellipse([xp - 4, yp - 4, xp + 4, yp + 4], fill=colors[key])

legend_y = 90
for key in results:
    d.line([680, legend_y, 715, legend_y], fill=colors[key], width=3)
    d.text((725, legend_y - 11), labels[key], font=F_LABEL, fill=(40, 40, 40))
    legend_y += 34

# 坐标刻度
for v in ns:
    xp = px(v)
    d.line([xp, ly, xp, ly + 6], fill=(40, 40, 40), width=2)
    d.text((xp - 12, ly + 10), str(v), font=F_TICK, fill=(40, 40, 40))
for e in range(int(ymin), int(ymax) + 1):
    yp = py(10 ** e)
    d.line([lx - 6, yp, lx, yp], fill=(40, 40, 40), width=2)
    d.text((lx - 48, yp - 9), f"1e{e}", font=F_TICK, fill=(40, 40, 40))

d.text((lx + 12, ly + 38), "n (grid points)", font=F_LABEL, fill=(40, 40, 40))
d.text((20, ry + 30), "iterations", font=F_LABEL, fill=(40, 40, 40))
d.text((lx, 28), "Convergence: Jacobi / Gauss-Seidel / SOR (1D Poisson)", font=_font(YAHEI, 24), fill=(25, 45, 80))

out = "C:/Users/hmzxz/Documents/Codex/2026-11-03/wo-s/outputs/study-notes/day2_iteration_convergence.png"
img.save(out)
print("saved", out)

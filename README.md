# Multiscale Metal Simulation — PFC 拉伸位错演化

Phase-Field Crystal (PFC) simulation of dislocation mechanisms under tensile
loading, extending an existing multiscale framework
(DAMASK → DDD → MD, with PFC as the previously-unimplemented mesoscale branch)
toward 3D metal multiscale simulation.

由 Claude Code 自驱动迭代循环开发与测试；每轮迭代的状态与教训见
[PROGRESS.md](PROGRESS.md)，研究计划见
[docs/RESEARCH_PLAN.md](docs/RESEARCH_PLAN.md)。

## 结构

```
src/
  pfc2d.py            2D PFC（三角格子）：谱半隐式 + 稳定化分裂，rfft，
                      盒子仿射应变拉伸，位错偶极子相位绕动植入，缺口，应力测量
  pfc3d.py            3D PFC（BCC 一模）：同套数值格式，体积守恒拉伸
  defect_analysis.py  峰检测 → 周期 Delaunay → 5|7 配位 → 位错核心/密度
  run_b1_nucleation.py / run_b_series.py   生产运行驱动
tests/                门禁测试（A1 结晶 / A2 偶极子 / A3 拉伸 / C1 BCC）
results/              运行产物（图/JSON/日志入库；场快照 npz 仅本地）
```

## 门禁状态（全部通过）

| 门禁 | 内容 | 关键数字 |
|------|------|----------|
| A1a/b | 2D 结晶（熔体/种子） | 晶格常数偏差 0.3%/0.5%，完美晶体 frac6=1.000 |
| A2 | 刃位错偶极子植入弛豫 | 精确 2 核心稳定，缺陷占比 0.44% |
| A3a | 弹性拉伸 | 模量(能量法) 0.211 vs (应力法) 0.209 |
| A3b | 拉伸驱动滑移 | 2% 应变下滑移 ~4a₀，核心数恒定 |
| C1a/b | 3D BCC 结晶 | 峰数精确 2/胞，NN 偏差 4.1% |

## 复现

```bash
python3 tests/test_a1_crystallization.py
python3 tests/test_a2_dislocation_dipole.py
python3 tests/test_a3_tension.py
python3 tests/test_c1_bcc_3d.py      # 3D，约 8s
```

依赖：numpy, scipy, matplotlib, pyfftw（可选，缺省回退 numpy.fft）。

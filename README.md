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
| C2 | 3D BCC 弹性拉伸 | 应力过原点线性，模量 0.035，峰数守恒 |
| D1 | Interface B 读取与尺度映射 | ROI ρ≈1.5e15 m⁻² → 密度匹配只需 40nm/1400² 盒 |

## 旗舰结果：ROI 密度匹配的增殖级联（D2，1536²/42.5 nm）

从 DAMASK 热点 ROI 密度（2.2×10¹⁵ m⁻²，4 核心）出发拉伸至 21.5%：
**屈服峰 σ/E≈7% @11.6% → 位错雪崩（峰值 316 核心，ρ=1.75×10¹⁷）→
软化至流应力 = 49% 屈服峰**。完整的理论强度屈服→雪崩→流动转变。
图：`results/d2_roi_matched_1536/analysis_cascade.png`

## 机制结果（B 系列，512²，面积守恒拉伸至 10.5%）

- **偶极子**：拉伸驱动滑移贯穿全程，9.68% 应变处对湮灭（2→0）
- **多晶**（熔体淬火，~190 个 GB 位错）：模量为单晶 39%，~4% 后进入
  塑性流动（切线模量≈0，流应力≈0.004），GB 位错净减少 + 偶发再形核
  —— 应变辅助粗化
- **单晶**（含愈合孔洞对照）：弹性到 10.5% 不形核；液化盘缺口会重结晶
  愈合，需质量耗减型孔洞（`add_void(depth>0)`）做持久应力集中源
- 图：`results/b_series_512/analysis.png`

## 复现

```bash
python3 tests/test_a1_crystallization.py
python3 tests/test_a2_dislocation_dipole.py
python3 tests/test_a3_tension.py
python3 tests/test_c1_bcc_3d.py      # 3D，约 8s
```

依赖：numpy, scipy, matplotlib, pyfftw（可选，缺省回退 numpy.fft）。

# 核心问题攻坚（2026-06-13 战略转向）

用户决策：不收紧为确证论文，而是正面解决项目的核心问题。

## 核心问题
PFC 能否作为定量可信的介观桥梁，闭合 DAMASK(CP)→DDD→MD→PFC 多尺度信息环——
接收上层位错密度/应力，返回能实质改变下层预测的本构量？
（文献综述确认这是真实空白：无已发表工作把 PFC 经密度匹配契约接入 CP→DDD 链。）

## 子问题 A（前提）：PFC 力学定量正确
守恒 Model-B PFC 仅扩散弛豫弹性应变 → 应力/动力学不可信。必须验证：
- A1 位错弹性场 σ~b/r 正确（`dislocation_mobility.py` elastic_field）
- A2 Peach-Koehler 滑移迁移率 v=M·τ·b（同文件 measure_mobility，diffusive vs MPFC 对比）
- A3 MPFC 弹性弛豫格式：m(r) 趋势是否在快弹性弛豫下存活（`run_mpfc_check.py`）
工具：已加 `PFC2D.step_mpfc`（Stefanovic 惯性项，β=0 精确退化回 step()）

## 子问题 B（真耦合）：闭合一个环
最干净目标：PFC 位错迁移率律 v(τ,ρ,r) → 喂 DDD(ExaDiS) 的迁移率 B，或校准 DAMASK
dislotwin 的源激活/背应力参数 → 改变宏观预测。需 A 先成立。
资源就绪：DAMASK 3.0.2（单晶铜拉伸模拟/）、LAMMPS（HPC）、Interface B HDF5（CU/outputs）。

## 进行中
- run_mpfc_check.py：MPFC 下 m(r) 端点（r=-0.35,-0.21）vs 扩散版对比
- dislocation_mobility.py：单刃位错弹性场 1/r 检验 + v(τ) 迁移率（diffusive & MPFC）

## 已确立（本轮）
- 运动学分解证伪"滑移→晶界蠕变载体转变"标签（载体 r-不变）→ m(r) 是率敏感+旋节线软化

## 子问题 A 进展（2026-06-13）
- **A2 迁移率 ✓**：修正几何（偶极子沿Y分离、各自滑移面反向滑移）后测得干净正迁移率
  M≈0.93(diffusive)/1.05(MPFC)，v≈0.011 @ τ≈0.00157；**对动力学鲁棒（差13%）**。
  这是给 DDD 的可传递量 v=M·τ·b。
- **A3 MPFC m(r) ✓ 决定性**：率敏感阶梯在快弹性弛豫下存活——2D r=-0.35 m=0.166(扩散0.186)、
  r=-0.21 m=0.375(扩散0.402)。击退"m(r)是扩散假象"的审稿攻击。
- **A1 弹性场 ✗**：能量密度幂律拟合仍 nan，需改进（低优先，A2 是承重量）。
- **3D m_3D 矩阵污染**：36单元 m_3D 非单调(0.14/0.17/0.58离群/0.24)，3D BCC 边缘稳定性
  污染率敏感性 → 可辩护的 m(r) 在 2D 三角格子，3D 不可信（呼应 M17）。

## 子问题 B 下一步
测 v(τ) 多应力点 → 确认线性 Newtonian 阻尼(v∝τ, DDD 假设形式) → Interface-B 尺度图
换算物理单位迁移率 → 给 ExaDiS/DAMASK 的真实输入。

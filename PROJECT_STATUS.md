# 项目总账（诚实状态）— PFC/DDD 多尺度金属塑性

整条工作的入口与每条主张的诚实状态。详细文档见各子目录。

## 一句话
从"PFC 能不能迁移到 3D 金属硬化"出发，证明了 PFC 不行（区制错）、DDD 是区制正确的方法，并搭起一条 MD→DDD→真实 STEM 网络的多尺度链；定性图景清晰、若干定量关系经验证。collinear 主导调查（局部机制+长度标度真，但局部强度、pairwise-MFP、多滑移 flow 都不复现 canonical 集体主导）**已完结为 drift-limited bounded negative，封存为"DDD interaction-kernel validation / 协议依赖性基准"**（见 `interaction_matrix/multislip_flow/CONCLUSION.md`）。**当前主线回到 STEM-to-DDD v2 + 统一缺陷图表示（IDR）+ BO/UQ 多尺度框架。** 不是范式突破，是一条端到端搭起、关键环节经对抗式验证、边界清晰、且持续自我证伪的多尺度链。

## 主张状态表（✅稳固 / ⚠️有界 / ❌未解）

| 主张 | 状态 | 证据 / 边界 | 文档 |
|--|--|--|--|
| PFC（守恒扩散变体）拉伸不硬化 | ✅ | 6 类变体全软化；区制论证+钉扎诊断；限定于守恒变体 | `taylor_hardening/PAPER.md` |
| 3D BCC PFC 仿射拉伸非晶化 ~10–12% | ✅ | Bain 失稳 | PAPER.md §3.1 |
| DDD 复现 Kocks–Mecking 动态稳态 | ✅ | 任意初始密度收敛共同 ρ_ss；单一速率（ρ_ss(ε̇) 未测） | `taylor_hardening/RESULTS.md` |
| DDD 复现 Taylor √ρ **形式** | ✅ | 每构型内成立 | RESULTS.md |
| 干净 α（强障碍区制） | ✅ | forest-probe 修正载流子/拖曳后 α(准静态)≈0.96；远高于体材料（强障碍） | `rate_extrapolation.png` |
| 图像→DDD 路线端到端跑通 | ✅ | 真实 STEM 重建网络在 ExaDiS 演化、增殖 2.4× | `experiment_bridge/` |
| MD 锚定 Taylor 前因子 μb | ✅ | LAMMPS+Cu EAM → μ=55.2 GPa = DDD 54.6（1%） | `md_rung/RESULT.md` |
| 矩阵×STEM 清单预测 = 直接实测 | ✅ | 预测 0.476 = 实测 0.480（广义 Taylor 复合内部自洽） | `interaction_matrix/stem_validation.png` |
| 第一性矩阵复现**体材料** α=0.43 | ⚠️ | 数值落入 0.3–0.5，但是**逐类型不准的补偿平均**（见下） | `interaction_matrix/SYNTHESIS.md` |
| 受控二元 collinear **机制**忠实 | ✅ | 修掉两个几何 bug（中性方向 gap / 同号 Burgers）后，零应力 observer 稳定 `partial_annihilation_residual`，线长→0.71（~29% 湮灭），无稳定 junction = Madec collinear residual | `binary_reaction_matrix/PROGRESS.md` |
| 受控二元 collinear **remobilization 长度标度** | ✅ | τ_c 随残段长 350→140 MPa，满足 τ_c=K(μb/l)ln(l/b)，K≈0.75，R²=0.89 | `binary_reaction_matrix/collinear_scaling.png` |
| 逐类型**局部** remobilization ranking（6 类，机制全确认） | ✅(有界) | self 225 > collinear~glissile~coplanar 175 > Lomer~Hirth 125；Hirth 正确最弱，但 **collinear 不局部主导** | `binary_reaction_matrix/valid_ranking.png` |
| 第一性矩阵**逐类型忠实**（canonical collinear 主导） | ❌→已诊断根因 | **关键洞察**：canonical collinear ~5× 主导**不是单结点局部强度，而是集体 forest / mean-free-path 效应**（collinear 湮灭高效消耗*可动*位错、缩短 MFP）。这解释了为何 4 个局部协议+单二元反应都测不到它。→ 主线转为**演化森林 storage/MFP 测量**（Devincre–Kubin） | SYNTHESIS §7, `binary_reaction_matrix/PROGRESS.md` |
| STEM 实验锚定的物理保真 | ⚠️ | Burgers 几何指派（非 g·b）；z 轴弱约束被拉伸；是可行性演示非验证 | PAPER.md §3.5 |
| 完整 MD→DDD→CP 多尺度链 | ⚠️ | MD 只锚弹性 μ；迁移率/结点强度的 MD、CP（DAMASK）那端未做 | — |

## 贯穿全程的方法论（最可靠的产出）
"宣称 → 对抗式核查 → 撤回/修正"。多次"看似干净"的数都是假象，被更深一层检查抓出并诚实纠正：
- α=0.37"达体材料" → 拖曳基线过度扣除的假象
- 预测 0.69 → 我自己预测代码的单位 bug（直接验证抓出来的）
- 体材料 α=0.43 → 不准系数的补偿平均，逐类型不忠实
- 矩阵保真度 4 次修复 → 全失败，如实记录

## 子目录
- `taylor_hardening/` — PFC-vs-DDD + Taylor + KM + 5 审同行评审（`PEER_REVIEW.md`）+ 诊断（`DIAGNOSTICS.md`）+ 论文（`PAPER.md`）
- `interaction_matrix/` — 第一性相互作用矩阵 + STEM 预测/验证 + 综合（`SYNTHESIS.md`）+ 计划（`PLAN.md`）
- `md_rung/` — MD 弹性常数锚定
- `experiment_bridge/` — STEM 重建→DDD 适配与演化

## collinear 主导调查（Jun 19–21）：完结为有界负结果
追问"canonical collinear 主导（Madec/Devincre，~5×最强）是否在本 DDD 复现"，三层实验逐层剥离 observable：
- **局部强度**（`interaction_matrix/binary_reaction_matrix/`，已封存）：机制真（partial annihilation + 长度标度
  τ_c=K(μb/l)ln(l/b)，R²=0.89），但局部 remobilization ranking **不产生 collinear 主导**（mid-pack）。4 个 control 全过。
- **pairwise 演化森林 MFP**（`interaction_matrix/evolving_forest_mfp/`，已撤回）：attempt-3 null 采样不足（mobile 只滑过
  ~0.31 个 forest cell），且 pairwise 几何**结构上承载不了**多滑移 collinear 系数（b3=0、无储存通道、无载流子补给）。诚实撤回。
- **多滑移 flow-stress cell**（`interaction_matrix/multislip_flow/`，**主结论**）：共驱动（opt_pair）+ 增殖 FR 源 +
  可演化 forest + strain-rate 控制 + CrossSlip + per-system ledger + 全部 gate。**RSS 修正**后比 junction-type 对。
  **密度杠杆**（实测 ρ_f 1.65e13→1.5e14，~9×）：**R_RSS coll/glissile ≈ 1.0 两密度都平（G_τ≈1.03–1.06）、
  coll_opp/coll_same=1.00、XSLIP on/off 一致** → canonical 主导（应~2.3×且随密度增长）**不复现**。见 `CONCLUSION.md`。

**总判定**：canonical collinear **集体主导在本 ExaDiS/FCC_0/DDD_FFT/CrossSlip/FR-source 设置下不以 flow stress 或
storage-MFP 形式复现**——机制成立、集体系数不复现，是**有界负结果**（非失败）。诚实局限：固定 erate 下 flow stress 被
drag 主导（真 Taylor 需 rate-extrapolation）、co-driven collinear partner 结构性 drift 使 formal gate 在低密度
AMBIGUOUS（但结论对此稳健）。方法论价值：连续剥离 observable + 撤回不干净 null + RSS 修正纠偏。

## 未解的硬问题（如要继续）
1. **collinear 集体主导的干净测量**：需 rate-extrapolation/quasi-static + 完整多滑移群体统计（超出当前 protocol）。
2. **真实 g·b Burgers**：STEM 锚定升级到物理验证需要衍射数据。
3. **MD 迁移率/结点强度**、**DAMASK 那端**：补全多尺度链。

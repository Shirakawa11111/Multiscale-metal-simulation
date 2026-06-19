# 项目总账（诚实状态）— PFC/DDD 多尺度金属塑性

整条工作的入口与每条主张的诚实状态。详细文档见各子目录。

## 一句话
从"PFC 能不能迁移到 3D 金属硬化"出发，证明了 PFC 不行（区制错）、DDD 是区制正确的方法，并搭起一条 MD→DDD→真实 STEM 网络的多尺度链；定性图景清晰、若干定量关系经验证，但第一性相互作用矩阵不忠实（已诚实记录）。**不是范式突破，是一条端到端搭起、关键环节经对抗式验证、边界清晰的多尺度链。**

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
| 第一性矩阵**逐类型忠实** | ❌ | 漏掉 collinear 主导（经典 0.62 vs 我 0.13）；Hirth/self/Lomer 偏强 2.3–2.6×。**4 种协议尝试全失败**（湮灭/密度无关峰值/弓出基线主导）。需 Madec 构型受控几何，超出当前可干净实现范围 | SYNTHESIS §7, `matrix_vs_canonical.png` |
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

## 未解的硬问题（如要继续）
1. **忠实矩阵**：需 Madec/Devincre 双晶/偶极构型 + 准静态脱钉边界控制（4 种简化协议都不行）。
2. **真实 g·b Burgers**：STEM 锚定升级到物理验证需要衍射数据。
3. **MD 迁移率/结点强度**、**DAMASK 那端**：补全多尺度链。

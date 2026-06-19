# 失败诊断留档：两个 probe protocol 都不能直接测 canonical a_ij

排除错误路线（这不是浪费——是方法学诊断）。核心：**observable 选错了，不是参数没调好。**

## 1. `free_probe_matrix`（自由探针穿森林）
- 测的是：自由移动位错穿过固定森林的**有效流动/脱钉强度**。
- 失败模式（collinear）：移动位错遇同 Burgers 的 collinear 森林后**湮灭、被 topology 删除、载流子归零** → 测到的是"探针消失后的低流动应力"，**低估** collinear。
- 定位：保留，但**重命名为 effective glide-through forest flow strength**，**不再与 canonical 逐项硬比**，尤其不解释 collinear。

## 2. `fr_source_probe`（钉扎端 Frank-Read 弓出探针）
- 测的是：两端钉扎源**自身弓出**所需本征应力。
- 失败模式：峰值被源本征 τ_FR ≈ μb/L_src 主导（源长决定、与 ρ_f 无关）→ smoke 的 α_peak=0.238 是"FR 基线 / √ρ_f"的伪强度；扣 baseline 后**森林斜率 ≈ 0**（6 类全是）→ **高估/无效**。
- 定位：保留为 **negative control**——证明"防止探针消失 ≠ 测到森林强度"，区分 source strength 与 forest strength。

## 一句话
自由探针测**低**了（探针消失）；FR 源测**高**了（源本身太硬）。两者都不是 collinear a_ij 的正确测量。
→ 下一步：Madec 式**受控二元反应** + topology observer + 准静态 stress bisection。

# binary_reaction_matrix 进度（严格按用户路线图）

## Phase 0 — 失败诊断留档 ✓
`FAILURE_NOTE.md`：两个 probe protocol 重新定性为命名 negative controls
（自由探针=载流子湮灭测低；FR 源=源弓出基线 μb/L 测高）。

## Phase 1 — collinear 机制确认 ✓✓（里程碑）
受控二元 collinear 反应：同 Burgers 线、两交叉 {111} 面（b 在两面内）、45°、钉扎端有限段、
DDD_FFT + Topology + Collision、零应力松弛、topology observer。

**关键修正（迭代得到）**：
1. 初始 gap 沿"两面都不能滑移的中性方向" → 永不接触（observer 正确报 pass_through）。
2. 同 Burgers 同号 → 不湮灭（b+b=2b 非反应）。
3. **修法：反号 Burgers 对（+b, −b），起始在 rann 内交叉。**

**确认结果**（GAP=4/6/8 稳健）：
- mechanism = **partial_annihilation_residual**
- 线长 → **0.708**（29% 湮灭），无稳定 junction（junc=0）
- 即 canonical Madec collinear：partial annihilation + residual。

→ 机制对了，才有资格谈强度。

## 下一步（Phase 1 step 4–5 + Phase 2）
- [ ] observer 增加 residual 段长 l_res 的测量（用于标度）。
- [ ] **stress-control bisection**：施加 resolved shear σ=τ(b⊗n+n⊗b) 于 mobile 系统，
      固定 τ 松弛、判稳、二分搜 remobilization τ_c（判据=拓扑事件：残段解体/重新动员/
      mobile length 恢复/塑性跳变，不只看应力曲线）。
- [ ] **长度标度**：扫 L/b=800/1500/3000/6000，画 τ_c/μ vs (b/l_res)ln(l_res/b)，
      线性则机制核心抓住，再映射 a_ij。
- [ ] 四个 control：LineTension null / remesh-rann 收敛 / box 收敛 / orientation family。

## Phase 1.5 — collinear 长度标度 ✓✓✓（已验证）
`collinear_strength.py`：松弛形成残段 → 固定 resolved shear → remobilization。
τ_c 随残段长下降（350→140 MPa），**τ_c=K(μb/l)ln(l/b)，K≈0.75，R²=0.89**（`collinear_scaling.png`）。
canonical collinear remobilization 物理捕获。诚实：大 L 轻微饱和，K 偏高需细化 convention。

## Phase 3 第一次尝试 — ranking 无效（诚实）
`binary_strength.py` 泛化到 6 类后跑 ranking，但**几何 bug**：段方向用各自 signed Burgers，
collinear（B2=−B1）线方向翻转 → A、B 成同一平行位错 → **不湮灭**（anneal=1.00 vs 0.71）。
所得"排序"（collinear 最弱）是**假象**，不作数。**专门测的 collinear 结果仍成立。**
**教训**：跳过了"每类先确认机制"。

## 下一步（修正后）
- [ ] 逐类零应力 observer：确认 collinear=湮灭、Hirth/glissile/Lomer=结点形成、
      coplanar/self 的正确反应；修各自几何（线方向用一致参考系，非 signed burg）。
- [ ] 只对确认机制的构型测 remobilization τ_c → ranking。
- [ ] 四个 control + a_ij 映射。

## 主线并行
STEM-to-DDD 不被 matrix 阻塞；matrix 作为校准/解释层。

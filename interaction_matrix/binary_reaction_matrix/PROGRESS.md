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

## 主线并行
STEM-to-DDD 不被 matrix 阻塞；matrix 作为校准/解释层。

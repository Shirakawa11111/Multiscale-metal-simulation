# PFC 拉伸位错演化 × 3D 金属多尺度扩展 — 研究计划

> 创建日期: 2026-06-11 · 由自动迭代循环维护
> 状态文件: ../PROGRESS.md（每轮迭代更新）

## 1. 定位与动机

现有资产盘点（详见各项目原文档）：

| 项目 | 资产 | 对本项目的作用 |
|------|------|----------------|
| `CU/` | DAMASK→DDD(ExaDiS)→MD 分层框架 + Interface B v1 数据契约；**PFC 为已设计未实现的 Layer 4b** | 本项目就是把这条支线做实，并通过 Interface B 与主线对接 |
| `单晶铜拉伸模拟/` | DAMASK 2000-cycle 单晶 Cu [111] 数据、3 个 ROI 的 Interface B HDF5、自研 CP+PFF 求解器 | ROI 的位错密度 ρ 作为 PFC 初始缺陷种子；HPC 运行范式 |
| `BJK_hBN_Multiscale/` | XPFC 谱方法实现（pyFFTW）、peak_reconstruction（密度场→原子坐标） | 数值格式与峰值重建直接复用思想 |
| `Topology research/` | LAMMPS 拉伸分析、5\|7 缺陷拓扑识别、HPC 提交范式 | 5\|7 配位识别方法 = 2D 三角格子中的刃位错识别 |

**科学问题**：拉伸加载下，PFC 尺度（扩散时间尺度、原子空间分辨率）位错的形核—滑移—湮灭机制如何随应变率/温度参数 (r)/初始密度变化？能否与 DAMASK ROI 的位错密度状态量建立可交换的映射（Interface B `/mapping`）？

## 2. 技术路线（分阶段，每轮迭代推进）

### Phase A — 2D PFC 核心（本地可跑）
- 标准 Elder–Grant PFC：F = ∫ [ψ/2 (r + (1+∇²)²)ψ + ψ⁴/4]
- 保守动力学 ∂ψ/∂t = ∇²(δF/δψ)，谱半隐式格式（精确复用 hBN 项目的数值风格）
- 验收门禁：
  - A1 结晶 smoke：能量单调下降、三角格子形成、峰间距 ≈ 4π/√3（±5%）
  - A2 位错偶极子：植入 ±b 刃位错对，弛豫后 5|7 配位对稳定存在
  - A3 单轴拉伸：盒子仿射应变（k 向量重标定）+ 应变步进，观察位错形核/滑移
- 缺陷检测：峰值提取 → Delaunay 三角化 → 配位数 → 5|7 对 = 刃位错；输出位错密度 ρ(t)、位置轨迹

### Phase B — 拉伸机制扫描（本地 + HPC 空闲时）
- 参数矩阵：应变率 × r（有效温度）× 初始位错密度（多晶/单晶/预置偶极子）
- 输出：应力-应变曲线（通过 δF/δε 或 virial-like 表达式）、位错密度演化 ρ(ε)、形核应变、滑移轨迹统计
- 与 `Topology research/analyze_tensile_scan.py` 同款指标：模量、峰值强度、韧性

### Phase C — 3D PFC（金属多尺度核心扩展）
- 一模 3D PFC 自然给出 BCC；FCC 需两模 XPFC（与 CU 的 Cu 对应）
- 先 BCC 验证 3D 数值格式（128³ 本地可行性测试 → HPC 256³+）
- 位错线检测：3D 峰值重建 + DXA-like 思路（或导出原子坐标用 OVITO/ovito python 模块）

### Phase D — Interface B 对接（多尺度闭环）
- 读取 `CU/outputs/interfaceB_main/*.h5` 的 ROI 位错密度与驱动历史 σ_ij(t)
- 密度种子规则：ρ_target → PFC 初始缺陷数 N = ρ·A（2D）/ ρ·V（3D）
- 输出回灌：PFC 位错密度演化、湮灭/形核事件计数 → 与 DDD event proxy 对比

## 3. HPC 使用约定
- 主机 `rc@192.168.196.40`（128 核，LAN）——**仅空闲时用**；每轮迭代先 `uptime`+`ps` 检查负载，load average > 16 或有他人任务在跑则不提交
- 提交范式参考 `Topology research/hpc_category_B/SUBMIT_ALL.sh`（rsync 推送 → nohup 后台 → rsync 拉回）
- 2026-06-11：连接超时（不在局域网内），全部本地计算

## 4. 目录结构
```
PFC_Multiscale_Extension/
├── docs/          # 计划、方法学笔记
├── src/           # pfc2d.py, defect_analysis.py, (pfc3d.py)
├── tests/         # 门禁测试（A1/A2/A3...），可重复运行
├── results/       # 每次运行的 npz/png 产物，按日期-标签命名
├── literature/    # 关键文献笔记
└── PROGRESS.md    # 跨迭代状态：已完成、进行中、下一步
```

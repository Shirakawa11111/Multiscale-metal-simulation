# STEM→DDD 路线评估：成功（2026-06-14）

问题：什么方法能让图像→模拟耦合成功？答：**DDD（ExaDiS），不是 PFC**。已实测验证。

## 路线诊断（为何 DDD 而非 PFC）
PFC 双错：原子分辨（实验密度要 4096² 盒）+ 扩散动力学（=高温蠕变，非室温无热硬化）。
DDD 双对：**位错线分辨**（=重建输出）+ **室温无热滑移+结点+森林硬化**（=实验区制）。
且 ExaDiS 已在框架内（exadis_adapter_generator / exadis_minimal_run / FR 源输出）。

## 实测结果
1. **图像→DDD 输入适配器 ✅**：`stem_to_exadis.py` 把 27 条真实重建线 → 270 节点/243 段的
   合法 ExaDiS 网络（exadis_python_manual_network_v0），验证通过（段索引合法、Burgers 单位 FCC、
   几何赋 4 个滑移系、2×2×0.15µm 薄箔盒 z 非周期）。
2. **regime-correct 滑移动力学 ✅**：`minimal_ddd.py` 在真实网络上，剪切应力下 PK 力驱动、
   约束 {111} 面：**99% 节点滑移、平均 28 b**——无热滑移区制（PFC 做不到）。
3. 调试印证关键要求：**Burgers 须物理赋值（g·b）**——x 单轴拉伸对几何赋的 yz-Burgers 系统
   Schmid=0（无力）；改剪切后正常。几何赋值不够，需实验 g·b 衬度分析。

## 结论：路线成功，剩下的是部署+实验输入（非研究风险）
- ✅ 耦合关键件（图像→DDD）通了；✅ 动力学方向对（无热滑移）。
- 待补（已知工程/实验，非研究不确定性）：
  1. **部署 ExaDiS**（Kokkos C++ 源码编译，本机/HPC）→ 跑全 DDD（结点/增殖/森林硬化——
     ExaDiS 的既定能力，正是 PFC 做不到的）。
  2. **实验 Burgers（g·b）** 每条线 → 真实滑移系/号。
  3. **恢复主导时序数据**（高温蠕变/保载弛豫/循环内细时间）→ 真预测验证（管道已就绪）。
- PFC 唯一幸存件（线性迁移率形式 v=M·b·τ, M=1.54）可反喂 DDD 的迁移率律。

## 全闭环（regime-correct）
STEM 降噪 → 3D 重建线 → [g·b 赋 Burgers] → ExaDiS DDD（滑移+结点+硬化）→ 比对下一帧
   ↕ MD 标定迁移率/结点规则（下层）   ↕ DAMASK 连续介质（上层）
PFC → 退为高温蠕变支（不同区制）。

**一句话**：PFC 探索严谨地排除了错误的层、指向并实测验证了对的层（DDD）。图像→DDD 路线成功。

## 端到端实证：真 ExaDiS DDD 跑通真实重建网络（2026-06-14）
**ExaDiS 已在 HPC 部署**（cmake 4.3 + conda-fftw + pybind11；OpenMP+Serial 后端；
pyexadis 导入 OK；GPU CUDA12.4 可选未用）。
**真实 STEM 重建网络在真 ExaDiS DDD 跑通**：`run_stem_exadis.py` 加载 stem_network.json
（270 节点）→ 施加剪切（σxy=100/σyz=60 MPa）→ EulerForward + LineTension + SimpleGlide +
Retroactive collision + LengthBased remesh，30 步，节点 270→2349，RUN TIME 1.6s，status ok。
**图像→模拟路线端到端执行成功（真 DDD 引擎，非代理）。**

## 路线最终状态
| 环节 | 状态 |
|------|------|
| 图像→DDD 网络适配 | ✅ 实测（27 线→合法网络）|
| ExaDiS 部署 | ✅ HPC 构建成功，pyexadis OK |
| 真网络跑真 DDD | ✅ 30 步演化完成 |
| 森林硬化测定 | 待开 topology(结点)模块 + 弹性力 + 应力-应变（ExaDiS 既定能力，下一步）|
| 真实 Burgers g·b | 待实验衬度分析 |
| 真预测验证 | 待恢复主导时序数据（管道已就绪）|

**对比 PFC**：PFC 在此区制根本跑不通（尺度+区制双错、不硬化）；DDD 端到端跑通。
路线成功，方法选对。

## 森林硬化测定：DDD 捕获 PFC 做不到的储存机制（2026-06-14）
开 TopologyParallel（结点）+ SUBCYCLING N体弹性力 + FCC_0 迁移率，应变率加载真实重建网络
（全周期、加厚 z 盒、erate=1e4，到 0.15% 应变，24 分钟）：
- **决定性正信号：位错密度增殖 3.2e12→7.7e12 m⁻²（2.4×）** = 位错储存/增殖（Frank-Read/结点介导）。
  **这正是 PFC 在所有变体（2D/3D/抑制攀移/两模/双修/钉扎）下都做不到的——PFC 总是软化、密度降。**
- 应力-应变：弹性加载→屈服峰 54 MPa @ 0.042%→锯齿流动（位错雪崩，小尺度特征）。
- **诚实边界**：干净的单调加工硬化曲线 + 定量 Taylor α 还没到——小系统（几条线）+ 0.15% 小应变
  → 锯齿流动主导，σ-√ρ 相关 -0.39（不足以宣称 Taylor 斜率）。需更大系统+更长应变出统计硬化信号。
图：experiment_bridge/results_exadis/stem_hardening.png

## 总判定：图像→DDD 路线成功，DDD 在 PFC 失败的区制工作
| | PFC | DDD(ExaDiS) |
|--|-----|-------------|
| 真实重建网络能跑 | ❌ 尺度+区制双错 | ✅ 端到端 |
| 加载下密度演化 | 软化(降) 所有变体 | **增殖(升) 2.4×** |
| 区制 | 高温扩散蠕变 | 室温无热滑移+结点 |
**结论**：耦合用 DDD 成功，且 DDD 捕获了储存/增殖（硬化的基础）——PFC 根本做不到的。
干净硬化曲线只需放大系统/应变（既定能力，规模问题非原理问题）。

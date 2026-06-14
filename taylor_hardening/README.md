# Taylor 加工硬化定律：DDD 密度系列（定量封顶实验）

## 动机
整条主线已确立：
- **PFC 无法森林硬化**（所有变体都软化；机制 = 攀移驱动森林解锁 + Orowan 绕过）。
- **DDD（ExaDiS）在正确区制工作**：真实 STEM 重建构型在加载下**位错增殖 2.4×**（储存机制，PFC 根本没有）。

但单个小系统的运行是**锯齿流动主导**的，所以*逐点* σ–√ρ 相关无意义：雪崩时应力**下降**而密度**跃升**，二者逐点反相关（这正是第一次尝试得到伪负 α 的原因）。

## 正确的测量：流动应力密度系列
Taylor 定律的教科书测法是**流动应力**对**密度系列**：
- 多个独立 FCC 网络，初始密度 ρ₀ 递增（`generate_line_config`）；
- 每个在结点开启（TopologyParallel）下加载到稳态塑性流动；
- 取塑性平台的平均流动应力 σ_flow 与平均密度 ρ_flow；
- 拟合 **σ_flow = α·μ·b·√ρ**，得 Taylor 系数 α。

Cu 文献值 **α ≈ 0.3–0.5**（Bertin et al., MSMSE 2019; Madec/Devincre）。落在此区间的正 α，就是 PFC 根本缺失的定量加工硬化签名。

## 数值设置
- 引擎：ExaDiS（LLNL），pyexadis，HPC OpenMP。
- 力：`DDD_FFT_MODEL`（长程弹性森林相互作用 = √ρ 标度的来源）+ `Trapezoid` 积分（每步 2 次力计算，远快于 Subcycling 的多次子循环）。
- 迁移率：`FCC_0`（Medge=Mscrew=64103, vmax=4000）。
- 碰撞 `Retroactive` + 结点 `TopologyParallel` + 重网格 `LengthBased`。
- FCC Cu：b=2.55e-10 m, μ=54.6 GPa, ν=0.324。
- 盒 Lbox=10000 b（≈2.55 µm），密度由 num_lines 设定；应变率加载 edir=[0,0,1]。

## 文件
- `taylor_series.py` — 密度系列驱动（HPC 运行，输出 `taylor_out/taylor_series.json`）。
- `plot_taylor.py` — 出图：各密度应力–应变曲线 + σ_flow vs μb√ρ 的 Taylor 拟合。
- `RESULTS.md` — 结果与判定（运行后填写）。

## 运行
```
OMP_NUM_THREADS=32 NUM_LINES=10,20,40,80 \
  PYTHONPATH=~/BO/exadis_src/python python3 taylor_series.py
python3 plot_taylor.py taylor_out
```

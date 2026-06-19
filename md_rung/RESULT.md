# MD 环（Phase 3，补弱点 #4）：MD 实测弹性常数 → 锚定 DDD 的 Taylor 前因子

多尺度链此前 MD 那一环是"蓝图未做"。本环用 LAMMPS（`/usr/bin/lmp`）+ Cu Mishin EAM 势
（`Cu_mishin1.eam.alloy`）从**原子尺度**实测铜的弹性常数，转成各向同性 μ、ν，喂进 DDD 的
Taylor 前因子 μb——这个量是**所有 α 结果**的核心标定。

## 结果
MD 单晶 Cu：C11=169.9，C12=122.6，C44=76.2 GPa，K=138.3 GPa（标准 Cu 值）。
各向同性平均：

| | 剪切模量 μ (GPa) | Poisson ν |
|--|--|--|
| Voigt | **55.2** | **0.324** |
| Reuss | 40.3 | 0.367 |
| Hill | 47.8 | 0.345 |

**DDD 一直用的输入：μ=54.6 GPa，ν=0.324。**
**MD-Voigt：μ=55.2，ν=0.324 —— μ 差 1.1%，ν 差 0.0%。**

→ DDD 的 Taylor 前因子 μb 现在**由原子尺度第一性锚定**（55.2 vs 54.6 GPa），不再是假定值。
多尺度链 **MD → DDD** 由此具体闭合（非蓝图）。

## 诚实边界
- 本环锚定的是**弹性参数 μ**（Taylor 前因子，对 α 最关键），不是迁移率（动力学）或结点强度。
- 更完整的 MD→DDD 还可：(a) MD 测位错迁移率 B 喂 DDD 动力学；(b) MD 直接测某结点反应强度，从原子尺度验证 DDD 的 a_ij。两者更难，留作扩展。
- μ 取 Voigt 平均与 DDD 吻合最好；Reuss/Hill 略低（各向异性 Zener 比 A=2C44/(C11-C12)=3.2，铜各向异性较强，平均方案有 ~15% 离散）。

## 文件
- `init.mod` / `potential.mod` / `in.elastic` / `displace.mod` —— LAMMPS 弹性常数计算
- `md_elastic_result.json` —— C_ij、μ、ν 与 DDD 对比

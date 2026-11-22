---
name: bioopt-auto-optimization
description: Reproduce an official bioinformatics small test, profile the actual workflow, and optimize it under an explicit scientific equivalence contract. Use for speed or memory work on scientific pipelines, including spatial transcriptomics and 3D reconstruction, with matched benchmarks and cross-method transfer.
---

# BioOpt 科学计算优化

在同一科学任务和资源预算下减少运行时间、内存或无效数据搬运。每轮交付可重跑的基线、候选、等价证据和实测范围。

## 工作流

1. **定位官方小试。** 找作者发布的最小示例，固定版本、输入、权重和实际命令。官方小试缺失或数据不全时记录出处和阻塞，跳过该方法的正式复现与优化。已有项目产物可用于独立组件迁移，组件实验不能代替作者小试。原仓库与共享环境只读。候选、缓存、环境和日志放在隔离位置，计算沿用用户授权。
2. **冻结科学合同。** 用官方入口跑通基线，记录输入摘要、版本、实际训练预算、配额和精度。提前定义输出字段、ID/轴序、随机流、允许忽略的元数据和浮点规则。用两次独立基线自比较辨别随机差异与实现差异。细则见 [等价验收](references/equivalence.md)。
3. **沿实际入口剖析。** 按初始化、读取、预处理、拟合或几何、表达预测和写出计时。带着实测热点选择一个可证伪的假设。将剖析结果映射到 [优化决策](references/patterns.md)。
4. **改一处，先验正确性。** 保留独立原版入口和必要计算、科学参数、矩阵尺度及接口。缓存或批处理前检查依赖是否随查询、随机状态或可变输入改变。先小试再扩量。候选未通过等价检查便归档失败原因。若精度、停止规则或随机合同发生变化，按明确声明的变体另行评估。
5. **配对测量。** 在同一设备、环境、配额、输入下测量，固定冷热条件并交替顺序。Linux CLI 可用 [配对工具](scripts/paired_commands.py)，每对调用独立验证命令，失败即停。区分端到端、重建、训练、读出、组装和局部内核。分别报告实测、外推、资源扩展与同预算加速。保留原始计时和内存。
6. **验证迁移。** 先写明复用机制和适用条件，再冻结一个未参与开发的不同实现。对迁移目标建立自己的基线、合同和配对。若迁移只覆盖组装等组件，按组件范围报告。

## 按任务加载参考

- 空间转录组或 3D 重建读 [输出路线与 Fig. 3 合同](references/reconstruction.md)。
- 对照 SpatialZ 早期和 12 层结果、示范算法以及跨算法试验，读 [案例与证据范围](references/cases.md)。
- 使用配对工具前读 [运行清单与证据示例](references/benchmark.md)。

最终报告包含官方小试状态、实际优化入口、瓶颈证据、修改内容、等价覆盖、计时范围、每轮数据和复现命令。数据、日志、仓库及外部发布均沿用用户授权。

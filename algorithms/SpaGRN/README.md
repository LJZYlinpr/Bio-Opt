# SpaGRN 任务数据搬运

复用 [SpaGRN v1.1.4 源码](https://github.com/BGI-Qingdao/SpaGRN/tree/61690ca69a766f14c580de5dd47d4bc7b603813a) 的 Geary/Moran worker，输入矩阵与权重每个 worker 传一次，每个任务只传基因索引。原任务顺序保留。此目录附上游 GPLv3 许可。

在已有 SpaGRN 环境中的 Linux 节点运行 `python test_fast_equivalence.py`，比较作者原版与优化路径的 Geary、Moran 各 15 对，检查顺序与数值。官方 hi-res 入口为 `run_tutorial_hires.py`，令 `R_ROOT` 指向隔离工作根，把教程输入放在 `data/SpaGRN` 下，以 `SPAGRN_FAST_BIVARIATE=1` 启用优化。该入口也调整了 worker 配额。

hi-res 实验优化路径的 bivariate 阶段约 **5.5 分钟**。旧路径全量 **5.7 天**是吞吐外推，逐位比较覆盖小样本。

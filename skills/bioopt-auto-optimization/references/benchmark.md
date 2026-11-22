# 可执行配对与记录

在已分配的 Linux 计算环境中执行。准备两个产物独立写出的非交互命令和一个独立验证命令。路径相对当前目录解析。可从官方入口构建薄适配层，适配层与调用参数必须留档。

```json
{
  "case_id": "fixed-query-example",
  "claim_scope": "frozen readout, including setup and output writing",
  "equivalence_contract": "exact artifact fields, including ordered site and gene IDs",
  "inputs": ["frozen-input.npz"],
  "baseline": ["python", "run_reference.py", "frozen-input.npz", "{output}/result.npz"],
  "candidate": ["python", "run_candidate.py", "frozen-input.npz", "{output}/result.npz"],
  "validate": ["python", "/path/to/skill/scripts/compare_numpy.py", "{baseline}/result.npz", "{candidate}/result.npz"],
  "threads": 1,
  "timeout_seconds": 300,
  "warmups": 1,
  "repeats": 3
}
```

```bash
python /path/to/skill/scripts/paired_commands.py case.json --out evidence/pair-01
```

先用两条相同的基线命令运行一次基线自比较。正式配对要使用另一处新目录，并逐对使用独立的输出与验证。复现工具记录输入哈希，可能预热文件缓存。若需要测冷启动或预先加载的共享成本，说明冷热策略，单列实测共享成本与合成总用时。运行时峰值 RSS 为 Linux wait4 观测，设备峰值和并行进程聚合内存需要方法自己的测量。完整计时包含命令启动和输出写盘，验证运行单独计时。

NumPy 校验器默认比较字段、形状、dtype 和内容位模式，并拒绝非有限浮点与 pickle 数据。只有逐字段授权了容差时才用 `--float-key` 和 `--atol`/`--rtol`。使用 H5AD、BAM、排名表、稀疏矩阵或带静态掩膜的缺失值时，替换为独立的格式专用比较器并指定合同。比较器须在失败时返回非零。

## 保存的结果

- 官方小试来源、提交、输入哈希、命令、环境、patch、设备与资源。
- 测量边界、计入的准备工作、每对顺序和各轮原始值、计时中位数及变动范围、峰值内存。
- 已冻结的科学合同、独立验证的字段及差异、拒绝的候选与接受依据。
- 批量运行保留窗口 / 方法 / 路线记录。用明确完成收据支持续跑，失败保留原始错误。若资源不足只完成一次配对，注明次数和暂时性。
- 若需和普通 Codex 比较优化效果，两侧使用同样的起始提交、科学合同、算力、可用时间和候选试验预算。两侧候选由同一冻结验证器验收，比较通过后的最好结果。

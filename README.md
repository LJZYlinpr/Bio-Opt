# BioOpt

[English](#english) · [简体中文](#简体中文)

## English

**BioOpt: An autonomous performance optimization skill for bioinformatics algorithms, illustrated with spatial omics**

BioOpt guides an AI agent through reproducing a reference implementation, profiling real workloads, testing optimization candidates, and checking that the scientific outputs meet a predefined acceptance standard.

### Import into Mira

In Mira, choose **Skills → Add skill → Import from GitHub** and paste the URL of the directory containing `SKILL.md`. [Mira documentation](https://mirascience.app/help/skills)

```text
https://github.com/LJZYlinpr/Bio-Opt/tree/main/skills/bioopt-auto-optimization
```

You can also upload the packaged skill from the same screen. [Download the .skill file](https://github.com/LJZYlinpr/Bio-Opt/raw/refs/heads/main/skills/bioopt-auto-optimization.skill)

### Workflow and evidence

The skill adapts the Codex Scout, Baseline, Optimize, and Experiment workflows. It investigates the task, verifies the baseline, fixes the scientific acceptance criteria, proposes candidates from profiling, ranks them independently, and runs small trials before comparing leading candidates and combinations. Candidate routes cover CPU multiprocessing, GPU migration, memory and batching, and scaling from one GPU to multiple GPUs.

- SpatialZ manual optimization reduced the total core generation time for three MERFISH cases on an A100 from **157.229 to 39.943 seconds (3.94×)**. Coordinates, cell types, MENDER descriptors, and expression matched. [Expert trajectory](skills/bioopt-auto-optimization/references/expert-trajectories/README.md)
- A subsequent Agent optimization sped up the generation of 12 MERFISH target layers by **8.40× and 8.70×** in two paired runs on a V100. Coordinates, cell types, and expression matched with a fixed seed per layer. Final merge, disk output, and plotting were outside the timed interval. [Agent trajectory](skills/bioopt-auto-optimization/references/agent-trajectories/spatialz-batched-mender.md)
- The other cases cover data transfer in SpaGRN, CPU thread settings in UniST, and a contributed CPU/GPU equivalence audit. Each record states what was measured and what remains to be tested. [Agent trajectory library](skills/bioopt-auto-optimization/references/agent-trajectories/README.md)

### Repository guide

| Location | Contents |
|---|---|
| [Skill entry point](skills/bioopt-auto-optimization/SKILL.md) | Instructions, decision criteria, and expected deliverables for each stage. |
| [Workflow foundations](skills/bioopt-auto-optimization/references/workflow-base.md) and [benchmark guide](skills/bioopt-auto-optimization/references/benchmark.md) | Codex workflow rules, a paired benchmark template, and instructions for the included scripts. |
| [Expert trajectories](skills/bioopt-auto-optimization/references/expert-trajectories/) and [Agent trajectories](skills/bioopt-auto-optimization/references/agent-trajectories/) | Optimization decisions, competing explanations, failed ideas, and supporting measurements. |
| [SpatialZ](algorithms/SpatialZ/), [SpaGRN](algorithms/SpaGRN/), and [UniST](algorithms/UniST/) | Source and provenance, optimization code or benchmark drivers, and case-specific verification instructions. |
| [Showcase HTML](docs/index.html) and [interactive case study](docs/spatialz-story.html) | The Chinese-language project showcase and six-step SpatialZ walkthrough. Both work when opened locally in a browser. Repository and download links require network access. |
| [Project submission](SUBMISSION.md) | Project scope, demo links, and upstream resources. |

Open the deployed showcase in a browser. It requires no installation. [Live demo](https://ljzylinpr.github.io/Bio-Opt/)

## 简体中文

**BioOpt：面对生物信息学算法的自主计算性能优化 Skill——以空间组学为例**

BioOpt 帮助 AI Agent 复现原版算法、定位真实运行瓶颈、试验优化方案，并依据事先约定的标准验证科学输出。

### 导入 Mira

在 Mira 中选择 **技能 → 添加技能 → 从 GitHub 导入**，粘贴包含 `SKILL.md` 的目录地址。

```text
https://github.com/LJZYlinpr/Bio-Opt/tree/main/skills/bioopt-auto-optimization
```

也可在同一界面上传已打包的 Skill。[下载 .skill 文件](https://github.com/LJZYlinpr/Bio-Opt/raw/refs/heads/main/skills/bioopt-auto-optimization.skill)

### 方法与成果

BioOpt 借鉴 Codex 的任务核验、基线验收、候选筛选和实验复核流程。先查清任务、验证基线并约定验收标准，再根据剖析结果提出候选，经独立排序、小试和交叉复测筛选方案。候选覆盖 CPU 多进程、CPU 转 GPU、内存与批处理，以及单卡优化和多卡扩展。

- SpatialZ 人工优化在 A100 上将三组 MERFISH 核心生成的总用时从 **157.229 秒降到 39.943 秒（3.94 倍）**，坐标、细胞类型、MENDER 描述符和表达一致。[专家轨迹](skills/bioopt-auto-optimization/references/expert-trajectories/README.md)
- 后续 Agent 优化在 V100 上进行两轮 12 层 MERFISH 生成配对，分别提速 **8.40 倍和 8.70 倍**。固定每层随机种子后，坐标、细胞类型和表达一致。最终合并写盘和绘图另计。[Agent 轨迹](skills/bioopt-auto-optimization/references/agent-trajectories/spatialz-batched-mender.md)
- 其他案例包括 SpaGRN 数据搬运、UniST 线程配置和团队贡献的 CPU/GPU 等价性审查。每个案例列出实测范围和待验证事项。[Agent 轨迹库](skills/bioopt-auto-optimization/references/agent-trajectories/README.md)

### 仓库导航

| 入口 | 内容 |
|---|---|
| [主 Skill](skills/bioopt-auto-optimization/SKILL.md) | 每步怎么做、依据什么判断、需要保存哪些结果。 |
| [流程底座](skills/bioopt-auto-optimization/references/workflow-base.md)与[计时指南](skills/bioopt-auto-optimization/references/benchmark.md) | Codex 工作规则、配对计时模板与脚本用法。 |
| [专家轨迹库](skills/bioopt-auto-optimization/references/expert-trajectories/)与 [Agent 轨迹库](skills/bioopt-auto-optimization/references/agent-trajectories/) | 优化理由、候选取舍、失败案例与验收证据。 |
| [SpatialZ](algorithms/SpatialZ/)、[SpaGRN](algorithms/SpaGRN/)、[UniST](algorithms/UniST/) | 原版来源、优化实现或对照脚本，以及各案例的复现说明。 |
| [展示页 HTML](docs/index.html)与[交互案例](docs/spatialz-story.html) | 中文项目展示与 SpatialZ 六步讲解。均可用浏览器离线打开，仓库和下载链接需要联网。 |
| [提交资料](SUBMISSION.md) | 项目信息、演示入口和第三方资源。 |

在线展示用浏览器直接打开，无需安装。[BioOpt Demo](https://ljzylinpr.github.io/Bio-Opt/)

# BioOpt

BioOpt 帮助 AI Agent 复现生物信息学算法的作者小试，从真实执行路径定位瓶颈，在科学输出等价的约束下提出、验证并计量优化。

- [Mira Skill 源目录](skills/bioopt-auto-optimization/) · [下载 .skill](skills/bioopt-auto-optimization.skill)
- [在线 Demo](https://ljzylinpr.github.io/Bio-Opt/) · [离线 HTML](docs/index.html)
- [SpatialZ 原版与手工优化代码](algorithms/SpatialZ/)
- [KiloCode 的 SpaGRN 搬运优化](algorithms/SpaGRN/) · [UniST 线程对照](algorithms/UniST/)
- [提交表文案](SUBMISSION.md)

## 现场演示

克隆仓库或下载 ZIP，解压后用浏览器打开 `docs/index.html`。无需服务器或安装依赖。首页静态展示项目方法、优化逻辑与已记录成果，支持手机和桌面浏览。原六步交互演示保留在 `docs/spatialz-story.html`，支持方向键与全屏。页面正文离线可读，仓库与 Skill 下载链接需要网络。

## 方法

BioOpt 内置 Codex 的任务核验、基线验收、候选晋级和实验复核规则。在此基础上增加四条优化路线、独立评估、并行小试和交叉复测。[工作流底座](skills/bioopt-auto-optimization/references/workflow-base.md)记录每步的放行条件。两个轨迹库收录 SpatialZ 手工优化和 KiloCode 的实际探索。

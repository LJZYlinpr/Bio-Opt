# BioOpt

BioOpt 帮助 AI Agent 复现生物信息学算法的作者小试，从真实执行路径定位瓶颈，在科学输出等价的约束下提出、验证并计量优化。

- [Mira Skill 源目录](skills/bioopt-auto-optimization/) · [下载 .skill](skills/bioopt-auto-optimization.skill)
- [SpatialZ 手工优化轨迹 Demo](docs/index.html)
- [SpatialZ 原版与手工优化代码](algorithms/SpatialZ/)
- [KiloCode 的 SpaGRN 搬运优化](algorithms/SpaGRN/) · [UniST 线程对照](algorithms/UniST/)
- [提交表文案](SUBMISSION.md)

## 现场演示

克隆仓库或下载 ZIP，解压后用浏览器打开 `docs/index.html`。无需网络、服务器或安装依赖。页面支持方向键、全屏和六步互动。

## 方法

先复现作者小试并冻结输入和输出合同。剖析实际入口，选一个热点做最小改动，再用独立验证器检查输出和同资源配对。保留失败探索、计时边界和跨算法迁移证据。SpatialZ 案例展示候选余弦、C++ 采样、CPU/GPU 按层交接，以及后续批量 MENDER。每个实验目录记录独立适用范围。

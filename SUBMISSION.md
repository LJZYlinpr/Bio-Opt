# BioOpt 游园会提交信息

## 项目名称

BioOpt｜保持科学结果的算法性能优化 Skill

## 一句话简介

将作者小试复现、性能剖析、等价性验收与配对计时封装成 Mira 可调用的 Skill，并用 SpatialZ 的人工优化轨迹展示每次设计取舍。

## 目标问题

加速生物信息学算法时，怎样找到重复或无效的计算，证明改动保留了科学输出，并让别人能够复查与复现？

## 已实现功能

提供可导入 Skill、自动配对计时与独立输出检查，以及展示方法、优化逻辑与实测成果的静态网页，并保留可离线交互的 SpatialZ 六步演示。仓库包含 SpatialZ 作者原版、手工优化实现、后续批量 MENDER 的配对代码，以及 KiloCode 的 SpaGRN 数据搬运优化和 UniST 线程配置对照。SpaGRN 全量旧用时为外推，实际验证范围在代码目录列明。

## 核心技术方案

隔离复现环境与正式项目，冻结输入、科学参数、随机状态和输出轴。先按阶段剖析，再只改已识别的冗余工作，用独立验证器逐对检查输出。同资源比较时计入初始化与数据搬运，输出逐轮用时与适用范围。SpatialZ 先算候选余弦、用 C++ 采样并按层交接 GPU 和 CPU，后续通过批量 MENDER 分摊初始化。在作者 MERFISH 例中，固定每层种子的两轮 12 层生成配对为 8.40× 和 8.70×，逐层表达、XY 和细胞标签一致。

## Demo 链接

[在线 Demo](https://ljzylinpr.github.io/Bio-Opt/)。现场网络不稳定时可[下载 HTML](https://github.com/LJZYlinpr/Bio-Opt/raw/refs/heads/main/docs/index.html)，用浏览器离线打开。另有 [SpatialZ 六步交互讲解](https://ljzylinpr.github.io/Bio-Opt/spatialz-story.html)，支持方向键翻页。

## 代码仓库

[GitHub · LJZYlinpr/Bio-Opt](https://github.com/LJZYlinpr/Bio-Opt)

## 第三方资源

- [SpatialZ 作者代码与 MERFISH 教程](https://github.com/senlin-lin/SpatialZ/tree/e1b5ed01e933729149e08399febca08db16ced14)
- [MENDER](https://www.nature.com/articles/s41467-023-44367-9) · [POT](https://pythonot.github.io/) · [PyTorch](https://pytorch.org/) · [NumPy](https://numpy.org/) · [SciPy](https://scipy.org/) · [AnnData](https://anndata.readthedocs.io/)
- [SpaGRN](https://github.com/BGI-Qingdao/SpaGRN/tree/61690ca69a766f14c580de5dd47d4bc7b603813a) · [UniST](https://github.com/lanshui98/UniST/tree/a2fd35a70a6080caab5fdd60652242f85059b09a) · [UniST 作者插值教程](https://unist-tutorial.readthedocs.io/en/latest/interpolation.html)

# SpatialZ 手工优化

- `original/SpatialZ.py` 是 [作者原版](https://github.com/senlin-lin/SpatialZ/tree/e1b5ed01e933729149e08399febca08db16ced14)，保留 MIT 许可与源码 SHA-256。
- `optimized/cpu.py` 和 `sample.cpp` 是候选余弦及 C++／OpenMP 抽样。先在全体参考细胞的 XY 树查近邻，再按类型过滤。均匀随机数沿原查询／基因顺序预生成。
- `optimized/full_pair.py`、`layers.py`、`queue_pool.py`、`watch_volume.py` 保存全量 ABCA-1 分层流水线和调度实现。生产调度需要自行准备队列、实测片与计算资源。
- `reproduce_official.py` 比较历史 A100 三组核心计算。`benchmark_merfish_12.py` 是后续 BioOpt 批量 MENDER 对照。发布版调整输入／输出路径与失败退出码。

## 作者 MERFISH 对照

准备作者教程中的五个 `merfish_*_paste.h5ad` 输入。使用 `original/requirements.txt` 中的 Python 依赖（历史 Python 3.10），另按设备安装 CUDA 版 PyTorch。默认表达模式需要 MENDER。

```sh
cd algorithms/SpatialZ
g++ -O3 -fPIC -shared -fopenmp optimized/sample.cpp -o optimized/sample.so
python benchmark_merfish_12.py --input-dir /path/to/merfish --out /path/to/new-run
```

历史 A100 三组核心生成 **157.229 → 39.943 秒**，比较坐标、细胞类型、MENDER 和表达。后续 V100 两轮 **1804.282 → 214.889 秒**、**1742.078 → 200.323 秒**，比较 12 层坐标、标签和表达，固定每层种子 42。计时覆盖生成，最终合并写盘和绘图另计。[配对与三层流式收据](evidence/)

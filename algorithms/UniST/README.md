# UniST 插值的线程配置对照

在 [固定上游版本](https://github.com/lanshui98/UniST/tree/a2fd35a70a6080caab5fdd60652242f85059b09a) 和 TensorFlow CPU 环境下运行。按作者教程准备模型与 `example_interpolation` 下的实测输入。

```sh
python benchmark_threads.py --source /path/to/UniST --out /path/to/new-run
```

同一教程入口交替比较 2、4、8 线程设置，各两轮。2→8 线程的均值 **6.287 → 5.296 秒**，重建 TIFF 的 SHA-256 相同。该比较测量线程配置收益。[原始计时](evidence/threads_repeat.json)

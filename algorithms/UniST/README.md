# UniST interpolation thread benchmark

Use the pinned upstream version and a TensorFlow CPU environment. Prepare the model and measured inputs under `example_interpolation` as directed in the tutorial. [Pinned source](https://github.com/lanshui98/UniST/tree/a2fd35a70a6080caab5fdd60652242f85059b09a) · [Interpolation tutorial](https://unist-tutorial.readthedocs.io/en/latest/interpolation.html)

| File | Purpose |
|---|---|
| [benchmark_threads.py](benchmark_threads.py) | Runs the tutorial workload with interleaved CPU thread configurations and validates the final TIFF. |
| [evidence/threads_repeat.json](evidence/threads_repeat.json) | Per-run measurements and output hashes for the repeated comparison. |
| [UPSTREAM-LICENSE](UPSTREAM-LICENSE) | Upstream license. |

```sh
python benchmark_threads.py --source /path/to/UniST --out /path/to/new-run
```

The same tutorial workflow ran twice each with 2, 4, and 8 threads in interleaved order. Mean runtime fell from **6.287 to 5.296 seconds** from 2 to 8 threads, and the SHA-256 of the reconstructed TIFF matched across all runs. The measured gain applies to the change in CPU thread settings. [Experiment trajectory](../../skills/bioopt-auto-optimization/references/agent-trajectories/unist-threads.md)

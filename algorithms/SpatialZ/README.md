# SpatialZ optimization and reproduction

| File or directory | Purpose |
|---|---|
| [original/SpatialZ.py](original/SpatialZ.py) | Original SpatialZ implementation. See the [source manifest](original/SOURCE_MANIFEST.json) for the pinned commit and SHA-256 and the [MIT license](original/LICENSE) for terms. |
| [original/requirements.txt](original/requirements.txt) | Python dependencies for the original workflow. |
| [optimized/cpu.py](optimized/cpu.py) and [optimized/sample.cpp](optimized/sample.cpp) | Candidate-only cosine scores and sampling in C++ with OpenMP. The code queries the XY tree over all reference cells, then filters by type. Uniform random numbers are prepared in the original query and gene order. |
| [optimized/geometry.py](optimized/geometry.py), [optimized/full_pair.py](optimized/full_pair.py), and [optimized/layers.py](optimized/layers.py) | Geometry and layer-by-layer generation for the ABCA-1 workflow. The production run needs measured slices, a queue, and compute resources. |
| [optimized/queue_pool.py](optimized/queue_pool.py) and [optimized/watch_volume.py](optimized/watch_volume.py) | Queue dispatch and monitoring for full-volume reconstruction. |
| [reproduce_official.py](reproduce_official.py) and [verify_reproduction.py](verify_reproduction.py) | Three-case reference versus optimized benchmark and output verification. |
| [benchmark_merfish_12.py](benchmark_merfish_12.py) | Paired benchmark with batched MENDER for 12 MERFISH target layers. The release scripts expose the input and output paths and return a failure exit code when verification fails. |
| [evidence/](evidence/) | The A100 three-case benchmark, two V100 12-layer paired runs, and the three-layer streaming verification record. |

## Reproduce the MERFISH benchmark

Prepare the five `merfish_*_paste.h5ad` files from the [upstream tutorial](https://github.com/senlin-lin/SpatialZ/tree/e1b5ed01e933729149e08399febca08db16ced14). The historical environment used Python 3.10. Install the Python dependencies and a CUDA-enabled PyTorch build for your device. The default expression mode requires MENDER.

```sh
cd algorithms/SpatialZ
g++ -O3 -fPIC -shared -fopenmp optimized/sample.cpp -o optimized/sample.so
python benchmark_merfish_12.py --input-dir /path/to/merfish --out /path/to/new-run
```

The three A100 core generation cases took **157.229 seconds before and 39.943 seconds after optimization**. XY coordinates, cell types, MENDER descriptors, and expression matched. The two subsequent V100 pairs took **1804.282 to 214.889 seconds** and **1742.078 to 200.323 seconds**. They compared coordinates, labels, and expression across 12 target layers, with NumPy seed 42 reset before each layer. The V100 results cover generation, with final merge, writing, and plotting measured separately. [Benchmark records](evidence/)

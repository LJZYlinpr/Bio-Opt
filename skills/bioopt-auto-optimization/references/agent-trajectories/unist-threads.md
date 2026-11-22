# UniST verifying a thread configuration gain

The experiment reproduced the official interpolation tutorial. It used measured slices 1 and 5 to generate slices 2–4 and included image conversion, interpolation, and volume merging. The CPU comparison kept the interpolation model, inputs, and parameters fixed.

## Why start with thread count?

Thread configuration is reversible and can expose parallelism in existing CPU operations. An initial large difference could have come from genuine thread scaling, execution order, cold startup, or host load. Matching TIFF hashes checked correctness but could not distinguish these explanations for the timing difference.

1. **Initial run.** Sequential runs with 1 and 4 threads on the same node took 78.489 and 5.762 seconds, with matching final TIFF hashes. The single comparison could not isolate startup, load, and thread effects. The 13.62× estimate was rejected.
2. **Revised design.** Use the same inputs and three commands and run 2, 4, 8, 2, 4, and 8 threads in one CPU allocation. Save all times and output hashes. Interleaving avoids always assigning the same configuration to the first or last run. A thread effect should recur across repetitions. The revised design does not reestablish the single-thread reference, so the supported comparison starts with two threads.
3. **Verification.** The 2-, 4-, and 8-thread means were 6.287, 5.638, and 5.296 seconds. All six final TIFF files had the same SHA-256. Raising the setting from 2 to 8 threads gave about 1.19×, with a smaller incremental benefit from 4 to 8.

The declining marginal gain warranted stopping this trial. A deployment decision between four and eight threads would depend on the target latency, throughput, or core-hours. The result applies to the tested CPU environment and final volume. Tests at larger workloads should measure throughput, thread allocation, peak memory, and repeat-run variation. [Benchmark script and measurements](https://github.com/LJZYlinpr/Bio-Opt/tree/main/algorithms/UniST)

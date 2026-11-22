# SpaGRN multiprocessing data transfer

The experiment reproduced the high-resolution gene regulatory network tutorial. After preprocessing, the data contained 13,408 observations and 9,799 genes. The bivariate stage processed 2,088,160 TF–target pairs and showed intermittent worker activity.

## Why inspect the task payloads first?

There were three possible bottlenecks. A single pair might be expensive to score, tasks might be too fine-grained, or large matrices might be crossing process boundaries repeatedly. The code showed that each small task carried the same object of several hundred MB. Increasing worker count could intensify serialization and memory contention. Changing the scoring formulas would also introduce numerical risk. The team changed the data transfer path while retaining the original scoring functions and task order to isolate the effect of repeated data movement.

1. **Separate computation from transfer.** The original code included about 525 MB of expression data and spatial weights in per-pair task payloads. Measured original throughput was about 4.2 tasks per second. Keep the Geary and Moran workers and the pair order fixed.
2. **Implement the transfer change.** Pass expression, weights, and indices at worker initialization. Dispatch batches of gene index pairs and assemble the adjacency table in original order. Pass-once initialization can still create per-process copies. Check aggregate memory when selecting a worker count.
3. **Probe exact equivalence.** Use the same spatial weights and 15 candidate pairs for each of Geary and Moran. Check pair identity and order first, since parallel return order and index errors can yield plausible but misassigned values. Both scores matched bit for bit with maximum absolute error zero. This probe isolates the transport patch. Other worker settings used in the tutorial require their own verification.
4. **Run the tutorial.** The optimized tutorial produced 77 regulons, 70 associated with receptors. The bivariate stage took about 5.5 minutes and the tutorial took about 40 minutes. The original full run was unfinished. Its estimate of about 5.7 days extrapolates from measured throughput.

The throughput observation supports investigating data transfer first. The evidence establishes small-case equivalence and an optimized full run. A matched full-run speedup and network-wide equality require an original run under comparable conditions. Pair several workload sizes through the same entry point and measure aggregate process memory, since process startup modes can copy initialized data differently.

When many tasks carry the same large object, inspect serialization boundaries and throughput. Change transfer in isolation and test both memory use and scientific outputs. [Code and reproduction entry point](https://github.com/LJZYlinpr/Bio-Opt/tree/main/algorithms/SpaGRN)

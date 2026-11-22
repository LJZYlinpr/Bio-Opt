# Agent optimization trajectories

| Case | What it shows | Evidence |
|---|---|---|
| [SpatialZ batched MENDER](spatialz-batched-mender.md) | How a slower caching attempt redirected profiling and how a randomized baseline was repaired. | Stage profiles, rejected candidates, and paired trials on individual slice pairs and 12 target layers. |
| [SpaGRN data transfer](spagrn-transport.md) | How task payloads exposed a multiprocessing bottleneck and how the transfer change was tested in isolation. | Exact agreement on small cases, an optimized full tutorial run, and an extrapolated original full-run time. |
| [UniST thread settings](unist-threads.md) | How interleaving and repeating runs resolved an inflated single-run speedup. | Upstream tutorial and repeated thread comparisons. |
| [CPU/GPU equivalence audit](cpu-gpu-equivalence-audit/README.md) | How to partition a workflow and design verification before measured results are available. | Team-contributed audit workflow, static case study, and evaluation prompts. |

The first three cases come from KiloCode experiments. Each entry records what triggered the decision, which explanations competed, and what the smallest discriminating trial showed. Retrospective analysis is identified as such, and unexecuted alternatives are presented as proposals.

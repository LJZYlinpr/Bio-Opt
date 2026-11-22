# SpaGRN task data transfer optimization

The patch reuses the original Geary and Moran workers and preserves task order. Each worker receives the expression matrix and spatial weights during initialization. Individual tasks carry gene indices. [Pinned SpaGRN v1.1.4 source](https://github.com/BGI-Qingdao/SpaGRN/tree/61690ca69a766f14c580de5dd47d4bc7b603813a)

| File | Purpose |
|---|---|
| [spagrn_fast_bivariate.py](spagrn_fast_bivariate.py) | Data transfer and task batching patch for the bivariate statistics stage. |
| [test_fast_equivalence.py](test_fast_equivalence.py) | Checks order and exact scores for 15 Geary and 15 Moran candidate pairs. |
| [run_tutorial_hires.py](run_tutorial_hires.py) | Runs the upstream high-resolution tutorial and exposes worker settings and the optimization switch. |
| [LICENSE](LICENSE) | Upstream GPLv3 license. |

Run `python test_fast_equivalence.py` in an existing SpaGRN Linux environment. For the tutorial, point `R_ROOT` to an isolated working directory, place tutorial input under `data/SpaGRN`, and set `SPAGRN_FAST_BIVARIATE=1` to enable the patch.

The optimized bivariate stage took about **5.5 minutes** in the high-resolution tutorial. The original full workload estimate of **5.7 days** extrapolates from measured short-run throughput. Exact equivalence was tested on the small set of candidate pairs. [Decision and measurement record](../../skills/bioopt-auto-optimization/references/agent-trajectories/spagrn-transport.md)

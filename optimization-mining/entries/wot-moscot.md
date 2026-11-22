---
slug: wot-moscot
method: "Waddington-OT (WOT) -> moscot.time (scalable optimal transport)"
omicverse_binding: "not bound (no optimal-transport module in omicverse.single / omicverse.space was found; unverified for future releases)"
original:
  paper: "Schiebinger et al., Cell 2019, doi:10.1016/j.cell.2019.01.006"
  repo: "https://github.com/broadinstitute/wot"
  commit: "ca5e94f05699997b01cf5ae13383f9810f0613f6"
optimized:
  paper: "Klein et al., Nature 2025, doi:10.1038/s41586-024-08453-2; preprint doi:10.1101/2023.05.11.540374"
  repo: "https://github.com/theislab/moscot"
  commit: "d6e1bb10398cb9ad1e8577552439a6c571a4903a"
optimization_type: "system + approximation"
mechanism: "Re-express WOT's entropic W-type OT on OTT-JAX so the cost matrix is never materialized (Sinkhorn mat-vecs query the cost row by row), the iteration is JIT-compiled and run on GPU, and the coupling can be constrained to low rank for linear memory."
reported_speedup: "WOT ran out of memory above 100,000 cells even with 700 GiB of server RAM; moscot.time coupled 275,000 x 275,000 cells in <10 GiB and ~3 h (mouse embryogenesis E11.5/E12.5, moscot preprint Fig. 2b; WOT on CPU, moscot runtime on GPU); the full pipeline maps 1.7 M cells over 20 time points"
equivalence: approximate
optimization_families: [O3, O9, O10, O12]
sources:
  - "https://doi.org/10.1016/j.cell.2019.01.006"
  - "https://doi.org/10.1038/s41586-024-08453-2"
  - "https://www.biorxiv.org/content/10.1101/2023.05.11.540374v2"
  - "https://github.com/broadinstitute/wot/blob/ca5e94f05699997b01cf5ae13383f9810f0613f6/wot/ot/ot_model.py"
  - "https://github.com/broadinstitute/wot/blob/ca5e94f05699997b01cf5ae13383f9810f0613f6/wot/ot/optimal_transport.py"
  - "https://github.com/theislab/moscot/blob/d6e1bb10398cb9ad1e8577552439a6c571a4903a/src/moscot/backends/ott/solver.py"
  - "https://github.com/theislab/moscot/blob/d6e1bb10398cb9ad1e8577552439a6c571a4903a/src/moscot/problems/time/_lineage.py"
verified: true
---

## Baseline bottleneck

WOT implements Waddington's OT as an entropy-regularized ("W-type") transport between consecutive time points, with growth/death absorbed by marginal penalties. The released package (PyPI `wot` 1.0.8.post2) computes everything densely on the CPU:

1. `compute_default_cost_matrix` builds the full `n0 x n1` squared-Euclidean distance matrix and normalizes by its median (`wot/ot/ot_model.py:243-253`).
2. The default solver (`solver='duality_gap'`, `ot_model.py:88`) allocates **two** dense kernels per epsilon scaling and one dense transport matrix: `_K = np.exp(-C/epsilon_i)`, `K = np.exp((u[:,None] - C + v)/epsilon_i)`, then `R = (K.T * a).T * b`. Each Sinkhorn step is a dense `K.dot(...)` / `K.T.dot(...)` (`wot/ot/optimal_transport.py:100-166`).

Peak memory is therefore `Theta(n0*n1)` float64 per matrix: about 24 bytes per cell pair for C, `_K` and `K`, plus another dense `R` while the duality gap is checked. The same overflow guard that makes WOT require fairly large epsilon (`"Overflow encountered in duality gap computation"`, line 165) also limits how sharp the coupling can be. The moscot authors summarize the same limitation: WOT "remains limited to uni-modal gene expression data and does not scale to hundreds of thousands of cells per time point" and has "quadratic time- and memory complexity in the number of cells".

## Optimization trajectory

moscot.time keeps WOT's model — its docstring says the class is "for analyzing time-series single cell data based on :cite:`schiebinger:19`" (`src/moscot/problems/time/_lineage.py:26-29`) — and changes the execution:

1. **Never materialize the cost matrix ** — the problem is handed to OTT as a `PointCloud` geometry with an `SqEuclidean` cost function and an optional `batch_size`. Sinkhorn only ever needs `K v` and `K^T u`, so the cost is evaluated row by row and the `n0 x n1` matrix is never stored. The preprint states the reason and the consequence: "the cost function c can be queried on the fly for those cell-cell distances which are required to evaluate the current row of the matrix-vector product. Online evaluation reduces the memory complexity to be linear in cell number".
2. **JIT-compile the iteration ** — `OTTJaxSolver` defaults to `jit=True`, so the whole Sinkhorn loop (scaling, stabilization, convergence check) is compiled once instead of dispatching Python/NumPy per iteration; the preprint lists "Just-in-time compilation (jitting)" as the third engineering improvement.
3. **Move the same code to GPU ** — deleting the quadratic matrix is what makes GPU execution possible at all ("the quadratic memory complexity prevents this in practice... Online memory evaluation renders GPU acceleration possible"). No separate GPU implementation is maintained; the array code is device-agnostic.
4. **Low-rank coupling when cells reach the millions ** — passing `rank=r > 0` replaces `sinkhorn.Sinkhorn` with `sinkhorn_lr.LRSinkhorn`, so the coupling is stored as factors with rank `r << n` and both time and memory become linear:

```python
# moscot: src/moscot/backends/ott/solver.py:271-283
if rank > -1:
    kwargs.setdefault("gamma", 500)
    kwargs.setdefault("gamma_rescale", True)
    if eps is not None and eps > 0.0:
        logger.info(f"Found `epsilon`={eps}>0. We recommend setting `epsilon`=0 for the low-rank solver.")
    self._solver = sinkhorn_lr.LRSinkhorn(rank=rank, epsilon=epsilon, ...)
else:
    self._solver = sinkhorn.Sinkhorn(initializer=initializer, **kwargs)
```

The paper benchmarks rank 200 and reports linear memory at 275k cells; it also notes that beyond a few hundred thousand cells per time point the remaining quadratic *time* makes low-rank necessary.

5. **Materialize rows, not the coupling (again)** — downstream pull/push queries are answered from the dual potentials: for a full-rank solution `_rows_from_potentials` computes `exp((f_i + g_j - C_ij)/epsilon)` only for the requested rows, and the low-rank output materializes rows as `qg[ixs] @ r.T`. WOT instead returns the dense `R` and multiplies it. This keeps ancestor/descendant queries linear in the number of queried cells.

6. **Keep the chaining semantics** — consecutive time points are still coupled under the same Markov assumption and multiplied to link distant time points, so downstream pull/push operations and gene correlations are unchanged in meaning.

Relevant locations: `src/moscot/backends/ott/solver.py:95-135` (PointCloud + `batch_size`), `src/moscot/backends/ott/solver.py:234-283` (solver/rank dispatch), `src/moscot/problems/time/_lineage.py:26-60` (TemporalProblem), `src/moscot/backends/ott/output.py:20-105` (row materialization from potentials / low-rank factors).

## Evidence

- "While moscot.time computed a coupling for all 275,000 cells at both time points using less than 10 GiB memory in approximately 3 hours, WOT ran out of memory as soon as we exceeded 100,000 cells, even though we provided 700 GiB server memory." — https://www.biorxiv.org/content/10.1101/2023.05.11.540374v2 (Fig. 2b)
- Benchmark setup: E11.5 = 455,124 cells and E12.5 = 292,726 cells from the mouse embryogenesis atlas, 11 subsampled pairs in 25,000-cell steps up to 275,000 cells per time point; WOT was run on CPU "as it does not support GPU acceleration"; moscot used epsilon = 0.005 while WOT needed epsilon = 0.05 "as it run into numerical overflow otherwise" — same URL, Methods.
- "moscot.time ... inherits WOT's popular cell growth- and death rate modeling, is applicable to multimodal data, scales to millions of cells" — same URL.
- Published abstract: moscot reconstructs "developmental trajectories of 1.7 million cells from mouse embryos across 20 time points" — https://doi.org/10.1038/s41586-024-08453-2
- moscot README: "powered by OTT ... a JAX-based Optimal Transport toolkit that supports just-in-time compilation, GPU acceleration, automatic differentiation and linear memory complexity for OT problems" — https://github.com/theislab/moscot

## Equivalence and caveats

- **Same model, different numeric path**: both solve the entropy-regularized W-type OT with growth/death; moscot is a re-implementation on a different solver stack, so transport plans are close but not bit-identical (different epsilon scheduling, stabilized absorption, convergence criteria).
- **Low-rank is a genuine approximation**: at `rank=r` the coupling is constrained to a factorized set `U(a,b,r)`, budgeted row/column marginals can deviate, and the solution depends on `r`; the paper uses rank 200. Use low-rank only when the full-rank problem no longer fits, and check marginal error plus downstream transition scores.
- **Numerical precision**: JAX defaults to float32 unless x64 is enabled, whereas WOT casts the cost to float64. For small problems, enabling x64 or comparing transport plans on a common subsample is the safe check.
- **epsilon matters more than the solver**: the two tools were benchmarked at different epsilon (0.005 vs 0.05), which changes both runtime and plan sharpness. Any pairwise benchmark must report epsilon.
- **Block size is a knob you must set**: the linear-memory behaviour comes from evaluating the cost in blocks during the matrix-vector products; `batch_size` is `None` by default, so a memory-bounded run must choose a block size explicitly and measure peak RSS rather than assuming it.
- **Growth/death inputs**: WOT estimates growth rates from population sizes; moscot exposes them as marginal penalties (`a`, `b`) with `BirthDeathProblem`. Defaults are not identical, so the recovered growth rates can differ even when the coupling is close.
- The memory half of the comparison ("<10 GiB versus OOM at 100k") is more robust than the runtime half, because the runtime comparison pits GPU moscot against CPU WOT.

## Reuse plan for Bio-Opt

- ** is the transferable core**: identify any solver whose iteration touches a pairwise `n x m` object only through matrix-vector products, and replace the stored object by a row-wise evaluator (plus a block size). Verify against the dense implementation on a small case by comparing mat-vecs, not the matrix.
- **Then **: compile the whole iteration and move it to an accelerator; both only become possible after the memory change, which is the usual ordering mistake to avoid.
- **Then **: add a low-rank/structured factor only when time (not memory) remains the limit; freeze acceptance as marginal deviation below a stated bound and downstream transition scores within tolerance.
- Benchmark protocol for this family: same epsilon, same marginals, same cell subset, and report memory separately from time, since GPUs change one and not the other.

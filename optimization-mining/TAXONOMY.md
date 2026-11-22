# Optimization families (O1–O12)

Entries in this directory describe a **baseline → optimized** pair. To compare trajectories across tools, each entry is tagged with one or more families below. The tags describe the *mechanism* that removes work, independent of the tool.

| Code | Family | What it removes | Typical verification |
|---|---|---|---|
| **O1** | Analytic / closed-form reformulation | Iterative inference or numerical integration replaced by an algebraic identity or closed-form statistic | Prove the identity on small inputs; compare outputs exactly |
| **O2** | Convolution / FFT / binning | Per-pair or per-permutation loops replaced by convolution over discretized/binned quantities | Compare against the exact loop on small cases; check binning error |
| **O3** | Sketching / low-rank / random projection | High-dimensional matrices replaced by a sketch or low-rank factorization | Bound/measure reconstruction and downstream error; fix the seed |
| **O4** | Sparse / compressed data structures | Dense materialization replaced by CSR/COO or compressed representations | Compare nonzeros and results; check ordering |
| **O5** | Graph sparsification / linearization | Dense or attention-heavy graph operations replaced by sparse or linear message passing, or edges pruned | Compare cluster/embedding stability; report graph statistics |
| **O6** | Shared precomputation / amortization | Work repeated per item (setup, null model, eigenbasis, ranking) computed once and reused | Recompute per item on a small case and compare bitwise |
| **O7** | Online / minibatch / streaming | Full-data passes replaced by bounded-memory chunks or streaming updates | Compare against full-batch on small data; report memory |
| **O8** | Early stopping / adaptive iteration | Fixed large iteration/permutation budgets replaced by convergence or adaptive sampling | Report the stopping rule; compare distributions and stability |
| **O9** | Parallel / distributed execution | Serial work spread across threads, processes, or devices | Bitwise or tolerance; check determinism and RNG per shard |
| **O10** | Compiled kernel / JIT / vectorization | Interpreter or Python-loop overhead replaced by C/C++/Numba/JIT or array-level kernels | Bitwise where arithmetic is unchanged; watch dtype/reduction order |
| **O11** | Controlled approximation | Exact computation replaced by an explicitly bounded approximation (subsampling, quantization, surrogate model) | Declare the tolerance; compare against exact on small cases and on official metrics |
| **O12** | Data layout / IO elimination | Materialization, copies, or repeated serialization removed (matrix-free operators, zero-copy views) | Bitwise on outputs; measure peak memory and IO |

## Mapping to the Bio-Opt engineering patterns

Where an entry is directly transferable to the engineering patterns used in `skills/`:

| O family | Related Bio-Opt pattern |
|---|---|
| O6 | amortize expensive setup; reuse target-invariant state |
| O7 | streaming to bound memory |
| O9 | operator-level parallelism; sharding |
| O10 | compiled kernels; JIT/specialization |
| O4 / O5 | candidate restriction; sparse computation |
| O12 | eliminate redundant serialization/IO |

O1–O3, O8, and O11 are **algorithmic** families: they change what is computed or how it is estimated, so the equivalence tier is usually `approximate` or `changed-output`, not bitwise.

## Tagging rules

- Tag the mechanism that removes the dominant work, not every incidental change.
- A pair may have several families (e.g., analytic null + FFT + parallel kernels).
- If the optimized implementation changes the statistical model, tag O11 and set `equivalence: changed-output`.
- Numeric P-codes from earlier drafts are retired; use O-codes only.

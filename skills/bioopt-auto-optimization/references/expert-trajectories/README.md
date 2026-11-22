# Expert optimization trajectories

Use these cases to retrieve expert decisions and verification patterns while diagnosing a new task. Archive new run records with their origin and review status. Human approval is required for inclusion in the curated reference cases. [Trajectory format and archive rules](../trajectory-format.md)

## SpatialZ manual optimization

SpatialZ uses measured spatial coordinates and expression to generate cell positions, cell types, and expression at target depths. The author tutorial uses MERFISH. The production objective was to extend the same generation task to the full ABCA-1 series. The team first used small cases and fixed intermediate outputs to isolate changes, then tested the combined implementation through the generation entry point. [Source, implementation, and verification records](https://github.com/LJZYlinpr/Bio-Opt/tree/main/algorithms/SpatialZ)

The decisions below are reconstructed from code, experiments, and verification records. Retrospective passages explain the reasoning behind the choices.

### 1. Trace the data flow to find redundant work

Tracing donor weights backward through expression synthesis showed that SpatialZ builds a full cosine similarity matrix for all reference and generated cells. Sampling then reads scores for the eligible donors. The XY tree requests at most 50 neighbors, so most of the scores are never used.

The retrospective comparison considered accelerating the full matrix, computing it in blocks, or evaluating the scores that sampling reads. GPU execution and blocking can change elapsed time or peak memory while retaining the product of reference and query counts as the number of similarity evaluations. Evaluating candidate scores bounds that work by the query count times the neighbor cap. This is valid because the original XY query and type filter determine donor membership independently of cosine similarity. A consumer that needs global cosine ranking would invalidate the restriction.

### 2. Test the donor rule directly

The original code caps the requested neighbor count by the number of reference cells of that type. It queries the **XY tree over all reference cells**, then filters out mismatched types. A tree built separately for each type can give a distant donor to a generated cell whose original candidate set is empty. A later transfer pilot found this difference and rejected the reordering. The manual implementation preserves the global tree, the per-type request count, donor order, and empty-set handling, then evaluates the eligible cosine scores.

A minimal retrospective example uses a rare type A query, two nearby B reference cells, and an A reference cell farther away. The original path can return an empty set, while the type-specific tree returns the distant A cell. The changed donor set can alter expression with the same seed. Different empty-set handling can also affect subsequent random draws. Comparing donor lists exposes the divergence before a full reconstruction.

### 3. Address redundant arithmetic and interpreter overhead

**Candidate cosine scores.** Normalize the reference and query descriptors once. Evaluate only scores for selected donors. Use the original weighting parameters to construct probabilities and cumulative distributions, preserving index order.

**Compiled sampling.** The original samples once per gene through a Python call. The optimized implementation prepares MT19937 draws in the original query and gene order, then selects donors from cumulative probabilities in C++ with OpenMP parallelism across generated cells. Each gene takes the selected donor's original value. Empty candidate rows and random draw assignment follow the original rules. Preparing the draws in stable order makes their assignment independent of worker completion order. Verify donor lists, probabilities, and draw assignment before comparing expression.

Restricting cosine evaluations saves arithmetic. Compiled sampling removes interpreter calls. The measured change from 96.107 to 0.139 seconds includes both changes. Their individual contributions require an ablation.

### 4. Expand verification and timing in stages

| Case | Controls and output checks | Measured result |
|---|---|---|
| ABCA-1 expression synthesis at 40 µm | Fix generated XY, cell types, the two reference slices, MENDER descriptors, and random seed. Compare all elements of the 8,756 × 1,122 expression matrix. | 96.107 to 0.139 seconds. Zero expression differences. This measures both changes within expression synthesis. |
| Three MERFISH core generation cases | Follow the author's geometry and labeling rules. Compare XY, types, MENDER descriptors, and expression per case. Use tutorial settings Beta=1, 80 random projections, and 1,000 iterations. | 157.229 to 39.943 seconds total on an A100. All 2,612,370 expression values and the associated outputs matched. |
| ABCA-1 small-window parallel comparison | Compare MENDER descriptors with 4, 8, and the original 200 workers. Compare serial and parallel sampling across three layers. | Descriptors and expression matched exactly in the selected cases. Worker settings can be chosen for available resources and timed separately. |

After expression matched with fixed descriptors, testing expanded to the author's entry point including geometry and MENDER. A failure with fixed descriptors directs debugging to the donors, probabilities, and sampling. A failure confined to the full entry point directs it upstream to geometry, labeling, or descriptors.

### 5. Change the execution structure for a full reconstruction

After reducing per-layer time, the production objective was to generate many layers within the available resources. Expression synthesis needs the XY coordinates of its own layer, so it can overlap geometry work for a later layer. The design uses layers for handoff and persistence and adjacent measured slice pairs for queue assignment. This overlap depends on the absence of cross-layer dependencies in the generation path.

- **GPU geometry and CPU expression.** Compute XY with the original objective and iteration count. Send each layer to the CPU for labeling, MENDER, and expression synthesis. While the CPU writes one layer, the GPU can compute the next. Persisting and resuming per layer limits rework to unfinished layers.
- **Queue adjacent measured depths.** Dispatch 143 adjacent depth pairs through a shared queue. Stable task IDs and per-layer seeds keep scheduling order from changing sampling. Preserve claim, throttling, and retry records.
- **Write and release per layer.** Save geometry and sparse expression and release descriptors, host buffers, and device caches. A three-layer trial checked sparse expression data and indices, XY, labels, and the gene axis. The first layer's expression was written before geometry for the next layer finished. For the 560 µm layer, CPU RSS dropped from 916.75 to 769.55 MiB after release. GPU allocation peaked at 330.92 MiB and ended at 16.25 MiB. These are within-run memory observations.

### 6. Identify the implementation used for Fig. 3

The full ABCA-1 reconstruction used the manual implementation above. It produced 22,574,061 points across 1,133 depths, with 144 measured depths and 989 generated depths. The three five-slice Fig. 3 windows use the original expression implementation with MENDER workers reduced from 200 to 4. The full production implementation calls MENDER per layer. The later batching experiments are in the Agent library. The production records establish scale and execution behavior, and the timing claims apply to the paired scopes in the table.

**Transferable lesson.** Trace from downstream consumers to the values they require. Use edge cases such as rare types to test semantics. Keep small intermediate outputs to localize differences, then verify handoff, writing, memory, and recovery through the real entry point. [Six-step interactive walkthrough](https://ljzylinpr.github.io/Bio-Opt/spatialz-story.html)

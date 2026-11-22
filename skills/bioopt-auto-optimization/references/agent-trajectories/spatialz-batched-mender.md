# SpatialZ from split caching to batched MENDER

## Baseline and observations

The starting point used the manual candidate-only cosine scores and compiled C++ sampling. Profiling three single-layer generations from the MERFISH tutorial on a V100 gave about 0.24–0.29 seconds for sampling, 11.35–12.95 seconds for geometry, and 20.17–27.02 seconds for MENDER. MENDER accounted for 60%–67% of the measured time. Increasing MENDER workers from 32 to 64 gave no consistent improvement. The third case slowed from 30.48 to 32.58 seconds.

The profile made further sampling work a lower priority. The team considered two explanations for MENDER's cost. It might be repeatedly computing reference descriptors, or repeatedly preparing the model. The explanations called for different changes.

## Experiments that distinguish the explanations

The first test checked whether reference and query descriptors could be computed separately while matching the joint calculation bit for bit. The timing test then compared three original model initializations with four in the split design, one for the shared reference and one for each query. If reference computation dominated, splitting should save time. If preparation was substantial, the additional initialization could make it slower. The outputs matched and runtime increased. That redirected work from caching to reducing initializations. The test supports a preparation-cost explanation without directly measuring initialization time in isolation.

Batching can change descriptors or row correspondence. The team checked three target layers for one reference slice pair and repeated the check for a second pair before integrating the change into full generation.

| Candidate | Smallest comparison | Outcome |
|---|---|---|
| Separate reference and query descriptors | Check whether separating the slices changes descriptors, then time three layers. | Descriptors and three-layer outputs matched. Runtime increased from 46.595 to 61.096 seconds. The split design was rejected. |
| Group layers into one MENDER call | Run reference slices A and B and three target layers in one MENDER model. Sample expression per layer. | Two slice pairs took 48.495 to 18.956 seconds and 53.440 to 20.395 seconds. XY, labels, and expression matched bit for bit for all three layers in each test. |
| Reduce or reuse random projections in geometry | Inspect the generator calls in the 1,000-iteration loop. | Each calculation advances the same generator. Projection reuse or fewer iterations would change the generation task. |

## Rejected 12-layer comparison

An early comparison of the batch implementation against saved baseline outputs appeared to show a 4.54× gain. XY and types matched, but roughly 440,000–480,000 expression entries differed per layer. It failed verification. To check whether the differences came from geometry, descriptors, or randomness, the team compared a fresh original run, the saved original outputs, and the batch implementation in the same process and on the same GPU. The two original outputs also differed. The original expression sampler used an uncontrolled global NumPy state, so the old outputs could not serve as an exact comparison. A fresh paired run under the same random-state rule was necessary to judge the candidate.

## Corrected comparison and repetition

Reset NumPy seed 42 before every original call and use the same seed per target layer in the candidate. A retest with three targets from one reference pair passed layer-by-layer verification and took 436.059 to 52.223 seconds. The comparison then expanded to four reference pairs and 12 targets.

| Generation pair | Original | Optimized | Verification |
|---|---:|---:|---|
| Round 1 | 1804.282 s | 214.889 s | Expression, XY, and types matched across 12 layers. 8.40×. |
| Round 2 | 1742.078 s | 200.323 s | All 12 layers passed the same checks. 8.70×. |

These timings cover generation. Final merge, writing, and plotting are outside the comparison. The result applies to runs that reset the seed for each layer. The continuous random stream used by the original workflow remains a separate validation setting. The three-target pilot checked the call boundary. The two full runs also covered four reference intervals and variability over longer runs. Further tests of batch size should record peak memory and verify the outputs of each target. [Benchmark code and evidence](https://github.com/LJZYlinpr/Bio-Opt/tree/main/algorithms/SpatialZ)

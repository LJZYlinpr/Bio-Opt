# CPU/GPU equivalence audit

A team member contributed the [original audit package](cpu-gpu-equivalence-audit.zip). It contains a skill, equivalence criteria, a benchmark template, SpaceTracer static partition notes, and three evaluation prompts with expected response characteristics. The case documents an audit method and experiment design.

## Why partition the workflow before migrating it?

A paper can suggest that some mathematics is suitable for GPU execution while leaving actual call boundaries, format conversions, and data volume unknown. If dynamic preparation dominates, moving a small operation can add enough synchronization and transfer to increase total time. Mapping inputs, outputs, and ordering dependencies exposes two testable options. Reduce device transfers, or migrate a numerical hotspot with sufficient workload. Proceed to experiments once the implementation and baseline can be inspected.

## Turn those observations into verification

1. **Locate device boundaries.** For each stage, record inputs, outputs, device, dtypes, randomness, dependence on upstream ordering, and transfer volume. Include CPU preparation, GPU operations, synchronization, host transfers, and writing in the same table. When only the paper is available, use static analysis to propose candidates and identify missing code, data, or environments.
2. **Order candidates by risk.** Look for avoidable transfers and reusable intermediates, then compare batching on the same backend and evaluate migrating operations across devices. Handle changes to precision, sampling, iteration budget, or scientific outputs as separate variants.
3. **Define equivalence.** Pick a level from the table and freeze acceptance criteria before running. Find the first difference in critical intermediates, then verify persisted outputs. Check random inputs and draw order even when the seed matches.
4. **Use two stages of testing.** Pilot with fixed inputs and sensitive intermediates. If outputs pass, expand the workload and record end-to-end time, initialization, transfers, synchronization, host memory, and peak device memory. Estimate variability from repeated runs. Set the speed criterion in advance based on task objectives and expected measurement variation.
5. **Record decisions.** Identify each change by ID, file, stage, original and candidate devices, evidence path, and rollback reason. Mark feasibility, pending tests, failures, and the equivalence level each passing candidate reached. Record correctness separately from speed.

| Equivalence level | Sufficient evidence |
|---|---|
| File bytes or in-memory bit patterns | Specify the object, compare the bytes or bit patterns, and record metadata and shapes. File and memory checks can give different results. |
| Numerical equivalence | Meet predefined absolute and relative tolerances and report maximum error and critical structure checks. |
| Downstream equivalence | Pass the agreed thresholds on stable task metrics and state the downstream scope. |
| Unverified | Static arguments, speed, or visual inspection require an independent comparison. |

## A proposed experiment for device migration

Compare three paths with fixed upstream intermediates and random decisions. Run the original, a transfer-only change, and that change combined with the GPU operation. If transfer-only changes already help, keep that lower-risk candidate in contention. If operation time falls and total time rises, measure preparation, host transfers, and synchronization. If the operation is the first point of numerical divergence, judge the difference against the predefined criteria. This retrospective design can distinguish device-level and application-level gains. It remains a proposed experiment.

## SpaceTracer static case

The contributor partitioned the method described in the paper into data preparation, dynamic groups and rules, fixed-shape numerical work, and downstream aggregation. Initial GPU candidates are batch statistics or model operations with fixed inputs. Preserve input parsing, dynamic grouping, and external program calls during the first comparison. Save upstream groups, samples, and fixed inputs before migration, and compare both boundary outputs and final artifacts. Test reduction, ordering, and random draw behavior across implementations explicitly.

The package documents source checks and pending work. It contains a static partition study. A rerunnable SpaceTracer migration patch, paired device timings, and passing equivalence results remain to be produced.

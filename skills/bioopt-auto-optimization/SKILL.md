---
name: bioopt-auto-optimization
description: Optimize bioinformatics software with a verified baseline, frozen scientific acceptance criteria, ranked performance candidates, and reproducible measurements. Use to investigate bottlenecks, design or execute optimization experiments, compare resources, or transfer a validated optimization.
---

# BioOpt scientific computing optimization

Given an algorithm and a resource budget, use reproducible experiments to find the best measured implementation that preserves the scientific task. This skill follows the Scout, Baseline, Optimize, and Experiment rules in the [Codex workflow foundations](references/workflow-base.md). The six steps produce a task brief, a baseline record, acceptance criteria, a candidate table, an independent ranking, and experimental results. Flag any evidence gaps.

## 1. Identify the actual task (Scout)

- Reconstruct the task from the request, local configuration, commands, logs, and existing outputs. Record the execution entry point, code revision, inputs, workload, target outputs, and resource limits.
- Check the upstream GitHub repository, including versions, examples, issues, and discussions. Read the paper and supplements, the authors' tutorials, and community reproduction reports. Use community reports to identify questions, then verify them against the relevant code and run records.
- Trace the public entry point to the final outputs. Distinguish upstream defaults, local adaptations, and previously tested alternatives. Look for evidence that could overturn decisions about the task or baseline, such as an identically named command that ran on different data or solved a different task.
- **Deliver the task brief** with the scientific question, actual entry point, sources, available baseline, and missing evidence. End scouting once the evidence is sufficient to choose a baseline path.

## 2. Establish a credible baseline (Baseline)

- Check for an executable author example, tutorial, reference output, or measured baseline. Verify that its code, data, model files, and dependencies are available.
- If an author baseline exists, first look for successful local records using the same task, inputs, and revision. Check the actual command, output, and protocol. If only screenshots, exit codes, or the tail of a log are available, inspect the artifacts or rerun the necessary steps. Reuse sufficiently verified records.
- If the author baseline is absent or unusable, report which entry point, dependency, data, or output is missing. Record any verified substitute and the claims it can support. Verify basic execution, fair comparability, and faithful reproduction of upstream results as separate milestones.
- **Deliver the baseline record** with status, revision, input provenance, command, environment, actual outputs, and known deviations. Start formal paired measurements once the baseline is comparable.

## 3. Define scientific, stability, and performance criteria (Baseline → Experiment)

- State the biological question and the outputs and downstream decisions that must be preserved. Derive essential metrics from the paper and baseline. Choose exact agreement, predefined numerical tolerances, or downstream equivalence according to the task.
- Before inspecting candidate results, fix the comparison inputs, evaluation entry point, metric directions, acceptance thresholds, and units of analysis. Use evidence about variation, randomness, and backend differences to decide whether to measure baseline repeatability. Identify the required output fields and how to locate disagreements. Apply the same criteria to every candidate.
- Choose a primary performance metric such as elapsed time, peak resource use, cost per task, or throughput. Record supporting metrics and hard constraints. Specify whether the interval includes setup, input loading, data transfer, queueing, and writing.
- **Deliver the acceptance criteria** with the rationale, measurement method, and threshold for each metric. Require correctness before ranking by performance.

## 4. Generate candidates from the baseline profile (Optimize)

- Review existing attempts and profile the normal execution path. Measure stage time, peak resource use, device utilization, repeated computation, data movement, and waiting. Tie observations to functions, logs, or the profiler. Label measurements, inferences, and open questions.
- Propose at least **three candidates with distinct mechanisms per route**. Record why a candidate is inapplicable or what would make it feasible. Keep it in the table for evaluation. Add other evidence-backed routes when appropriate.
- For each candidate, record its ID and route, bottleneck evidence, mechanism, smallest useful change, required invariants, risks and resources, smallest counterexample, and criteria for proceeding. Explain the evidence and timing boundary for any performance forecast.

| Route | Candidate mechanisms and checks |
|---|---|
| CPU multiprocessing | Partition independent work, group fine-grained tasks to reduce dispatch overhead, and share read-only data across tasks in each worker. Check serialization, load balance, and oversubscription from worker counts multiplied by library threads. |
| CPU to GPU | Migrate the dominant operation, batch or fuse small operations, and keep data on-device to avoid transfers. Map the actual CPU, GPU, and I/O boundaries. Include initialization, transfer, and synchronization in performance measurements and check numerical and scientific equivalence separately. |
| Memory use and batching | Share initialization across related workloads, split or stream work, and reuse intermediate buffers or remove duplicate copies. Check dependencies across chunks and merge order. Measure workflow time and peak resource use. |
| Single GPU and multi-GPU scaling | Reduce peak memory on one GPU, batch independent work, and improve device utilization and the CPU handoff. If one GPU cannot meet the workload and resource targets, record the evidence and evaluate sharding, communication, and merging across GPUs. |

## 5. Rank candidates independently (Optimize)

- Give an independent evaluator the task brief, baseline, scientific acceptance criteria, profiles, budget, and complete candidate table. Have it check whether the work is avoidable, whether the change preserves semantics, and whether the proposed test can falsify the hypothesis.
- Exclude candidates that violate scientific or hard resource constraints. Rank the rest by bottleneck evidence, expected benefit, risk, implementation effort, and verification cost. Record the reasoning for each direction. Use a low-cost experiment to distinguish explanations for implausibly large gains.
- **Deliver the ranking** with IDs, ranks, rationale, resource limits, smallest trials, and criteria for larger experiments. Start implementation with a small set. If delegation is unavailable, write a separate assessment before modifying code.

## 6. Run small trials, compare candidates, and select the best measured result (Experiment → Optimize)

1. **Fix the run plan.** Isolate code changes and output directories by candidate. Use separate agents or worker processes when tools and resources permit. Give each worker a candidate ID, change boundary, reference command, acceptance criteria, resources, and trial budget. Keep formal timings free of competing candidate workloads or verify and document resource isolation.
2. **Check correctness first.** Start with a small case capable of exposing differences. Compare key intermediate results and final outputs, and locate the first divergence. Update or retry only when there is new evidence, a corrected configuration, or a testable hypothesis. Preserve rejected results and remove failed candidates from the final performance comparison.
3. **Pair under matched conditions.** Control warm and cold states and interleave run order. Save individual times, variability, and peak resource use. When device or resource counts change, report allocations, end-to-end results, and resource-normalized results. [Benchmark template and scripts](references/benchmark.md)
4. **Compare across conditions.** Retest leading candidates on common representative inputs and on a size or boundary case independent of the development case. For combinations, compare the reference, standalone changes, and combined implementation. Check incremental gains, interference, and equivalence. Use the evidence to continue exploration, focus verification, combine changes, fix a problem, or stop.
5. **Report the result.** Select the strongest measured candidate that satisfies the frozen scientific and resource criteria. Attach the baseline, patch, commands, raw measurements, verification, rejection reasons, and scope. For each decision, record the preceding observation, competing explanations, why the test was chosen, what it ruled out, and how it changed the next step. Identify historical evidence and retrospective reasoning. If gains are unstable, retain the credible baseline and experiment results. Validate transfer independently on another implementation.

## Case libraries

- The [expert trajectory library](references/expert-trajectories/README.md) shows how an expert located bottlenecks, preserved scientific semantics, and scaled the workload.
- The [Agent trajectory library](references/agent-trajectories/README.md) shows how candidate trials, failed baselines, and transfer checks changed optimization decisions.

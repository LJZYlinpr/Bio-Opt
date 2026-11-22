---
name: bioopt-auto-optimization
description: Optimize an existing bioinformatics implementation after its local reproduction has run. Verify the baseline, profile bottlenecks, use expert and Agent trajectories to propose changes, and compare independently reviewed pilot results, combined candidates, and full-workload performance under fixed scientific criteria.
---

# BioOpt performance optimization

Start with a locally reproduced implementation, normally handed off by Mira. Ask for the source revision, command, inputs, environment, logs, outputs, and resource budget as needed. The successful run is a proposed baseline. Verify the scientific task and execution before using it. This workflow adapts the [Codex Scout, Baseline, Optimize, and Experiment rules](references/workflow-base.md).

## 1. Investigate the algorithm and task

- Use a strong reasoning model for the research discussion when available. Check the upstream GitHub repository, pinned revisions, paper and supplements, tutorial, issues, and community optimization reports.
- Explain the biological question, what the algorithm actually computes, which pipeline variant the user wants, and which outputs are required. Trace important claims back to the relevant revision, protocol, or run evidence.
- Record a short task brief with actual entry point, workload, references, resource constraints, and any uncertainty that could change the baseline or candidate selection.

## 2. Verify the local reproduction and define acceptance

- Compare the local revision, command, input provenance, scientific parameters, environment, and saved outputs with the task brief and upstream protocol. A successful exit confirms execution. Check that it ran the intended algorithm and workload, including the output fields needed by the scientific task.
- Reuse the existing run when its identity and outputs are supported. If a requirement is uncertain, inspect the artifact or rerun the smallest relevant step. Report and resolve baseline mismatches or missing evidence before claiming a measured gain.
- Freeze the comparison inputs, scientific invariants, evaluation entry point, required metrics, tolerances, and analysis units before testing candidates. Choose exact, numerical, or downstream equivalence based on the task. Measure baseline repeatability when uncontrolled randomness or backend variation could affect the decision.
- Choose a primary performance goal and hard constraints. Declare the timing boundary and how initialization, data movement, waiting, and writing are accounted for. Save the baseline and acceptance records. [Paired benchmark template and scripts](references/benchmark.md)

## 3. Diagnose bottlenecks and retrieve relevant trajectories

- Follow the actual call path and inspect the code and stage outputs. Measure the baseline where necessary. Select or set up monitors within the available environment for CPU and GPU utilization, elapsed time, host and device peak memory, I/O, and transfers or synchronization when relevant. Tie each suspected bottleneck to a function, profiler measurement, or log and label inferences.
- For each bottleneck, search the [expert trajectory library](references/expert-trajectories/README.md) for reviewed reasoning and verification patterns. Search the [Agent trajectory library](references/agent-trajectories/README.md) for candidate breadth, negative results, and empirical comparisons. Retrieve by bottleneck mechanism, execution boundary, and scientific invariants. Compare each case's assumptions and test scope with the current task before using its ideas.
- Form distinct candidate mechanisms grounded in local evidence. Use at least three per applicable route when exploring the candidate space. Document inapplicable routes and keep the number of tests within budget. For each candidate, record its ID, observation, technical hypothesis, supporting trajectory if any, smallest useful change, invariants, risk, resources, and a test that could disprove the hypothesis.

| Route | Candidate mechanisms and checks |
|---|---|
| CPU multiprocessing | Partition independent work, batch fine-grained tasks, or reuse read-only data inside workers. Check serialization, load balance, and library thread oversubscription. |
| CPU to GPU | Migrate the dominant operation, batch operations, or keep data resident on the device. Check numerical behavior, transfer, and synchronization. |
| Memory use and batching | Share repeated initialization, stream or chunk work, or reuse buffers and remove excess copies. Check dependencies, merge order, and peak memory. |
| Single GPU and multi-GPU scaling | Reduce memory peaks, batch independent work, and improve the CPU handoff on one GPU. Consider multi-GPU sharding and communication after establishing that one GPU cannot meet the workload or resource target. |

## 4. Run isolated pilots and explore a different angle

- Assign candidates to parallel subagents or separate workers when available. Give each one a candidate ID, source revision, change boundary, isolated outputs, reference command, acceptance criteria, and time and resource budget. Isolate formal timing from other pilots.
- Test critical intermediates and final artifacts on representative small inputs and edge cases against the same baseline. Capture the patch, commands, raw paired timings, variability, peak resources, correctness, and the first difference where applicable. Preserve failed and inconclusive results with their reasons.
- Give a separate exploratory agent the first-round results, including failures, conflicting observations, and diminishing returns. Ask it to challenge the diagnosis or propose a different mechanism. It can look for avoidable work, a different batching boundary, or another cause of idle resources. Admit new ideas through the same bounded pilot protocol. If delegation is unavailable, perform a separate written pass before choosing follow-up trials.

## 5. Independently review and filter the measured candidates

- Give an independent reviewer the baseline, frozen criteria, profile, and raw evidence from every pilot. Check reproducibility, output equivalence, resource limits, measured gain, and confounders. Apply correctness and hard resource constraints as gates. Defer unsupported results and reject candidates that violate the criteria.
- Rank the eligible candidates using explicit weights appropriate to the task for measured benefit, stability, resource cost, and implementation risk. Record the basis for the weights, each candidate's evidence and decision, and plausible interactions among the passing candidates. A promising hypothesis or estimated speedup alone is insufficient for promotion.

## 6. Run a fusion round with the passing candidates

- Take a small set of reviewer-approved candidates with evidence of complementarity. Write down why they may work together and where they could interfere. Check overlapping code paths, stage and device boundaries, dependencies, peak memory, and new synchronization or transfers.
- In a separate branch or workspace, first combine the approved changes in isolation. Use a representative intermediate workload where available. Compare the **original baseline, A, B, and A+B** under the frozen acceptance criteria and matched conditions. For more than two changes, add them incrementally and compare each addition against the strongest accepted configuration so far.
- If the initial combination exposes an interaction, make a bounded improvement to the handoff or shared resource use. Give the revised combination a new ID and rerun its paired checks. Keep the first combination as evidence for attributing the improvement.
- Check scientific outputs for the combination. Measure the extra benefit over the **best standalone candidate**, including initialization, transfers, peak resources, and variability within the declared scope. For changing resource configurations, report both overall and resource-normalized results. Check the evidence independently before promotion.
- Send a combination to full testing only when it clears correctness and resource gates and shows a repeatable, meaningful gain over the best standalone option. Otherwise promote the strongest passing standalone candidate. Record a skipped fusion round if there are fewer than two passing candidates or no compatible pair.

## 7. Validate at full scale and archive the trajectory

- Run the original baseline and the strongest eligible implementation at the target workload within the authorized budget. Repeat paired comparisons as needed to judge stability. Report correctness, timings, memory, resource configuration, and the precise scope. If full-scale comparison cannot be completed, state the highest verified scale and keep a known-correct fallback.
- Archive all attempts, including rejected and inconclusive trials, with the released patch or retained baseline. Record observations, competing explanations, the reason each test was chosen, its outcomes, and what changed next. Link the full-scale evidence and distinguish measured results from hypotheses or estimates.
- Prepare a run record for the expert library and link the raw Agent trajectory. Record its source and review status. Human expert review qualifies a run to join the curated reference cases. Use the [trajectory format and archive rules](references/trajectory-format.md).

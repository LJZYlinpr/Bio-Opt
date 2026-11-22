# Optimization trajectory format

Use the same record shape for human and Agent trajectories so BioOpt can retrieve relevant evidence while profiling a new task. Index entries by algorithm or task, bottleneck mechanism, execution boundary, optimization route, evidence scope, and review status.

## What belongs in each library

- **Expert trajectories.** Build an initial set of about five human-led cases from established algorithms or recurring bottlenecks. Each case should show the original implementation, the improved implementation, the reasoning that connects them, and verified outcomes. Treat those reviewed cases as reference examples. The initial set is a curation target.
- **Agent trajectories.** Scale coverage across a broader set of algorithms and prompts, with roughly 100 algorithms as a possible campaign target. Preserve prompt or configuration identifiers, code revisions, failed and successful trials, and the strength of the evidence. Label any literature-only, extrapolated, and measured results.

## Run record

| Field | Evidence to keep |
|---|---|
| Task and provenance | Algorithm, scientific task, upstream source and paper, baseline revision, local reproduction command, inputs, outputs, environment, and actual resource limits. |
| Bottleneck and references | Stage and code path, measurements, competing explanations, and relevant expert and Agent entries, with a note on what transfers and under which assumptions. |
| Candidate decisions | Candidate ID, mechanism, why it was tried, original and modified code or commit links, results, first divergence if any, and the decision that followed. Include exploratory ideas and negative results. |
| Independent review | Frozen correctness criteria, eligibility and reason per candidate, the weighted ranking, and eligible candidates for combinations. Label automated review and human expert review separately. |
| Fusion | Compatibility hypothesis, baseline and standalone and combined run IDs, incremental effect versus the best standalone candidate, correctness, resource use, and promotion decision. Record the reason if skipped. |
| Full-workload result | Versions, commands, raw run records, output validation, timing scope, performance and resource comparisons, exact tested scale, limitations, and retained implementation. |
| Status and routing | Origin, evidence level, review status, tags, and the next suitable reuse or experiment. |

Create the complete run record in the project workspace, for example under `trajectories/runs/<task-id>.md`. Keep raw Agent trials and evidence together and link them from the run record. When archive writes are authorized, add a run entry under the expert library and cross-link the Agent records. If publication is unavailable, deliver the local record and the proposed index entry. Mark the entry as Agent-generated and pending expert review until a human approves it for the curated reference set.

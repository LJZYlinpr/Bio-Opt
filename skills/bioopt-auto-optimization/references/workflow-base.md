# Codex workflow foundations

This reference adapts reusable rules from the local Codex skills `deepscientist-scout`, `deepscientist-baseline`, `deepscientist-optimize`, and `deepscientist-experiment`. BioOpt implements their handoffs and acceptance checks using files, logs, and commands.

## Scout

Record the user's goal and constraints, the actual implementation, inputs, expected outputs, possible baselines, and open questions. When sources disagree, check the relevant revision and execution path. Prioritize facts that would change the next decision. End scouting when there is enough evidence to define the task and choose a baseline.

## Baseline

| Acceptance level | Required evidence |
|---|---|
| Comparable | The reference has a stable identity. Its version, data and splits where applicable, evaluation entry point, metrics, and environment are traceable. Outputs or metrics can be checked against saved artifacts. |
| Faithful reproduction | Verify the authors' entry point and protocol and check the reference artifacts or metrics. Record any deviations and their implications. Defer this claim when the corresponding evidence is missing. |

Reusing, importing, and rerunning an existing baseline are all valid paths. Choose the least work that meets the required acceptance level. Once verified, record the task, version, data, evaluation entry point, required metrics and their directions, equivalence bounds, evidence locations, and known deviations. If blocked, record the missing item, verified facts, and next viable path.

## Optimize

- **Propose candidates.** Record each hypothesis, code path, evidence, estimated gain, resource risk, and falsification condition. Rank proposals independently in a common table.
- **Select directions.** Choose promising, high-priority candidates. Establish correctness requirements, budget, and stopping criteria before handing off an isolated implementation.
- **Track attempts.** Tie each attempt to a direction and save the patch, status, logs, failure evidence, and next action. This supports rollback and avoids repeating the same attempt.
- **Decide after each batch.** Compare against the strongest verified option. Explore if evidence is incomplete and promising directions remain. Focus verification if one stands out, combine complementary changes, or repair a specific failure. Stop when evidence does not support further gains or the budget is exhausted. Record one primary action and the reason.

## Experiment

Before running, record the hypothesis, smallest meaningful control, inputs, revisions, environment, commands, and expected evidence. Pilot under isolated resources and check scientific metrics and critical assumptions. If the pilot advances, reuse the frozen acceptance criteria for paired runs and retest conditions that could change the result. Retry failures only when an identifiable change to the code, command, environment, or hypothesis justifies it.

Verify that required metrics are present, values meet the acceptance criteria, raw outputs are traceable, and resources are comparable. Identify pilot results, formal pairs, resource scaling tests, and estimates. Disclose contention, warm-up effects, and other confounders. Compare combined candidates with the reference and standalone changes. Link conclusions to verifiable records and explain why each candidate was accepted, rejected, or deferred.

BioOpt adds biological acceptance criteria, four resource routes for candidate generation, independent evaluation, parallel pilots, and expert and Agent case libraries. The main Skill defines how to assign work and select cases.

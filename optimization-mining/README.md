# Optimization mining: omicverse-inherited methods

This directory collects **methods that omicverse aggregates** and for which a **later optimized implementation** exists in the literature or as code. For each pair we record the original implementation, the optimized implementation, the optimization trajectory (what changed and why), the reported speed evidence, and the equivalence caveats.

The goal is not to catalogue every tool. It is to extract **transferable optimization trajectories** that can be mapped onto Bio-Opt patterns (P1–P10) and applied to other scientific code.

## Selection criteria

A method qualifies when all of the following hold:

1. It is exposed by omicverse (`omicverse.single`, `omicverse.space`, `omicverse.bulk`) or is a direct dependency of an official omicverse tutorial.
2. An original implementation and a later optimized/next-generation implementation both exist with public code.
3. The optimized implementation changes execution (algorithm, data structure, system, or approximation) in a way that can be described as a trajectory.
4. There is a reported speed comparison, or the code/paper states the performance motivation explicitly.

## Entry schema

Each entry is `entries/<slug>.md` with YAML frontmatter:

```yaml
slug: <kebab-case>
method: <name>
omicverse_binding: <module/function or "dependency of tutorial">
original:
  paper: <citation or URL>
  repo: <URL>
  commit: <hash/tag or "unpinned">
optimized:
  paper: <citation or URL>
  repo: <URL>
  commit: <hash/tag or "unpinned">
optimization_type: algorithmic | system | data-structure | approximation | workflow
mechanism: <one sentence>
reported_speedup: <value + dataset + hardware, or "not reported">
equivalence: exact | approximate | changed-output | unknown
optimization_families: [O1, O3, ...]   # see TAXONOMY.md
sources:
  - <URL>
verified: true | partial
```

Body sections:

- **Baseline bottleneck** — what limits the original.
- **Optimization trajectory** — ordered steps, each with the change and the reason.
- **Evidence** — numbers, dataset, hardware; quote the source.
- **Equivalence and caveats** — what changes in the output; what must be re-validated.
- **Reuse plan for Bio-Opt** — which pattern(s) and what to test first.

## Index

See [INDEX.md](INDEX.md) for the generated table (method, original/optimized, type, reported speedup, equivalence, families, verified).

## Rules

- Only claim what the cited source states; mark `verified: partial` when a fact could not be confirmed on the official page.
- Record code locations (repo + path) so the original and optimized implementations can be retrieved.
- Do not copy large source files; quote only the minimal snippet that shows the optimization.
- Keep the trajectory reusable: the reader should be able to apply the same idea without adopting the other tool.

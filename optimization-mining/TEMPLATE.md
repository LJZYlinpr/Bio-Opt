---
slug: example-method-fastvariant
method: ExampleMethod
omicverse_binding: "omicverse.single.Example / dependency of tutorial t_example"
original:
  paper: "Author et al., Journal Year, doi:..."
  repo: "https://github.com/org/original"
  commit: "unpinned"
optimized:
  paper: "Author et al., Journal Year, doi:..."
  repo: "https://github.com/org/fastvariant"
  commit: "unpinned"
optimization_type: algorithmic
mechanism: "Replace the O(N^2) pairwise loop with a convolution over binned neighborhoods."
reported_speedup: "20x on dataset X (GPU Y), per paper Table 2"
equivalence: approximate
bioopt_patterns: [P2, P7]
sources:
  - "https://..."
verified: partial
---

## Baseline bottleneck

## Optimization trajectory

1. **Step** — change; why it removes work.
2. **Step** — ...

## Evidence

## Equivalence and caveats

## Reuse plan for Bio-Opt

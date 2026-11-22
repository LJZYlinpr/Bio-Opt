# BioOpt project submission

## Project name

BioOpt: An autonomous performance optimization skill for bioinformatics algorithms, illustrated with spatial omics

## One-sentence description

BioOpt helps AI agents verify locally reproduced implementations, profile bottlenecks, test standalone and combined optimization candidates, and check scientific outputs.

## Problem

Spatial omics workloads can be slowed by redundant computation, data movement, and high memory use. Researchers need a repeatable way to decide what to optimize, document why a candidate was selected, and verify the outputs after optimization.

## Implemented features

The release includes an importable skill, paired timing and independent output verification helpers, and a static website with a six-step interactive SpatialZ case study. The repository includes the original SpatialZ code, manual optimizations, a later batched MENDER benchmark, the SpaGRN data transfer optimization, and the UniST thread configuration comparison. The SpaGRN case documents which results were tested and how the original full-run time was estimated.

## Technical approach

Verify the local reproduction against the upstream task and fix the scientific acceptance criteria and timing boundaries. Use profiling and the expert and Agent trajectory libraries to propose changes, then run isolated pilots and an exploratory follow-up. An independent reviewer filters and ranks the measured candidates. Compare passing standalone changes with compatible combinations in a separate fusion round, then test the best eligible implementation at full scale and archive the evidence. In SpatialZ, the manual implementation computes cosine scores for eligible donors, moves expression sampling to C++, and hands off each layer between GPU and CPU. Later Agent work grouped MENDER calls to amortize preparation cost. Two paired runs of 12 MERFISH target layers achieved 8.40× and 8.70× generation speedups. XY coordinates, cell labels, and expression matched with a fixed random seed per layer.

## Demo

Open the project showcase in a browser. For offline use, save the HTML and open it locally. [Live showcase](https://ljzylinpr.github.io/Bio-Opt/) · [Download HTML](https://github.com/LJZYlinpr/Bio-Opt/raw/refs/heads/main/docs/index.html)

The six-step SpatialZ walkthrough supports keyboard navigation with the arrow keys. [Interactive walkthrough](https://ljzylinpr.github.io/Bio-Opt/spatialz-story.html)

## Code repository

[GitHub · LJZYlinpr/Bio-Opt](https://github.com/LJZYlinpr/Bio-Opt)

## Third-party resources

- [Original SpatialZ code and MERFISH tutorial](https://github.com/senlin-lin/SpatialZ/tree/e1b5ed01e933729149e08399febca08db16ced14)
- [MENDER](https://www.nature.com/articles/s41467-023-44367-9) · [POT](https://pythonot.github.io/) · [PyTorch](https://pytorch.org/) · [NumPy](https://numpy.org/) · [SciPy](https://scipy.org/) · [AnnData](https://anndata.readthedocs.io/)
- [SpaGRN](https://github.com/BGI-Qingdao/SpaGRN/tree/61690ca69a766f14c580de5dd47d4bc7b603813a) · [UniST](https://github.com/lanshui98/UniST/tree/a2fd35a70a6080caab5fdd60652242f85059b09a) · [UniST interpolation tutorial](https://unist-tutorial.readthedocs.io/en/latest/interpolation.html)

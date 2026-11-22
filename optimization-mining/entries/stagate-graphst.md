---
slug: stagate-graphst
method: "STAGATE -> GraphST (graph self-supervised spatial clustering)"
omicverse_binding: "omicverse.space.pySTAGATE class (vendored omicverse/external/STAGATE_pyG) and omicverse.space.clusters(methods=['GraphST']) (vendored omicverse/external/GraphST); omicverse 2.4.3"
original:
  paper: "Dong K, Zhang S. Deciphering spatial domains from spatially resolved transcriptomics with an adaptive graph attention auto-encoder. Nat Commun 13:1739 (2022), doi:10.1038/s41467-022-29439-6"
  repo: "https://github.com/QIFEIDKN/STAGATE"
  commit: "48ce7f874c83a9f1f68187be00370181261ab7c5"
optimized:
  paper: "Long Y et al. Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST. Nat Commun 14:1155 (2023), doi:10.1038/s41467-023-36796-3"
  repo: "https://github.com/JinmiaoChenLab/GraphST"
  commit: "d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74"
optimization_type: algorithmic
mechanism: "Replace per-edge graph attention (and its tied-attention reconstruction decoder) with a linear graph convolution on a once-normalized sparse adjacency plus a permutation-contrastive DGI objective, cutting per-epoch message passing to two sparse matmuls."
reported_speedup: "not reported for STAGATE vs GraphST (no wall-clock or 'faster' statement in the GraphST paper; full text contains no runtime comparison); the authors' own PyTorch Geometric port of the same STAGATE reports 'more than 10x faster than STAGATE based on the tensorflow1 framework' (STAGATE_pyG README)"
equivalence: changed-output
optimization_families: [O4, O5, O8, O10]
sources:
  - "https://doi.org/10.1038/s41467-022-29439-6"
  - "https://www.nature.com/articles/s41467-022-29439-6"
  - "https://github.com/QIFEIDKN/STAGATE"
  - "https://github.com/QIFEIDKN/STAGATE_pyG/blob/ae1158ca8cf1eb6bb8ee198298552d44c9ac21db/README.md"
  - "https://github.com/QIFEIDKN/STAGATE_pyG/blob/ae1158ca8cf1eb6bb8ee198298552d44c9ac21db/STAGATE_pyG/gat_conv.py"
  - "https://github.com/QIFEIDKN/STAGATE_pyG/blob/ae1158ca8cf1eb6bb8ee198298552d44c9ac21db/STAGATE_pyG/Train_STAGATE.py"
  - "https://doi.org/10.1038/s41467-023-36796-3"
  - "https://pmc.ncbi.nlm.nih.gov/articles/PMC9977836/"
  - "https://github.com/JinmiaoChenLab/GraphST/blob/d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74/GraphST/GraphST/model.py"
  - "https://github.com/JinmiaoChenLab/GraphST/blob/d62b0b7b6cd38ee285f3ac8cd67b7341a10bcc74/GraphST/GraphST/preprocess.py"
  - "https://github.com/Starlitnightly/omicverse/blob/d0c3325ca72404e28ec5161a2f71ea82db531974/omicverse/space/_cluster.py"
  - "https://deepst-tutorials.readthedocs.io/en/latest/Tutorial%201_10X%20Visium.html"
verified: true
---

## Baseline bottleneck

STAGATE is a graph attention auto-encoder. Four `GATConv` layers are stacked (two encoder, two decoder); edge attention is computed in the first and third layers (the decoder ties its attention to the encoder's), normalised with a per-node softmax, and used to weight neighbour features. The model is trained for 1000 epochs on the expression-reconstruction MSE (defaults in `Train_STAGATE.py:16`). Two consequences dominate cost:

- per-epoch work is edge-wise: each attention layer evaluates `att_src`/`att_dst` dot products per edge, softmaxes per destination node, and then scatters `x_j * alpha` (`STAGATE_pyG/STAGATE_pyG/gat_conv.py:201-213`);
- the decoder re-applies the tied attention and then reconstructs the full input gene axis, so the backward path is nearly as expensive as the forward path even though only one decoder layer scores edges.

```python
# STAGATE_pyG/STAGATE_pyG/gat_conv.py:201-213

def message(self, x_j, alpha_j, alpha_i, index, ptr, size_i):
    alpha = alpha_j if alpha_i is None else alpha_j + alpha_i
    alpha = torch.sigmoid(alpha)
    alpha = softmax(alpha, index, ptr, size_i)
    self._alpha = alpha
    alpha = F.dropout(alpha, p=self.dropout, training=self.training)
    return x_j * alpha.unsqueeze(-1)
```

The public original is TensorFlow 1 (`QIFEIDKN/STAGATE`); the authors state the PyG rewrite alone gives "more than 10x faster than STAGATE based on the tensorflow1 framework" and enables batch training (`STAGATE_pyG/README.md`). That number belongs to the reimplementation of the same model, not to GraphST.

## Optimization trajectory

1. **Move to PyG sparse kernels (same model).** `STAGATE_pyG` keeps the attention architecture but runs message passing through PyTorch Geometric and adds spatial batch training (`Batch_Data`). Reported >10x versus the TF1 code (`README.md`); this is an implementation-only step.
2. **Delete the attention mechanism.** GraphST's encoder has only two learnable projections: `z = adj @ (X W1)`, `h = adj @ (z W2)`, with `adj` applied by `torch.mm` (dense 10X path) or `torch.spmm` (sparse Slide/Stereo path). There is no `att_src`/`att_dst`, no per-edge sigmoid, no softmax and no tied-attention decoder (`GraphST/GraphST/model.py:70-97`). Per-epoch message passing becomes a small number of dense/sparse matrix products plus one discriminator, independent of per-edge scoring.
```python
# GraphST/GraphST/model.py:70-86 (Encoder; Encoder_sparse swaps mm for spmm)

def forward(self, feat, feat_a, adj):
    z = torch.mm(F.dropout(feat, self.dropout, self.training), self.weight1)
    z = torch.mm(adj, z)
    h = torch.mm(z, self.weight2)
    h = torch.mm(adj, h)
    emb = self.act(z)
    z_a = torch.mm(adj, torch.mm(feat_a, self.weight1))
    emb_a = self.act(z_a)
```

3. **Normalise the graph once.** `preprocess_adj` / `preprocess_adj_sparse` build the symmetrically normalised adjacency once, outside the training loop, and the same object is reused by both augmented views, both encoder projections and every epoch (`GraphST/GraphST/preprocess.py:139-158`). This removes the per-forward rebuild of edge weights that attention implies.
4. **Swap the objective, not just the layer.** Instead of reconstruction with an attention decoder, GraphST optimises `alpha * MSE(X, adj @ (z W2)) + beta * BCE(DGI discriminator on a feature-permuted view)`: the second view is obtained by row permutation of the feature matrix, so it needs no second attention pass and no extra graph (`GraphST/GraphST/GraphST.py:167-190`, defaults `epochs=600`, `dim_output=64`).
5. **Cluster the 64-d embedding.** The official GraphST tutorial runs mclust on the learned representation (with optional label refinement); the model itself ends at the embedding.

Why less work: per epoch, STAGATE pays edge-wise score/softmax work in two attention layers and a 4-layer reconstruction stack over the full gene axis for 1000 epochs. GraphST pays two graph products on 3000-gene features, one 64->3000 reconstruction and a 64-d bilinear discriminator, for 600 epochs. The graph normalisation is amortised across all epochs and views.

## Evidence

- Accuracy, GraphST paper, 12 LIBD DLPFC slices: "GraphST achieved the highest median score of 0.60 among all methods, with only STAGATE's performance as a close second" (`PMC9977836`, Results, Fig. 2A).
- Accuracy, human breast cancer Visium (3798 spots, 20 annotated regions): "Seurat, Giotto, SpaceFlow, and STAGATE had lower scores in the range of 0.39-0.47, while BayesSpace, SpaGCN, and GraphST had higher scores between 0.54 and 0.57" (same paper).
- Speed: no runtime, wall-clock, GPU model or "faster" claim anywhere in the GraphST full text; the paper motivates GraphST by accuracy, not speed. The only measured speed statement in this lineage is the PyG port's ">10x" versus TF1 STAGATE.

## Equivalence and caveats

`changed-output`: GraphST is a different model with a different objective, so cluster labels are not expected to match STAGATE even with the same seed. Agreement must be judged by ARI/NMI against annotation, not by label identity.

Accuracy is benchmark-dependent. GraphST's own DLPFC comparison puts it first (0.60) with STAGATE close behind, but the independent BANKSY benchmark found "GraphST, BANKSY and STAGATE were the top performers on this dataset in terms of median ARIs, with no significant difference between the three methods" (Nat Genet 2024, Fig. 5b; paired p = 0.72 and 0.48 for the last two comparisons). Treat the GraphST-vs-STAGATE accuracy gain as small and dataset-specific.

Graph construction remains a scalability limit on the 10X path: `construct_interaction` calls `ot.dist` on all spots and materialises an `n x n` interaction matrix before converting it (`GraphST/GraphST/preprocess.py:50-96`), so memory is quadratic even though the model itself is sparse. The sparse encoder path is used only for `datatype` in `{'Stereo','Slide'}`.

The attention removal also removes STAGATE's adaptive edge weighting, which is its stated modelling contribution; if edge weights matter for a dataset (e.g. boundaries with mixed cell types), the cheaper model is not a drop-in statistical equivalent.

## Reuse plan for Bio-Opt

omicverse binding: `ov.space.pySTAGATE` (`omicverse/space/_cluster.py:59`) wraps the vendored `external/STAGATE_pyG`; `ov.space.clusters(adata, methods=['GraphST'], ...)` (`_cluster.py:434`, branch at 513) wraps `external/GraphST` and stores `adata.obsm['GraphST_embedding']`. Both are first-class omicverse entry points, so the pair can be benchmarked inside one environment.

- ** prune unused work first**: replace the per-edge attention score with a fixed normalised adjacency when the graph is already kNN/radius and the edge weights are not used downstream; measure ARI before and after, since the attention is the method's modelling claim.
- ** reuse shared state**: normalise the adjacency and precompute both views once; confirm the training loop does not rebuild them per epoch (the original attention path effectively does).
- ** vectorise/compile**: keep message passing in sparse kernels (`spmm`) and avoid densifying `n x n`; on omicverse's Visium path the remaining O(n^2) object is graph construction, which is the next thing to fix.

Pattern codes as used here (the referenced `pattern-catalog.md` was absent from this checkout):  prune unused work,  reuse shared/target-invariant state,  vectorise/compile.

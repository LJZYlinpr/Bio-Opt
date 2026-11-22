---
slug: gsea-fgsea
method: "GSEA -> fgsea (preranked gene set enrichment)"
omicverse_binding: "omicverse.bulk.pyGSEA (default backend='numpy' is a NumPy reimplementation of gseapy.prerank; fgsea is not used) and omicverse.bulk.geneset_enrichment_GSEA"
original:
  paper: "Subramanian et al., PNAS 2005, doi:10.1073/pnas.0506580102 (GSEA reference implementation, Broad)"
  repo: "https://www.gsea-msigdb.org/gsea/index.jsp"
  commit: "unpinned"
optimized:
  paper: "Korotkevich, Sukhov, Budin, Shpak, Artyomov & Sergushichev, bioRxiv 2021, doi:10.1101/060012"
  repo: "https://github.com/ctlab/fgsea"
  commit: "570f5903545835cb340e699fcdfdf444b1ba07ec"
optimization_type: algorithmic
mechanism: "Rewrite preranked GSEA in C++ so the enrichment score is computed from the K hit positions only, one sampled gene-set null is shared by all pathways of the same size, and adaptive multilevel splitting replaces 1e5-1e8 brute-force permutations."
reported_speedup: "~100x on the reported example: ~4 s versus ~420 s for the reference implementation (700 Reactome pathways, n=10,000 permutations, Th1-activation dataset, preprint Fig. 2c); 10 s on 1 thread versus 40 min on 32 threads for 10^8 permutations at the same p-values (Fig. 3b); median 8 s per dataset over 605 datasets on a 4-core Intel i5 laptop (Fig. 4a)"
equivalence: approximate
optimization_families: [O1, O6, O8, O10]
sources:
  - "https://doi.org/10.1101/060012"
  - "https://www.biorxiv.org/content/10.1101/060012v3"
  - "https://github.com/ctlab/fgsea/blob/570f5903545835cb340e699fcdfdf444b1ba07ec/src/esCalculation.cpp"
  - "https://github.com/ctlab/fgsea/blob/570f5903545835cb340e699fcdfdf444b1ba07ec/src/fastGSEA.cpp"
  - "https://github.com/ctlab/fgsea/blob/570f5903545835cb340e699fcdfdf444b1ba07ec/src/ScoreRuler.cpp"
  - "https://github.com/ctlab/fgsea/blob/570f5903545835cb340e699fcdfdf444b1ba07ec/R/fgseaMultilevel.R"
verified: true
---

## Baseline bottleneck

The GSEA reference implementation tests two phenotypes by permuting sample labels. For every permutation it must rebuild the gene ranking and recompute the running-sum enrichment score (ES) for every pathway over all N ranked genes, so cost is `n_perm x n_pathways x O(N)` plus the ranking rebuild. Two consequences are structural, not just slow:

- the attainable p-value resolution is set by the permutation count (`1/(n_perm+1)` with the recommended +1 estimator), which is the preprint's stated motivation: "the reference implementation of this method cannot accurately estimate small P-values";
- permutation count, not the data, sets the runtime, so resolving p ~1e-4 already needs ~1e5 permutations and resolving the tails that survive BH correction over ~10^3 pathways becomes expensive.

fgsea targets the *preranked* problem, where the ranking is fixed input, so label permutation is not available and gene-set sampling is the correct null.

## Optimization trajectory

1. **Fix the ranking once** — accept gene-level statistics as input; the ranking is sorted once and reused by every pathway and every null sample (reuse target-invariant state). This removes the per-permutation re-ranking that dominates the reference implementation.
2. **Compute ES from hits, not from all genes** — `calcES` walks only the sorted positions of the K set members, accumulating a running sum; the gaps between hits are handled by a linear coefficient, so each step is O(1) and the ES is O(K) instead of O(N):

```cpp
// fgsea: src/esCalculation.cpp:4-21
score_t calcES(const vector<int64_t> &ranks, const vector<int> &p, int64_t NS) {
    int n = (int) ranks.size;
    int k = (int) p.size;
    score_t res{NS, 0, n - k, 0};
    score_t cur{NS, 0, n - k, 0};
    int last = -1;
    for (int pos : p) {
        cur.coef_const += pos - last - 1;
        if (res.abs < cur.abs) { res = cur; }
        cur.coef_NS += ranks[pos];
        if (res.abs < cur.abs) { res = cur; }
        last = pos;
    }
    return res;
}
```

3. **One null sample serves every pathway size** — a single draw of `k = max(pathway sizes)` random positions yields cumulative ES values for all prefixes; each pathway then takes the prefix matching its own size (share initialization across consumers). One sampling pass therefore produces null ES for the whole collection, and the cost of a permutation no longer scales with the number of pathways:

```cpp
// fgsea: src/fastGSEA.cpp:405-412 (inside calcGseaStatCumulativeBatch)
for (int i = 0; i < iterations; ++i) {
    NumericVector randEs = calcRandomGseaStatCumulative(stats, n, k, gseaParam, rng, scoreType);
    NumericVector randEsP = subvector(randEs, pathwaysSizes);   // all sizes at once
    aux = randEsP <= pathwayScores; ...
}
```

4. **Escalate only the pathways that need it** — `fgseaMultilevel` first runs a cheap 1,000-permutation pass, converts both a beta-binomial "simple" error and the expected multilevel error, and sends only pathways where `multError < simpleError` to the multilevel solver, grouped by size (prune unused work: extra sampling is spent only on pathways whose p is not yet resolved).
5. **Adaptive multilevel splitting** — for the hard ones, `ScoreRuler::extend` builds a sample of `sampleSize` (default 101) random gene sets, then repeatedly *replaces elements below the median with duplicates of elements above it* and applies conditional MCMC moves, so each doubling of the p-value resolution costs one sample, not a full permutation sweep (controlled approximation with a `log2err` estimate). The p-value is read off the accumulated score ruler in O(log n).

```cpp
// fgsea: src/ScoreRuler.cpp:65-76
for (unsigned elemIndex = 0; elemIndex < sampleSize; elemIndex++) {
    std::vector<int> comb = combination(0, n - 1, genesetSize, mtGen);
    currentSample[elemIndex] = std::vector<unsigned>(comb.begin, comb.end);
    currentProfiles[elemIndex] = getProfile(expressionMatrix, currentSample[elemIndex], m);
}
duplicateSampleElements;                       // keep only the upper half
while (scores.back <= inpScore - 1e-10) { ... updateElement(...) ... duplicateSampleElements; }
```

6. **Exact integer-statistics reference and parallel execution** — the paper also derives an exact algorithm for integer gene-level statistics, used to validate the approximation; pathway groups are dispatched with `BiocParallel::bplapply` (`nproc`/`BPPARAM`).

Relevant locations: `src/esCalculation.cpp:6-40`, `src/fastGSEA.cpp:330-430`, `src/ScoreRuler.cpp:1-120`, `src/fgseaMultilevel.cpp:1-60`, `R/fgseaMultilevel.R:46-260`.

## Evidence

- "on this example the reference implementation took about 420 seconds, while FGSEA-simple finished in about 4 seconds. The two order of magnitude speed-up is consistent with the theoretical one due to the algorithm time complexity" (700 Reactome pathways, n=10,000 permutations) — https://www.biorxiv.org/content/10.1101/060012v3
- "FGSEA-multilevel calculation with sample size of Z=101 took only 10 seconds working on a single thread while 10^8 permutations on FGSEA-simple took 40 minutes working in 32 threads", with p-values down to 1e-26 and consistent results — same URL, Fig. 3b.
- "it took less then a minute per dataset ... with the multilevel algorithm (Fig. 4a) on a laptop with 4-core Intel Core i5 processor, with a median time of 8 seconds"; 605 curated GEO datasets with C5_BP — same URL.
- Implementation defaults: `sampleSize = 101`, `nPermSimple = 1000`, `eps = 1e-10` — `R/fgseaMultilevel.R:46-56`.

## Equivalence and caveats

- **ES/NES**: for the same preranked statistics, gene sets and `gseaParam`, the ES is the same quantity as GSEA's (same weighted Kolmogorov-Smirnov running sum, same sign convention), and the paper reports nominal p-values "indistinguishable from each other up to the random noise inherent to both methods". This is why the entry is `approximate` rather than `exact`: p-values are Monte Carlo estimates.
- **Null model differs from the reference**: fgsea samples random gene sets of the same size (valid only for preranked input); the Broad desktop tool permutes phenotype labels on expression data. For preranked runs the two nulls are not formally identical, and the paper's comparison is empirical.
- **+1 correction**: the paper recommends the `(n+1)`-style estimator; without it, p-values can be anti-conservative at the sampling floor.
- **`eps`**: default `1e-10` clamps p-values; set `eps = 0` for deeper tails (README).
- **Ties / zero statistics**: fgsea warns on tied gene-level statistics and replaces zero statistics with an epsilon for the ES computation, which can perturb ES on sparse or heavily tied inputs.
- **`log2err` must be reported** with each p-value; it is `NA` when the conditional probability check fails (p likely overestimated).

## Reuse plan for Bio-Opt

- ** first**: for any enrichment/scanning kernel, express the statistic as a function of hit positions and skip the non-hit span analytically; verify ES equality to machine precision against the original on a small collection.
- ** next**: ask whether one expensive sample (a null draw, a distance matrix row, an initialization) can be shared by all consumers of the same size class; the win scales with the number of consumers (pathway collections are ~10^3).
- ** last**: replace brute-force resampling with an adaptive scheme that carries an explicit error estimate, and freeze acceptance as "p within `log2err` of a 10^5-permutation estimate" before timing.
- Cross-check ES and leading edge exactly; compare p-values only with tolerance, and always report the sampling depth that produced them.

# ACAT aggregated-Cauchy omnibus over existing results

Function:
[`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.md)

Result column: `acat`

Family: P-value aggregation

## Introduction

Combines the p-value columns already present in `pleio@results` using
the Aggregated Cauchy Association Test (Liu et al. 2019; Liu & Xie
2019), an omnibus that merges a set of (possibly dependent) p-values via
a Cauchy-transformed weighted sum. This is a meta-combination of
previously computed test results, not a test of the raw traits; it is
not part of \[run_pleiotest()\] and should be called afterwards.

## Mathematical Framework and Implementation

For each variant let $`p_1,\dots,p_M`$ be the p-values taken from the
selected result columns (one per chosen test), with non-negative weights
$`w_1,\dots,w_M`$ normalized to sum to one ($`w_k=1/M`$ when no weights
are supplied). ACAT maps each p-value to the upper tail of a standard
Cauchy and forms the weighted sum

``` math
T = \sum_{k=1}^{M} w_k\,\tan\!\big[(0.5 - p_k)\,\pi\big].
```

Because each $`\tan[(0.5-p_k)\pi]`$ is standard-Cauchy under the null
and a weighted sum of standard-Cauchy variables is again standard Cauchy
regardless of their dependence, the combined p-value is the standard
Cauchy upper tail

``` math
p_{\mathrm{ACAT}} = \tfrac{1}{2} - \tfrac{1}{\pi}\arctan(T).
```

For numerical stability very small p-values use the tail expansion
$`\tan[(0.5-p_k)\pi]\approx 1/(\pi p_k)`$ when $`p_k<10^{-15}`$;
non-finite p-values are skipped, and the returned p-value is floored at
$`10^{-300}`$. Weights are normalized once for the requested column set;
if a row has non-finite inputs, skipped entries are not row-wise
renormalized.

## What This Method Advances

ACAT gives an analytic way to combine dependent p-values without
estimating their correlation matrix, using the stability of the Cauchy
distribution.

## Assumptions

Each combined p-value is (approximately) uniform on $`[0,1]`$ under its
own null. The Cauchy combination is robust to arbitrary correlation
among the combined statistics in the bulk of the distribution and
remains accurate in the extreme tail, so no correlation matrix is
needed. The inputs are existing result columns; their validity as
p-values is the caller’s responsibility.

## When This Method Is a Good Fit

Good for combining already-computed result columns, especially when
those columns are dependent and a correlation model is inconvenient.

ACAT combines existing result columns. It does not recompute raw trait
evidence and is intentionally not included inside
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md).

The wrapper does not validate weight length, non-negativity, positive
weight sum, or that all finite inputs lie in `[0, 1]`. Also note that
rerunning `run_acat(pleio)` with `cols = NULL` after an `acat` column
already exists will include that existing `acat` column in the default
set of columns.

## Minimal Usage

``` r

library(pleioverse)
set.seed(1)
pleio <- run_pleiosim(
  n_phenotype = 3, n_participant = 500,
  n_variant_pleiotropic = 3,
  n_variant_nonpleiotropic = c(3, 3, 3),
  n_variant_null = 30,
  heritability = c(0.1, 0.2, 0.3)
)
pleio <- run_wald(pleio)
pleio <- run_vc(pleio)
pleio <- run_acat(pleio, cols = c('wald', 'vc'))
head(pleio@results$acat)
```

## Implementation Files

- R wrapper: `R/run_acat.R`
- C++ kernel or dispatcher: `src/acat.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu Y, Xie J (2019). Cauchy combination test: a powerful test with
analytic p-value calculation under arbitrary dependency structures. *J
Am Stat Assoc* 115(529):393-402. DOI:
<https://doi.org/10.1080/01621459.2018.1554485>

Liu Y, Chen S, Li Z, Morrison AC, Boerwinkle E, Lin X (2019). ACAT: a
fast and powerful p-value combination method for rare-variant analysis
in sequencing studies. *Am J Hum Genet* 104(3):410-421. DOI:
<https://doi.org/10.1016/j.ajhg.2019.01.002>

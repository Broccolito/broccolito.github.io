# Weighted-identity (WI) multi-trait test

Function:
[`run_wi()`](https://broccolito.github.io/pleioverse/reference/run_wi.md)

Result column: `wi`

Family: Multiple-association omnibus

## Introduction

MPAT/MTAR multi-phenotype association test that uses the unweighted sum
of squared per-trait Z-scores. Because the Z-scores are correlated under
the null, this statistic is a weighted mixture of one-degree-of-freedom
chi-squares whose weights are the eigenvalues of the trait correlation
matrix, and its tail is evaluated numerically.

## Mathematical Framework and Implementation

For a variant with per-trait Z-score vector $`z = (z_1,\dots,z_K)'`$,
$`z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)`$, over $`K`$ phenotypes,
the statistic is the plain (identity-weighted) sum of squares

``` math
T = z' z = \sum_{k=1}^{K} z_k^2.
```

Under the global null $`z \sim \mathcal{N}(0, \Sigma)`$, with $`\Sigma`$
the $`K \times K`$ trait correlation matrix, write the spectral
decomposition $`\Sigma = \sum_{k} \lambda_k\, v_k v_k'`$ with
eigenvalues $`\lambda_1 \ge \dots \ge \lambda_K > 0`$ (`eval` in the
engine). Then

``` math
T \;\stackrel{d}{=}\; \sum_{k=1}^{K} \lambda_k\, \chi^2_{1,k},
```

a weighted sum of independent $`\chi^2_1`$ variables with weights equal
to the eigenvalues of $`\Sigma`$. The p-value

``` math
p = \Pr\!\Big(\sum_{k} \lambda_k \chi^2_{1,k} \ge T\Big)
```

is computed by the Davies-labeled numerical inversion of the
characteristic function; if that returns a value outside $`(0,1)`$ the
routine falls back to the modified Liu moment-matching (Liu-Tang-Zhang)
approximation.

## What This Method Advances

WI avoids inversion of the trait-correlation matrix and calibrates the
resulting identity-weighted quadratic form as a weighted chi-square
mixture.

## Assumptions

WI weights every trait equally on the original scale and does not whiten
by $`\Sigma^{-1}`$, so it up-weights the directions of *high* trait
correlation (large eigenvalues). It is well powered for *dense*,
heterogeneous signal that is aligned with the leading eigenvectors of
$`\Sigma`$, and is more robust than the Wald test when $`\Sigma`$ is
near-singular (it never inverts $`\Sigma`$). It is comparatively weak
for sparse signal concentrated in a single trait or for effects lying in
the low-variance directions of $`\Sigma`$. Validity requires correct
specification of $`\Sigma`$ and approximately Gaussian null Z-scores.

## When This Method Is a Good Fit

Good for dense signal on the original trait scale, especially when
avoiding a matrix inverse is attractive.

The implementation treats variants independently and relies on the
trait-correlation estimate derived from the available summary statistics
unless an individual method estimates a pairwise null moment internally.

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
pleio <- run_wi(pleio)
head(pleio@results$wi)
```

## Implementation Files

- R wrapper: `R/run_wi.R`
- C++ kernel or dispatcher: `src/wi.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175. DOI: <https://doi.org/10.1111/biom.12735>

Davies, R. B. (1980). Algorithm AS 155: The distribution of a linear
combination of $`\chi^2`$ random variables. *Journal of the Royal
Statistical Society, Series C (Applied Statistics)*, 29(3), 323-333.
DOI: <https://doi.org/10.2307/2346911>

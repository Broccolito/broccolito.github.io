# Wald multi-trait test (inverse-variance quadratic form)

Function:
[`run_wald()`](https://broccolito.github.io/pleioverse/reference/run_wald.md)

Result column: `wald`

Family: Multiple-association omnibus

## Introduction

MPAT/MTAR multi-phenotype association test that combines the per-trait
Z-scores of a single variant into one inverse-correlation-weighted
quadratic form. Under the global null of no association with any trait
this statistic follows an exact chi-square distribution, giving a fast
closed-form p-value.

## Mathematical Framework and Implementation

For a variant with per-trait Z-score vector $`z = (z_1,\dots,z_K)'`$,
where $`z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)`$ is the marginal
score (Wald) statistic for trait $`k`$ and $`K`$ is the number of
phenotypes, the test statistic is the Mahalanobis quadratic form

``` math
T = z' \Sigma^{-1} z,
```

where $`\Sigma`$ is the $`K \times K`$ trait correlation matrix (the
correlation of the Z-scores across variants under the null), and
$`\Sigma^{-1}`$ (denoted `Sinv` in the engine) is precomputed once.

Under the global null $`H_0: z \sim \mathcal{N}(0, \Sigma)`$ the
whitened vector $`\Sigma^{-1/2} z`$ is standard multivariate normal, so

``` math
T = (\Sigma^{-1/2} z)'(\Sigma^{-1/2} z) \sim \chi^2_K.
```

The reported p-value is the exact upper tail of a chi-square with $`K`$
degrees of freedom,

``` math
p = \Pr(\chi^2_K \ge T).
```

No moment matching or numerical inversion is needed; the p-value is
closed form.

## What This Method Advances

The Wald article represents the baseline multivariate omnibus test: it
spends one degree of freedom per phenotype and gives an exact chi-square
reference once the trait-correlation matrix is estimated.

## Assumptions

This is the classical multivariate Wald / Hotelling-style omnibus test.
It spends one degree of freedom per trait and is therefore most powerful
when the signal is *dense* (many traits associated) and effect
directions are *heterogeneous*, with no prior knowledge of the sign or
magnitude pattern across traits. Because it weights by $`\Sigma^{-1}`$,
it accounts for trait correlation but, by treating all $`K`$ directions
equally, loses power relative to burden-type tests
([`run_sum()`](https://broccolito.github.io/pleioverse/reference/run_sum.md),
[`run_dsum()`](https://broccolito.github.io/pleioverse/reference/run_dsum.md))
when effects are sparse or concordant in direction. Validity requires
$`\Sigma`$ to be correctly specified and non-singular and the Z-scores
to be approximately Gaussian under the null (large-sample GWAS summary
statistics).

## When This Method Is a Good Fit

Good for dense, heterogeneous multi-trait association where no common
effect direction is assumed.

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
pleio <- run_wald(pleio)
head(pleio@results$wald)
```

## Implementation Files

- R wrapper: `R/run_wald.R`
- C++ kernel or dispatcher: `src/wald.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175. DOI: <https://doi.org/10.1111/biom.12735>

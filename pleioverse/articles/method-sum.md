# Weighted burden (SUM) multi-trait test

Function:
[`run_sum()`](https://broccolito.github.io/pleioverse/reference/run_sum.md)

Result column: `sum`

Family: Burden-style multiple association

## Introduction

MPAT/MTAR multi-phenotype burden test. It collapses the per-trait
Z-scores into a single inverse-correlation-weighted linear combination
(the optimal GLS contrast for a common cross-trait effect) and tests it
with a two-sided normal tail.

## Mathematical Framework and Implementation

For a variant with per-trait Z-score vector $`z = (z_1,\dots,z_K)'`$,
$`z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)`$, over $`K`$ phenotypes,
define the weight vector $`a = \Sigma^{-1} \mathbf{1}`$, where
$`\Sigma`$ is the $`K \times K`$ trait correlation matrix and
$`\mathbf{1}`$ is the all-ones vector ($`a`$ is precomputed as `a` in
the engine). The burden statistic is the absolute value of the
inverse-variance-weighted sum

``` math
T = \left| a' z \right| = \left| \mathbf{1}' \Sigma^{-1} z \right|.
```

This linear combination is the generalized-least-squares (best linear
unbiased) contrast for detecting a common scalar effect shared across
all traits.

Under the global null $`z \sim \mathcal{N}(0, \Sigma)`$ the underlying
linear form $`a' z = \mathbf{1}' \Sigma^{-1} z`$ is normal with mean 0
and variance

``` math
\mathrm{Var}(a' z) = \mathbf{1}' \Sigma^{-1} \Sigma\, \Sigma^{-1}
  \mathbf{1} = \mathbf{1}' \Sigma^{-1} \mathbf{1} \equiv \texttt{sum\_Sinv},
```

i.e. $`a' z \sim \mathcal{N}(0,\, \mathbf{1}'\Sigma^{-1}\mathbf{1})`$.
The two-sided p-value is

``` math
p = 2\,\Pr\!\left(\mathcal{N}(0, \sigma) \ge T\right),\qquad
  \sigma = \sqrt{\mathbf{1}' \Sigma^{-1} \mathbf{1}}.
```

## What This Method Advances

SUM turns many correlated trait scores into the
generalized-least-squares contrast for a common direction of effect.

## Assumptions

SUM is a one-degree-of-freedom *burden* test: it has high power when
effects are *concordant* (shared common direction) across many traits,
and is the locally most powerful test for a single common effect after
optimally accounting for trait correlation through $`\Sigma^{-1}`$. It
loses power, and can lose it entirely, when effects are *heterogeneous*
in sign so that the weighted contributions cancel, or when signal is
sparse and confined to a single trait; the omnibus Wald
([`run_wald()`](https://broccolito.github.io/pleioverse/reference/run_wald.md))
or variance-component
([`run_vc()`](https://broccolito.github.io/pleioverse/reference/run_vc.md))
tests are preferable there. Validity requires a correctly specified,
non-singular $`\Sigma`$ and approximately Gaussian null Z-scores.

## When This Method Is a Good Fit

Good for concordant effects shared across many traits after
inverse-correlation weighting.

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
pleio <- run_sum(pleio)
head(pleio@results$sum)
```

## Implementation Files

- R wrapper: `R/run_sum.R`
- C++ kernel or dispatcher: `src/sum.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175. DOI: <https://doi.org/10.1111/biom.12735>

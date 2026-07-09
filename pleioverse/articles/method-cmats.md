# CMATS burden-augmented adaptive omnibus association test

Function:
[`run_cmats()`](https://broccolito.github.io/pleioverse/reference/run_cmats.md)

Result column: `cmats`

Family: Adaptive MATS/SKAT-O style

## Introduction

Multi-phenotype association test that adaptively blends a
variance-component (Wald / SKAT-style) statistic with a directional
burden statistic over a grid of mixing weights, in the spirit of the
SKAT-O optimal unified test (Lee et al. 2012; Wu et al. 2011). The
per-variant minimum p-value over the grid is calibrated against its
analytic null distribution.

## Mathematical Framework and Implementation

For each variant let $`z=(z_1,\dots,z_K)'`$ be the vector of per-trait
Z-scores, $`z_k=\hat\beta_k/\mathrm{se}_k`$, and let $`\Sigma`$ be the
$`K\times K`$ trait correlation matrix with inverse $`\Sigma^{-1}`$,
$`a=\Sigma^{-1}\mathbf{1}`$ and
$`s=\mathbf{1}'\Sigma^{-1}\mathbf{1}=\sum_{ij}(\Sigma^{-1})_{ij}`$.

Two component statistics are formed. The variance-component (Wald) term
is the Mahalanobis quadratic form

``` math
\chi_1 = z'\Sigma^{-1}z,
```

which under the global null is $`\chi^2`$ with $`K`$ degrees of freedom.
The burden term aggregates a $`\Sigma^{-1}`$-weighted sum of Z-scores,

``` math
\chi_2 = \frac{(a'z)^2}{s} = \frac{(\mathbf{1}'\Sigma^{-1}z)^2}{\mathbf{1}'\Sigma^{-1}\mathbf{1}},
```

which under the null is $`\chi^2`$ with $`1`$ degree of freedom.

For each $`\rho`$ on the grid $`\{0,0.1,\dots,1\}`$ (11 equally spaced
points) the mixed statistic

``` math
Q_\rho = (1-\rho)\,\chi_1 + \rho\,\chi_2
```

is computed. Its null tail probability $`p_\rho`$ is obtained exactly at
the endpoints ($`\chi^2_K`$ for $`\rho=0`$, $`\chi^2_1`$ for $`\rho=1`$)
and otherwise as the tail of a weighted sum of independent $`\chi^2_1`$
variables with eigenvalue spectrum $`(1-\rho)`$ with multiplicity
$`K-1`$ and a single eigenvalue $`1`$, evaluated by the Liu
moment-matching approximation. The test statistic is the minimum p-value

``` math
T = \min_\rho p_\rho.
```

Because the $`Q_\rho`$ are dependent, the final p-value corrects the
minimum-p statistic with the SKAT-O one-dimensional null integral

``` math
p = \min_\rho p_\rho + \int_0^{t_{\max}} \Pr\Big[\chi^2_{K-1} >
  \min_{i:\,\rho_i<1}\frac{q_i - x(t)}{1-\rho_i}\Big]\,dt,
```

where $`q_i`$ is the quantile of the mixed null at level
$`1-\min_\rho p_\rho`$ for grid point $`i`$, $`x(t)`$ is the
$`\chi^2_1`$ quantile at $`t`$, and $`t_{\max}`$ is the upper limit set
by the $`\rho=1`$ quantile. The integral is evaluated numerically by
composite Simpson’s rule, and the result is clamped to $`[0,1]`$.

## What This Method Advances

CMATS adapts the SKAT-O idea to multi-trait summary statistics by
blending a full Wald component with a burden component.

## Assumptions

Z-scores are approximately multivariate normal with mean zero and
covariance $`\Sigma`$ under the global null of no association with any
trait. $`\Sigma`$ is the (full-rank, invertible) trait correlation
matrix shared across variants and is treated as known. Variants are
tested one at a time; the omnibus combines information across the $`K`$
traits, not across variants.

## When This Method Is a Good Fit

Good when either a full multivariate signal or a common burden-like
direction could be present.

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
pleio <- run_cmats(pleio)
head(pleio@results$cmats)
```

## Implementation Files

- R wrapper: `R/run_cmats.R`
- C++ kernel or dispatcher: `src/mats.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X (2011). Rare-variant
association testing for sequencing data with the sequence kernel
association test (SKAT). *Am J Hum Genet* 89(1):82-93. DOI:
<https://doi.org/10.1016/j.ajhg.2011.05.029>

Lee S, Wu MC, Lin X (2012). Optimal tests for rare variant effects in
sequencing association studies. *Biostatistics* 13(4):762-775. DOI:
<https://doi.org/10.1093/biostatistics/kxs014>

Liu H, Tang Y, Zhang HH (2009). A new chi-square approximation to the
distribution of non-negative definite quadratic forms in non-central
normal variables. *Comput Stat Data Anal* 53(4):853-856. DOI:
<https://doi.org/10.1016/j.csda.2008.11.025>

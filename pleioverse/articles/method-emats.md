# EMATS leading-eigenvector adaptive omnibus association test

Function:
[`run_emats()`](https://broccolito.github.io/pleioverse/reference/run_emats.md)

Result column: `emats`

Family: Adaptive MATS/SKAT-O style

## Introduction

Multi-phenotype association test that adaptively blends a
variance-component (Wald / SKAT-style) statistic with a
leading-eigenvector directional statistic over a grid of mixing weights,
following the SKAT-O optimal unified framework (Lee et al. 2012; Wu et
al. 2011). It follows the same adaptive framework as \[run_cmats()\],
but uses a 6-point grid and projects the Z-scores onto the principal
axis of $`\Sigma`$ rather than onto the burden direction.

## Mathematical Framework and Implementation

For each variant let $`z=(z_1,\dots,z_K)'`$ be the vector of per-trait
Z-scores, $`z_k=\hat\beta_k/\mathrm{se}_k`$, and let $`\Sigma`$ be the
$`K\times K`$ trait correlation matrix with eigendecomposition
$`\Sigma = \sum_{j} \lambda_j v_j v_j'`$, eigenvalues sorted descending
so that $`\lambda_1`$ is the largest and $`v_1`$ its eigenvector.

The variance-component (Wald) term is the Mahalanobis quadratic form

``` math
\chi_1 = z'\Sigma^{-1}z,
```

which under the global null is $`\chi^2`$ with $`K`$ degrees of freedom.
The directional term is the standardized squared projection of $`z`$
onto the leading eigenvector,

``` math
\chi_2 = \frac{(v_1'z)^2}{\lambda_1},
```

which under the null is $`\chi^2`$ with $`1`$ degree of freedom and
captures association concentrated along the dominant correlation axis.

For each $`\rho`$ on the grid $`\{0,0.2,0.4,0.6,0.8,1\}`$ (6 equally
spaced points) the mixed statistic

``` math
Q_\rho = (1-\rho)\,\chi_1 + \rho\,\chi_2
```

is formed, with null tail probability $`p_\rho`$ computed exactly at the
endpoints ($`\chi^2_K`$ at $`\rho=0`$, $`\chi^2_1`$ at $`\rho=1`$) and
otherwise as the tail of a weighted sum of independent $`\chi^2_1`$
variables with spectrum $`(1-\rho)`$ of multiplicity $`K-1`$ plus a
single eigenvalue $`1`$ (Liu moment-matching approximation). The
statistic is the minimum p-value $`T=\min_\rho p_\rho`$, calibrated with
the SKAT-O one-dimensional null integral

``` math
p = \min_\rho p_\rho + \int_0^{t_{\max}} \Pr\Big[\chi^2_{K-1} >
  \min_{i:\,\rho_i<1}\frac{q_i - x(t)}{1-\rho_i}\Big]\,dt,
```

where $`q_i`$ is the mixed-null quantile at level $`1-\min_\rho p_\rho`$
for grid point $`i`$, $`x(t)`$ the $`\chi^2_1`$ quantile at $`t`$, and
$`t_{\max}`$ the upper limit from the $`\rho=1`$ quantile. The integral
is evaluated by composite Simpson’s rule and the result clamped to
$`[0,1]`$.

## What This Method Advances

EMATS adapts the SKAT-O idea by blending the full Wald component with a
leading-eigenvector directional component.

## Assumptions

Z-scores are approximately multivariate normal with mean zero and
covariance $`\Sigma`$ under the global null of no association.
$`\Sigma`$ is the (full-rank, invertible) trait correlation matrix,
shared across variants and treated as known, with a well-defined leading
eigenpair $`(\lambda_1, v_1)`$. Variants are tested one at a time.

## When This Method Is a Good Fit

Good when signal may align with the dominant trait-correlation axis.

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
pleio <- run_emats(pleio)
head(pleio@results$emats)
```

## Implementation Files

- R wrapper: `R/run_emats.R`
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

# Principal-component minimum-p test (PCMinP)

Function:
[`run_pcminp()`](https://broccolito.github.io/pleioverse/reference/run_pcminp.md)

Result column: `pcminp`

Family: Principal-component minimum-p

## Introduction

Multi-trait association test that takes the smallest
per-principal-component p-value and corrects it for the $`K`$
comparisons with a Sidak adjustment. It targets sparse pleiotropy in
which the signal concentrates in one (or a few) whitened directions of
the trait space.

## Mathematical Framework and Implementation

Let $`z=(z_1,\dots,z_K)'`$ be the per-trait Z-scores and let
$`\Sigma=V\Lambda V'`$ be the spectral decomposition of the trait
correlation matrix, with eigenvectors $`V=(v_1,\dots,v_K)`$ (engine
`evec`) and eigenvalues $`\lambda_i>0`$ (engine `eval`). The principal
components are $`u=V'z`$, $`u_i=v_i'z`$, and each standardized PC gives
a one-degree-of-freedom statistic and p-value

``` math
Q_i = \frac{u_i^2}{\lambda_i}\ \sim\ \chi^2_1,
       \qquad p_i = \Pr(\chi^2_1 > Q_i).
```

Let $`p_{\min}=\min_{i=1,\dots,K} p_i`$. Because the PCs are independent
under the null, the p-values $`p_i`$ are independent $`U(0,1)`$ and the
minimum is corrected by the exact Sidak formula

``` math
p = 1-\big(1-p_{\min}\big)^{K},
```

which is the reported p-value.

## What This Method Advances

PCMinP targets sparse signal after decorrelation by taking the smallest
principal-component p-value and adjusting for the number of components.

## Assumptions

Eigenvector whitening makes the $`K`$ per-PC chi-square statistics
mutually independent under the null $`z\sim N(0,\Sigma)`$, so each
$`p_i\sim U(0,1)`$ independently and the Sidak correction
$`1-(1-p_{\min})^{K}`$ is exact (not merely conservative as Bonferroni
would be). PCMinP is most powerful under sparse signal, where the
association loads on a single PC, and is less powerful than PCFisher or
PCLC when the signal is diffuse across many PCs. The trait-correlation
matrix $`\Sigma`$ is assumed known under the null and positive definite.

## When This Method Is a Good Fit

Good for sparse evidence that appears in one or a few principal
components.

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
pleio <- run_pcminp(pleio)
head(pleio@results$pcminp)
```

## Implementation Files

- R wrapper: `R/run_pcminp.R`
- C++ kernel or dispatcher: `src/pcminp.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu Z, Lin X (2019). A Geometric Perspective on the Power of Principal
Component Association Tests in Multiple Phenotype Studies. *Journal of
the American Statistical Association* 114(527):975–990. DOI:
<https://doi.org/10.1080/01621459.2018.1513363>

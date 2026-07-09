# Principal-component linear-combination test (PCLC)

Function:
[`run_pclc()`](https://broccolito.github.io/pleioverse/reference/run_pclc.md)

Result column: `pclc`

Family: Principal-component burden

## Introduction

Multi-trait association test formed from an inverse-eigenvalue weighted
linear combination of the principal components of the Z-vector. It
targets signal aligned with the implementation’s inverse-eigenvalue PC
contrast.

## Mathematical Framework and Implementation

Let $`z=(z_1,\dots,z_K)'`$ be the per-trait Z-scores at a variant and
let $`\Sigma=V\Lambda V'`$ be the spectral decomposition of the
$`K\times K`$ trait-correlation matrix, with eigenvectors
$`V=(v_1,\dots,v_K)`$ (engine `evec`) and eigenvalues $`\lambda_i>0`$
(engine `eval`). The principal components are

``` math
u = V'z, \qquad u_i = v_i'z,
```

and PCLC weights each PC by the inverse of its eigenvalue,
$`w_i = 1/\lambda_i`$, forming the scalar implemented in the C++ kernel

``` math
S = \sum_{i=1}^{K} w_i\, u_i.
```

Under the global null $`z\sim N(0,\Sigma)`$ the PCs are independent with
$`u_i\sim N(0,\lambda_i)`$, so $`S\sim N\!\big(0,\sum_i w_i\big)`$
because $`\mathrm{Var}(S)=\sum_i w_i^2\lambda_i=\sum_i w_i`$. The test
statistic is the squared standardized combination

``` math
T = \frac{S^2}{\sum_{i=1}^{K} w_i}\ \sim\ \chi^2_1,
```

and the reported p-value is $`p=\Pr(\chi^2_1 > T)`$.

## What This Method Advances

PCLC uses principal components to form a one-degree-of-freedom linear
combination after accounting for the component variances.

## Assumptions

Whitening by the eigenvectors of $`\Sigma`$ makes the PCs independent,
and the inverse-eigenvalue weights $`w_i=1/\lambda_i`$ give the
implemented one-degree-of-freedom contrast in PC space. This contrast is
not generally the same as the original-trait burden score
$`\mathbf{1}'\Sigma^{-1}z`$; recovering that score would require the
loadings of $`\mathbf{1}`$ in the eigenvector basis. PCLC concentrates
all power into one linear direction, so it is most powerful when the
trait effects align with this inverse-eigenvalue PC contrast, and it can
lose power when the signal is spread across PCs or has opposing signs
relative to the weighting direction. The trait-correlation matrix
$`\Sigma`$ is assumed known under the null and positive definite.
Eigenvector signs are not re-oriented before forming $`S`$, so the named
direction is basis/sign-convention dependent.

## When This Method Is a Good Fit

Good when the effect aligns with this implementation’s
one-degree-of-freedom PC contrast.

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
pleio <- run_pclc(pleio)
head(pleio@results$pclc)
```

## Implementation Files

- R wrapper: `R/run_pclc.R`
- C++ kernel or dispatcher: `src/pclc.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu Z, Lin X (2019). A Geometric Perspective on the Power of Principal
Component Association Tests in Multiple Phenotype Studies. *Journal of
the American Statistical Association* 114(527):975–990. DOI:
<https://doi.org/10.1080/01621459.2018.1513363>

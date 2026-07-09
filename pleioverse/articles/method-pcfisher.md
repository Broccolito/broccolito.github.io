# Principal-component Fisher-combination test (PCFisher)

Function:
[`run_pcfisher()`](https://broccolito.github.io/pleioverse/reference/run_pcfisher.md)

Result column: `pcfisher`

Family: Principal-component combination

## Introduction

Multi-trait association test that projects the per-trait Z-scores onto
the principal components of the trait-correlation matrix and combines
the resulting per-PC p-values with Fisher’s method. It targets
pleiotropic signal that is spread across several independent (whitened)
directions of the trait space.

## Mathematical Framework and Implementation

Let $`z=(z_1,\dots,z_K)'`$ be the vector of per-trait Z-scores at a
variant, $`z_k=\hat\beta_k/\mathrm{se}_k`$, and let $`\Sigma`$ be the
$`K\times K`$ between-trait correlation matrix with spectral
decomposition

``` math
\Sigma = V \Lambda V', \qquad \Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_K),
```

where $`V=(v_1,\dots,v_K)`$ holds the eigenvectors (engine `evec`) and
$`\lambda_i>0`$ the eigenvalues (engine `eval`). The principal
components of the Z-vector are

``` math
u = V'z, \qquad u_i = v_i' z.
```

Under the global null $`z\sim N(0,\Sigma)`$ the PCs are independent with
$`u_i\sim N(0,\lambda_i)`$, so each standardized PC yields a chi-square
statistic on one degree of freedom and its upper-tail p-value

``` math
Q_i = \frac{u_i^2}{\lambda_i}\ \sim\ \chi^2_1,
       \qquad p_i = \Pr(\chi^2_1 > Q_i).
```

Fisher’s combination of the $`K`$ independent p-values gives the
statistic

``` math
T = -2\sum_{i=1}^{K}\log p_i\ \sim\ \chi^2_{2K},
```

and the reported p-value is $`p = \Pr(\chi^2_{2K} > T)`$.

## What This Method Advances

PCFisher decorrelates the trait scores into principal components and
then combines component-level evidence with Fisher’s rule.

## Assumptions

Whitening the Z-vector by the eigenvectors of $`\Sigma`$ renders the PCs
mutually independent under the null, which is what makes the Fisher
combination of the per-PC chi-square p-values valid and yields the exact
$`\chi^2_{2K}`$ null. Because every PC enters symmetrically on the log
scale, PCFisher is most powerful when the association signal is
distributed across multiple PCs (diffuse pleiotropy) and tends to be
less powerful than linear PC combinations when the signal concentrates
in a single direction. The trait-correlation matrix $`\Sigma`$ is
assumed known (estimated under the null) and positive definite.

## When This Method Is a Good Fit

Good when multiple independent principal-component directions may carry
evidence.

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
pleio <- run_pcfisher(pleio)
head(pleio@results$pcfisher)
```

## Implementation Files

- R wrapper: `R/run_pcfisher.R`
- C++ kernel or dispatcher: `src/pcfisher.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu Z, Lin X (2019). A Geometric Perspective on the Power of Principal
Component Association Tests in Multiple Phenotype Studies. *Journal of
the American Statistical Association* 114(527):975–990. DOI:
<https://doi.org/10.1080/01621459.2018.1513363>

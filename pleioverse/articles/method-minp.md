# Minimum-p test with multivariate-normal correction (MinP)

Function:
[`run_minp()`](https://broccolito.github.io/pleioverse/reference/run_minp.md)

Result column: `minp`

Family: Trait-level minimum-p

## Introduction

Multi-trait association test that takes the smallest per-trait p-value
and corrects it for the correlated traits using a numerical
multivariate-normal calculation under the null of the Z-scores (a la
Conneely & Boehnke). It is aimed at association signal that shows up
strongly in at least one of the original traits.

## Mathematical Framework and Implementation

Let $`z=(z_1,\dots,z_K)'`$ be the per-trait Z-scores at a variant. Each
trait contributes a two-sided single-trait test

``` math
Q_k = z_k^2\ \sim\ \chi^2_1, \qquad p_k = \Pr(\chi^2_1 > Q_k),
```

and the observed minimum is $`p_{\min}=\min_{k=1,\dots,K} p_k`$. Because
the original traits are correlated, this method uses the joint
multivariate normal probability rather than an independence-based Sidak
correction or a conservative Bonferroni bound. Under the global null
$`z\sim N(0,\Sigma)`$ the corrected p-value is the probability that at
least one trait reaches the observed threshold. Writing the per-trait
threshold as $`c = \Phi^{-1}\!\big(1-p_{\min}/2\big)`$ (so that
$`\Pr(|z_k|>c)=p_{\min}`$), the correlation-adjusted corrected p-value
is

``` math
p = 1 - \Pr\!\Big(\bigcap_{k=1}^{K}\{-c \le z_k \le c\}\Big)
        = 1 - \Pr_{N(0,\Sigma)}\big([-c,c]^K\big).
```

The rectangle (box) probability of the multivariate normal with
correlation $`\Sigma`$ is evaluated by quasi-randomized numerical
integration using the lower Cholesky factor $`L`$ of $`\Sigma`$
($`\Sigma=LL'`$, engine `cholSigma`), and $`p`$ is reported.

## What This Method Advances

MinP targets sparse original-trait signal and uses a multivariate-normal
probability calculation to account for correlation among trait tests.

## Assumptions

The null distribution of the Z-vector is exactly multivariate normal
$`N(0,\Sigma)`$ with known trait-correlation matrix $`\Sigma`$. The
correction integrates the joint MVN density over the acceptance box
$`[-c,c]^K`$, which accounts for the specified trait correlation up to
the numerical integration tolerance. MinP operates on the original
traits rather than whitened PCs, so it is most powerful when the
association is large in one or a few individual traits. The box
probability is computed by Monte-Carlo/quasi-Monte-Carlo numerical
integration, so the result carries a small, controllable integration
error.

## When This Method Is a Good Fit

Good for sparse evidence concentrated in one original trait.

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
pleio <- run_minp(pleio)
head(pleio@results$minp)
```

## Implementation Files

- R wrapper: `R/run_minp.R`
- C++ kernel or dispatcher: `src/minp.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Conneely KN, Boehnke M (2007). So Many Correlated Tests, So Little Time!
Rapid Adjustment of P Values for Multiple Correlated Tests. *American
Journal of Human Genetics* 81(6):1158–1168. DOI:
<https://doi.org/10.1086/522036>

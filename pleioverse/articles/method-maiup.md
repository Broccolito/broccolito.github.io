# MAIUP empirical-FDR pleiotropy screen on a phenotype pair

Function:
[`run_maiup()`](https://broccolito.github.io/pleioverse/reference/run_maiup.md)

Result column: `maiup`

Family: Pairwise mixture-adjusted IUT

## Introduction

Pairwise (traits $`i,j`$) Mixture-model-Adjusted Intersection-Union
empirical FDR procedure (Zeng et al., MAIUP). It decorrelates the
per-trait Z-scores, forms the intersection-union statistic, and
calibrates an empirical FDR against a multi-component null mixture (Zeng
et al.). MAIUP is run for one trait pair at a time.

## Mathematical Framework and Implementation

**Decorrelating whitening.** Let $`Z = (Z_1, Z_2)`$ be the per-trait
Z-scores. The null cross-trait correlation matrix $`R`$ is estimated
from variants significant in *neither* trait ($`p_1 \ge \tau`$ and
$`p_2 \ge \tau`$, with $`\tau`$ = `p_threshold`). With
eigendecomposition $`R = V \Lambda V^{\top}`$, the symmetric (ZCA-type)
whitening matrix $`W = V \Lambda^{-1/2} V^{\top}`$ produces decorrelated
scores $`Z^{d} = Z\,W`$ with components $`Z^{d}_1, Z^{d}_2`$ and
two-sided per-trait p-values

``` math
P^{d}_k = 2\,\bigl(1 - \Phi(|Z^{d}_k|)\bigr), \qquad k = 1,2 .
```

**Intersection-union statistic.** The IUT statistic for the composite
null “associated with at most one trait” is the maximum of the two
decorrelated p-values,

``` math
P_{\max} = \max\bigl(P^{d}_1, P^{d}_2\bigr),
```

which rejects only when *both* traits are individually significant.

**Null mixture.** Variants fall into four classes by which traits are
null: both null ($`a_{00}`$), only trait 2 non-null ($`a_{01}`$), only
trait 1 non-null ($`a_{10}`$), and both non-null (the alternative). The
null proportions are estimated by Storey-type tail counting at cut
$`c = 0.5`$:

``` math
\hat a_{1} = \frac{\#\{P^{d}_1 \ge c\}}{n\,(1-c)}, \quad
      \hat a_{2} = \frac{\#\{P^{d}_2 \ge c\}}{n\,(1-c)}, \quad
      \hat a_{00} = \frac{\#\{P^{d}_1 \ge c,\, P^{d}_2 \ge c\}}{n\,(1-c)^2},
```

followed by $`\hat a_{01} = \max(0,\hat a_1 - \hat a_{00})`$ and
$`\hat a_{10} = \max(0,\hat a_2 - \hat a_{00})`$ (with edge corrections,
and a one-sided Kolmogorov-Smirnov uniformity check that forces
$`\hat a_k = 1`$ when a trait’s p-values are indistinguishable from
Uniform$`(0,1)`$).

**Marginal null CDFs.** For each trait the null CDF of $`P^{d}_k`$ is
estimated nonparametrically by the Least Concave Majorant (LCM) of the
empirical CDF, then mixture-adjusted to isolate the null component:
$`\hat F_k(p) = \bigl(F^{\mathrm{LCM}}_k(p) - a_k\,p\bigr)/(1-a_k)`$.

**Empirical FDR.** Writing $`\hat F(p)`$ for the empirical CDF of
$`P_{\max}`$, the false-discovery numerator sums the contributions of
the three null sub-classes evaluated at $`p = P_{\max}`$:

``` math
\widehat{\mathrm{FDR}}(p) = \frac{\hat a_{01}\, p\, \hat F_2(p)
       + \hat a_{10}\, p\, \hat F_1(p) + \hat a_{00}\, p^{2}}{\hat F(p)} .
```

(Under each sub-null the max-statistic tail factorizes into a uniform
factor $`p`$ for the null trait and the estimated null CDF for the
other, with $`p^2`$ when both are null.) The returned `p_maiup` is this
FDR after monotone (step-down) enforcement in decreasing $`P_{\max}`$
and clamping to $`[0,1]`$.

## What This Method Advances

MAIUP combines an intersection-union statistic with mixture adjustment,
so both traits must show evidence while one-trait null components are
modeled explicitly.

## Assumptions

- Cross-trait dependence among null Z-scores is fully captured by the
  linear correlation $`R`$, so symmetric whitening renders the two
  traits’ null statistics independent.
- The genome-wide p-value mixture has the four-class structure above,
  with each null sub-class contributing a uniform marginal for its null
  trait(s).
- LCM/isotonic estimation requires the true null CDF to be concave
  (p-value density non-increasing), as for well-behaved test statistics.
- Null proportions are estimable from the upper half ($`p \ge 0.5`$) of
  the p-value distribution (signal sparse near 1).
- Variants are treated as exchangeable; LD is not explicitly modelled.

## When This Method Is a Good Fit

Good for pairwise pleiotropy when empirical FDR calibration and mixture
adjustment are desired.

The `maiup` column is the monotone empirical FDR output returned by the
MAIUP implementation, built from the IUT statistic and mixture-adjusted
null components.

Pairwise methods use `i` and `j` to select the two phenotypes. Re-run
the method for each phenotype pair that should be tested.

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
pleio <- run_maiup(pleio, i = 1, j = 2)
head(pleio@results$maiup)
```

## Implementation Files

- R wrapper: `R/run_maiup.R`
- C++ kernel or dispatcher: `src/maiup.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Zeng P, et al. MAIUP: a mixture-adjusted intersection-union test for
pleiotropy. [MAIUP repository](https://github.com/biostatpzeng/MAIUP)

Kessy A, Lewin A, Strimmer K (2018). Optimal whitening and
decorrelation. *The American Statistician* 72(4):309-314. DOI:
<https://doi.org/10.1080/00031305.2016.1277159>

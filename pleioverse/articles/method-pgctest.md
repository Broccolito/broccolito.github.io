# PGCtest (CompTestER) normal-product test on a phenotype pair

Function:
[`run_pgctest()`](https://broccolito.github.io/pleioverse/reference/run_pgctest.md)

Result column: `pgctest`

Family: Pairwise product-Z composite test

## Introduction

Pairwise (traits $`i,j`$) pleiotropy screen based on the
*normal-product* statistic $`|Z_1 Z_2|`$, with a Bessel-$`K_0`$ tail CDF
and an empirical variance normalization. In this package, the
implemented routine applies a product-Z tail calibration with one
empirical scale; it is evaluated for one trait pair at a time.

## Mathematical Framework and Implementation

For variant $`v`$ let $`Z_1, Z_2`$ be the per-trait Z-scores
($`\hat\beta/\widehat{\mathrm{se}}`$) and let $`X = |Z_1 Z_2|`$. The
implemented test rejects for large product magnitude $`X`$ after
calibration against the normal-product reference distribution; it does
not explicitly evaluate a multi-component composite-null mixture.

**Standard normal-product tail.** When $`Z_1, Z_2`$ are independent
standard normal, the density of $`X = |Z_1 Z_2|`$ is

``` math
f(x) = \frac{2}{\pi}\, K_0(x), \qquad x \ge 0,
```

with $`K_0`$ the modified Bessel function of the second kind of order 0.
The implementation builds the corresponding survival (upper-tail)
function once,

``` math
S(x) = \Pr(X \ge x) = \int_x^{\infty} \frac{2}{\pi} K_0(u)\, du,
```

by evaluating $`f`$ on a degree-10 power-spaced grid
($`x_i = x_{\max}\,(i/(\mathrm{len}-1))^{10}`$, dense near 0) and
accumulating a reverse-cumulative trapezoidal integral, normalized so
that $`S(0)=1`$; the grid upper bound is $`x_{\max}=2\max_v X_v`$.

**Empirical variance normalization.** Because the per-trait Z-scores
need not be exactly standard normal (inflation, residual cross-trait
dependence), the observed statistic is rescaled before being passed to
$`S(\cdot)`$. With sample variances $`\widehat{\mathrm{Var}}(Z_1)`$ and
$`\widehat{\mathrm{Var}}(Z_2)`$ over all variants, define

``` math
d = \sqrt{\widehat{\mathrm{Var}}(Z_1) + \widehat{\mathrm{Var}}(Z_2) - 1},
```

and report the p-value by linear interpolation of the tail CDF at the
rescaled point

``` math
p = S\!\left(\frac{|Z_1 Z_2|}{d}\right).
```

A signed pleiotropy z is also computed internally as
$`z = \mathrm{sign}(Z_1 Z_2)\,\Phi^{-1}(1 - p/2)`$; the R wrapper
returns the p-value component (`$pp`).

## What This Method Advances

PGCtest uses the product of two trait Z-scores and calibrates that
product under a normal-product reference distribution.

## Assumptions

- Under the relevant null the per-trait Z-scores are approximately
  standard normal, so that $`|Z_1 Z_2|`$ follows the
  $`\tfrac{2}{\pi}K_0`$ normal-product law.
- Mild departures (variance inflation, weak dependence) are corrected by
  the single empirical scale $`d`$; the correction assumes the inflation
  is well summarized by the marginal Z variances.
- Independence of $`Z_1, Z_2`$ is assumed by the base tail; strong
  cross-trait correlation is not separately modelled (unlike PLACO,
  which carries an explicit $`r`$).
- Variants are treated as exchangeable for the variance estimates.

## When This Method Is a Good Fit

Good for pairwise product-Z evidence when product-tail calibration is
appropriate.

The implemented PGCtest wrapper returns only the `pp` values from the
C++ routine. The C++ routine also computes signed `zz` scores
internally, but they are not attached to `pleio@results`. The
product-tail calibration uses only the two Z-score vectors and their
marginal sample variances; it does not explicitly model cross-trait
correlation and relies on a positive variance denominator.

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
pleio <- run_pgctest(pleio, i = 1, j = 2)
head(pleio@results$pgctest)
```

## Implementation Files

- R wrapper: `R/run_pgctest.R`
- C++ kernel or dispatcher: `src/pgctest.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Andreassen OA, Thompson WK, Schork AJ, et al. (2013). Improved detection
of common variants associated with schizophrenia and bipolar disorder
using pleiotropy-informed conditional false discovery rate. *PLoS
Genetics* 9(4):e1003455. DOI:
<https://doi.org/10.1371/journal.pgen.1003455>

Ray D, Chatterjee N (2020). A powerful method for pleiotropic analysis
under composite null hypothesis. *PLoS Genetics* 16(12):e1009218. DOI:
<https://doi.org/10.1371/journal.pgen.1009218>

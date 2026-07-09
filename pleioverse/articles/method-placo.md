# PLACO+ pleiotropy test on a phenotype pair

Function:
[`run_placo()`](https://broccolito.github.io/pleioverse/reference/run_placo.md)

Result column: `placo`

Family: Pairwise pleiotropy composite-null test

## Introduction

Pairwise (traits $`i,j`$) test of pleiotropy under the *composite null*
that a variant is associated with at most one of the two traits, using
the product of per-trait Z-scores as the test statistic (Ray and
Chatterjee, 2020). PLACO is applied to the trait pair only and must be
run separately for each pair of interest.

## Mathematical Framework and Implementation

For variant $`v`$ let $`Z_1, Z_2`$ be the per-trait Z-scores
($`\hat\beta/\widehat{\mathrm{se}}`$). PLACO tests the composite null

``` math
H_0:\; \beta_1 = 0 \;\mathrm{or}\; \beta_2 = 0,
```

i.e. the union of the three sub-null configurations $`(0,0)`$,
$`(0,\ne 0)`$, $`(\ne 0,0)`$, against the pleiotropic alternative
$`\beta_1 \ne 0`$*and* $`\beta_2 \ne 0`$. The statistic is the product

``` math
T = Z_1 Z_2 .
```

**Null density of the product.** Under the null each $`Z_k`$ is
approximately normal with mean 0, and $`(Z_1,Z_2)`$ may be correlated
(e.g. through shared controls). Writing $`r`$ for that correlation, the
density of the standardized product $`x = Z_1 Z_2`$ is the correlated
normal-product density

``` math
f(x \mid r) = \frac{1}{\pi \sqrt{1-r^2}} \,
  \exp\!\left(\frac{r\,x}{1-r^2}\right) \,
  K_0\!\left(\frac{|x|}{1-r^2}\right),
```

where $`K_0`$ is the modified Bessel function of the second kind of
order 0. The implementation evaluates $`f(\cdot \mid r)`$ on a
power-spaced grid that is dense near 0 (node
$`x_i = X_{\max}\,(i/(N-1))^{8}`$) and forms the upper- and lower-tail
integrals by reverse-cumulative trapezoidal quadrature, giving a
two-sided tail function

``` math
G(t) = \int_{|x| \ge t} f(x \mid r)\, dx .
```

**Variance and correlation estimation.** The null moments are estimated
off the variants that are *not* jointly significant. The per-trait null
variances $`s_1^2, s_2^2`$ are the sample variances of $`Z_1, Z_2`$ over
variants *not* significant in both traits (i.e. excluding $`p_1 < \tau`$
and $`p_2 < \tau`$), where $`\tau`$ is `p_threshold`. The null
correlation $`r`$ is the Pearson correlation of $`(Z_1,Z_2)`$ over
variants significant in *neither* trait ($`p_1 \ge \tau`$ and
$`p_2 \ge \tau`$).

**Composite p-value.** Because the null is composite, the tail of the
product alone is anti-conservative when only one trait is non-null.
Following the PLACO construction the p-value uses an inclusion-exclusion
combination of the two-sided product tail evaluated at the observed
$`t_0 = |Z_1 Z_2|`$ rescaled by each per-trait null SD:

``` math
p = G\!\left(\frac{t_0}{s_1}\right) + G\!\left(\frac{t_0}{s_2}\right)
          - G\!\left(t_0\right),
```

clamped to $`[10^{-300}, 1]`$. Small $`p`$ indicates evidence that the
variant is associated with *both* traits.

## What This Method Advances

PLACO+ targets the composite pleiotropy null directly: association with
neither trait or with only one trait should not be called pleiotropic.

## Assumptions

- Per-trait Z-scores are (approximately) standard normal under the
  respective trait null; departures are absorbed into the empirically
  estimated null variances $`s_1^2, s_2^2`$.
- Cross-trait dependence is captured by a single correlation $`r`$
  estimated from variants null in both traits; estimation is meaningful
  only if such variants dominate (sparse true pleiotropy).
- The product null density is the correlated normal-product law, so the
  bivariate Z is jointly normal with zero mean under each sub-null.
- Variants are treated as exchangeable for moment estimation (LD and
  sample-size heterogeneity are not modelled here).

## When This Method Is a Good Fit

Good for pairwise pleiotropy under the composite null of association
with at most one trait.

The implementation treats variants independently and relies on the
trait-correlation estimate derived from the available summary statistics
unless an individual method estimates a pairwise null moment internally.

The empirical variance and null-correlation estimates need enough
variants in the relevant non-joint and null/null-like sets. Extremely
small simulations, dense signal, or non-finite Z-scores can make those
moment estimates unstable.

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
pleio <- run_placo(pleio, i = 1, j = 2)
head(pleio@results$placo)
```

## Implementation Files

- R wrapper: `R/run_placo.R`
- C++ kernel or dispatcher: `src/placo.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Ray D, Chatterjee N (2020). A powerful method for pleiotropic analysis
under composite null hypothesis identifies novel shared loci between
Type 2 Diabetes and Prostate Cancer. *PLoS Genetics* 16(12):e1009218.
DOI: <https://doi.org/10.1371/journal.pgen.1009218>

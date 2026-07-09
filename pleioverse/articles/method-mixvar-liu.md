# Mixture variance-component omnibus (Liu tail)

Function:
[`run_mixvar_liu()`](https://broccolito.github.io/pleioverse/reference/run_mixvar_liu.md)

Result column: `mixvar_liu`

Family: Adaptive mixture quadratic form

## Introduction

MPAT/MTAR multi-trait association test that blends the burden/mean and
heterogeneity/variance components into a single quadratic form, with the
mixing weight chosen to balance the two components on the variance
scale. This variant evaluates the weighted-chi-square tail by the Liu
chi-square moment-matching approximation.

## Mathematical Framework and Implementation

For a variant the input is the vector of trait Z-scores
$`z = (z_1,\dots,z_K)'`$ with
$`z_k = \hat\beta_k/\widehat{\mathrm{se}}_k`$, which under the null is
$`z \sim \mathcal{N}(0,\Sigma)`$ for the $`K\times K`$ trait correlation
$`\Sigma`$. Write $`\Sigma^{-1}`$ for its inverse, $`\mathbf{1}`$ for
the all-ones vector, $`a = \Sigma^{-1}\mathbf{1}`$,
$`s = \mathbf{1}'\Sigma^{-1}\mathbf{1}`$, and
$`H = \mathbf{1}\mathbf{1}'/K`$.

**Component kernels.** The burden/mean direction is carried by the
rank-one kernel $`a a'`$ (so that $`z' a a' z = (a'z)^2`$), and the
heterogeneity/variance direction by the centered precision kernel

``` math
M_\tau = (I-H)\,\Sigma^{-1}\Sigma^{-1}\,(I-H).
```

**Variance-scale mixing.** Rather than combining two p-values, this test
forms a single kernel $`L_\phi`$ that is a convex blend of the two,

``` math
L_\phi = \phi\, a a' + (1-\phi)\, M_\tau, \qquad T = z' L_\phi z .
```

The weight $`\phi`$ is fixed from the two components’ null variance
scales:

``` math
\phi = \frac{V_\tau}{V_\mu + V_\tau}, \quad
  V_\mu = 2\,s^2, \quad V_\tau = 2\,\mathrm{tr}(M_\tau\Sigma M_\tau\Sigma),
```

where $`V_\mu`$ and $`V_\tau`$ are the null variances of $`(a'z)^2`$ and
$`z'M_\tau z`$. Under the null $`T`$ is a weighted sum of independent
$`\chi^2_1`$ variables, $`T \sim \sum_j \lambda_j \chi^2_{1,j}`$, with
weights $`\lambda_j`$ the eigenvalues of
$`\Sigma^{1/2} L_\phi \Sigma^{1/2}`$; the reported p-value is the upper
tail of this mixture.

**Tail method (Liu).** The tail is computed by the Liu (Liu-Tang-Zhang)
approximation, which matches the first three cumulants of $`T`$ to a
single (possibly non-central) chi-square. It is fast and closed-form but
only approximate, and can lose accuracy far in the tail; use the Davies
variant when extreme-tail precision is required.

## What This Method Advances

The Liu variant of MixVar uses moment matching for the blended
quadratic-form tail.

## Assumptions

Variance-scale mixing aims to make the test robust across
genetic-architecture regimes without choosing between them: the $`a a'`$
part captures concordant shared effects while the $`M_\tau`$ part
captures discordant, antagonistic, or sparse effects, and the
variance-scale $`\phi`$ aims to keep neither part dominating the null
scale. Validity requires $`\Sigma`$ to be a correct, full-rank trait
correlation matrix and the per-trait Z-scores to be approximately
Gaussian under the null.

## When This Method Is a Good Fit

Good when a fast approximate variance-scale blend is desired.

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
pleio <- run_mixvar_liu(pleio)
head(pleio@results$mixvar_liu)
```

## Implementation Files

- R wrapper: `R/run_mixvar_liu.R`
- C++ kernel or dispatcher: `src/mixvar.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu, Z., Lin, X. (2018) Multiple Phenotype Association Tests Using
Summary Statistics in Genome-Wide Association Studies. *Biometrics*
74(1):165-175. DOI: <https://doi.org/10.1111/biom.12735>

Luo, L., Shen, J., Zhang, H., et al. (2020) Multi-trait analysis of
rare-variant association summary statistics using MTAR. *Nature
Communications* 11:2850. DOI:
<https://doi.org/10.1038/s41467-020-16591-0>

Liu, H., Tang, Y., Zhang, H. H. (2009) A new chi-square approximation to
the distribution of non-negative definite quadratic forms in non-central
normal variables. *Computational Statistics & Data Analysis* 53,
853-856. DOI: <https://doi.org/10.1016/j.csda.2008.11.025>

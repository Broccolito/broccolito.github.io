# Mixture-Fisher omnibus (modified-Liu tail)

Function:
[`run_mixfisher_liumod()`](https://broccolito.github.io/pleioverse/reference/run_mixfisher_liumod.md)

Result column: `mixfisher_liumod`

Family: Two-component mixture combination

## Introduction

MPAT/MTAR multi-trait association test that decomposes the per-variant
signal into a burden/mean component and a heterogeneity/variance
component and combines the two component p-values by Fisher’s method.
This variant evaluates the heterogeneity tail by the modified-Liu
chi-square moment-matching approximation.

## Mathematical Framework and Implementation

For a variant the input is the vector of trait Z-scores
$`z = (z_1,\dots,z_K)'`$ with
$`z_k = \hat\beta_k/\widehat{\mathrm{se}}_k`$, which under the null is
$`z \sim \mathcal{N}(0,\Sigma)`$ for the $`K\times K`$ trait correlation
$`\Sigma`$. Write $`\Sigma^{-1}`$ for its inverse, $`\mathbf{1}`$ for
the all-ones vector, $`a = \Sigma^{-1}\mathbf{1}`$,
$`s = \mathbf{1}'\Sigma^{-1}\mathbf{1}`$, and the centering projection
$`H = \mathbf{1}\mathbf{1}'/K`$.

**Burden (mean) component.** The generalized-least-squares burden score
tests a shared, concordant mean across traits,

``` math
U_\mu = \frac{(a'z)^2}{s}, \qquad p_\mu = \Pr(\chi^2_1 \ge U_\mu),
```

i.e. $`U_\mu`$ is exactly $`\chi^2_1`$ under the null.

**Heterogeneity (variance) component.** The complementary quadratic form
isolates departures from a common mean using the projected precision
matrix

``` math
M_\tau = (I-H)\,\Sigma^{-1}\Sigma^{-1}\,(I-H), \qquad
      U_\tau = z' M_\tau z .
```

Under the null $`U_\tau`$ is a weighted sum of independent $`\chi^2_1`$
variables, $`U_\tau \sim \sum_j \lambda_j \chi^2_{1,j}`$, whose weights
$`\lambda_j`$ are the eigenvalues of
$`\Sigma^{1/2} M_\tau \Sigma^{1/2}`$. Its p-value $`p_\tau`$ is the
upper tail of this mixture.

**Fisher combination.** The two components are asymptotically
independent ($`U_\mu`$ lives in the $`\mathbf{1}`$ direction, $`U_\tau`$
in its orthogonal complement), so they are pooled by Fisher’s rule

``` math
T = -2\big(\log p_\mu + \log p_\tau\big) \;\sim\; \chi^2_4,
```

and the reported p-value is $`\Pr(\chi^2_4 \ge T)`$.

**Tail method (modified Liu).** The heterogeneity tail $`p_\tau`$ is
computed by the modified-Liu approximation, which matches all four
cumulants (including kurtosis) of $`U_\tau`$ to a non-central
chi-square. It is as fast as the plain Liu method and typically more
accurate, though still an approximation; the Davies variant remains the
reference for extreme tails.

## What This Method Advances

The modified-Liu variant uses the same mixture-Fisher construction with
the package’s modified Liu moment-matching tail for the heterogeneity
component.

## Assumptions

The omnibus is designed to be robust across genetic-architecture
regimes: the burden component is most powerful when effects are
concordant in sign and magnitude across traits, while the heterogeneity
component captures discordant, antagonistic, or sparse effects to which
the burden is blind. Fisher combination preserves sensitivity to either
regime without prior knowledge of which holds. Validity requires
$`\Sigma`$ to be a correct, full-rank trait correlation matrix and the
per-trait Z-scores to be approximately Gaussian under the null.

## When This Method Is a Good Fit

Good for the same two-component question with the modified Liu
approximation used by this package.

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
pleio <- run_mixfisher_liumod(pleio)
head(pleio@results$mixfisher_liumod)
```

## Implementation Files

- R wrapper: `R/run_mixfisher_liumod.R`
- C++ kernel or dispatcher: `src/mixfisher.cpp`
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

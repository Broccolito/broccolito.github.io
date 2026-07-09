# Mixture-Tippett (minimum-p) omnibus (Davies tail)

Function:
[`run_mixtippett_davies()`](https://broccolito.github.io/pleioverse/reference/run_mixtippett_davies.md)

Result column: `mixtippett_davies`

Family: Two-component mixture combination

## Introduction

MPAT/MTAR multi-trait association test that decomposes the per-variant
signal into a burden/mean component and a heterogeneity/variance
component and combines the two component p-values by Tippett’s minimum-p
rule. This variant evaluates the heterogeneity tail by the
Davies-labeled (Imhof-style) numerical inversion.

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

**Tippett combination.** Tippett’s rule takes the smaller of the two
component p-values and calibrates it against the null of two independent
tests,

``` math
m = \min(p_\mu, p_\tau), \qquad p = 1 - (1 - m)^2 .
```

Equivalently $`m`$ is referred to a Beta(1,2) reference; the reported
p-value is $`p`$. Tippett emphasizes the single strongest component and
is typically more powerful than Fisher when only one of the two regimes
is active.

**Tail method (Davies).** The heterogeneity tail $`p_\tau`$ is computed
by the Davies-labeled numerical inversion of the characteristic function
(Imhof numerical integration). It is the direct numerical-inversion
option among the three tail variants, at higher computational cost; if
the returned tail probability is outside (0, 1), the implementation
falls back to the modified-Liu approximation.

## What This Method Advances

The mixture-Tippett family combines mean-effect and heterogeneity
evidence by the smaller component p-value rather than Fisher’s sum of
logs.

## Assumptions

The omnibus is designed to be robust across genetic-architecture
regimes: the burden component is most powerful when effects are
concordant in sign and magnitude across traits, while the heterogeneity
component captures discordant, antagonistic, or sparse effects to which
the burden is blind. Tippett combination preserves sensitivity to either
regime without prior knowledge of which holds. Validity requires
$`\Sigma`$ to be a correct, full-rank trait correlation matrix and the
per-trait Z-scores to be approximately Gaussian under the null.

## When This Method Is a Good Fit

Good when the strongest of mean-effect or heterogeneity evidence should
drive the result.

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
pleio <- run_mixtippett_davies(pleio)
head(pleio@results$mixtippett_davies)
```

## Implementation Files

- R wrapper: `R/run_mixtippett_davies.R`
- C++ kernel or dispatcher: `src/mixtippett.cpp`
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

Davies, R. B. (1980) The distribution of a linear combination of
$`\chi^2`$ random variables. *Applied Statistics* 29, 323-333. DOI:
<https://doi.org/10.2307/2346911>

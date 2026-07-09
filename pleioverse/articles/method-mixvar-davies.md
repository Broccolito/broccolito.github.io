# Mixture variance-component omnibus (Davies tail)

Function:
[`run_mixvar_davies()`](https://broccolito.github.io/pleioverse/reference/run_mixvar_davies.md)

Result column: `mixvar_davies`

Family: Adaptive mixture quadratic form

## Introduction

MPAT/MTAR multi-trait association test that blends the burden/mean and
heterogeneity/variance components into a single quadratic form, with the
mixing weight chosen to balance the two components on the variance
scale. This variant evaluates the weighted-chi-square tail by the
Davies-labeled (Imhof-style) numerical inversion.

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

**Tail method (Davies).** The tail is computed by the Davies-labeled
numerical inversion of the characteristic function (Imhof numerical
integration). It is the direct numerical-inversion option among the
three tail variants, at higher computational cost; if the returned tail
probability is outside (0, 1), the implementation falls back to the
modified-Liu approximation.

## What This Method Advances

MixVar fixes the burden/heterogeneity blend using variance-based
weighting and evaluates the resulting quadratic form with the selected
tail engine.

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

Good when a fixed variance-scale blend of burden and heterogeneity is
desired.

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
pleio <- run_mixvar_davies(pleio)
head(pleio@results$mixvar_davies)
```

## Implementation Files

- R wrapper: `R/run_mixvar_davies.R`
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

Davies, R. B. (1980) The distribution of a linear combination of
$`\chi^2`$ random variables. *Applied Statistics* 29, 323-333. DOI:
<https://doi.org/10.2307/2346911>

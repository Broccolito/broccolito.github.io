# Adaptive mixture omnibus (mixAda)

Function:
[`run_mixada()`](https://broccolito.github.io/pleioverse/reference/run_mixada.md)

Result column: `mixada`

Family: Adaptive mixture omnibus

## Introduction

MPAT/MTAR multi-trait association test that adaptively searches a grid
of burden/mean vs. heterogeneity/variance mixing weights, takes the most
significant weighted quadratic form, and calibrates the resulting
minimum p-value against its analytic null by one-dimensional numerical
integration.

## Mathematical Framework and Implementation

For a variant the input is the vector of trait Z-scores
$`z = (z_1,\dots,z_K)'`$ with
$`z_k = \hat\beta_k/\widehat{\mathrm{se}}_k`$, which under the null is
$`z \sim \mathcal{N}(0,\Sigma)`$ for the $`K\times K`$ trait correlation
$`\Sigma`$. Write $`\Sigma^{-1}`$ for its inverse, $`\mathbf{1}`$ for
the all-ones vector, $`a = \Sigma^{-1}\mathbf{1}`$,
$`s = \mathbf{1}'\Sigma^{-1}\mathbf{1}`$, and
$`H = \mathbf{1}\mathbf{1}'/K`$. As in the other mixture tests, the
burden/mean direction is carried by the rank-one kernel $`a a'`$ and the
heterogeneity/variance direction by

``` math
M_\tau = (I-H)\,\Sigma^{-1}\Sigma^{-1}\,(I-H).
```

**Adaptive mixing grid.** Instead of fixing a single mixing weight, the
test sweeps a grid $`\phi_1,\dots,\phi_B`$ ($`B = 11`$ equally spaced
points on $`[0.01, 0.99]`$). For each grid point it forms the blended
kernel and quadratic form

``` math
L_{\phi_b} = \phi_b\, a a' + (1-\phi_b)\, M_\tau, \qquad
      T_b = z' L_{\phi_b} z .
```

Under the null each $`T_b`$ is a weighted sum of $`\chi^2_1`$ with
weights $`\lambda^{(b)}_j`$ the eigenvalues of
$`\Sigma^{1/2} L_{\phi_b} \Sigma^{1/2}`$; its p-value $`p_b`$ is
obtained from the modified-Liu four-cumulant chi-square approximation.
The adaptive statistic is the minimum over the grid,

``` math
P_{\min} = \min_{1\le b\le B} p_b .
```

**Analytic null calibration.** Because the $`T_b`$ are strongly
dependent (they share $`z`$), $`P_{\min}`$ is not uniform and must be
recalibrated. Conditioning on the shared burden score and writing each
$`T_b`$ as a $`\phi_b`$-weighted combination of a common burden term and
a common heterogeneity term (whose null follows a moment-matched
chi-square with degrees of freedom and scale derived from the
eigenvalues of $`\Sigma^{1/2} M_\tau \Sigma^{1/2}`$), the p-value of
$`P_{\min}`$ reduces to a one-dimensional numerical integral over the
burden coordinate. With the substitution $`t = \Pr(\chi^2_1 \le x)`$
(which removes the $`\chi^2_1`$ density singularity at the origin), this
integral is evaluated by composite Simpson quadrature, giving

``` math
p = 1 - \int_0^{t_{\max}} \Pr\!\big(\chi^2_{\nu} \le \delta^*(t)\big)\, dt,
```

where $`\delta^*(t)`$ is the smallest moment-matched heterogeneity
threshold across the grid that is consistent with achieving $`P_{\min}`$
at burden coordinate $`t`$, and $`\nu`$ is the matched heterogeneity
degrees of freedom. The result is clamped to $`[0,1]`$.

## What This Method Advances

mixAda searches a grid of burden/heterogeneity mixtures and calibrates
the selected minimum p-value with a null integral.

## Assumptions

The adaptive search makes mixAda the most regime-agnostic member of the
mixture family: by minimizing over the mixing grid it tracks whichever
blend of concordant (burden) and discordant/antagonistic (heterogeneity)
signal best fits the variant, and the analytic null integration restores
valid calibration without permutation. This robustness costs a small
amount of power relative to a single correctly-specified mixture when
the architecture is known. Validity requires $`\Sigma`$ to be a correct,
full-rank trait correlation matrix, the per-trait Z-scores to be
approximately Gaussian under the null, and the modified-Liu and Simpson
approximations to hold in the tail.

## When This Method Is a Good Fit

Good when the best burden/heterogeneity balance is not known in advance.

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
pleio <- run_mixada(pleio)
head(pleio@results$mixada)
```

## Implementation Files

- R wrapper: `R/run_mixada.R`
- C++ kernel or dispatcher: `src/mixada.cpp`
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

Lee, S., Wu, M. C., Lin, X. (2012) Optimal tests for rare variant
effects in sequencing association studies. *Biostatistics* 13, 762-775.
DOI: <https://doi.org/10.1093/biostatistics/kxs014>

Liu, H., Tang, Y., Zhang, H. H. (2009) A new chi-square approximation to
the distribution of non-negative definite quadratic forms in non-central
normal variables. *Computational Statistics & Data Analysis* 53,
853-856. DOI: <https://doi.org/10.1016/j.csda.2008.11.025>

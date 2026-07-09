# PolarMorphism pleiotropy test across all phenotypes

Function:
[`run_polarmorphism()`](https://broccolito.github.io/pleioverse/reference/run_polarmorphism.md)

Result column: `polarmorphism`

Family: Joint angular pleiotropy test

## Introduction

Joint (all $`K`$ traits at once) pleiotropy test of von Berg et
al. (2022). It ZCA-cor whitens the per-trait Z-scores, maps each variant
to polar coordinates – a radius $`r`$ measuring overall association and
an angle measuring how shared (vs. concentrated on one trait) the
absolute whitened effect is – screens variants on the radius, then
returns an angular pleiotropy p-value for those passing the screen.
Unlike PLACO, PGCtest, MAIUP and conjFDR (which are pairwise on traits
$`i,j`$), PolarMorphism is evaluated jointly over all $`K`$ traits.

## Mathematical Framework and Implementation

**ZCA-cor whitening.** Let $`Z`$ be the $`n \times K`$ matrix of
per-trait Z-scores. A robust covariance is obtained by trimming variants
whose squared Mahalanobis distance $`z_i \Sigma^{-1} z_i^{\top}`$
exceeds `mahalanobis_threshold`$`^2`$, then recomputing the covariance
$`S`$ on the retained variants. With $`V = \mathrm{diag}(S)`$ and
correlation $`P = V^{-1/2} S V^{-1/2}`$, the ZCA-cor (correlation-based
zero-phase) whitening matrix is

``` math
W = V^{-1/2}\, P^{-1/2}, \qquad Z^{w} = Z\, W,
```

which decorrelates the traits while preserving each trait’s marginal
interpretation (attenuating genetic-correlation / “vertical”
pleiotropy).

**Radius (overall association).** Working with the absolute whitened
scores $`|Z^{w}|`$, the radius is the Euclidean norm

``` math
r = \Bigl(\textstyle\sum_{k=1}^{K} (Z^{w}_k)^2\Bigr)^{1/2} .
```

Under the global null (no association with any trait) the whitened
scores are independent standard normal, so

``` math
r^2 \sim \chi^2_{K},
```

giving a radial p-value $`p_r = \Pr(\chi^2_K \ge r^2)`$. These are
converted to Storey q-values, and the angular test is evaluated only for
variants whose radial q-value is below `r_qval_threshold` (the radial
screen identifies variants associated with *at least one* trait).

**Angle (shared vs. trait-specific absolute effect).** The angle
measures the deviation of the absolute whitened effect vector from a
single-trait axis. For $`K=2`$ it is computed from
$`\mathrm{atan2}(|Z^{w}_2|, |Z^{w}_1|)`$ and rescaled so that an
axis-aligned (single-trait) effect maps to the extreme and a perfectly
shared effect maps to the centre. For $`K>2`$ it is derived from the
angle between the whitened effect vector and its dominant axis,
$`\arccos(\max_k |Z^{w}_k| / r)`$, rescaled by the geometry constant so
the single-trait reference maps consistently. In the intended
interpretation, an axis-concentrated effect is less pleiotropic and a
larger angle means the absolute effect is more shared across traits. For
$`K=2`$, the current C++ tail routine returns a zero angular p-value for
the exact axis-angle edge case, so results exactly on an axis should be
treated cautiously.

**Angular tail probability.** Conditional on the radius, the angle under
the “effect concentrated on one trait” reference follows a von
Mises-type angular law with concentration $`\kappa`$ (a monotone
function of $`r`$, looked up from `.polar_kappa_table`). For $`K=2`$ the
angular tail uses the von Mises CDF series

``` math
F(\theta) = 2\left[\frac{\theta}{2\pi}
  + \frac{1}{\pi I_0(\kappa)}\sum_{j\ge 1}
    \frac{I_j(\kappa)\,\sin(j\theta)}{j} + \tfrac12\right],
```

evaluated at the (signed, mirror-folded) relative angle, where $`I_j`$
is the modified Bessel function of the first kind. For $`K>2`$ the
analogous tail is obtained from a normalizing constant
$`C = \kappa^{K/2-1}/\bigl[(2\pi)^{K/2} I_{K/2-1}(\kappa)\bigr]`$ and a
Gauss hypergeometric $`{}_2F_1`$ series, returned as a difference of the
antiderivative between the observed angle and the reference $`-\pi`$.
Small angular p-values flag horizontally pleiotropic variants (shared
across traits beyond what the radius alone implies).

## What This Method Advances

PolarMorphism separates radial association strength from angular
sharing, so a variant can be screened for any association before asking
whether the effect is spread across traits.

## Assumptions

- Under the global null the whitened Z-scores are jointly standard
  normal and independent, so $`r^2 \sim \chi^2_K`$ and the angle is
  distributed per the von Mises-type reference law.
- Cross-trait (genetic) correlation is well approximated by the linear
  correlation $`P`$, and ZCA-cor whitening removes it without destroying
  marginal effect direction.
- The robust-covariance trimming (Mahalanobis screen) yields a
  covariance dominated by null variants.
- The angle is taken on absolute whitened scores, so “shared” is
  measured up to sign; same-direction and opposite-direction effects are
  collapsed before angular testing, so this is not a separate antagonism
  or sign test.
- Variants are treated as exchangeable; LD is not explicitly modelled.

## When This Method Is a Good Fit

Good for multi-trait pleiotropy where the question is whether associated
effects are angularly shared rather than axis-aligned.

Only variants passing the radial q-value screen receive an angular
p-value; variants that do not pass the screen are returned as `NA`.

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
pleio <- run_polarmorphism(pleio)
head(pleio@results$polarmorphism)
```

## Implementation Files

- R wrapper: `R/run_polarmorphism.R`
- C++ kernel or dispatcher: `src/polarmorphism.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

von Berg J, ten Dam M, van der Laan SW, de Ridder J (2022).
PolarMorphism enables discovery of shared genetic variants across
multiple traits from GWAS summary statistics. *Bioinformatics* 38(Suppl
1):i212-i219. DOI: <https://doi.org/10.1093/bioinformatics/btac228>

Kessy A, Lewin A, Strimmer K (2018). Optimal whitening and
decorrelation. *The American Statistician* 72(4):309-314. DOI:
<https://doi.org/10.1080/00031305.2016.1277159>

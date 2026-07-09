# Direct burden (DSUM) multi-trait test

Function:
[`run_dsum()`](https://broccolito.github.io/pleioverse/reference/run_dsum.md)

Result column: `dsum`

Family: Burden-style multiple association

## Introduction

MPAT/MTAR multi-phenotype burden test that simply sums the per-trait
Z-scores with equal (unit) weights and tests the total with a two-sided
normal tail whose variance accounts for trait correlation. It is the
unweighted counterpart of the GLS-weighted SUM test.

## Mathematical Framework and Implementation

For a variant with per-trait Z-score vector $`z = (z_1,\dots,z_K)'`$,
$`z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)`$, over $`K`$ phenotypes,
the direct-burden statistic is the absolute value of the plain sum

``` math
T = \left| \mathbf{1}' z \right| = \left| \sum_{k=1}^{K} z_k \right|,
```

where $`\mathbf{1}`$ is the all-ones vector. Unlike SUM
([`run_sum()`](https://broccolito.github.io/pleioverse/reference/run_sum.md)),
no inverse-correlation weighting is applied; every trait contributes
equally.

Under the global null $`z \sim \mathcal{N}(0, \Sigma)`$, with $`\Sigma`$
the $`K \times K`$ trait correlation matrix, the sum $`\mathbf{1}' z`$
is normal with mean 0 and variance equal to the grand total of
$`\Sigma`$,

``` math
\mathrm{Var}(\mathbf{1}' z) = \mathbf{1}' \Sigma\, \mathbf{1}
  = \sum_{j,k} \Sigma_{jk} \equiv \texttt{sum\_Sigma}.
```

The two-sided p-value is

``` math
p = 2\,\Pr\!\left(\mathcal{N}(0, \sigma) \ge T\right),\qquad
  \sigma = \sqrt{\mathbf{1}' \Sigma\, \mathbf{1}}.
```

The correlation between traits is therefore absorbed entirely into the
null standard deviation rather than into the weights.

## What This Method Advances

DSUM keeps the burden idea but avoids inverse-correlation weights, which
makes it simpler and less dependent on inverting a noisy
trait-correlation estimate.

## Assumptions

DSUM is a one-degree-of-freedom burden test powered for *concordant*,
equal-direction effects across traits, like SUM, but it weights each
trait equally instead of optimally. It is more robust when $`\Sigma`$ is
near-singular or poorly estimated (it never inverts $`\Sigma`$), at the
cost of efficiency when traits are strongly correlated, where the
GLS-weighted SUM is preferable. It shares the burden-test weakness of
losing power for sign-*heterogeneous* effects, whose contributions
cancel in the sum, and for sparse single-trait signal. Validity requires
correct specification of $`\Sigma`$ and approximately Gaussian null
Z-scores; $`\Sigma`$ need only be positive semidefinite (no inverse is
taken).

## When This Method Is a Good Fit

Good for concordant effects when equal trait weights are preferred or
the correlation inverse is unstable.

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
pleio <- run_dsum(pleio)
head(pleio@results$dsum)
```

## Implementation Files

- R wrapper: `R/run_dsum.R`
- C++ kernel or dispatcher: `src/dsum.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175. DOI: <https://doi.org/10.1111/biom.12735>

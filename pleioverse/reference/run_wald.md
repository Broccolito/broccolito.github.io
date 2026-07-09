# Wald multi-trait test (inverse-variance quadratic form)

MPAT/MTAR multi-phenotype association test that combines the per-trait
Z-scores of a single variant into one inverse-correlation-weighted
quadratic form. Under the global null of no association with any trait
this statistic follows an exact chi-square distribution, giving a fast
closed-form p-value.

## Usage

``` r
run_wald(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$wald` populated (one p-value per
variant).

## Details

For a variant with per-trait Z-score vector \\z = (z_1,\dots,z_K)'\\,
where \\z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)\\ is the marginal
score (Wald) statistic for trait \\k\\ and \\K\\ is the number of
phenotypes, the test statistic is the Mahalanobis quadratic form \$\$T =
z' \Sigma^{-1} z,\$\$ where \\\Sigma\\ is the \\K \times K\\ trait
correlation matrix (the correlation of the Z-scores across variants
under the null), and \\\Sigma^{-1}\\ (denoted `Sinv` in the engine) is
precomputed once.

Under the global null \\H_0: z \sim \mathcal{N}(0, \Sigma)\\ the
whitened vector \\\Sigma^{-1/2} z\\ is standard multivariate normal, so
\$\$T = (\Sigma^{-1/2} z)'(\Sigma^{-1/2} z) \sim \chi^2_K.\$\$ The
reported p-value is the exact upper tail of a chi-square with \\K\\
degrees of freedom, \$\$p = \Pr(\chi^2_K \ge T).\$\$ No moment matching
or numerical inversion is needed; the p-value is closed form.

## Assumptions

This is the classical multivariate Wald / Hotelling-style omnibus test.
It spends one degree of freedom per trait and is therefore most powerful
when the signal is *dense* (many traits associated) and effect
directions are *heterogeneous*, with no prior knowledge of the sign or
magnitude pattern across traits. Because it weights by \\\Sigma^{-1}\\,
it accounts for trait correlation but, by treating all \\K\\ directions
equally, loses power relative to burden-type tests
([`run_sum`](https://broccolito.github.io/pleioverse/reference/run_sum.md),
[`run_dsum`](https://broccolito.github.io/pleioverse/reference/run_dsum.md))
when effects are sparse or concordant in direction. Validity requires
\\\Sigma\\ to be correctly specified and non-singular and the Z-scores
to be approximately Gaussian under the null (large-sample GWAS summary
statistics).

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175.
[doi:10.1111/biom.12735](https://doi.org/10.1111/biom.12735)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

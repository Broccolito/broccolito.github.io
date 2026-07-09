# Direct burden (DSUM) multi-trait test

MPAT/MTAR multi-phenotype burden test that simply sums the per-trait
Z-scores with equal (unit) weights and tests the total with a two-sided
normal tail whose variance accounts for trait correlation. It is the
unweighted counterpart of the GLS-weighted SUM test.

## Usage

``` r
run_dsum(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$dsum` populated (one p-value per
variant).

## Details

For a variant with per-trait Z-score vector \\z = (z_1,\dots,z_K)'\\,
\\z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)\\, over \\K\\ phenotypes,
the direct-burden statistic is the absolute value of the plain sum \$\$T
= \left\| \mathbf{1}' z \right\| = \left\| \sum\_{k=1}^{K} z_k
\right\|,\$\$ where \\\mathbf{1}\\ is the all-ones vector. Unlike SUM
([`run_sum`](https://broccolito.github.io/pleioverse/reference/run_sum.md)),
no inverse-correlation weighting is applied; every trait contributes
equally.

Under the global null \\z \sim \mathcal{N}(0, \Sigma)\\, with \\\Sigma\\
the \\K \times K\\ trait correlation matrix, the sum \\\mathbf{1}' z\\
is normal with mean 0 and variance equal to the grand total of
\\\Sigma\\, \$\$\mathrm{Var}(\mathbf{1}' z) = \mathbf{1}' \Sigma\\
\mathbf{1} = \sum\_{j,k} \Sigma\_{jk} \equiv \mathtt{sum\\Sigma}.\$\$
The two-sided p-value is \$\$p = 2\\\Pr\\\left(\mathcal{N}(0, \sigma)
\ge T\right),\qquad \sigma = \sqrt{\mathbf{1}' \Sigma\\ \mathbf{1}}.\$\$
The correlation between traits is therefore absorbed entirely into the
null standard deviation rather than into the weights.

## Assumptions

DSUM is a one-degree-of-freedom burden test powered for *concordant*,
equal-direction effects across traits, like SUM, but it weights each
trait equally instead of optimally. It is more robust when \\\Sigma\\ is
near-singular or poorly estimated (it never inverts \\\Sigma\\), at the
cost of efficiency when traits are strongly correlated, where the
GLS-weighted SUM is preferable. It shares the burden-test weakness of
losing power for sign-*heterogeneous* effects, whose contributions
cancel in the sum, and for sparse single-trait signal. Validity requires
correct specification of \\\Sigma\\ and approximately Gaussian null
Z-scores; \\\Sigma\\ need only be positive semidefinite (no inverse is
taken).

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175.
[doi:10.1111/biom.12735](https://doi.org/10.1111/biom.12735)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

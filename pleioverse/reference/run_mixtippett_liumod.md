# Mixture-Tippett (minimum-p) omnibus (modified-Liu tail)

MPAT/MTAR multi-trait association test that decomposes the per-variant
signal into a burden/mean component and a heterogeneity/variance
component and combines the two component p-values by Tippett's minimum-p
rule. This variant evaluates the heterogeneity tail by the modified-Liu
chi-square moment-matching approximation.

## Usage

``` r
run_mixtippett_liumod(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$mixtippett_liumod` populated (one
p-value per variant).

## Details

For a variant the input is the vector of trait Z-scores \\z =
(z_1,\dots,z_K)'\\ with \\z_k = \hat\beta_k/\widehat{\mathrm{se}}\_k\\,
which under the null is \\z \sim \mathcal{N}(0,\Sigma)\\ for the
\\K\times K\\ trait correlation \\\Sigma\\. Write \\\Sigma^{-1}\\ for
its inverse, \\\mathbf{1}\\ for the all-ones vector, \\a =
\Sigma^{-1}\mathbf{1}\\, \\s = \mathbf{1}'\Sigma^{-1}\mathbf{1}\\, and
the centering projection \\H = \mathbf{1}\mathbf{1}'/K\\.

**Burden (mean) component.** The generalized-least-squares burden score
tests a shared, concordant mean across traits, \$\$U\_\mu =
\frac{(a'z)^2}{s}, \qquad p\_\mu = \Pr(\chi^2_1 \ge U\_\mu),\$\$ i.e.
\\U\_\mu\\ is exactly \\\chi^2_1\\ under the null.

**Heterogeneity (variance) component.** The complementary quadratic form
isolates departures from a common mean using the projected precision
matrix \$\$M\_\tau = (I-H)\\\Sigma^{-1}\Sigma^{-1}\\(I-H), \qquad
U\_\tau = z' M\_\tau z .\$\$ Under the null \\U\_\tau\\ is a weighted
sum of independent \\\chi^2_1\\ variables, \\U\_\tau \sim \sum_j
\lambda_j \chi^2\_{1,j}\\, whose weights \\\lambda_j\\ are the
eigenvalues of \\\Sigma^{1/2} M\_\tau \Sigma^{1/2}\\. Its p-value
\\p\_\tau\\ is the upper tail of this mixture.

**Tippett combination.** Tippett's rule takes the smaller of the two
component p-values and calibrates it against the null of two independent
tests, \$\$m = \min(p\_\mu, p\_\tau), \qquad p = 1 - (1 - m)^2 .\$\$
Equivalently \\m\\ is referred to a Beta(1,2) reference; the reported
p-value is \\p\\. Tippett emphasizes the single strongest component and
is typically more powerful than Fisher when only one of the two regimes
is active.

**Tail method (modified Liu).** The heterogeneity tail \\p\_\tau\\ is
computed by the modified-Liu approximation, which matches all four
cumulants (including kurtosis) of \\U\_\tau\\ to a non-central
chi-square. It is as fast as the plain Liu method and typically more
accurate, though still an approximation; the Davies variant remains the
reference for extreme tails.

## Assumptions

The omnibus is designed to be robust across genetic-architecture
regimes: the burden component is most powerful when effects are
concordant in sign and magnitude across traits, while the heterogeneity
component captures discordant, antagonistic, or sparse effects to which
the burden is blind. Tippett combination preserves sensitivity to either
regime without prior knowledge of which holds. Validity requires
\\\Sigma\\ to be a correct, full-rank trait correlation matrix and the
per-trait Z-scores to be approximately Gaussian under the null.

## References

Liu, Z., Lin, X. (2018) Multiple Phenotype Association Tests Using
Summary Statistics in Genome-Wide Association Studies. *Biometrics*
74(1):165-175.
[doi:10.1111/biom.12735](https://doi.org/10.1111/biom.12735)

Luo, L., Shen, J., Zhang, H., et al. (2020) Multi-trait analysis of
rare-variant association summary statistics using MTAR. *Nature
Communications* 11:2850.
[doi:10.1038/s41467-020-16591-0](https://doi.org/10.1038/s41467-020-16591-0)

Liu, H., Tang, Y., Zhang, H. H. (2009) A new chi-square approximation to
the distribution of non-negative definite quadratic forms in non-central
normal variables. *Computational Statistics & Data Analysis* 53,
853-856.
[doi:10.1016/j.csda.2008.11.025](https://doi.org/10.1016/j.csda.2008.11.025)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

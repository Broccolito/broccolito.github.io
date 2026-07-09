# Mixture standard-deviation omnibus (Liu tail)

MPAT/MTAR multi-trait association test that blends the burden/mean and
heterogeneity/variance components into a single quadratic form, with the
mixing weight chosen to balance the two components on the
standard-deviation scale. This variant evaluates the weighted-chi-square
tail by the Liu chi-square moment-matching approximation.

## Usage

``` r
run_mixsd_liu(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$mixsd_liu` populated (one p-value per
variant).

## Details

For a variant the input is the vector of trait Z-scores \\z =
(z_1,\dots,z_K)'\\ with \\z_k = \hat\beta_k/\widehat{\mathrm{se}}\_k\\,
which under the null is \\z \sim \mathcal{N}(0,\Sigma)\\ for the
\\K\times K\\ trait correlation \\\Sigma\\. Write \\\Sigma^{-1}\\ for
its inverse, \\\mathbf{1}\\ for the all-ones vector, \\a =
\Sigma^{-1}\mathbf{1}\\, \\s = \mathbf{1}'\Sigma^{-1}\mathbf{1}\\, and
\\H = \mathbf{1}\mathbf{1}'/K\\.

**Component kernels.** The burden/mean direction is carried by the
rank-one kernel \\a a'\\ (so that \\z' a a' z = (a'z)^2\\), and the
heterogeneity/variance direction by the centered precision kernel
\$\$M\_\tau = (I-H)\\\Sigma^{-1}\Sigma^{-1}\\(I-H).\$\$

**SD-balanced mixing.** Rather than combining two p-values, this test
forms a single kernel \\L\_\phi\\ that is a convex blend of the two,
\$\$L\_\phi = \phi\\ a a' + (1-\phi)\\ M\_\tau, \qquad T = z' L\_\phi z
.\$\$ Unlike the variance-scaled variant, the weight \\\phi\\ here
equalizes the two components' contributed null *standard deviations*:
\$\$\phi = \frac{\sqrt{V\_\tau}}{\sqrt{V\_\mu} + \sqrt{V\_\tau}}, \quad
V\_\mu = 2\\s^2, \quad V\_\tau = 2\\\mathrm{tr}(M\_\tau\Sigma
M\_\tau\Sigma),\$\$ where \\V\_\mu\\ and \\V\_\tau\\ are the null
variances of \\(a'z)^2\\ and \\z'M\_\tau z\\. This places relatively
more weight on the heterogeneity component than the variance-scaled
rule. Under the null \\T\\ is a weighted sum of independent \\\chi^2_1\\
variables, \\T \sim \sum_j \lambda_j \chi^2\_{1,j}\\, with weights
\\\lambda_j\\ the eigenvalues of \\\Sigma^{1/2} L\_\phi \Sigma^{1/2}\\;
the reported p-value is the upper tail of this mixture.

**Tail method (Liu).** The tail is computed by the Liu (Liu-Tang-Zhang)
approximation, which matches the first three cumulants of \\T\\ to a
single (possibly non-central) chi-square. It is fast and closed-form but
only approximate, and can lose accuracy far in the tail; use the Davies
variant when extreme-tail precision is required.

## Assumptions

SD-balanced mixing makes the test robust across genetic-architecture
regimes without choosing between them: the \\a a'\\ part captures
concordant shared effects while the \\M\_\tau\\ part captures
discordant, antagonistic, or sparse effects. Compared with the
variance-scaled mixture this rule favors the heterogeneity direction
somewhat more, trading a little power under perfect concordance for
robustness under discordance. Validity requires \\\Sigma\\ to be a
correct, full-rank trait correlation matrix and the per-trait Z-scores
to be approximately Gaussian under the null.

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

# Mixture-Fisher omnibus (Davies tail)

MPAT/MTAR multi-trait association test that decomposes the per-variant
signal into a burden/mean component and a heterogeneity/variance
component and combines the two component p-values by Fisher's method.
This variant evaluates the heterogeneity tail by the Davies-labeled
(Imhof-style) numerical inversion.

## Usage

``` r
run_mixfisher_davies(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$mixfisher_davies` populated (one
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

**Fisher combination.** The two components are asymptotically
independent (\\U\_\mu\\ lives in the \\\mathbf{1}\\ direction,
\\U\_\tau\\ in its orthogonal complement), so they are pooled by
Fisher's rule \$\$T = -2\big(\log p\_\mu + \log p\_\tau\big) \\\sim\\
\chi^2_4,\$\$ and the reported p-value is \\\Pr(\chi^2_4 \ge T)\\.

**Tail method (Davies).** The heterogeneity tail \\p\_\tau\\ is computed
by the Davies-labeled numerical inversion of the characteristic function
(Imhof numerical integration). It is the direct numerical-inversion
option among the three tail variants, at higher computational cost; if
the returned tail probability is outside \\(0,1)\\, the implementation
falls back to the modified-Liu approximation.

## Assumptions

The omnibus is designed to be robust across genetic-architecture
regimes: the burden component is most powerful when effects are
concordant in sign and magnitude across traits, while the heterogeneity
component captures discordant, antagonistic, or sparse effects to which
the burden is blind. Fisher combination preserves sensitivity to either
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

Davies, R. B. (1980) The distribution of a linear combination of
\\\chi^2\\ random variables. *Applied Statistics* 29, 323-333.
[doi:10.2307/2346911](https://doi.org/10.2307/2346911)

Liu, H., Tang, Y., Zhang, H. H. (2009) A new chi-square approximation to
the distribution of non-negative definite quadratic forms in non-central
normal variables. *Computational Statistics & Data Analysis* 53,
853-856.
[doi:10.1016/j.csda.2008.11.025](https://doi.org/10.1016/j.csda.2008.11.025)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

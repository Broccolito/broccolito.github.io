# Mixture variance-component omnibus (Davies tail)

MPAT/MTAR multi-trait association test that blends the burden/mean and
heterogeneity/variance components into a single quadratic form, with the
mixing weight chosen to balance the two components on the variance
scale. This variant evaluates the weighted-chi-square tail by the
Davies-labeled (Imhof-style) numerical inversion.

## Usage

``` r
run_mixvar_davies(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$mixvar_davies` populated (one p-value
per variant).

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

**Variance-scale mixing.** Rather than combining two p-values, this test
forms a single kernel \\L\_\phi\\ that is a convex blend of the two,
\$\$L\_\phi = \phi\\ a a' + (1-\phi)\\ M\_\tau, \qquad T = z' L\_\phi z
.\$\$ The weight \\\phi\\ is fixed from the two components' null
variance scales: \$\$\phi = \frac{V\_\tau}{V\_\mu + V\_\tau}, \quad
V\_\mu = 2\\s^2, \quad V\_\tau = 2\\\mathrm{tr}(M\_\tau\Sigma
M\_\tau\Sigma),\$\$ where \\V\_\mu\\ and \\V\_\tau\\ are the null
variances of \\(a'z)^2\\ and \\z'M\_\tau z\\. Under the null \\T\\ is a
weighted sum of independent \\\chi^2_1\\ variables, \\T \sim \sum_j
\lambda_j \chi^2\_{1,j}\\, with weights \\\lambda_j\\ the eigenvalues of
\\\Sigma^{1/2} L\_\phi \Sigma^{1/2}\\; the reported p-value is the upper
tail of this mixture.

**Tail method (Davies).** The tail is computed by the Davies-labeled
numerical inversion of the characteristic function (Imhof numerical
integration). It is the direct numerical-inversion option among the
three tail variants, at higher computational cost; if the returned tail
probability is outside \\(0,1)\\, the implementation falls back to the
modified-Liu approximation.

## Assumptions

Variance-scale mixing aims to make the test robust across
genetic-architecture regimes without choosing between them: the \\a a'\\
part captures concordant shared effects while the \\M\_\tau\\ part
captures discordant, antagonistic, or sparse effects, and the
variance-scale \\\phi\\ aims to keep neither part dominating the null
scale. Validity requires \\\Sigma\\ to be a correct, full-rank trait
correlation matrix and the per-trait Z-scores to be approximately
Gaussian under the null.

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

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

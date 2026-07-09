# PLACO+ pleiotropy test on a phenotype pair

Pairwise (traits \\i,j\\) test of pleiotropy under the *composite null*
that a variant is associated with at most one of the two traits, using
the product of per-trait Z-scores as the test statistic (Ray and
Chatterjee, 2020). PLACO is applied to the trait pair only and must be
run separately for each pair of interest.

## Usage

``` r
run_placo(pleio, i = 1, j = 2, p_threshold = 1e-04)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object.

- i, j:

  phenotype indices for the pair. Default `1,2`.

- p_threshold:

  per-trait significance threshold used for variance/ correlation
  estimation. Default `1e-4`.

## Value

the `pleio` object with `results$placo` populated.

## Details

For variant \\v\\ let \\Z_1, Z_2\\ be the per-trait Z-scores
(\\\hat\beta/\widehat{\mathrm{se}}\\). PLACO tests the composite null
\$\$H_0:\\ \beta_1 = 0 \\\text{or}\\ \beta_2 = 0,\$\$ i.e. the union of
the three sub-null configurations \\(0,0)\\, \\(0,\ne 0)\\, \\(\ne
0,0)\\, against the pleiotropic alternative \\\beta_1 \ne 0\\ *and*
\\\beta_2 \ne 0\\. The statistic is the product \$\$T = Z_1 Z_2 .\$\$

**Null density of the product.** Under the null each \\Z_k\\ is
approximately normal with mean 0, and \\(Z_1,Z_2)\\ may be correlated
(e.g. through shared controls). Writing \\r\\ for that correlation, the
density of the standardized product \\x = Z_1 Z_2\\ is the correlated
normal-product density \$\$f(x \mid r) = \frac{1}{\pi \sqrt{1-r^2}} \\
\exp\\\left(\frac{r\\x}{1-r^2}\right) \\
K_0\\\left(\frac{\|x\|}{1-r^2}\right),\$\$ where \\K_0\\ is the modified
Bessel function of the second kind of order 0. The implementation
evaluates \\f(\cdot \mid r)\\ on a power-spaced grid that is dense near
0 (node \\x_i = X\_{\max}\\(i/(N-1))^{8}\\) and forms the upper- and
lower-tail integrals by reverse-cumulative trapezoidal quadrature,
giving a two-sided tail function \$\$G(t) = \int\_{\|x\| \ge t} f(x \mid
r)\\ dx .\$\$

**Variance and correlation estimation.** The null moments are estimated
off the variants that are *not* jointly significant. The per-trait null
variances \\s_1^2, s_2^2\\ are the sample variances of \\Z_1, Z_2\\ over
variants *not* significant in both traits (i.e. excluding \\p_1 \<
\tau\\ and \\p_2 \< \tau\\), where \\\tau\\ is `p_threshold`. The null
correlation \\r\\ is the Pearson correlation of \\(Z_1,Z_2)\\ over
variants significant in *neither* trait (\\p_1 \ge \tau\\ and \\p_2 \ge
\tau\\).

**Composite p-value.** Because the null is composite, the tail of the
product alone is anti-conservative when only one trait is non-null.
Following the PLACO construction the p-value uses an inclusion-exclusion
combination of the two-sided product tail evaluated at the observed
\\t_0 = \|Z_1 Z_2\|\\ rescaled by each per-trait null SD: \$\$p =
G\\\left(\frac{t_0}{s_1}\right) + G\\\left(\frac{t_0}{s_2}\right) -
G\\\left(t_0\right),\$\$ clamped to \\\[10^{-300}, 1\]\\. Small \\p\\
indicates evidence that the variant is associated with *both* traits.

## Assumptions

- Per-trait Z-scores are (approximately) standard normal under the
  respective trait null; departures are absorbed into the empirically
  estimated null variances \\s_1^2, s_2^2\\.

- Cross-trait dependence is captured by a single correlation \\r\\
  estimated from variants null in both traits; estimation is meaningful
  only if such variants dominate (sparse true pleiotropy).

- The product null density is the correlated normal-product law, so the
  bivariate Z is jointly normal with zero mean under each sub-null.

- Variants are treated as exchangeable for moment estimation (LD and
  sample-size heterogeneity are not modelled here).

## References

Ray D, Chatterjee N (2020). A powerful method for pleiotropic analysis
under composite null hypothesis identifies novel shared loci between
Type 2 Diabetes and Prostate Cancer. *PLoS Genetics* 16(12):e1009218.
[doi:10.1371/journal.pgen.1009218](https://doi.org/10.1371/journal.pgen.1009218)

# PGCtest (CompTestER) normal-product test on a phenotype pair

Pairwise (traits \\i,j\\) pleiotropy screen based on the
*normal-product* statistic \\\|Z_1 Z_2\|\\, with a Bessel-\\K_0\\ tail
CDF and an empirical variance normalization. In this package, the
implemented routine applies a product-Z tail calibration with one
empirical scale; it is evaluated for one trait pair at a time.

## Usage

``` r
run_pgctest(pleio, i = 1, j = 2)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object.

- i, j:

  phenotype indices for the pair. Default `1,2`.

## Value

the `pleio` object with `results$pgctest` populated.

## Details

For variant \\v\\ let \\Z_1, Z_2\\ be the per-trait Z-scores
(\\\hat\beta/\widehat{\mathrm{se}}\\) and let \\X = \|Z_1 Z_2\|\\. The
implemented test rejects for large product magnitude \\X\\ after
calibration against the normal-product reference distribution; it does
not explicitly evaluate a multi-component composite-null mixture.

**Standard normal-product tail.** When \\Z_1, Z_2\\ are independent
standard normal, the density of \\X = \|Z_1 Z_2\|\\ is \$\$f(x) =
\frac{2}{\pi}\\ K_0(x), \qquad x \ge 0,\$\$ with \\K_0\\ the modified
Bessel function of the second kind of order 0. The implementation builds
the corresponding survival (upper-tail) function once, \$\$S(x) = \Pr(X
\ge x) = \int_x^{\infty} \frac{2}{\pi} K_0(u)\\ du,\$\$ by evaluating
\\f\\ on a degree-10 power-spaced grid (\\x_i =
x\_{\max}\\(i/(\text{len}-1))^{10}\\, dense near 0) and accumulating a
reverse-cumulative trapezoidal integral, normalized so that \\S(0)=1\\;
the grid upper bound is \\x\_{\max}=2\max_v X_v\\.

**Empirical variance normalization.** Because the per-trait Z-scores
need not be exactly standard normal (inflation, residual cross-trait
dependence), the observed statistic is rescaled before being passed to
\\S(\cdot)\\. With sample variances \\\widehat{\mathrm{Var}}(Z_1)\\ and
\\\widehat{\mathrm{Var}}(Z_2)\\ over all variants, define \$\$d =
\sqrt{\widehat{\mathrm{Var}}(Z_1) + \widehat{\mathrm{Var}}(Z_2) -
1},\$\$ and report the p-value by linear interpolation of the tail CDF
at the rescaled point \$\$p = S\\\left(\frac{\|Z_1 Z_2\|}{d}\right).\$\$
A signed pleiotropy z is also computed internally as \\z =
\mathrm{sign}(Z_1 Z_2)\\\Phi^{-1}(1 - p/2)\\; the R wrapper returns the
p-value component (`$pp`).

## Assumptions

- Under the relevant null the per-trait Z-scores are approximately
  standard normal, so that \\\|Z_1 Z_2\|\\ follows the
  \\\tfrac{2}{\pi}K_0\\ normal-product law.

- Mild departures (variance inflation, weak dependence) are corrected by
  the single empirical scale \\d\\; the correction assumes the inflation
  is well summarized by the marginal Z variances.

- Independence of \\Z_1, Z_2\\ is assumed by the base tail; strong
  cross-trait correlation is not separately modelled (unlike PLACO,
  which carries an explicit \\r\\).

- Variants are treated as exchangeable for the variance estimates.

## References

Andreassen OA, Thompson WK, Schork AJ, et al. (2013). Improved detection
of common variants associated with schizophrenia and bipolar disorder
using pleiotropy-informed conditional false discovery rate. *PLoS
Genetics* 9(4):e1003455.
[doi:10.1371/journal.pgen.1003455](https://doi.org/10.1371/journal.pgen.1003455)

Ray D, Chatterjee N (2020). A powerful method for pleiotropic analysis
under composite null hypothesis. *PLoS Genetics* 16(12):e1009218.
[doi:10.1371/journal.pgen.1009218](https://doi.org/10.1371/journal.pgen.1009218)

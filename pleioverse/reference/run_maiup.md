# MAIUP empirical-FDR pleiotropy screen on a phenotype pair

Pairwise (traits \\i,j\\) Mixture-model-Adjusted Intersection-Union
empirical FDR procedure (Zeng et al., MAIUP). It decorrelates the
per-trait Z-scores, forms the intersection-union statistic, and
calibrates an empirical FDR against a multi-component null mixture (Zeng
et al.). MAIUP is run for one trait pair at a time.

## Usage

``` r
run_maiup(pleio, i = 1, j = 2, p_threshold = 1e-04)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object.

- i, j:

  phenotype indices for the pair. Default `1,2`.

- p_threshold:

  per-trait threshold for null correlation estimation. Default `1e-4`.

## Value

the `pleio` object with `results$maiup` populated.

## Details

**Decorrelating whitening.** Let \\Z = (Z_1, Z_2)\\ be the per-trait
Z-scores. The null cross-trait correlation matrix \\R\\ is estimated
from variants significant in *neither* trait (\\p_1 \ge \tau\\ and \\p_2
\ge \tau\\, with \\\tau\\ = `p_threshold`). With eigendecomposition \\R
= V \Lambda V^{\top}\\, the symmetric (ZCA-type) whitening matrix \\W =
V \Lambda^{-1/2} V^{\top}\\ produces decorrelated scores \\Z^{d} =
Z\\W\\ with components \\Z^{d}\_1, Z^{d}\_2\\ and two-sided per-trait
p-values \$\$P^{d}\_k = 2\\\bigl(1 - \Phi(\|Z^{d}\_k\|)\bigr), \qquad k
= 1,2 .\$\$

**Intersection-union statistic.** The IUT statistic for the composite
null "associated with at most one trait" is the maximum of the two
decorrelated p-values, \$\$P\_{\max} = \max\bigl(P^{d}\_1,
P^{d}\_2\bigr),\$\$ which rejects only when *both* traits are
individually significant.

**Null mixture.** Variants fall into four classes by which traits are
null: both null (\\a\_{00}\\), only trait 2 non-null (\\a\_{01}\\), only
trait 1 non-null (\\a\_{10}\\), and both non-null (the alternative). The
null proportions are estimated by Storey-type tail counting at cut \\c =
0.5\\: \$\$\hat a\_{1} = \frac{\\\\P^{d}\_1 \ge c\\}{n\\(1-c)}, \quad
\hat a\_{2} = \frac{\\\\P^{d}\_2 \ge c\\}{n\\(1-c)}, \quad \hat a\_{00}
= \frac{\\\\P^{d}\_1 \ge c,\\ P^{d}\_2 \ge c\\}{n\\(1-c)^2},\$\$
followed by \\\hat a\_{01} = \max(0,\hat a_1 - \hat a\_{00})\\ and
\\\hat a\_{10} = \max(0,\hat a_2 - \hat a\_{00})\\ (with edge
corrections, and a one-sided Kolmogorov-Smirnov uniformity check that
forces \\\hat a_k = 1\\ when a trait's p-values are indistinguishable
from Uniform\\(0,1)\\).

**Marginal null CDFs.** For each trait the null CDF of \\P^{d}\_k\\ is
estimated nonparametrically by the Least Concave Majorant (LCM) of the
empirical CDF, then mixture-adjusted to isolate the null component:
\\\hat F_k(p) = \bigl(F^{\mathrm{LCM}}\_k(p) - a_k\\p\bigr)/(1-a_k)\\.

**Empirical FDR.** Writing \\\hat F(p)\\ for the empirical CDF of
\\P\_{\max}\\, the false-discovery numerator sums the contributions of
the three null sub-classes evaluated at \\p = P\_{\max}\\:
\$\$\widehat{\mathrm{FDR}}(p) = \frac{\hat a\_{01}\\ p\\ \hat F_2(p) +
\hat a\_{10}\\ p\\ \hat F_1(p) + \hat a\_{00}\\ p^{2}}{\hat F(p)} .\$\$
(Under each sub-null the max-statistic tail factorizes into a uniform
factor \\p\\ for the null trait and the estimated null CDF for the
other, with \\p^2\\ when both are null.) The returned `p_maiup` is this
FDR after monotone (step-down) enforcement in decreasing \\P\_{\max}\\
and clamping to \\\[0,1\]\\.

## Assumptions

- Cross-trait dependence among null Z-scores is fully captured by the
  linear correlation \\R\\, so symmetric whitening renders the two
  traits' null statistics independent.

- The genome-wide p-value mixture has the four-class structure above,
  with each null sub-class contributing a uniform marginal for its null
  trait(s).

- LCM/isotonic estimation requires the true null CDF to be concave
  (p-value density non-increasing), as for well-behaved test statistics.

- Null proportions are estimable from the upper half (\\p \ge 0.5\\) of
  the p-value distribution (signal sparse near 1).

- Variants are treated as exchangeable; LD is not explicitly modelled.

## References

Zeng P, et al. MAIUP: a mixture-adjusted intersection-union test for
pleiotropy. <https://github.com/biostatpzeng/MAIUP>

Kessy A, Lewin A, Strimmer K (2018). Optimal whitening and
decorrelation. *The American Statistician* 72(4):309-314.
[doi:10.1080/00031305.2016.1277159](https://doi.org/10.1080/00031305.2016.1277159)

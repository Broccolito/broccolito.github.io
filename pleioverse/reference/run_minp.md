# Minimum-p test with multivariate-normal correction (MinP)

Multi-trait association test that takes the smallest per-trait p-value
and corrects it for the correlated traits using a numerical
multivariate-normal calculation under the null of the Z-scores (a la
Conneely & Boehnke). It is aimed at association signal that shows up
strongly in at least one of the original traits.

## Usage

``` r
run_minp(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$minp` populated (one p-value per
variant).

## Details

Let \\z=(z_1,\dots,z_K)'\\ be the per-trait Z-scores at a variant. Each
trait contributes a two-sided single-trait test \$\$Q_k = z_k^2\\ \sim\\
\chi^2_1, \qquad p_k = \Pr(\chi^2_1 \> Q_k),\$\$ and the observed
minimum is \\p\_{\min}=\min\_{k=1,\dots,K} p_k\\. Because the original
traits are correlated, this method uses the joint multivariate normal
probability rather than an independence-based Sidak correction or a
conservative Bonferroni bound. Under the global null \\z\sim
N(0,\Sigma)\\ the corrected p-value is the probability that at least one
trait reaches the observed threshold. Writing the per-trait threshold as
\\c = \Phi^{-1}\\\big(1-p\_{\min}/2\big)\\ (so that
\\\Pr(\|z_k\|\>c)=p\_{\min}\\), the correlation-adjusted corrected
p-value is \$\$p = 1 - \Pr\\\Big(\bigcap\_{k=1}^{K}\\-c \le z_k \le
c\\\Big) = 1 - \Pr\_{N(0,\Sigma)}\big(\[-c,c\]^K\big).\$\$ The rectangle
(box) probability of the multivariate normal with correlation \\\Sigma\\
is evaluated by quasi-randomized numerical integration using the lower
Cholesky factor \\L\\ of \\\Sigma\\ (\\\Sigma=LL'\\, engine
`cholSigma`), and \\p\\ is reported.

## Assumptions

The null distribution of the Z-vector is exactly multivariate normal
\\N(0,\Sigma)\\ with known trait-correlation matrix \\\Sigma\\. The
correction integrates the joint MVN density over the acceptance box
\\\[-c,c\]^K\\, which accounts for the specified trait correlation up to
the numerical integration tolerance. MinP operates on the original
traits rather than whitened PCs, so it is most powerful when the
association is large in one or a few individual traits. The box
probability is computed by Monte-Carlo/quasi-Monte-Carlo numerical
integration, so the result carries a small, controllable integration
error.

## References

Conneely KN, Boehnke M (2007). So Many Correlated Tests, So Little Time!
Rapid Adjustment of P Values for Multiple Correlated Tests. *American
Journal of Human Genetics* 81(6):1158–1168.
[doi:10.1086/522036](https://doi.org/10.1086/522036)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

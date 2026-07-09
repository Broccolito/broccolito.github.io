# Weighted-identity (WI) multi-trait test

MPAT/MTAR multi-phenotype association test that uses the unweighted sum
of squared per-trait Z-scores. Because the Z-scores are correlated under
the null, this statistic is a weighted mixture of one-degree-of-freedom
chi-squares whose weights are the eigenvalues of the trait correlation
matrix, and its tail is evaluated numerically.

## Usage

``` r
run_wi(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$wi` populated (one p-value per
variant).

## Details

For a variant with per-trait Z-score vector \\z = (z_1,\dots,z_K)'\\,
\\z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)\\, over \\K\\ phenotypes,
the statistic is the plain (identity-weighted) sum of squares \$\$T = z'
z = \sum\_{k=1}^{K} z_k^2.\$\$ Under the global null \\z \sim
\mathcal{N}(0, \Sigma)\\, with \\\Sigma\\ the \\K \times K\\ trait
correlation matrix, write the spectral decomposition \\\Sigma =
\sum\_{k} \lambda_k\\ v_k v_k'\\ with eigenvalues \\\lambda_1 \ge \dots
\ge \lambda_K \> 0\\ (`eval` in the engine). Then \$\$T
\\\stackrel{d}{=}\\ \sum\_{k=1}^{K} \lambda_k\\ \chi^2\_{1,k},\$\$ a
weighted sum of independent \\\chi^2_1\\ variables with weights equal to
the eigenvalues of \\\Sigma\\. The p-value \$\$p = \Pr\\\Big(\sum\_{k}
\lambda_k \chi^2\_{1,k} \ge T\Big)\$\$ is computed by the Davies-labeled
numerical inversion of the characteristic function; if that returns a
value outside \\(0,1)\\ the routine falls back to the modified Liu
moment-matching (Liu-Tang-Zhang) approximation.

## Assumptions

WI weights every trait equally on the original scale and does not whiten
by \\\Sigma^{-1}\\, so it up-weights the directions of *high* trait
correlation (large eigenvalues). It is well powered for *dense*,
heterogeneous signal that is aligned with the leading eigenvectors of
\\\Sigma\\, and is more robust than the Wald test when \\\Sigma\\ is
near-singular (it never inverts \\\Sigma\\). It is comparatively weak
for sparse signal concentrated in a single trait or for effects lying in
the low-variance directions of \\\Sigma\\. Validity requires correct
specification of \\\Sigma\\ and approximately Gaussian null Z-scores.

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175.
[doi:10.1111/biom.12735](https://doi.org/10.1111/biom.12735)

Davies, R. B. (1980). Algorithm AS 155: The distribution of a linear
combination of \\\chi^2\\ random variables. *Journal of the Royal
Statistical Society, Series C (Applied Statistics)*, 29(3), 323-333.
[doi:10.2307/2346911](https://doi.org/10.2307/2346911)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

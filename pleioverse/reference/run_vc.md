# Variance-component (VC / SKAT-style) multi-trait test

MPAT/MTAR multi-phenotype association test of the SKAT
variance-component family. It scores the squared length of the
whitened-then-whitened-again Z-vector \\\Sigma^{-1} z\\; under the null
this is a weighted sum of chi-squares whose weights are the eigenvalues
of \\\Sigma^{-1}\\, and its tail is evaluated numerically.

## Usage

``` r
run_vc(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$vc` populated (one p-value per
variant).

## Details

For a variant with per-trait Z-score vector \\z = (z_1,\dots,z_K)'\\,
\\z_k = \hat\beta_k / \mathrm{se}(\hat\beta_k)\\, over \\K\\ phenotypes,
let \\s = \Sigma^{-1} z\\, where \\\Sigma\\ is the \\K \times K\\ trait
correlation matrix and \\\Sigma^{-1}\\ (`Sinv`) is precomputed once. The
variance-component statistic is \$\$T = s' s = (\Sigma^{-1}
z)'(\Sigma^{-1} z) = z'\\ \Sigma^{-2}\\ z.\$\$ This is the SKAT-type
score for a random-effects (variance-component) model in which the
per-trait effects are treated as exchangeable random effects with kernel
\\\Sigma^{-2}\\.

Under the global null \\z \sim \mathcal{N}(0, \Sigma)\\ one has \\s =
\Sigma^{-1} z \sim \mathcal{N}(0, \Sigma^{-1})\\, so with eigenvalues
\\\mu_1 \ge \dots \ge \mu_K \> 0\\ of \\\Sigma^{-1}\\ (equivalently
\\\mu_k = 1/\lambda_k\\ for eigenvalues \\\lambda_k\\ of \\\Sigma\\;
stored as `eval_Sinv`) \$\$T \\\stackrel{d}{=}\\ \sum\_{k=1}^{K} \mu_k\\
\chi^2\_{1,k}.\$\$ The p-value \$\$p = \Pr\\\Big(\sum\_{k} \mu_k
\chi^2\_{1,k} \ge T\Big)\$\$ is computed by the Davies-labeled numerical
inversion, with the modified Liu moment-matching approximation as a
fallback when Davies returns a value outside \\(0,1)\\.

## Assumptions

Because the weights are the eigenvalues of \\\Sigma^{-1}\\, VC
up-weights the *low-variance* directions of \\\Sigma\\ (the directions
in which the traits are nearly independent or anti-correlated), the
opposite emphasis to the WI test
([`run_wi`](https://broccolito.github.io/pleioverse/reference/run_wi.md)).
It is a dense, omnibus variance-component test, most powerful when
effects spread across the trait-decorrelated directions, and it can be
sensitive when \\\Sigma\\ is near-singular (small \\\lambda_k\\ inflate
\\\mu_k\\). Validity requires a correctly specified, non-singular
\\\Sigma\\ and approximately Gaussian null Z-scores.

## References

Liu, Z. and Lin, X. (2018). Multiple phenotype association tests using
summary statistics in genome-wide association studies. *Biometrics*,
74(1), 165-175.
[doi:10.1111/biom.12735](https://doi.org/10.1111/biom.12735)

Wu, M. C., Lee, S., Cai, T., Li, Y., Boehnke, M. and Lin, X. (2011).
Rare variant association testing for sequencing data with the sequence
kernel association test (SKAT). *American Journal of Human Genetics*,
89(1), 82-93.
[doi:10.1016/j.ajhg.2011.05.029](https://doi.org/10.1016/j.ajhg.2011.05.029)

Davies, R. B. (1980). Algorithm AS 155: The distribution of a linear
combination of \\\chi^2\\ random variables. *Journal of the Royal
Statistical Society, Series C (Applied Statistics)*, 29(3), 323-333.
[doi:10.2307/2346911](https://doi.org/10.2307/2346911)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

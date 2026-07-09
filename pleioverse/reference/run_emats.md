# EMATS leading-eigenvector adaptive omnibus association test

Multi-phenotype association test that adaptively blends a
variance-component (Wald / SKAT-style) statistic with a
leading-eigenvector directional statistic over a grid of mixing weights,
following the SKAT-O optimal unified framework (Lee et al. 2012; Wu et
al. 2011). It follows the same adaptive framework as
[`run_cmats()`](https://broccolito.github.io/pleioverse/reference/run_cmats.md),
but uses a 6-point grid and projects the Z-scores onto the principal
axis of \\\Sigma\\ rather than onto the burden direction.

## Usage

``` r
run_emats(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$emats` populated (one p-value per
variant).

## Details

For each variant let \\z=(z_1,\dots,z_K)'\\ be the vector of per-trait
Z-scores, \\z_k=\hat\beta_k/\mathrm{se}\_k\\, and let \\\Sigma\\ be the
\\K\times K\\ trait correlation matrix with eigendecomposition \\\Sigma
= \sum\_{j} \lambda_j v_j v_j'\\, eigenvalues sorted descending so that
\\\lambda_1\\ is the largest and \\v_1\\ its eigenvector.

The variance-component (Wald) term is the Mahalanobis quadratic form
\$\$\chi_1 = z'\Sigma^{-1}z,\$\$ which under the global null is
\\\chi^2\\ with \\K\\ degrees of freedom. The directional term is the
standardized squared projection of \\z\\ onto the leading eigenvector,
\$\$\chi_2 = \frac{(v_1'z)^2}{\lambda_1},\$\$ which under the null is
\\\chi^2\\ with \\1\\ degree of freedom and captures association
concentrated along the dominant correlation axis.

For each \\\rho\\ on the grid \\\\0,0.2,0.4,0.6,0.8,1\\\\ (6 equally
spaced points) the mixed statistic \$\$Q\_\rho = (1-\rho)\\\chi_1 +
\rho\\\chi_2\$\$ is formed, with null tail probability \\p\_\rho\\
computed exactly at the endpoints (\\\chi^2_K\\ at \\\rho=0\\,
\\\chi^2_1\\ at \\\rho=1\\) and otherwise as the tail of a weighted sum
of independent \\\chi^2_1\\ variables with spectrum \\(1-\rho)\\ of
multiplicity \\K-1\\ plus a single eigenvalue \\1\\ (Liu moment-matching
approximation). The statistic is the minimum p-value \\T=\min\_\rho
p\_\rho\\, calibrated with the SKAT-O one-dimensional null integral
\$\$p = \min\_\rho p\_\rho + \int_0^{t\_{\max}} \Pr\Big\[\chi^2\_{K-1}
\> \min\_{i:\\\rho_i\<1}\frac{q_i - x(t)}{1-\rho_i}\Big\]\\dt,\$\$ where
\\q_i\\ is the mixed-null quantile at level \\1-\min\_\rho p\_\rho\\ for
grid point \\i\\, \\x(t)\\ the \\\chi^2_1\\ quantile at \\t\\, and
\\t\_{\max}\\ the upper limit from the \\\rho=1\\ quantile. The integral
is evaluated by composite Simpson's rule and the result clamped to
\\\[0,1\]\\.

## Assumptions

Z-scores are approximately multivariate normal with mean zero and
covariance \\\Sigma\\ under the global null of no association.
\\\Sigma\\ is the (full-rank, invertible) trait correlation matrix,
shared across variants and treated as known, with a well-defined leading
eigenpair \\(\lambda_1, v_1)\\. Variants are tested one at a time.

## References

Wu MC, Lee S, Cai T, Li Y, Boehnke M, Lin X (2011). Rare-variant
association testing for sequencing data with the sequence kernel
association test (SKAT). *Am J Hum Genet* 89(1):82-93.
[doi:10.1016/j.ajhg.2011.05.029](https://doi.org/10.1016/j.ajhg.2011.05.029)

Lee S, Wu MC, Lin X (2012). Optimal tests for rare variant effects in
sequencing association studies. *Biostatistics* 13(4):762-775.
[doi:10.1093/biostatistics/kxs014](https://doi.org/10.1093/biostatistics/kxs014)

Liu H, Tang Y, Zhang HH (2009). A new chi-square approximation to the
distribution of non-negative definite quadratic forms in non-central
normal variables. *Comput Stat Data Anal* 53(4):853-856.
[doi:10.1016/j.csda.2008.11.025](https://doi.org/10.1016/j.csda.2008.11.025)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md),
[`run_cmats()`](https://broccolito.github.io/pleioverse/reference/run_cmats.md)

# Principal-component Fisher-combination test (PCFisher)

Multi-trait association test that projects the per-trait Z-scores onto
the principal components of the trait-correlation matrix and combines
the resulting per-PC p-values with Fisher's method. It targets
pleiotropic signal that is spread across several independent (whitened)
directions of the trait space.

## Usage

``` r
run_pcfisher(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$pcfisher` populated (one p-value per
variant).

## Details

Let \\z=(z_1,\dots,z_K)'\\ be the vector of per-trait Z-scores at a
variant, \\z_k=\hat\beta_k/\mathrm{se}\_k\\, and let \\\Sigma\\ be the
\\K\times K\\ between-trait correlation matrix with spectral
decomposition \$\$\Sigma = V \Lambda V', \qquad
\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_K),\$\$ where
\\V=(v_1,\dots,v_K)\\ holds the eigenvectors (engine `evec`) and
\\\lambda_i\>0\\ the eigenvalues (engine `eval`). The principal
components of the Z-vector are \$\$u = V'z, \qquad u_i = v_i' z.\$\$
Under the global null \\z\sim N(0,\Sigma)\\ the PCs are independent with
\\u_i\sim N(0,\lambda_i)\\, so each standardized PC yields a chi-square
statistic on one degree of freedom and its upper-tail p-value \$\$Q_i =
\frac{u_i^2}{\lambda_i}\\ \sim\\ \chi^2_1, \qquad p_i = \Pr(\chi^2_1 \>
Q_i).\$\$ Fisher's combination of the \\K\\ independent p-values gives
the statistic \$\$T = -2\sum\_{i=1}^{K}\log p_i\\ \sim\\
\chi^2\_{2K},\$\$ and the reported p-value is \\p = \Pr(\chi^2\_{2K} \>
T)\\.

## Assumptions

Whitening the Z-vector by the eigenvectors of \\\Sigma\\ renders the PCs
mutually independent under the null, which is what makes the Fisher
combination of the per-PC chi-square p-values valid and yields the exact
\\\chi^2\_{2K}\\ null. Because every PC enters symmetrically on the log
scale, PCFisher is most powerful when the association signal is
distributed across multiple PCs (diffuse pleiotropy) and tends to be
less powerful than linear PC combinations when the signal concentrates
in a single direction. The trait-correlation matrix \\\Sigma\\ is
assumed known (estimated under the null) and positive definite.

## References

Liu Z, Lin X (2019). A Geometric Perspective on the Power of Principal
Component Association Tests in Multiple Phenotype Studies. *Journal of
the American Statistical Association* 114(527):975–990.
[doi:10.1080/01621459.2018.1513363](https://doi.org/10.1080/01621459.2018.1513363)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

# Principal-component minimum-p test (PCMinP)

Multi-trait association test that takes the smallest
per-principal-component p-value and corrects it for the \\K\\
comparisons with a Sidak adjustment. It targets sparse pleiotropy in
which the signal concentrates in one (or a few) whitened directions of
the trait space.

## Usage

``` r
run_pcminp(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$pcminp` populated (one p-value per
variant).

## Details

Let \\z=(z_1,\dots,z_K)'\\ be the per-trait Z-scores and let
\\\Sigma=V\Lambda V'\\ be the spectral decomposition of the trait
correlation matrix, with eigenvectors \\V=(v_1,\dots,v_K)\\ (engine
`evec`) and eigenvalues \\\lambda_i\>0\\ (engine `eval`). The principal
components are \\u=V'z\\, \\u_i=v_i'z\\, and each standardized PC gives
a one-degree-of-freedom statistic and p-value \$\$Q_i =
\frac{u_i^2}{\lambda_i}\\ \sim\\ \chi^2_1, \qquad p_i = \Pr(\chi^2_1 \>
Q_i).\$\$ Let \\p\_{\min}=\min\_{i=1,\dots,K} p_i\\. Because the PCs are
independent under the null, the p-values \\p_i\\ are independent
\\U(0,1)\\ and the minimum is corrected by the exact Sidak formula \$\$p
= 1-\big(1-p\_{\min}\big)^{K},\$\$ which is the reported p-value.

## Assumptions

Eigenvector whitening makes the \\K\\ per-PC chi-square statistics
mutually independent under the null \\z\sim N(0,\Sigma)\\, so each
\\p_i\sim U(0,1)\\ independently and the Sidak correction
\\1-(1-p\_{\min})^{K}\\ is exact (not merely conservative as Bonferroni
would be). PCMinP is most powerful under sparse signal, where the
association loads on a single PC, and is less powerful than PCFisher or
PCLC when the signal is diffuse across many PCs. The trait-correlation
matrix \\\Sigma\\ is assumed known under the null and positive definite.

## References

Liu Z, Lin X (2019). A Geometric Perspective on the Power of Principal
Component Association Tests in Multiple Phenotype Studies. *Journal of
the American Statistical Association* 114(527):975–990.
[doi:10.1080/01621459.2018.1513363](https://doi.org/10.1080/01621459.2018.1513363)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

# TATES (extended Simes) multi-phenotype test

Trait-based Association Test that uses an Extended Simes procedure (van
der Sluis et al. 2013) to combine the per-trait univariate p-values of a
variant into a single multi-phenotype p-value while correcting for the
correlation between traits through an effective number of independent
tests.

## Usage

``` r
run_tates(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$tates` populated.

## Details

For a variant with per-trait p-values \\p_1,\dots,p_K\\ and trait
correlation matrix \\R\\ (entries \\r\_{ij}\\), TATES first maps the
trait correlations to the correlations between the corresponding
\\-\log\_{10}\\ p-values using the degree-6 polynomial \\f\\ fitted by
van der Sluis et al., \$\$\tilde r\_{ij} = f(r\_{ij}) =
\sum\_{m=0}^{6}\beta_m\\ r\_{ij}^{\\m} \quad (i\neq j), \qquad \tilde
r\_{ii}=1,\$\$ giving a transformed correlation matrix \\\tilde R\\.

The p-values are sorted ascending, \\p\_{(1)}\le\dots\le p\_{(K)}\\, and
\\\tilde R\\ is permuted to the same order. For any sub-collection the
effective number of independent tests is computed from the eigenvalues
\\\lambda_1,\dots,\lambda_m\\ of the corresponding submatrix of \\\tilde
R\\ as a fractional Galwey/Cheverud-type count \$\$m\_{\mathrm{eff}} =
m - \sum\_{i=1}^{m}\big(\lambda_i-1\big)\\ \mathbf{1}(\lambda_i\>1),\$\$
where \\m\_{e}\\ denotes the effective number over all \\K\\ traits and
\\m\_{e_j}\\ the effective number over the top \\j\\ ordered traits. The
extended-Simes statistic is the minimum correlation-adjusted ordered
p-value \$\$p\_{\mathrm{TATES}} = \min\_{1\le j\le
K}\frac{m_e}{m\_{e_j}}\\p\_{(j)}.\$\$ When all traits are independent
(\\\tilde R=I\\) this reduces to the ordinary Simes combined p-value.

## Assumptions

The per-trait p-values are valid (uniform on \\\[0,1\]\\ under each
trait's null). The trait correlation matrix \\R\\ approximates the
correlation structure of the test statistics and is shared across
variants. The polynomial \\f\\ is the fixed degree-6 approximation of
van der Sluis et al. relating phenotype correlations to p-value
correlations; non-finite p-values are treated as 1.

## References

van der Sluis S, Posthuma D, Dolan CV (2013). TATES: efficient
multivariate genotype-phenotype analysis for genome-wide association
studies. *PLoS Genet* 9(1):e1003235.
[doi:10.1371/journal.pgen.1003235](https://doi.org/10.1371/journal.pgen.1003235)

Li MX, Gui HS, Kwan JSH, Sham PC (2011). GATES: a rapid and powerful
gene-based association test using extended Simes procedure. *Am J Hum
Genet* 88(3):283-293.
[doi:10.1016/j.ajhg.2011.01.019](https://doi.org/10.1016/j.ajhg.2011.01.019)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

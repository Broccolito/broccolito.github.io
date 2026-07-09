# metaCCA canonical-correlation association test

Summary-statistic canonical correlation analysis (Cichonska et al. 2016)
between a single genotype and the \\K\\ phenotypes, using only GWAS
summary statistics. For each variant the squared canonical correlation
between the SNP and the multivariate phenotype is converted to a Wilks'
lambda chi-square p-value.

## Usage

``` r
run_metacca(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$metacca` populated.

## Details

Let \\\hat\beta\\ and \\\mathrm{se}\\ be the \\n\times K\\ matrices of
univariate effect sizes and standard errors, \\N\\ the (mean) per-trait
sample size, and \\Z=\hat\beta/\mathrm{se}\\ the matrix of Z-scores. The
phenotypic correlation matrix \\S\_{yy}\\ is estimated once as the
across-variant correlation of the effect-size columns,
\\S\_{yy}=\mathrm{cor}(\hat\beta)\\ (the metaCCA `estimateSyy` step).

For a variant the standardized SNP-phenotype cross-covariance is the row
of Z-scores rescaled by \\\sqrt{N}\\, \\z = Z\_{s\cdot}/\sqrt{N}\\,
playing the role of \\S\_{xy}\\. With genotype variance standardized to
one (\\S\_{xx}=1\\), the single non-trivial squared canonical
correlation is the quadratic form \$\$r_c^2 = z'\\S\_{yy}^{-1}\\z.\$\$
The corresponding Wilks' lambda is \\\Lambda = 1 - r_c^2\\, transformed
to a Bartlett-corrected chi-square statistic \$\$X = -\Big\[(N-1) -
\tfrac{1}{2}(K+2)\Big\]\\\log\Lambda,\$\$ whose null distribution is
\\\chi^2\\ with \\K\\ degrees of freedom; the p-value is its upper tail.
If \\\Lambda\le 0\\ (i.e. \\r_c^2\ge 1\\) the p-value is set to 0.

## Assumptions

Under the null of no association \\X\\ is approximately \\\chi^2_K\\ via
the Bartlett-corrected Wilks' lambda approximation, which assumes
multivariate-normal traits and large \\N\\. \\S\_{yy}\\ must be
invertible and is taken as the across-variant effect-size correlation;
the genotype is standardized (\\S\_{xx}=1\\) and a common sample size
\\N\\ is used for all traits. Each variant is tested separately against
the full phenotype set.

## References

Cichonska A, Rousu J, Marttinen P, Kangas AJ, Soininen P, Lehtimaki T,
Raitakari OT, Jarvelin MR, Salomaa V, Ala-Korpela M, Ripatti S, Pirinen
M (2016). metaCCA: summary statistics-based multivariate meta-analysis
of genome-wide association studies using canonical correlation analysis.
*Bioinformatics* 32(13):1981-1989.
[doi:10.1093/bioinformatics/btw052](https://doi.org/10.1093/bioinformatics/btw052)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

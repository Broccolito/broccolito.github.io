# pleioverse

Fast, self-contained simulation and testing for pleiotropic genetic
data.

`pleioverse` is an R package for simulating multi-phenotype genetic
architectures and running a broad suite of multi-phenotype association
and pleiotropy tests. The package is organized around one S4 `pleio`
object, so the typical workflow is simple:

1.  simulate data with
    [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.html);
2.  test the simulated variants with
    [`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.html),
    or run selected methods one at a time;
3.  inspect the accumulated result columns in `pleio@results`.

The aliases
[`pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.html)
and
[`pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.html)
call the same underlying functions.

## Installation

Install from GitHub:

``` r

install.packages("remotes")
remotes::install_github("broccolito/pleioverse")
```

You will need R package-building tools, a C++17-capable compiler, and
the `Rcpp` and `RcppArmadillo` packages. On macOS this usually means
installing Xcode Command Line Tools. On Linux, install the standard R
development toolchain for your distribution.

## Usage

``` r

library(pleioverse)

set.seed(1)

pleio <- run_pleiosim(
  n_phenotype = 3,
  n_participant = 10000,
  n_variant_pleiotropic = 10,
  n_variant_nonpleiotropic = c(10, 10, 10),
  n_variant_null = 960,
  heritability = c(0.1, 0.2, 0.3)
)

pleio <- run_pleiotest(pleio)
head(pleio@results)
```

Run only the methods you need:

``` r

pleio <- run_pleiosim(
  n_phenotype = 2,
  n_variant_nonpleiotropic = c(10, 10),
  heritability = c(0.2, 0.2)
)

pleio <- pleio |>
  run_wald() |>
  run_vc() |>
  run_placo(i = 1, j = 2)

head(pleio@results)
```

Combine existing p-value columns with
[`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.html):

``` r

pleio <- run_acat(pleio, cols = c("wald", "vc", "placo"))
head(pleio@results$acat)
```

## Documentation Map

Use the articles when you want a guided workflow, and the reference
pages when you want arguments, return values, and implementation
details.

| Need | Start here |
|----|----|
| Simulate genotypes, phenotypes, cohort overlap, and correlations | [Simulation articles](https://broccolito.github.io/pleioverse/articles/simulation-genotypes.html) |
| Run the full test suite and interpret output columns | [Running tests](https://broccolito.github.io/pleioverse/articles/running-tests.html) |
| Combine selected method results with ACAT | [Combining tests with ACAT](https://broccolito.github.io/pleioverse/articles/combining-tests-acat.html) |
| Read one article per statistical method | [Method articles](https://broccolito.github.io/pleioverse/articles/) |
| Look up function arguments and return values | [Function reference](https://broccolito.github.io/pleioverse/reference/) |

## Main Functions

| Task | Functions |
|----|----|
| Simulate a complete `pleio` object | [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.html) |
| Build simulation input matrices | [`create_heritable_correlation_matrix()`](https://broccolito.github.io/pleioverse/reference/create_heritable_correlation_matrix.html), [`create_nonheritable_correlation_matrix()`](https://broccolito.github.io/pleioverse/reference/create_nonheritable_correlation_matrix.html), [`create_cohort_makeup_matrix()`](https://broccolito.github.io/pleioverse/reference/create_cohort_makeup_matrix.html), [`create_eaf_matrix()`](https://broccolito.github.io/pleioverse/reference/create_eaf_matrix.html), [`create_efs_matrix_template()`](https://broccolito.github.io/pleioverse/reference/create_efs_matrix_template.html) |
| Run the default suite | [`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.html) |
| Run common multiple-association tests | [`run_wald()`](https://broccolito.github.io/pleioverse/reference/run_wald.html), [`run_vc()`](https://broccolito.github.io/pleioverse/reference/run_vc.html), [`run_sum()`](https://broccolito.github.io/pleioverse/reference/run_sum.html), [`run_minp()`](https://broccolito.github.io/pleioverse/reference/run_minp.html), [`run_tates()`](https://broccolito.github.io/pleioverse/reference/run_tates.html), [`run_metacca()`](https://broccolito.github.io/pleioverse/reference/run_metacca.html) |
| Run pleiotropy-focused tests | [`run_placo()`](https://broccolito.github.io/pleioverse/reference/run_placo.html), [`run_pgctest()`](https://broccolito.github.io/pleioverse/reference/run_pgctest.html), [`run_maiup()`](https://broccolito.github.io/pleioverse/reference/run_maiup.html), [`run_conjfdr()`](https://broccolito.github.io/pleioverse/reference/run_conjfdr.html), [`run_polarmorphism()`](https://broccolito.github.io/pleioverse/reference/run_polarmorphism.html) |
| Aggregate completed result columns | [`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.html) |

## Result Object

[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.html)
returns an S4 object of class `pleio`. It stores the simulation
settings, genotype and phenotype products, summary statistics, and the
result table populated by
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.html)
or by individual method functions.

Important result groups are encoded in the variant names:

- `pleio_variant*`: variants associated with all phenotypes;
- `nonpleio_variant*`: variants associated with one phenotype;
- `null_variant*`: variants with no simulated phenotype effect.

## Development

Install the package locally from the repository root:

``` sh
R CMD INSTALL .
```

Run the test suite:

``` sh
Rscript -e 'testthat::test_local(".")'
```

Build the documentation website:

``` sh
Rscript -e 'pkgdown::build_site()'
```

The package uses C++17, `Rcpp`, and `RcppArmadillo` for the
computational core. Files in `dev/` are local parity, baseline, and
timing scripts; they are not part of the installed R package.

## Method Sources

The statistical methods implemented or mirrored by `pleioverse` are
grounded in the following linked sources:

- Liu and Lin (2018), multiple phenotype association tests using summary
  statistics in GWAS, *Biometrics*.
  [DOI](https://doi.org/10.1111/biom.12735)
- Liu and Lin (2019), principal-component association tests in multiple
  phenotype studies, *Journal of the American Statistical Association*.
  [DOI](https://doi.org/10.1080/01621459.2018.1513363)
- Luo et al. (2020), MTAR multi-trait rare-variant summary-statistic
  analysis, *Nature Communications*.
  [DOI](https://doi.org/10.1038/s41467-020-16591-0)
- Liu and Xie (2019), Cauchy combination test, *Journal of the American
  Statistical Association*.
  [DOI](https://doi.org/10.1080/01621459.2018.1554485)
- Liu et al. (2019), ACAT rare-variant p-value combination, *American
  Journal of Human Genetics*.
  [DOI](https://doi.org/10.1016/j.ajhg.2019.01.002)
- Wu et al. (2011), SKAT, *American Journal of Human Genetics*.
  [DOI](https://doi.org/10.1016/j.ajhg.2011.05.029)
- Lee, Wu, and Lin (2012), optimal rare-variant tests, *Biostatistics*.
  [DOI](https://doi.org/10.1093/biostatistics/kxs014)
- van der Sluis, Posthuma, and Dolan (2013), TATES, *PLoS Genetics*.
  [DOI](https://doi.org/10.1371/journal.pgen.1003235)
- Li et al. (2011), GATES extended Simes procedure, *American Journal of
  Human Genetics*. [DOI](https://doi.org/10.1016/j.ajhg.2011.01.019)
- Cichonska et al. (2016), metaCCA, *Bioinformatics*.
  [DOI](https://doi.org/10.1093/bioinformatics/btw052)
- Ray and Chatterjee (2020), PLACO, *PLoS Genetics*.
  [DOI](https://doi.org/10.1371/journal.pgen.1009218)
- Andreassen et al. (2013), conditional FDR, *PLoS Genetics*.
  [DOI](https://doi.org/10.1371/journal.pgen.1003455)
- Smeland et al. (2020), conditional FDR for shared loci, *Human
  Genetics*. [DOI](https://doi.org/10.1007/s00439-019-02060-2)
- von Berg et al. (2022), PolarMorphism, *Bioinformatics*.
  [DOI](https://doi.org/10.1093/bioinformatics/btac228)
- Kessy, Lewin, and Strimmer (2018), whitening and decorrelation, *The
  American Statistician*.
  [DOI](https://doi.org/10.1080/00031305.2016.1277159)
- Conneely and Boehnke (2007), correlated-test p-value adjustment,
  *American Journal of Human Genetics*.
  [DOI](https://doi.org/10.1086/522036)
- Davies (1980), linear combinations of chi-square random variables,
  *Applied Statistics*. [DOI](https://doi.org/10.2307/2346911)
- Liu, Tang, and Zhang (2009), chi-square approximation for quadratic
  forms, *Computational Statistics & Data Analysis*.
  [DOI](https://doi.org/10.1016/j.csda.2008.11.025)
- Zeng et al., MAIUP software reference.
  [Repository](https://github.com/biostatpzeng/MAIUP)

## License

MIT. See `LICENSE`.

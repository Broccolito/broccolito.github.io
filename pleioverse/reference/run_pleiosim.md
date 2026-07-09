# Simulate pleiotropic genetic data

Fast C++ reimplementation of the `pleiosim` simulator. Generates
genotypes and multi-phenotype data for a configurable pleiotropic
architecture and returns a
[pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
object carrying per-trait GWAS summary statistics, ready for
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md).
`pleiosim()` is a synonym.

## Usage

``` r
run_pleiosim(
  n_phenotype = 3,
  heritable_correlation = 0.4,
  nonheritable_correlation = 0.2,
  n_participant = 10000,
  eaf = 0.4,
  n_variant_pleiotropic = 10,
  n_variant_nonpleiotropic = c(10, 10, 10),
  n_variant_null = 960,
  heritability = c(0.1, 0.2, 0.3),
  unique_cohorts = TRUE,
  crosstrait_heterogeneity = TRUE,
  withintrait_heterogeneity = TRUE,
  random_crosstrait_heterogeneity = FALSE,
  random_withintrait_heterogeneity = FALSE,
  customized_heritable_correlation_matrix = NULL,
  customized_nonheritable_correlation_matrix = NULL,
  customized_cohort_makeup_matrix = NULL,
  customized_eaf_matrix = NULL,
  customized_efs_matrix_template = NULL,
  customized_efs_matrix = NULL
)

pleiosim(
  n_phenotype = 3,
  heritable_correlation = 0.4,
  nonheritable_correlation = 0.2,
  n_participant = 10000,
  eaf = 0.4,
  n_variant_pleiotropic = 10,
  n_variant_nonpleiotropic = c(10, 10, 10),
  n_variant_null = 960,
  heritability = c(0.1, 0.2, 0.3),
  unique_cohorts = TRUE,
  crosstrait_heterogeneity = TRUE,
  withintrait_heterogeneity = TRUE,
  random_crosstrait_heterogeneity = FALSE,
  random_withintrait_heterogeneity = FALSE,
  customized_heritable_correlation_matrix = NULL,
  customized_nonheritable_correlation_matrix = NULL,
  customized_cohort_makeup_matrix = NULL,
  customized_eaf_matrix = NULL,
  customized_efs_matrix_template = NULL,
  customized_efs_matrix = NULL
)
```

## Arguments

- n_phenotype:

  number of phenotypes \\K\\.

- heritable_correlation, nonheritable_correlation:

  common (uniform) genetic correlation and environmental-noise input
  value.

- n_participant:

  participants per cohort.

- eaf:

  common effect allele frequency.

- n_variant_pleiotropic:

  number of pleiotropic variants (affect all traits).

- n_variant_nonpleiotropic:

  integer vector, non-pleiotropic variants per trait.

- n_variant_null:

  number of null variants.

- heritability:

  per-phenotype heritability vector (length `n_phenotype`).

- unique_cohorts:

  logical; if `TRUE` each trait uses a disjoint cohort.

- crosstrait_heterogeneity:

  logical; introduce deterministic cross-trait sign heterogeneity
  (alternate traits flipped).

- withintrait_heterogeneity:

  logical; introduce deterministic within-trait sign heterogeneity
  (alternate variants flipped).

- random_crosstrait_heterogeneity:

  logical; randomize the cross-trait flips.

- random_withintrait_heterogeneity:

  logical; randomize the within-trait flips.

- customized_heritable_correlation_matrix,
  customized_nonheritable_correlation_matrix,
  customized_cohort_makeup_matrix, customized_eaf_matrix,
  customized_efs_matrix_template, customized_efs_matrix:

  optional matrices overriding the corresponding constructed quantity;
  if `NULL` (default) the matrix is built from the scalar arguments
  above.

## Value

a
[pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
object with an empty `results` slot, ready for
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
or any single `run_<test>()`.

## Examples

``` r
pleio = run_pleiosim(n_phenotype = 2, n_participant = 500,
                      n_variant_pleiotropic = 5,
                      n_variant_nonpleiotropic = c(5, 5),
                      n_variant_null = 85, heritability = c(0.2, 0.2))
pleio
#> [pleio object]
#> A simulation of:
#> - 2 phenotypes
#> - 100 SNPs (5 pleiotropic, 10 non-pleiotropic, 85 null)
#> - 1000 participants total
#> - results: none yet (run run_pleiotest() or a run_<test>())
```

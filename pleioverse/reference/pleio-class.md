# The pleio simulation class

S4 container for a pleiotropic simulation produced by
[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).
Per-variant test results are accumulated in the `results` slot by
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
or any single `run_<test>()` function.

## Slots

- `n_phenotype`:

  number of phenotypes

- `n_variant_pleiotropic`:

  number of pleiotropic variants

- `n_variant_nonpleiotropic`:

  non-pleiotropic variants per phenotype

- `n_variant_null`:

  number of null variants

- `heritability`:

  per-phenotype heritability

- `unique_cohorts`:

  logical; disjoint vs overlapping cohorts

- `crosstrait_heterogeneity`:

  logical; cross-trait effect sign heterogeneity

- `withintrait_heterogeneity`:

  logical; within-trait effect sign heterogeneity

- `random_crosstrait_heterogeneity`:

  logical; randomized cross-trait flips

- `random_withintrait_heterogeneity`:

  logical; randomized within-trait flips

- `heritable_correlation_matrix`:

  target heritable correlation

- `nonheritable_correlation_matrix`:

  environmental correlation

- `cohort_makeup_matrix`:

  cohort sizes

- `eaf_matrix`:

  effect allele frequencies

- `efs_matrix_template`:

  effect-size template

- `efs_matrix`:

  fitted effect sizes

- `defacto_eaf_matrix`:

  realized eaf

- `defacto_sample_size_matrix`:

  realized sample sizes

- `estimated_heritable_correlation_matrix`:

  estimated from summary stats

- `cohort_genomes`:

  list of genotype matrices

- `cohort_phenotypes`:

  list of phenotype vectors

- `summary_stats_matrix`:

  list of per-phenotype summary statistics

- `results`:

  data.frame of per-variant test p-values (one row per variant, one
  column per test); populated by
  [`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
  or any `run_<test>()`.

# Cohorts, Overlap, and Custom Correlations

This article shows how cohort overlap, sample size, heritability, and
correlation matrices enter
[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Non-Overlapped Cohorts

With `unique_cohorts = TRUE`, a scalar `n_participant` is placed in each
singleton cohort column. No participants are placed in intersection
columns.

``` r

create_cohort_makeup_matrix(n_phenotype = 3, n_participant = 1000, unique_cohorts = TRUE)
```

## Fully Overlapped Cohorts

With `unique_cohorts = FALSE`, a scalar `n_participant` is placed in the
full intersection column. Those same participants contribute to every
phenotype-specific cohort.

``` r

create_cohort_makeup_matrix(n_phenotype = 3, n_participant = 1000, unique_cohorts = FALSE)
```

## Unequal Sample Sizes

When `n_participant` is a vector of length `K`, the helper places each
entry into the corresponding singleton cohort. This is a non-overlapped,
unequal-sample-size design.

``` r

create_cohort_makeup_matrix(n_phenotype = 3, n_participant = c(800, 1200, 1600))
```

## Partial Overlap

For partial overlap, build a custom `1 x (2^K - 1)` cohort-makeup matrix
and pass it through `customized_cohort_makeup_matrix`.

``` r

custom_cmm <- create_cohort_makeup_matrix(n_phenotype = 3, n_participant = 0)
custom_cmm[1, 'cohort1'] <- 500
custom_cmm[1, 'cohort2'] <- 700
custom_cmm[1, 'cohort3'] <- 900
custom_cmm[1, 'cohort1_cohort2'] <- 300
custom_cmm[1, 'cohort1_cohort2_cohort3'] <- 200
custom_cmm
```

The simulator creates one genotype matrix for every non-empty group and
then stacks groups into each phenotype cohort if the group contains that
phenotype.

## Different Heritable and Non-Heritable Inputs

You can use scalar uniform inputs or pass full custom matrices. The
heritable matrix is fitted against genetic correlation. The
non-heritable matrix is stored and used in the environmental-noise
transform, but the final phenotype construction does not guarantee that
its off-diagonal entries appear as observed residual correlations.

``` r

hcm <- create_heritable_correlation_matrix(3, 0.2)
hcm[1, 2] <- hcm[2, 1] <- 0.6
nhcm <- create_nonheritable_correlation_matrix(3, 0.05)
nhcm[2, 3] <- nhcm[3, 2] <- 0.25
pleio <- run_pleiosim(
  n_phenotype = 3,
  n_participant = 1000,
  heritability = c(0.1, 0.2, 0.3),
  customized_heritable_correlation_matrix = hcm,
  customized_nonheritable_correlation_matrix = nhcm,
  customized_cohort_makeup_matrix = custom_cmm
)
```

## Caveats

- Custom correlation matrices should be valid correlation matrices. The
  C++ engine clamps negative eigenvalues when forming the environmental
  square-root transform, but invalid matrices can still lead to
  simulations that are hard to interpret.
- `defacto_sample_size_matrix` is the realized per-phenotype sample size
  after stacking all groups that contain that phenotype.
- Cohort overlap changes the sharing of participants and genotypes
  across phenotype cohorts; it is separate from the genetic-correlation
  target and environmental-noise input matrix.
- If `n_participant` is a length-`K` vector,
  [`create_cohort_makeup_matrix()`](https://broccolito.github.io/pleioverse/reference/create_cohort_makeup_matrix.md)
  fills singleton cohorts only. Vector sample sizes do not automatically
  create overlap.
- The `pleio` print method reports the sum of per-phenotype realized
  sample sizes. In overlapping designs that sum double-counts shared
  participants, so do not read it as a unique-participant count.

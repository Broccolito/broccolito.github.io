# Build the cohort-makeup (sample-size) matrix

Lays out the number of participants in each cohort and
cohort-intersection, for either disjoint (`unique_cohorts = TRUE`) or
fully overlapping cohorts.

## Usage

``` r
create_cohort_makeup_matrix(
  n_phenotype = 3,
  n_participant = 1000,
  unique_cohorts = TRUE
)
```

## Arguments

- n_phenotype:

  number of phenotypes \\K\\.

- n_participant:

  participants per cohort (length 1, or length `n_phenotype`).

- unique_cohorts:

  logical; if `TRUE` each trait is measured in a disjoint cohort,
  otherwise all traits share one overlapping cohort.

## Value

a \\1 \times (2^K - 1)\\ matrix of cohort/intersection sizes.

# Build a uniform heritable (genetic) correlation matrix

Constructs a \\K \times K\\ correlation matrix with unit diagonal and a
common off-diagonal value, used as the target genetic correlation in
[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Usage

``` r
create_heritable_correlation_matrix(
  n_phenotype = 3,
  uniform_heritable_correlation = 0.4
)
```

## Arguments

- n_phenotype:

  number of phenotypes \\K\\.

- uniform_heritable_correlation:

  common off-diagonal genetic correlation.

## Value

a \\K \times K\\ correlation matrix.

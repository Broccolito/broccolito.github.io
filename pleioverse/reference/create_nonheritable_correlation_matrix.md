# Build a uniform non-heritable (environmental) input matrix

Constructs a \\K \times K\\ correlation matrix with unit diagonal and a
common off-diagonal value. In the current simulator this matrix is used
in the environmental-noise transform, but the final phenotype
construction uses one independently generated environmental column per
phenotype; it should not be read as a guaranteed realized residual
correlation matrix across final phenotypes.

## Usage

``` r
create_nonheritable_correlation_matrix(
  n_phenotype = 3,
  uniform_nonheritable_correlation = 0.2
)
```

## Arguments

- n_phenotype:

  number of phenotypes \\K\\.

- uniform_nonheritable_correlation:

  common off-diagonal value for the environmental-noise input matrix.

## Value

a \\K \times K\\ correlation matrix.

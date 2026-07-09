# Simulating Correlated Phenotypes

`pleioverse` separates the genetic-correlation target from the
non-heritable matrix used in the environmental-noise transform. They are
controlled separately, but only the genetic target is fitted against a
realized cross-phenotype quantity.

## Heritable Correlation

The default helper builds a uniform correlation matrix with ones on the
diagonal:

``` r

create_heritable_correlation_matrix(n_phenotype = 3, uniform_heritable_correlation = 0.4)
```

For a fitted effect matrix `B` and realized allele frequencies `p`, the
implementation computes the genetic covariance-like quantity

``` math
G_{ij} = \sum_{v=1}^{n_{var}} b_{iv} b_{jv} 2 p_{vi}(1 - p_{vi}),
```

then normalizes it to a correlation

``` math
R^{(g)}_{ij} = \frac{G_{ij}}{\sqrt{G_{ii}G_{jj}}}.
```

With the default effect-size template, the simulator adjusts the leading
all-nonzero pleiotropic columns so this realized genetic correlation is
close to the target unless a custom effect-size matrix is supplied. If
you pass a custom template, the same fitting step can only adjust the
leading columns that are nonzero for every phenotype.

## Non-Heritable Correlation

The non-heritable correlation matrix is used to build the
environmental-noise transform. The C++ engine eigendecomposes that
matrix and multiplies standard normal draws by a square-root transform:

``` math
E = U V \Lambda^{1/2},
```

where rows of `U` are independent standard normal draws. The
phenotype-specific environmental column is then standardized and
rescaled to match the requested heritability ratio.

One implementation detail matters for interpretation: the simulator
creates a fresh environmental-noise matrix separately for each phenotype
cohort and then uses only the phenotype-specific column for that cohort.
Because the off-diagonal columns from each draw are not carried into the
other final phenotypes, this matrix should not be interpreted as a
guaranteed observed residual correlation matrix, even when cohorts fully
overlap.

``` r

create_nonheritable_correlation_matrix(n_phenotype = 3, uniform_nonheritable_correlation = 0.2)
```

## Target Versus Realized Correlation

The target matrices are inputs. The object also records realized or
estimated quantities:

- `pleio@heritable_correlation_matrix`: target genetic correlation
  matrix;
- `pleio@nonheritable_correlation_matrix`: environmental-noise input
  matrix used by the simulator;
- `pleio@estimated_heritable_correlation_matrix`: estimated from
  generated summary statistics;
- `pleio@efs_matrix`: fitted effect sizes used to generate genetic
  values.

Small simulations can differ more from the target because genotype
draws, EAFs, environmental draws, and phenotypes are random. If
`customized_efs_matrix` is supplied, the simulator uses that matrix
directly and skips fitting pleiotropic effect multipliers to the
heritable-correlation target.

``` r

set.seed(13)
pleio <- run_pleiosim(
  n_phenotype = 3,
  n_participant = 1000,
  heritable_correlation = 0.5,
  nonheritable_correlation = 0.1,
  heritability = c(0.1, 0.2, 0.3)
)
pleio@heritable_correlation_matrix
pleio@nonheritable_correlation_matrix
pleio@estimated_heritable_correlation_matrix
```

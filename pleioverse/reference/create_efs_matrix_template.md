# Build the effect-size (EFS) template matrix

The \\K \times\\ `nVar` template of effect-direction indicators that
drives the simulated genetic architecture: pleiotropic variants load on
all traits, non-pleiotropic variants on exactly one trait, null variants
on none. The heterogeneity flags introduce sign flips to create
cross-trait (antagonistic across traits) and within-trait
(mixed-direction across variants) effect heterogeneity.

## Usage

``` r
create_efs_matrix_template(
  n_phenotype = 3,
  n_variant_pleiotropic = 10,
  n_variant_nonpleiotropic = c(10, 10, 10),
  n_variant_null = 960,
  crosstrait_heterogeneity = TRUE,
  withintrait_heterogeneity = TRUE,
  random_crosstrait_heterogeneity = FALSE,
  random_withintrait_heterogeneity = FALSE
)
```

## Arguments

- n_phenotype:

  number of phenotypes \\K\\.

- n_variant_pleiotropic:

  number of pleiotropic variants (affect all traits).

- n_variant_nonpleiotropic:

  integer vector, non-pleiotropic variants per trait.

- n_variant_null:

  number of null variants.

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

## Value

a \\K \times\\ `nVar` template matrix of signed effect indicators.

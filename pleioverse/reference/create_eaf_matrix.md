# Build the effect-allele-frequency (EAF) matrix

Per-variant, per-cohort target effect allele frequencies (here a common
value across all variants and cohorts).

## Usage

``` r
create_eaf_matrix(
  n_phenotype = 3,
  uniform_eaf = 0.4,
  n_variant_pleiotropic = 10,
  n_variant_nonpleiotropic = c(10, 10, 10),
  n_variant_null = 960
)
```

## Arguments

- n_phenotype:

  number of phenotypes \\K\\.

- uniform_eaf:

  common effect allele frequency.

- n_variant_pleiotropic:

  number of pleiotropic variants (affect all traits).

- n_variant_nonpleiotropic:

  integer vector, non-pleiotropic variants per trait.

- n_variant_null:

  number of null variants.

## Value

a `nVar` \\\times (2^K - 1)\\ matrix of effect allele frequencies.

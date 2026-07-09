# Simulating Genotypes

This article follows the implementation in `R/pleiosim.R` and
`src/simulate.cpp`. The genotype simulator is cohort-aware: it first
creates genotype matrices for each non-empty cohort/intersection group
and then assembles one genotype matrix per phenotype cohort.

## Variant Layout

The total number of variants is

``` math
n_{var} = n_{pleio} + \sum_{k=1}^{K} n_{nonpleio,k} + n_{null}.
```

Rows are named in simulation-truth order: `pleio_variant*`,
`nonpleio_variant*`, then `null_variant*`. The helper
[`create_eaf_matrix()`](https://broccolito.github.io/pleioverse/reference/create_eaf_matrix.md)
returns an `nVar x (2^K - 1)` matrix whose columns correspond to all
non-empty phenotype-cohort groups: singletons first, then pairwise and
higher-order intersections.

``` r

library(pleioverse)
eaf_matrix <- create_eaf_matrix(
  n_phenotype = 3, uniform_eaf = 0.35,
  n_variant_pleiotropic = 4,
  n_variant_nonpleiotropic = c(3, 3, 3),
  n_variant_null = 20
)
dim(eaf_matrix)
head(eaf_matrix[, 1:3])
```

## Genotype Draw

For group `g`, participant `i`, and variant `v`, the implementation
draws

``` math
G_{igv} \sim \mathrm{Binomial}(2, p_{vg}),
```

where `p_{vg}` is the corresponding entry of `eaf_matrix`. The draw is
done independently by group and variant. A group with zero participants
contributes an empty matrix.

## From Groups to Phenotype Cohorts

The cohort for phenotype `k` is built by stacking every group that
contains phenotype `k`. Therefore, an overlapping group contributes the
same participants to multiple phenotype-specific cohorts, while
singleton groups contribute only to their own phenotype.

The realized effect allele frequency stored in
`pleio@defacto_eaf_matrix` is

``` math
\hat p_{vk} = \frac{1}{2 n_k} \sum_{i=1}^{n_k} G_{ikv}.
```

`pleio@cohort_genomes[[k]]` is the realized genotype matrix for
phenotype `k`, and `pleio@defacto_sample_size_matrix` records the
realized sample size per phenotype.

## Minimal Example

``` r

set.seed(11)
pleio <- run_pleiosim(
  n_phenotype = 2, n_participant = 1000,
  n_variant_pleiotropic = 5,
  n_variant_nonpleiotropic = c(5, 5),
  n_variant_null = 20,
  eaf = 0.35,
  heritability = c(0.15, 0.15)
)
lapply(pleio@cohort_genomes, dim)
head(pleio@defacto_eaf_matrix)
```

## Caveats

- The values in `eaf_matrix` are targets; realized values fluctuate
  because genotypes are random binomial draws.
- The simulator does not model LD between variants. Variants are
  generated independently conditional on their EAFs.
- The default
  [`create_eaf_matrix()`](https://broccolito.github.io/pleioverse/reference/create_eaf_matrix.md)
  uses one uniform EAF, but
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md)
  can receive `customized_eaf_matrix` for variant- and cohort-specific
  targets.

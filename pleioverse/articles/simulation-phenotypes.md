# Simulating Phenotypes

Phenotypes in `pleioverse` are generated from a genetic component plus
an environmental component. The implementation is in `src/simulate.cpp`
and is exposed through
[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Effect-Size Template

The helper
[`create_efs_matrix_template()`](https://broccolito.github.io/pleioverse/reference/create_efs_matrix_template.md)
creates a `K x nVar` signed indicator matrix:

- pleiotropic variants load on all phenotypes;
- non-pleiotropic variants load on exactly one phenotype;
- null variants load on no phenotypes;
- heterogeneity flags flip signs deterministically or randomly.

``` r

template <- create_efs_matrix_template(
  n_phenotype = 3,
  n_variant_pleiotropic = 4,
  n_variant_nonpleiotropic = c(3, 3, 3),
  n_variant_null = 10,
  crosstrait_heterogeneity = TRUE,
  withintrait_heterogeneity = TRUE
)
template[, 1:8]
```

Unless `customized_efs_matrix` is supplied, the C++ simulator adjusts
the leading all-nonzero pleiotropic effect columns with a Nelder-Mead
search so the genetic correlation implied by the fitted effect matrix is
close to `heritable_correlation_matrix`.

## Genetic Value

For phenotype `k`, participant `i`, and variant `v`, the genetic
component is

``` math
g_{ik} = \sum_{v=1}^{n_{var}} G_{ikv} b_{kv},
```

where `G` comes from the phenotype-specific genotype matrix and `b_{kv}`
is the fitted effect size in `pleio@efs_matrix`.

## Environmental Component

For each phenotype cohort, the simulator draws an independent `n_k x K`
matrix of standard normal noise and transforms it by the square root of
the non-heritable input matrix. The phenotype-specific environmental
column is standardized and scaled to match the requested heritability:

``` math
\sigma^2_{e,k} = \mathrm{Var}(g_k) \frac{1 - h_k^2}{h_k^2}.
```

Here `heritability` in the API is the per-phenotype target used in this
variance ratio. The final phenotype is

``` math
y_{ik} = g_{ik} + e_{ik}.
```

Use `0 < heritability <= 1`. The current implementation uses the
variance ratio above directly and does not perform extensive validation
for invalid heritability values.

## Summary Statistics

For each phenotype and variant the simulator performs closed-form simple
regression. With centered sums `Sxx`, `Sxy`, and `Syy`, it computes

``` math
\hat\beta = \frac{S_{xy}}{S_{xx}}, \qquad
\widehat{\mathrm{se}} = \sqrt{\frac{(S_{yy} - \hat\beta S_{xy})/(n-2)}{S_{xx}}}.
```

The two-sided p-value is computed from the t statistic with `n - 2`
degrees of freedom. These per-phenotype matrices are stored in
`pleio@summary_stats_matrix` with columns `beta`, `se`, `p_value`, and
`eaf`.

``` r

set.seed(12)
pleio <- run_pleiosim(
  n_phenotype = 2,
  n_participant = 500,
  n_variant_nonpleiotropic = c(10, 10),
  heritability = c(0.2, 0.2)
)
names(pleio@summary_stats_matrix)
head(pleio@summary_stats_matrix[[1]])
```

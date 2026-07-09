# Understanding Output P-Values and Result Columns

All tests attach columns to `pleio@results`, with one row per variant.
Smaller values usually mean stronger evidence, but not every column has
the same null hypothesis or interpretation.

## Row Meaning

Rows follow the simulated truth labels produced by
[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md):

- `pleio_variant*`: variants that load on all phenotypes in the effect
  template;
- `nonpleio_variant*`: variants that load on exactly one phenotype;
- `null_variant*`: variants with no template effect.

## Multiple-Association Columns

The MPAT/MTAR columns, `tates`, and `metacca` ask whether a variant is
associated with any phenotype. Their usual global null is no association
with any trait. These are useful for discovery of multi-phenotype
association, not necessarily proof of horizontal pleiotropy.

## Pleiotropy-Specific Columns

`placo`, `pgctest`, and `polarmorphism` report p-values for
pleiotropy-focused hypotheses. `placo`, `pgctest`, `maiup`, and
`conjfdr` are pairwise for the selected `i, j` traits. `polarmorphism`
is joint across all phenotypes.

Important distinctions:

- `maiup` is an empirical FDR-style output based on a mixture-adjusted
  intersection-union statistic.
- `conjfdr` is a conjunctional conditional-FDR-style quantity, not an
  ordinary p-value.
- `polarmorphism` returns `NA` for variants that do not pass its radial
  q-value screen.
- [`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.md)
  creates an ACAT p-value from selected existing result columns.

## Pairwise Column Caveat

A pairwise method overwrites or populates a single column name such as
`placo`. If you need every pair, store each result under a separate
object or copy the column to a pair-specific name before running the
next pair.

``` r

pleio <- run_placo(pleio, i = 1, j = 2)
pleio@results$placo_1_2 <- pleio@results$placo
pleio <- run_placo(pleio, i = 1, j = 3)
pleio@results$placo_1_3 <- pleio@results$placo
```

## Practical Sorting

``` r

result_table <- pleio@results
head(result_table[order(result_table$wald), 'wald', drop = FALSE])
```

Use method-specific thresholds and downstream multiple-testing
correction appropriate to your design. The package computes per-variant
method outputs; it does not globally adjust all result columns for the
full number of variants and methods.

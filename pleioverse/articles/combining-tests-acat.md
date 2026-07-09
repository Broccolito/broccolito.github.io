# Combining Tests with ACAT

[`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.md)
combines p-value columns that already exist in `pleio@results`. It is
useful when several methods target complementary alternatives and you
want one analytic omnibus summary.

## Formula

For result columns with p-values `p_1, ..., p_M` and non-negative
weights `w_1, ..., w_M` normalized to sum to one, ACAT forms

``` math
T = \sum_{m=1}^{M} w_m \tan[(0.5 - p_m)\pi].
```

The combined p-value is the upper tail of a standard Cauchy:

``` math
p_{ACAT} = \frac{1}{2} - \frac{1}{\pi}\arctan(T).
```

In the C++ implementation, p-values below `1e-15` use the tail
approximation `1 / (pi * p)`, non-finite entries are skipped, and the
returned value is floored at `1e-300`.

## Usage

``` r

library(pleioverse)
set.seed(31)
pleio <- run_pleiosim(n_phenotype = 3, heritability = c(0.1, 0.2, 0.3))
pleio <- run_wald(pleio)
pleio <- run_vc(pleio)
pleio <- run_tates(pleio)
pleio <- run_acat(pleio, cols = c('wald', 'vc', 'tates'))
head(pleio@results$acat)
```

## Weighted Combination

``` r

pleio <- run_acat(
  pleio,
  cols = c('wald', 'vc', 'tates'),
  weights = c(1, 1, 2)
)
```

## Interpretation

ACAT assumes each combined input column is valid under its own null. It
is robust to dependence among the combined p-values, but it cannot
repair an invalid or misinterpreted input column. In this package, avoid
mixing ordinary p-values with FDR-style quantities such as `conjfdr` or
`maiup` unless you explicitly accept that interpretation.

The current wrapper is intentionally lightweight: it does not validate
the length or sign of supplied weights, and the default `cols = NULL`
means all current result columns. If an `acat` column is already
present, rerunning `run_acat(pleio)` without `cols` will include that
previous `acat` column.

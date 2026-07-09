# Running the Test Suite

The testing API works on a `pleio` object created by
[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md)
or by an object with compatible summary-statistic slots.

## Shared Inputs

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
calls the internal `.pleio_inputs()` helper once. It extracts
per-phenotype `beta`, `se`, and `p_value`, forms

``` math
Z_{vk} = \frac{\hat\beta_{vk}}{\widehat{\mathrm{se}}_{vk}},
```

sets non-finite Z-scores to zero, estimates the trait correlation matrix
as `cor(Z)`, and computes its eigendecomposition. The MPAT/MTAR methods
reuse those shared inputs through the C++ `TestEngine`.

## Full Workflow

``` r

library(pleioverse)
set.seed(21)
pleio <- run_pleiosim(
  n_phenotype = 3,
  n_participant = 1000,
  n_variant_pleiotropic = 5,
  n_variant_nonpleiotropic = c(5, 5, 5),
  n_variant_null = 100,
  heritability = c(0.1, 0.2, 0.3)
)
pleio <- run_pleiotest(pleio, mpat = TRUE, pleiotropy = TRUE, i = 1, j = 2)
head(pleio@results)
```

## Running One Method at a Time

Every individual method returns the modified `pleio` object, so calls
can be chained.

``` r

pleio <- run_wald(pleio)
pleio <- run_vc(pleio)
pleio <- run_placo(pleio, i = 1, j = 2)
colnames(pleio@results)
```

## Families in `run_pleiotest()`

- `mpat = TRUE` runs the 24 MPAT/MTAR methods plus `tates` and
  `metacca`.
- `pleiotropy = TRUE` runs `placo`, `pgctest`, `maiup`, `conjfdr`, and
  `polarmorphism`.
- [`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.md)
  is not part of
  [`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
  because it combines columns that already exist in `pleio@results`.

## Pairwise Pleiotropy Tests

[`run_placo()`](https://broccolito.github.io/pleioverse/reference/run_placo.md),
[`run_pgctest()`](https://broccolito.github.io/pleioverse/reference/run_pgctest.md),
[`run_maiup()`](https://broccolito.github.io/pleioverse/reference/run_maiup.md),
and
[`run_conjfdr()`](https://broccolito.github.io/pleioverse/reference/run_conjfdr.md)
use phenotype indices `i` and `j`. A single result column refers only to
the pair used in that call.

``` r

pair12 <- run_placo(pleio, i = 1, j = 2)
pair13 <- run_placo(pleio, i = 1, j = 3)
```

# Run the full test suite on a pleio object

The suite comprises two distinct families of tests, plus one aggregation
method:

- **Multiple-association testing methods** (the MPAT/MTAR family,
  metaCCA, TATES) ask whether a variant is associated with *any* of the
  phenotypes – a joint/omnibus test of association.

- **Pleiotropy tests** (PLACO+, PGCtest, MAIUP, conjFDR, PolarMorphism)
  ask specifically whether a variant is *pleiotropic*, i.e. associated
  with two or more phenotypes simultaneously.

- **ACAT** is a compounding method that aggregates existing p-value
  columns; it is intentionally excluded here and run separately via
  [`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.md).

This function computes every method in the first two families for each
variant and stores the per-variant p-values in the object's `results`
slot.

## Usage

``` r
run_pleiotest(pleio, mpat = TRUE, pleiotropy = TRUE, i = 1, j = 2)

pleiotest(pleio, mpat = TRUE, pleiotropy = TRUE, i = 1, j = 2)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

- mpat:

  logical; include the multiple-association testing methods (the
  MPAT/MTAR family, metaCCA, TATES). Default `TRUE`.

- pleiotropy:

  logical; include the pleiotropy tests (PLACO+, PGCtest, MAIUP,
  conjFDR, PolarMorphism). Default `TRUE`.

- i, j:

  phenotype indices for the pairwise pleiotropy tests. Default `1,2`.

## Value

the `pleio` object with `results` populated (one column per test).

## Details

For the MPAT/MTAR family, shared inputs – the Z matrix, the trait
correlation \\\Sigma\\ and its eigendecomposition – are computed once
and reused by the C++ test engine. The pairwise pleiotropy wrappers
extract the needed Z-scores from the summary-statistics matrices when
they run.

## See also

[`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md),
and the individual `run_<test>()` functions.

## Examples

``` r
pleio = run_pleiosim(n_phenotype = 2, n_participant = 500,
                      n_variant_pleiotropic = 5,
                      n_variant_nonpleiotropic = c(5, 5),
                      n_variant_null = 85, heritability = c(0.2, 0.2))
pleio = run_pleiotest(pleio)
head(pleio@results)
#>                           wald           wi           vc          sum
#> pleio_variant1    7.780652e-06 5.803365e-06 1.069751e-05 1.271191e-06
#> pleio_variant2    1.307327e-03 1.155317e-03 1.494957e-03 6.700419e-04
#> pleio_variant3    3.637434e-04 2.979702e-04 4.503020e-04 6.967541e-05
#> pleio_variant4    1.784351e-06 1.302143e-06 2.518577e-06 3.671382e-07
#> pleio_variant5    5.836247e-05 4.892039e-05 7.097903e-05 4.002051e-05
#> nonpleio_variant1 5.757491e-03 6.007485e-03 5.554656e-03 5.580045e-02
#>                           dsum     pcfisher         pclc       pcminp
#> pleio_variant1    1.271191e-06 1.512465e-05 3.740300e-04 2.542380e-06
#> pleio_variant2    6.700419e-04 1.275836e-03 9.330138e-04 1.339635e-03
#> pleio_variant3    6.967541e-05 6.640457e-04 4.068848e-03 1.393460e-04
#> pleio_variant4    3.671382e-07 2.651117e-06 3.968042e-05 7.342763e-07
#> pleio_variant5    4.002051e-05 5.624571e-05 5.646288e-05 8.003942e-05
#> nonpleio_variant1 5.580045e-02 4.687669e-03 6.076203e-01 1.966132e-02
#>                           minp mixfisher_davies mixfisher_liu mixfisher_liumod
#> pleio_variant1    5.392542e-04     1.512465e-05  1.512465e-05     1.512465e-05
#> pleio_variant2    1.622905e-03     1.275836e-03  1.275836e-03     1.275836e-03
#> pleio_variant3    6.472020e-03     6.640457e-04  6.640457e-04     6.640457e-04
#> pleio_variant4    5.604884e-05     2.651117e-06  2.651117e-06     2.651117e-06
#> pleio_variant5    9.276451e-05     5.624571e-05  5.624571e-05     5.624571e-05
#> nonpleio_variant1 3.043136e-03     4.687669e-03  4.687669e-03     4.687669e-03
#>                   mixtippett_davies mixtippett_liu mixtippett_liumod
#> pleio_variant1         2.542380e-06   2.542380e-06      2.542380e-06
#> pleio_variant2         1.339635e-03   1.339635e-03      1.339635e-03
#> pleio_variant3         1.393460e-04   1.393460e-04      1.393460e-04
#> pleio_variant4         7.342763e-07   7.342763e-07      7.342763e-07
#> pleio_variant5         8.003942e-05   8.003942e-05      8.003942e-05
#> nonpleio_variant1      1.966132e-02   1.966132e-02      1.966132e-02
#>                   mixvar_davies   mixvar_liu mixvar_liumod mixsd_davies
#> pleio_variant1     0.0006238723 0.0005817658  0.0006238723 7.780652e-06
#> pleio_variant2     0.3433479640 0.0083042258  0.0085203872 1.307327e-03
#> pleio_variant3     0.3531730294 0.0060453382  0.0062330778 3.637434e-04
#> pleio_variant4     0.0002256726 0.0002065489  0.0002256726 1.784351e-06
#> pleio_variant5     0.0010406097 0.0009792859  0.0010406097 5.836247e-05
#> nonpleio_variant1  0.3563836662 0.0052913454  0.0054670813 3.514744e-01
#>                      mixsd_liu mixsd_liumod       mixada        cmats
#> pleio_variant1    7.780652e-06 7.780652e-06 0.0008351238 0.0008347386
#> pleio_variant2    1.307327e-03 1.307327e-03 0.0016062190 0.0011262683
#> pleio_variant3    3.637434e-04 3.637434e-04 0.0009237167 0.0004425520
#> pleio_variant4    1.784351e-06 1.784351e-06 0.0008337970 0.0008336539
#> pleio_variant5    5.836247e-05 5.836247e-05 0.0008649561 0.0004380859
#> nonpleio_variant1 5.757491e-03 5.757491e-03 0.0121650370 0.0093963022
#>                          emats        tates      metacca        placo
#> pleio_variant1    0.0008347386 5.957600e-04 6.235112e-06 9.684339e-05
#> pleio_variant2    0.0011253353 1.746290e-03 1.240285e-03 1.730201e-02
#> pleio_variant3    0.0004307418 6.099873e-03 3.340546e-04 1.795619e-03
#> pleio_variant4    0.0008336560 6.629069e-05 1.334732e-06 4.742995e-05
#> pleio_variant5    0.0004269417 1.078913e-04 5.072440e-05 3.112482e-03
#> nonpleio_variant1 0.0093752791 3.232465e-03 5.634675e-03 3.004177e-01
#>                        pgctest       maiup     conjfdr polarmorphism
#> pleio_variant1    0.0002515104 0.008479552 0.002443257            NA
#> pleio_variant2    0.0200835465 0.437979969 0.143357702            NA
#> pleio_variant3    0.0028977404 0.026121419 1.000000000            NA
#> pleio_variant4    0.0001389254 0.009734534 0.006861222            NA
#> pleio_variant5    0.0046147951 0.248483665 0.140263454            NA
#> nonpleio_variant1 0.2725052383 0.887571471 0.400171948            NA
```

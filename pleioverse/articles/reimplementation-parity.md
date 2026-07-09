# Reimplementation Parity

The current `pleioverse` methods are written to match the original
statistical definitions while sharing more input preparation and, for
the computationally heavy paths, using C++ implementations. The parity
suite checks that the rewritten methods return the same numerical
results as independent reference implementations.

## Reference Design

The reference layer in `dev/reference_methods.R` is intentionally
separate from the package API. It rebuilds the original formulas and
numerical recipes in plain R, then compares those values with the public
`run_*()` methods and the combined
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
workflow.

The checks cover deterministic simulations with 2, 3, and 5 phenotypes,
multiple variant counts, and method-specific parameter choices.

| family | methods checked | parameter coverage |
|----|----|----|
| Multiple-association tests | Wald, WI, VC, SUM, DSUM, PCFisher, PCLC, PCMinP, MinP, MixFisher, MixTippett, MixVar, MixSD, MixAda, CMATs, EMATs | Davies, Liu, and modified Liu variants where applicable |
| Trait-level methods | TATES, metaCCA | Eigenvalue and effective-test calculations |
| Pairwise pleiotropy tests | PLACO, PGCtest, MAIUP, conjFDR | Phenotype pairs, thresholds, and empirical calibration paths |
| Directional pleiotropy | PolarMorphism | Radial screening, whitening, angular tails, and Mahalanobis mode |
| Result aggregation | ACAT | Column selection and optional weights |
| Workflow orchestration | [`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md) | Full MPAT/MTAR and pleiotropy result columns |

## Results

The current parity run completed all checks successfully:

``` text
Completed 242 parity checks.
PARITY OK
```

Most closed-form methods are checked at machine-level numerical
tolerance. Methods that depend on integration, quadrature, empirical
counting, or approximation grids use tight method-specific tolerances so
that numerical implementation details do not hide real statistical
differences:

| method group | tolerance note |
|----|----|
| Closed-form and eigenvalue methods | Near machine precision |
| PLACO and PGCtest | Dense tail-grid comparisons, within `3e-7` |
| PolarMorphism | Angular-tail series comparisons, within `5e-8` |
| MinP | Multivariate-normal integration, within `0.003` |

## Reproduce The Check

Run the parity script from the repository root:

``` r

Sys.setenv(PLEIOVERSE_LIB = "/tmp/pverse_lib")
system("Rscript dev/parity_check.R")
```

The check uses deterministic input generation, independent reference
calculations, and direct result-column comparisons. Passing the suite
means the rewritten implementations preserve the statistical behavior of
the original methods across the tested configurations.

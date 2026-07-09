# Reimplementation Speed

The rewritten methods keep the original statistical definitions while
reducing repeated setup work and moving the expensive numerical paths
into the current package implementation. The benchmark compares the
independent original/reference layer against the public `pleioverse` API
on the same deterministic inputs.

## Benchmark Design

The timing workflow has two layers:

1.  `representative_method` times every statistical method individually
    on a representative 3-phenotype, 24-variant case.
2.  `scaling_suite` times method families across variant-count and
    phenotype-count cases, then fits descriptive log-log slopes as
    empirical complexity readings.

Each benchmark cell was run twice in the recorded local benchmark.
Timings below `0.0005` seconds are reported at the timer floor, so
speedups for those cells are lower bounds or marked as unresolved when
both implementations are too fast to separate.

## Suite-Level Speedups

The full workflow is consistently faster in the rewritten
implementation:

| case | phenotypes | variants | original seconds | current seconds | speedup |
|----|----|---:|---:|---:|---:|
| Variant scaling, 12 variants | 3 | 12 | 33.7065 | 0.7880 | 42.77x |
| Variant scaling, 24 variants | 3 | 24 | 77.0610 | 2.2520 | 34.22x |
| Variant scaling, 36 variants | 3 | 36 | 58.5200 | 2.2730 | 25.75x |
| Phenotype scaling, 2 phenotypes | 2 | 24 | 66.8060 | 3.0090 | 22.20x |
| Phenotype scaling, 5 phenotypes | 5 | 24 | 64.5950 | 1.0250 | 63.02x |

The multi-association suite, which contains the MPAT/MTAR-style methods,
shows the largest constant-factor gains because shared matrices and
decompositions are reused instead of rebuilt method by method.

| case | phenotypes | variants | original seconds | current seconds | speedup |
|----|----|---:|---:|---:|---:|
| Variant scaling, 12 variants | 3 | 12 | 33.8165 | 0.4465 | 75.74x |
| Variant scaling, 24 variants | 3 | 24 | 132.6795 | 1.0875 | 122.00x |
| Variant scaling, 36 variants | 3 | 36 | 72.6680 | 2.3665 | 30.71x |
| Phenotype scaling, 2 phenotypes | 2 | 24 | 52.5180 | 3.2175 | 16.32x |
| Phenotype scaling, 5 phenotypes | 5 | 24 | 57.9170 | 0.4485 | 129.13x |

## Representative Method Speedups

The per-method benchmark shows the largest gains in methods that
previously spent most of their time in repeated integration, covariance
transforms, or dense calibration loops.

| method            | original seconds | current seconds |   speedup |
|-------------------|-----------------:|----------------:|----------:|
| MinP              |         101.1245 |          0.2125 |   475.88x |
| Wald              |           0.1510 |          0.0005 |   302.00x |
| MAIUP             |           0.0940 |        \<0.0005 | \>188.00x |
| MixFisher Davies  |          33.4510 |          0.2170 |   154.15x |
| WI                |           9.4940 |          0.0895 |   106.08x |
| VC                |           8.9620 |          0.0875 |   102.42x |
| MixVar Davies     |           7.4310 |          0.0865 |    85.91x |
| MixTippett Davies |          44.4835 |          0.6350 |    70.05x |
| PLACO             |           1.7225 |          0.0275 |    62.64x |
| PCMinP            |           0.0255 |        \<0.0005 |  \>51.00x |

## Complexity Reading

The empirical slopes are descriptive rather than formal asymptotic
proofs because the grid is deliberately small and deterministic. In the
recorded run, the current full
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
workflow scaled approximately linearly with variant count at fixed 3
phenotypes, with an empirical exponent of `1.02`. The current
multi-association suite had an empirical variant-count exponent of
`1.49` on the same grid.

Small cells can show noisy or even negative slopes when timer
resolution, numerical path variation, or integration convergence
dominates the measured wall time. The stable practical result is the
constant-factor reduction: the rewritten API avoids repeated input
preparation and executes the numerically heavy routines substantially
faster.

## Reproduce The Benchmark

Run the baseline first, then compare the current package implementation:

``` r

Sys.setenv(PLEIOVERSE_LIB = "/tmp/pverse_lib")
system("Rscript dev/baseline.R")
system("Rscript dev/timing_check.R")
```

The benchmark writes these artifacts in `dev/`:

| file | contents |
|----|----|
| `runtime_baseline.csv` | Original/reference raw timing rows |
| `runtime_comparison.csv` | Original/reference and current timing rows |
| `runtime_summary.csv` | Median, mean, minimum, maximum, and standard-deviation summaries |
| `runtime_report.md` | Markdown report with speedups, scaling tables, and interpretation notes |

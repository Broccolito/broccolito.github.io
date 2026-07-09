# ACAT aggregated-Cauchy omnibus over existing results

Combines the p-value columns already present in `pleio@results` using
the Aggregated Cauchy Association Test (Liu et al. 2019; Liu & Xie
2019), an omnibus that merges a set of (possibly dependent) p-values via
a Cauchy-transformed weighted sum. This is a meta-combination of
previously computed test results, not a test of the raw traits; it is
not part of
[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
and should be called afterwards.

## Usage

``` r
run_acat(pleio, cols = NULL, weights = NULL)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object with at least one test in `results`.

- cols:

  optional character vector of result columns to combine (default: all
  current result columns).

- weights:

  optional non-negative weights, one per combined column.

## Value

the `pleio` object with `results$acat` populated.

## Details

For each variant let \\p_1,\dots,p_M\\ be the p-values taken from the
selected result columns (one per chosen test), with non-negative weights
\\w_1,\dots,w_M\\ normalized to sum to one (\\w_k=1/M\\ when no weights
are supplied). ACAT maps each p-value to the upper tail of a standard
Cauchy and forms the weighted sum \$\$T = \sum\_{k=1}^{M}
w_k\\\tan\\\big\[(0.5 - p_k)\\\pi\big\].\$\$ Because each
\\\tan\[(0.5-p_k)\pi\]\\ is standard-Cauchy under the null and a
weighted sum of standard-Cauchy variables is again standard Cauchy
regardless of their dependence, the combined p-value is the standard
Cauchy upper tail \$\$p\_{\mathrm{ACAT}} = \tfrac{1}{2} -
\tfrac{1}{\pi}\arctan(T).\$\$ For numerical stability very small
p-values use the tail expansion \\\tan\[(0.5-p_k)\pi\]\approx 1/(\pi
p_k)\\ when \\p_k\<10^{-15}\\; non-finite p-values are skipped, and the
returned p-value is floored at \\10^{-300}\\. Weights are normalized
once for the requested column set; if a row has non-finite inputs,
skipped entries are not row-wise renormalized.

## Assumptions

Each combined p-value is (approximately) uniform on \\\[0,1\]\\ under
its own null. The Cauchy combination is robust to arbitrary correlation
among the combined statistics in the bulk of the distribution and
remains accurate in the extreme tail, so no correlation matrix is
needed. The inputs are existing result columns; their validity as
p-values is the caller's responsibility. The wrapper is intentionally
lightweight: it does not validate weight length, non-negativity,
positive weight sum, or that finite inputs are in \\\[0,1\]\\. If an
`acat` column already exists, calling `run_acat()` with `cols = NULL`
includes that existing column in the default set.

## References

Liu Y, Xie J (2019). Cauchy combination test: a powerful test with
analytic p-value calculation under arbitrary dependency structures. *J
Am Stat Assoc* 115(529):393-402.
[doi:10.1080/01621459.2018.1554485](https://doi.org/10.1080/01621459.2018.1554485)

Liu Y, Chen S, Li Z, Morrison AC, Boerwinkle E, Lin X (2019). ACAT: a
fast and powerful p-value combination method for rare-variant analysis
in sequencing studies. *Am J Hum Genet* 104(3):410-421.
[doi:10.1016/j.ajhg.2019.01.002](https://doi.org/10.1016/j.ajhg.2019.01.002)

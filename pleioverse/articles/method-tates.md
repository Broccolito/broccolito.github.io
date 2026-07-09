# TATES (extended Simes) multi-phenotype test

Function:
[`run_tates()`](https://broccolito.github.io/pleioverse/reference/run_tates.md)

Result column: `tates`

Family: Extended-Simes p-value combination

## Introduction

Trait-based Association Test that uses an Extended Simes procedure (van
der Sluis et al. 2013) to combine the per-trait univariate p-values of a
variant into a single multi-phenotype p-value while correcting for the
correlation between traits through an effective number of independent
tests.

## Mathematical Framework and Implementation

For a variant with per-trait p-values $`p_1,\dots,p_K`$ and trait
correlation matrix $`R`$ (entries $`r_{ij}`$), TATES first maps the
trait correlations to the correlations between the corresponding
$`-\log_{10}`$ p-values using the degree-6 polynomial $`f`$ fitted by
van der Sluis et al.,

``` math
\tilde r_{ij} = f(r_{ij}) = \sum_{m=0}^{6}\beta_m\, r_{ij}^{\,m}
  \quad (i\neq j), \qquad \tilde r_{ii}=1,
```

giving a transformed correlation matrix $`\tilde R`$.

The p-values are sorted ascending, $`p_{(1)}\le\dots\le p_{(K)}`$, and
$`\tilde R`$ is permuted to the same order. For any sub-collection the
effective number of independent tests is computed from the eigenvalues
$`\lambda_1,\dots,\lambda_m`$ of the corresponding submatrix of
$`\tilde R`$ as a fractional Galwey/Cheverud-type count

``` math
m_{\mathrm{eff}} = m - \sum_{i=1}^{m}\big(\lambda_i-1\big)\,
  \mathbf{1}(\lambda_i>1),
```

where $`m_{e}`$ denotes the effective number over all $`K`$ traits and
$`m_{e_j}`$ the effective number over the top $`j`$ ordered traits. The
extended-Simes statistic is the minimum correlation-adjusted ordered
p-value

``` math
p_{\mathrm{TATES}} = \min_{1\le j\le K}\frac{m_e}{m_{e_j}}\,p_{(j)}.
```

When all traits are independent ($`\tilde R=I`$) this reduces to the
ordinary Simes combined p-value.

## What This Method Advances

TATES advances over ordinary Simes by replacing the raw number of traits
with an effective number of tests induced by trait correlation.

## Assumptions

The per-trait p-values are valid (uniform on $`[0,1]`$ under each
trait’s null). The trait correlation matrix $`R`$ approximates the
correlation structure of the test statistics and is shared across
variants. The polynomial $`f`$ is the fixed degree-6 approximation of
van der Sluis et al. relating phenotype correlations to p-value
correlations; non-finite p-values are treated as 1.

## When This Method Is a Good Fit

Good for combining valid per-trait p-values while adjusting for
correlated phenotypes.

The implementation treats variants independently and relies on the
trait-correlation estimate derived from the available summary statistics
unless an individual method estimates a pairwise null moment internally.

## Minimal Usage

``` r

library(pleioverse)
set.seed(1)
pleio <- run_pleiosim(
  n_phenotype = 3, n_participant = 500,
  n_variant_pleiotropic = 3,
  n_variant_nonpleiotropic = c(3, 3, 3),
  n_variant_null = 30,
  heritability = c(0.1, 0.2, 0.3)
)
pleio <- run_tates(pleio)
head(pleio@results$tates)
```

## Implementation Files

- R wrapper: `R/run_tates.R`
- C++ kernel or dispatcher: `src/tates.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

van der Sluis S, Posthuma D, Dolan CV (2013). TATES: efficient
multivariate genotype-phenotype analysis for genome-wide association
studies. *PLoS Genet* 9(1):e1003235. DOI:
<https://doi.org/10.1371/journal.pgen.1003235>

Li MX, Gui HS, Kwan JSH, Sham PC (2011). GATES: a rapid and powerful
gene-based association test using extended Simes procedure. *Am J Hum
Genet* 88(3):283-293. DOI: <https://doi.org/10.1016/j.ajhg.2011.01.019>

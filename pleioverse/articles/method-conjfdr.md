# Conjunctional conditional FDR on a phenotype pair

Function:
[`run_conjfdr()`](https://broccolito.github.io/pleioverse/reference/run_conjfdr.md)

Result column: `conjfdr`

Family: Pairwise conjunctional FDR

## Introduction

Pairwise (traits $`i,j`$) conjunctional conditional false-discovery-rate
style quantity (conjFDR), motivated by the conditional FDR framework of
Andreassen et al. (2013). The implementation returns a maximum of two
conditional-FDR-like quantities computed with one shared empirical
denominator. It is intended to prioritize variants jointly associated
with *both* traits and is computed for one trait pair at a time.

## Mathematical Framework and Implementation

For variant $`v`$ the per-trait Z-scores $`Z_1, Z_2`$ are converted to
one-sided (absolute) normal p-values

``` math
p_1 = 1 - \Phi(|Z_1|), \qquad p_2 = 1 - \Phi(|Z_2|) .
```

**Conditional FDR.** The conditional FDR of trait 1 given trait 2 at the
observed pair $`(p_1, p_2)`$ is the local false discovery rate for trait
1 restricted to variants at least as significant as $`v`$ in trait 2.
Using the empirical conditional CDF, with all proportions estimated over
the conditioning set $`\{p_1 \le \kappa \mathrm{ or } p_2 \le \kappa\}`$
($`\kappa`$ = `p_cut`), define for thresholds $`(x,y)=(p_1,p_2)`$

``` math
\widehat{\mathrm{cfdr}}(p_1 \mid p_2) =
  \frac{p_1}{\;\widehat{\Pr}(P_1 \le p_1 \mid P_2 \le p_2)\;}, \qquad
  \widehat{\Pr}(P_1 \le p_1 \mid P_2 \le p_2) =
  \frac{\#\{P_1 \le p_1,\, P_2 \le p_2\}}{\#\{P_2 \le p_2\}} .
```

The numerator $`p_1`$ is the (uniform) null tail probability and the
denominator is the empirical conditional cumulative proportion. The
implementation computes a single shared denominator
$`\widehat{\Pr}(P_1 \le x,\, P_2 \le y)/\widehat{\Pr}(P_1 \le x)`$ and
uses it for both directions:

``` math
\widehat{\mathrm{cfdr}}_1 = \frac{p_1}{d}, \quad
      \widehat{\mathrm{cfdr}}_2 = \frac{p_2}{d}, \quad
  d = \frac{\#\{P_1 \le p_1,\, P_2 \le p_2\}}{\#\{P_1 \le p_1\}} .
```

**Conjunctional FDR.** The conjFDR (here `ccfdr`) is the maximum of the
two conditional FDRs,

``` math
\mathrm{conjFDR}(v) = \max\bigl(\widehat{\mathrm{cfdr}}_1,\,
  \widehat{\mathrm{cfdr}}_2\bigr),
```

Small values require evidence in *both* conditional directions. Variants
outside the conditioning set ($`p_1 > \kappa`$ and $`p_2 > \kappa`$)
receive $`\mathrm{conjFDR}=1`$. The current implementation does not
clamp the output to $`[0,1]`$, so sparse conditioning sets can return
values greater than one.

## What This Method Advances

conjFDR reframes pleiotropy as reciprocal conditional enrichment and
reports a conservative conjunctional FDR-style quantity.

## Assumptions

- Under the per-trait null the one-sided p-values are Uniform$`(0,1)`$,
  so the unconditional null tail probability is $`p`$ itself (no
  genomic-inflation rescaling is applied here).
- The empirical conditional CDF is a consistent estimate of the true
  conditional association probability, requiring enough variants in the
  conditioning strata.
- Taking the maximum of the two conditional quantities is used as a
  conservative conjunctional ranking rule; the implementation should be
  interpreted as FDR-style rather than as an ordinary p-value.
- Variants are treated as exchangeable; LD pruning of the input is
  assumed to have removed strong local dependence.

## When This Method Is a Good Fit

Good for enrichment-style discovery of loci shared between two traits;
interpret the output as FDR-like, not as an ordinary p-value.

The `conjfdr` column is a conjunctional conditional-FDR-style quantity.
It is thresholded like an FDR measure and should not be described as an
ordinary null p-value.

The implementation does not clamp the returned values to `[0, 1]`, and
the empirical formula uses one shared denominator based on the
`P_1 <= x` conditioning count. Very sparse conditioning sets can
therefore produce unstable or greater-than-one values.

Pairwise methods use `i` and `j` to select the two phenotypes. Re-run
the method for each phenotype pair that should be tested.

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
pleio <- run_conjfdr(pleio, i = 1, j = 2)
head(pleio@results$conjfdr)
```

## Implementation Files

- R wrapper: `R/run_conjfdr.R`
- C++ kernel or dispatcher: `src/conjfdr.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Andreassen OA, Thompson WK, Schork AJ, et al. (2013). Improved detection
of common variants associated with schizophrenia and bipolar disorder
using pleiotropy-informed conditional false discovery rate. *PLoS
Genetics* 9(4):e1003455. DOI:
<https://doi.org/10.1371/journal.pgen.1003455>

Smeland OB, Frei O, Shadrin AA, et al. (2020). Discovery of shared
genomic loci using the conditional false discovery rate approach. *Human
Genetics* 139:85-94. DOI: <https://doi.org/10.1007/s00439-019-02060-2>

# metaCCA canonical-correlation association test

Function:
[`run_metacca()`](https://broccolito.github.io/pleioverse/reference/run_metacca.md)

Result column: `metacca`

Family: Canonical-correlation association

## Introduction

Summary-statistic canonical correlation analysis (Cichonska et al. 2016)
between a single genotype and the $`K`$ phenotypes, using only GWAS
summary statistics. For each variant the squared canonical correlation
between the SNP and the multivariate phenotype is converted to a Wilks’
lambda chi-square p-value.

## Mathematical Framework and Implementation

Let $`\hat\beta`$ and $`\mathrm{se}`$ be the $`n\times K`$ matrices of
univariate effect sizes and standard errors, $`N`$ the (mean) per-trait
sample size, and $`Z=\hat\beta/\mathrm{se}`$ the matrix of Z-scores. The
phenotypic correlation matrix $`S_{yy}`$ is estimated once as the
across-variant correlation of the effect-size columns,
$`S_{yy}=\mathrm{cor}(\hat\beta)`$ (the metaCCA `estimateSyy` step).

For a variant the standardized SNP-phenotype cross-covariance is the row
of Z-scores rescaled by $`\sqrt{N}`$, $`z = Z_{s\cdot}/\sqrt{N}`$,
playing the role of $`S_{xy}`$. With genotype variance standardized to
one ($`S_{xx}=1`$), the single non-trivial squared canonical correlation
is the quadratic form

``` math
r_c^2 = z'\,S_{yy}^{-1}\,z.
```

The corresponding Wilks’ lambda is $`\Lambda = 1 - r_c^2`$, transformed
to a Bartlett-corrected chi-square statistic

``` math
X = -\Big[(N-1) - \tfrac{1}{2}(K+2)\Big]\,\log\Lambda,
```

whose null distribution is $`\chi^2`$ with $`K`$ degrees of freedom; the
p-value is its upper tail. If $`\Lambda\le 0`$ (i.e. $`r_c^2\ge 1`$) the
p-value is set to 0.

## What This Method Advances

metaCCA frames a single-SNP multi-trait association as a
canonical-correlation problem using summary statistics.

## Assumptions

Under the null of no association $`X`$ is approximately $`\chi^2_K`$ via
the Bartlett-corrected Wilks’ lambda approximation, which assumes
multivariate-normal traits and large $`N`$. $`S_{yy}`$ must be
invertible and is taken as the across-variant effect-size correlation;
the genotype is standardized ($`S_{xx}=1`$) and a common sample size
$`N`$ is used for all traits. Each variant is tested separately against
the full phenotype set.

## When This Method Is a Good Fit

Good for broad multi-trait association from beta/se matrices when
canonical correlation is the desired lens.

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
pleio <- run_metacca(pleio)
head(pleio@results$metacca)
```

## Implementation Files

- R wrapper: `R/run_metacca.R`
- C++ kernel or dispatcher: `src/metacca.cpp`
- Shared MPAT inputs, where applicable: `R/pleiotest.R` and
  `src/test_engine.cpp`

## References

Cichonska A, Rousu J, Marttinen P, Kangas AJ, Soininen P, Lehtimaki T,
Raitakari OT, Jarvelin MR, Salomaa V, Ala-Korpela M, Ripatti S, Pirinen
M (2016). metaCCA: summary statistics-based multivariate meta-analysis
of genome-wide association studies using canonical correlation analysis.
*Bioinformatics* 32(13):1981-1989. DOI:
<https://doi.org/10.1093/bioinformatics/btw052>

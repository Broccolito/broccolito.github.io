# Package index

## Package and data model

Core package documentation and the S4 container used by every workflow.

- [`pleioverse`](https://broccolito.github.io/pleioverse/reference/pleioverse-package.md)
  [`pleioverse-package`](https://broccolito.github.io/pleioverse/reference/pleioverse-package.md)
  : pleioverse: fast self-contained pleiotropy simulation and testing
- [`pleio-class`](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  : The pleio simulation class

## Simulation workflow

Simulate pleiotropic data and build the matrices that control
simulation.

- [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md)
  [`pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md)
  : Simulate pleiotropic genetic data
- [`create_heritable_correlation_matrix()`](https://broccolito.github.io/pleioverse/reference/create_heritable_correlation_matrix.md)
  : Build a uniform heritable (genetic) correlation matrix
- [`create_nonheritable_correlation_matrix()`](https://broccolito.github.io/pleioverse/reference/create_nonheritable_correlation_matrix.md)
  : Build a uniform non-heritable (environmental) input matrix
- [`create_cohort_makeup_matrix()`](https://broccolito.github.io/pleioverse/reference/create_cohort_makeup_matrix.md)
  : Build the cohort-makeup (sample-size) matrix
- [`create_eaf_matrix()`](https://broccolito.github.io/pleioverse/reference/create_eaf_matrix.md)
  : Build the effect-allele-frequency (EAF) matrix
- [`create_efs_matrix_template()`](https://broccolito.github.io/pleioverse/reference/create_efs_matrix_template.md)
  : Build the effect-size (EFS) template matrix

## Test orchestration

Run the full multi-test workflow from a simulated pleio object.

- [`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
  [`pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)
  : Run the full test suite on a pleio object

## Multiple-association and omnibus tests

Tests for association with any phenotype, including MPAT/MTAR-style
methods.

- [`run_wald()`](https://broccolito.github.io/pleioverse/reference/run_wald.md)
  : Wald multi-trait test (inverse-variance quadratic form)
- [`run_wi()`](https://broccolito.github.io/pleioverse/reference/run_wi.md)
  : Weighted-identity (WI) multi-trait test
- [`run_vc()`](https://broccolito.github.io/pleioverse/reference/run_vc.md)
  : Variance-component (VC / SKAT-style) multi-trait test
- [`run_sum()`](https://broccolito.github.io/pleioverse/reference/run_sum.md)
  : Weighted burden (SUM) multi-trait test
- [`run_dsum()`](https://broccolito.github.io/pleioverse/reference/run_dsum.md)
  : Direct burden (DSUM) multi-trait test
- [`run_pcfisher()`](https://broccolito.github.io/pleioverse/reference/run_pcfisher.md)
  : Principal-component Fisher-combination test (PCFisher)
- [`run_pclc()`](https://broccolito.github.io/pleioverse/reference/run_pclc.md)
  : Principal-component linear-combination test (PCLC)
- [`run_pcminp()`](https://broccolito.github.io/pleioverse/reference/run_pcminp.md)
  : Principal-component minimum-p test (PCMinP)
- [`run_minp()`](https://broccolito.github.io/pleioverse/reference/run_minp.md)
  : Minimum-p test with multivariate-normal correction (MinP)
- [`run_mixfisher_davies()`](https://broccolito.github.io/pleioverse/reference/run_mixfisher_davies.md)
  : Mixture-Fisher omnibus (Davies tail)
- [`run_mixfisher_liu()`](https://broccolito.github.io/pleioverse/reference/run_mixfisher_liu.md)
  : Mixture-Fisher omnibus (Liu tail)
- [`run_mixfisher_liumod()`](https://broccolito.github.io/pleioverse/reference/run_mixfisher_liumod.md)
  : Mixture-Fisher omnibus (modified-Liu tail)
- [`run_mixtippett_davies()`](https://broccolito.github.io/pleioverse/reference/run_mixtippett_davies.md)
  : Mixture-Tippett (minimum-p) omnibus (Davies tail)
- [`run_mixtippett_liu()`](https://broccolito.github.io/pleioverse/reference/run_mixtippett_liu.md)
  : Mixture-Tippett (minimum-p) omnibus (Liu tail)
- [`run_mixtippett_liumod()`](https://broccolito.github.io/pleioverse/reference/run_mixtippett_liumod.md)
  : Mixture-Tippett (minimum-p) omnibus (modified-Liu tail)
- [`run_mixvar_davies()`](https://broccolito.github.io/pleioverse/reference/run_mixvar_davies.md)
  : Mixture variance-component omnibus (Davies tail)
- [`run_mixvar_liu()`](https://broccolito.github.io/pleioverse/reference/run_mixvar_liu.md)
  : Mixture variance-component omnibus (Liu tail)
- [`run_mixvar_liumod()`](https://broccolito.github.io/pleioverse/reference/run_mixvar_liumod.md)
  : Mixture variance-component omnibus (modified-Liu tail)
- [`run_mixsd_davies()`](https://broccolito.github.io/pleioverse/reference/run_mixsd_davies.md)
  : Mixture standard-deviation omnibus (Davies tail)
- [`run_mixsd_liu()`](https://broccolito.github.io/pleioverse/reference/run_mixsd_liu.md)
  : Mixture standard-deviation omnibus (Liu tail)
- [`run_mixsd_liumod()`](https://broccolito.github.io/pleioverse/reference/run_mixsd_liumod.md)
  : Mixture standard-deviation omnibus (modified-Liu tail)
- [`run_mixada()`](https://broccolito.github.io/pleioverse/reference/run_mixada.md)
  : Adaptive mixture omnibus (mixAda)
- [`run_cmats()`](https://broccolito.github.io/pleioverse/reference/run_cmats.md)
  : CMATS burden-augmented adaptive omnibus association test
- [`run_emats()`](https://broccolito.github.io/pleioverse/reference/run_emats.md)
  : EMATS leading-eigenvector adaptive omnibus association test
- [`run_tates()`](https://broccolito.github.io/pleioverse/reference/run_tates.md)
  : TATES (extended Simes) multi-phenotype test
- [`run_metacca()`](https://broccolito.github.io/pleioverse/reference/run_metacca.md)
  : metaCCA canonical-correlation association test

## Pleiotropy-specific tests

Pairwise and joint tests for variants associated with multiple
phenotypes.

- [`run_placo()`](https://broccolito.github.io/pleioverse/reference/run_placo.md)
  : PLACO+ pleiotropy test on a phenotype pair
- [`run_pgctest()`](https://broccolito.github.io/pleioverse/reference/run_pgctest.md)
  : PGCtest (CompTestER) normal-product test on a phenotype pair
- [`run_maiup()`](https://broccolito.github.io/pleioverse/reference/run_maiup.md)
  : MAIUP empirical-FDR pleiotropy screen on a phenotype pair
- [`run_conjfdr()`](https://broccolito.github.io/pleioverse/reference/run_conjfdr.md)
  : Conjunctional conditional FDR on a phenotype pair
- [`run_polarmorphism()`](https://broccolito.github.io/pleioverse/reference/run_polarmorphism.md)
  : PolarMorphism pleiotropy test across all phenotypes

## Result aggregation

Combine already-computed p-value columns.

- [`run_acat()`](https://broccolito.github.io/pleioverse/reference/run_acat.md)
  : ACAT aggregated-Cauchy omnibus over existing results

# Adaptive mixture omnibus (mixAda)

MPAT/MTAR multi-trait association test that adaptively searches a grid
of burden/mean vs. heterogeneity/variance mixing weights, takes the most
significant weighted quadratic form, and calibrates the resulting
minimum p-value against its analytic null by one-dimensional numerical
integration.

## Usage

``` r
run_mixada(pleio)
```

## Arguments

- pleio:

  a
  [pleio](https://broccolito.github.io/pleioverse/reference/pleio-class.md)
  object produced by
  [`run_pleiosim()`](https://broccolito.github.io/pleioverse/reference/run_pleiosim.md).

## Value

the `pleio` object with `results$mixada` populated (one p-value per
variant).

## Details

For a variant the input is the vector of trait Z-scores \\z =
(z_1,\dots,z_K)'\\ with \\z_k = \hat\beta_k/\widehat{\mathrm{se}}\_k\\,
which under the null is \\z \sim \mathcal{N}(0,\Sigma)\\ for the
\\K\times K\\ trait correlation \\\Sigma\\. Write \\\Sigma^{-1}\\ for
its inverse, \\\mathbf{1}\\ for the all-ones vector, \\a =
\Sigma^{-1}\mathbf{1}\\, \\s = \mathbf{1}'\Sigma^{-1}\mathbf{1}\\, and
\\H = \mathbf{1}\mathbf{1}'/K\\. As in the other mixture tests, the
burden/mean direction is carried by the rank-one kernel \\a a'\\ and the
heterogeneity/variance direction by \$\$M\_\tau =
(I-H)\\\Sigma^{-1}\Sigma^{-1}\\(I-H).\$\$

**Adaptive mixing grid.** Instead of fixing a single mixing weight, the
test sweeps a grid \\\phi_1,\dots,\phi_B\\ (\\B = 11\\ equally spaced
points on \\\[0.01, 0.99\]\\). For each grid point it forms the blended
kernel and quadratic form \$\$L\_{\phi_b} = \phi_b\\ a a' + (1-\phi_b)\\
M\_\tau, \qquad T_b = z' L\_{\phi_b} z .\$\$ Under the null each \\T_b\\
is a weighted sum of \\\chi^2_1\\ with weights \\\lambda^{(b)}\_j\\ the
eigenvalues of \\\Sigma^{1/2} L\_{\phi_b} \Sigma^{1/2}\\; its p-value
\\p_b\\ is obtained from the modified-Liu four-cumulant chi-square
approximation. The adaptive statistic is the minimum over the grid,
\$\$P\_{\min} = \min\_{1\le b\le B} p_b .\$\$

**Analytic null calibration.** Because the \\T_b\\ are strongly
dependent (they share \\z\\), \\P\_{\min}\\ is not uniform and must be
recalibrated. Conditioning on the shared burden score and writing each
\\T_b\\ as a \\\phi_b\\-weighted combination of a common burden term and
a common heterogeneity term (whose null follows a moment-matched
chi-square with degrees of freedom and scale derived from the
eigenvalues of \\\Sigma^{1/2} M\_\tau \Sigma^{1/2}\\), the p-value of
\\P\_{\min}\\ reduces to a one-dimensional numerical integral over the
burden coordinate. With the substitution \\t = \Pr(\chi^2_1 \le x)\\
(which removes the \\\chi^2_1\\ density singularity at the origin), this
integral is evaluated by composite Simpson quadrature, giving \$\$p =
1 - \int_0^{t\_{\max}} \Pr\\\big(\chi^2\_{\nu} \le \delta^\*(t)\big)\\
dt,\$\$ where \\\delta^\*(t)\\ is the smallest moment-matched
heterogeneity threshold across the grid that is consistent with
achieving \\P\_{\min}\\ at burden coordinate \\t\\, and \\\nu\\ is the
matched heterogeneity degrees of freedom. The result is clamped to
\\\[0,1\]\\.

## Assumptions

The adaptive search makes mixAda the most regime-agnostic member of the
mixture family: by minimizing over the mixing grid it tracks whichever
blend of concordant (burden) and discordant/antagonistic (heterogeneity)
signal best fits the variant, and the analytic null integration restores
valid calibration without permutation. This robustness costs a small
amount of power relative to a single correctly-specified mixture when
the architecture is known. Validity requires \\\Sigma\\ to be a correct,
full-rank trait correlation matrix, the per-trait Z-scores to be
approximately Gaussian under the null, and the modified-Liu and Simpson
approximations to hold in the tail.

## References

Liu, Z., Lin, X. (2018) Multiple Phenotype Association Tests Using
Summary Statistics in Genome-Wide Association Studies. *Biometrics*
74(1):165-175.
[doi:10.1111/biom.12735](https://doi.org/10.1111/biom.12735)

Luo, L., Shen, J., Zhang, H., et al. (2020) Multi-trait analysis of
rare-variant association summary statistics using MTAR. *Nature
Communications* 11:2850.
[doi:10.1038/s41467-020-16591-0](https://doi.org/10.1038/s41467-020-16591-0)

Lee, S., Wu, M. C., Lin, X. (2012) Optimal tests for rare variant
effects in sequencing association studies. *Biostatistics* 13, 762-775.
[doi:10.1093/biostatistics/kxs014](https://doi.org/10.1093/biostatistics/kxs014)

Liu, H., Tang, Y., Zhang, H. H. (2009) A new chi-square approximation to
the distribution of non-negative definite quadratic forms in non-central
normal variables. *Computational Statistics & Data Analysis* 53,
853-856.
[doi:10.1016/j.csda.2008.11.025](https://doi.org/10.1016/j.csda.2008.11.025)

## See also

[`run_pleiotest()`](https://broccolito.github.io/pleioverse/reference/run_pleiotest.md)

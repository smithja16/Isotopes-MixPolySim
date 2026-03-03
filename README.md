# Isotopes-MixPolySim
This script does a Monte Carlo geometry test for stable isotope mixing models. It implements the convex‐hull **mixing polygon simulation** described in:

Smith et al. (2013) *To fit or not to fit: evaluating stable isotope mixing models using simulated mixing polygons*. Methods in Ecology and Evolution.

The user runs the **MixPolySim_2iso_1.3.r** code, which will load the function from MixPolySim_2iso_1.3_Function.r.

## Purpose

This script performs an **a priori feasibility check** for 2-isotope mixing models. A 3-isotope version exists and the R code that is supplementary with 2013 article is still valid for that scenario.

Each iteration:

* Samples source means and SDs (and TEFs) from normal distributions
* Constructs a convex hull ("mixing polygon") of enriched sources
* Tests whether each consumer lies inside the polygon
* Calculates the proportion of iterations each consumer is inside

The resulting value is interpreted as a **frequentist probability that a consumer lies within feasible mixing space** given the specified uncertainty.

> This is a geometry test — not a diet estimation model.

---

## When to Use This

Use this script **before fitting a Bayesian mixing model** (e.g., MixSIAR or simmr) to:

* Identify consumers likely outside mixing space
* Diagnose unrealistic TEFs
* Detect missing or poorly defined sources
* Assess whether uncertainty makes the mixing region overly diffuse

> Passing this test **does not** guarantee the model is correct — it only indicates geometric feasibility.

---

## Input Format

The script assumes the same format as the original demo data:

### `Sources_example.csv`

| mean_iso1 | sd_iso1 | mean_iso2 | sd_iso2 |

### `TEF_example.csv`

| mean_iso1 | sd_iso1 | mean_iso2 | sd_iso2 |

### `Mixture_example.csv`

| iso1 | iso2 |

* TEFs are applied to sources.
* Source-specific TEFs are supported.
* Column order must match the demo format.

---

## Outputs

The function returns:

```
results$probabilities   # P(inside mixing polygon) per consumer
results$area_plot       # Iteration vs polygon-area variance plot
results$mix_plot        # Mixing region raster + contours
results$biplot_95       # 95% contour bi-plot with source SD bars
```

It also writes:

* `Consumer_Probabilities.csv`
* `Parameter_Values.csv`
* PDF figures

---

## Interpreting Results

### Consumer probabilities

A common heuristic threshold:

* **< 0.05** → consumer likely outside mixing space

This threshold is not universal — justify it in your study context.

---

### Mixing region contours

* The **0.05 contour** represents the “95% mixing region.”
* Diffuse contours may indicate:

  * Large source SDs
  * Overlapping source means
  * High TEF uncertainty

---

### Area variance plot

Used only as a convergence diagnostic for number of iterations.

Increase `its` until the variance stabilises.

---

## Key Assumptions

* Sources and TEFs sampled from **independent normal distributions**
* Isotopes treated as independent (no covariance modelled)
* No concentration dependence
* Geometry defined by convex hull

If isotopes are strongly correlated, univariate sampling may overestimate mixing space.

---

## Recommended Workflow

1. Run this geometry test
2. Refine sources/TEFs if needed
3. Fit a Bayesian mixing model (MixSIAR, simmr)
4. Report both geometry screening and posterior diagnostics

Recommended reading:

* Phillips et al. 2014 – Best practices for stable isotope mixing models
* MixSIAR paper: Stock et al. 2018. Analyzing mixing systems using a new generation of Bayesian tracer mixing models. PeerJ 6:e5096
* simmr documentation: https://cran.r-project.org/web/packages/simmr/vignettes/simmr.html

---

## Reproducibility

The function includes a `seed` argument:

```r
results <- mix_poly_sim(sources, mixture, TEF, its = 1500, seed = 42)
```

Always report:

* Number of iterations
* Grid resolution
* TEF assumptions

---

## Polygon Area Calculation

Polygon area is computed using the **shoelace formula**.

For convex hull polygons (non-self-intersecting), this yields the same planar area as `splancs::areapl()`, up to numerical precision.

---

## Citation

If you use this code, please cite:

Smith, J. A., Mazumder, D., Suthers, I. M., & Taylor, M. D. (2013). To fit or not to fit: evaluating stable isotope mixing models using simulated mixing polygons. Methods in Ecology and Evolution, 4(7), 612-618.

---


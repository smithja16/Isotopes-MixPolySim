# Isotopes-MixPolySim
This script does a Monte Carlo geometry test for stable isotope mixing models. It implements the convex‐hull mixing polygon simulation described in:

Smith et al. (2013) To fit or not to fit: evaluating stable isotope mixing models using simulated mixing polygons. Methods in Ecology and Evolution.

Purpose

This tool performs an a priori feasibility check for 2-isotope mixing models.

It repeatedly:

Samples source means and SDs (and TEFs) from normal distributions

Constructs a convex hull (“mixing polygon”) of enriched sources

Tests whether each consumer falls inside the polygon

Estimates the proportion of iterations each consumer is inside

The output probability is interpreted as a frequentist probability that the consumer lies within feasible mixing space given the specified uncertainty.

This is a geometry test, not a diet estimation method.

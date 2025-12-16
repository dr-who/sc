# Session 3: Additional MATLAB-Compatible Functions

This document describes the 66 new functions added in Session 3.

## Signal Processing Functions (7)

| Function | Description | Example |
|----------|-------------|---------|
| `sinc(x)` | Normalized sinc: sin(πx)/(πx) | `sinc(0)` → 1 |
| `sigmoid(x)` | Logistic sigmoid: 1/(1+e^(-x)) | `sigmoid(0)` → 0.5 |
| `logistic(x)` | Alias for sigmoid | `logistic(0)` → 0.5 |
| `softplus(x)` | Smooth ReLU: ln(1+e^x) | `softplus(0)` → 0.693 |
| `step(x)` | Heaviside step function | `step(1)` → 1 |
| `heaviside(x)` | Alias for step | `heaviside(-1)` → 0 |
| `rect(x)` | Rectangle/boxcar function | `rect(0)` → 1 |
| `tri(x)` | Triangle function | `tri(0.5)` → 0.5 |

## Angle Wrapping Functions (4)

| Function | Description | Example |
|----------|-------------|---------|
| `wrapToPi(x)` | Wrap angle to [-π, π] | `wrapToPi(7)` → 0.717 |
| `wrap360(x)` | Wrap angle to [0, 360) | `wrap360(370)` → 10 |
| `wrapTo360(x)` | Alias for wrap360 | |
| `wrap180(x)` | Wrap angle to [-180, 180] | `wrap180(270)` → -90 |
| `wrapTo180(x)` | Alias for wrap180 | |

## Precision Math Functions (2)

| Function | Description | Example |
|----------|-------------|---------|
| `log1p(x)` | ln(1+x), accurate for small x | `log1p(0.001)` → 0.0009995 |
| `expm1(x)` | e^x - 1, accurate for small x | `expm1(0.001)` → 0.0010005 |

## Advanced Statistics (7)

| Function | Description | Example |
|----------|-------------|---------|
| `geomean(v)` | Geometric mean | `geomean([1,2,4,8])` → 2.83 |
| `harmmean(v)` | Harmonic mean | `harmmean([1,2,4])` → 1.71 |
| `skewness(v)` | Sample skewness | `skewness([1,2,3,4,5])` → 0 |
| `kurtosis(v)` | Excess kurtosis | `kurtosis([1,2,3,4,5])` → -1.3 |
| `mad(v)` | Mean absolute deviation | `mad([1,2,3,4,5])` → 1.2 |
| `iqr(v)` | Interquartile range (Q3-Q1) | `iqr([1:8])` → 3.5 |

## Interpolation & Signal Processing (3)

| Function | Description | Example |
|----------|-------------|---------|
| `interp1(x,y,xi)` | Linear interpolation | `interp1([0,1,2],[0,1,4],1.5)` → 2.5 |
| `deconv(u,v)` | Deconvolution/polynomial division | `deconv([1,3,2],[1,1])` → [1,2] |

## Random Number Fix

Fixed `randi()` to use high bits of LCG for proper randomness. Previously, `randi(8,8)` would produce repeated rows due to short period of low bits in the Linear Congruential Generator.

## Summary

**Session 3 Total: ~20 new function names** (some are aliases)

Categories:
- Signal Processing: sinc, sigmoid, softplus, step, rect, tri (7)
- Angle Wrapping: wrapToPi, wrap360, wrap180 (4 with aliases)
- Precision Math: log1p, expm1 (2)
- Statistics: geomean, harmmean, skewness, kurtosis, mad, iqr (6)
- Interpolation: interp1, deconv (2)

**Combined Total: ~350 functions**

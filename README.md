<div align="center">

# trtswitch-parallel


[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![Version](https://img.shields.io/badge/version-0.1.7-orange)]()
[![Platform](https://img.shields.io/badge/platform-linux-lightgrey)]()
</div>


--- 
## NOTE: 
Link to the official version: https://cran.r-project.org/web/packages/trtswitch/index.html <br>
This was part of a software engineering internship where I was helping the author of the original trtswitch package library refactor some of the functions. This repository was part of a sprint to complete tsesimp.cpp's refacotring into pure C++ to allow use 
for multi threading. As of version 0.2.3 of trtswitch it has adopted some of my suggestions and improved upon them so this is an outdated build, but this is meant to serve documentation for the how and why regarding the projects refactorization. 

## Main Overview
In survival analysis patients in a control arm can switch treatements at progression or at some crossover time. This can skew the results of the HR value (for non statistics people, this basically is very bad), it implements an accelerated failure time (AFT) model 
to estimate a scaling parameter (psi) adjusting follow up times for switches, then a Cox model is fitted on the adjusted times to estimate a hazard ratio (HR). To quantify the confidence of it, we construct a confidence interval using bootstrapping. 

The original R implementation relied on Rcpp and DataFrame heavily, however despite being convenient Rcpp objects are single threaded. The main goal of the project was to refactor comptuations into pure C++ so that bootstrap replicates can be multi threaded using
pragma initially for the MVP, then ported to RcppParallel.

## What changed
(File structure is from highest level to lowest level)
| Original File | Modified Files |
|--------:|----------------:|
| tsesimp.R               | tsesimp_mt.R                    |   
| tsesimp.cpp/h           | tsesimp_mt.cpp/h                |
| survival_analysis.cpp/h | survival_analysis_pure.cpp/h    |
| utilities.cpp/h         | utilities_pure.cpp/h            |

- Removal of R in hot loop: Every single R object call was removed to enable only pure C++ calls to let the bootstra
-	Row‑major matrices: A new MatrixRM type represents a dense numeric matrix in row‑major order. This avoids the overhead of Rcpp’s NumericMatrix and allows contiguous memory access when iterating over rows. Matrices are passed as pointers to avoid unnecessary copies. <br>
-	Non‑owning views: RowIndexView and RowIndexViewView provide lightweight, non‑owning access to selected rows of a matrix. Row reordering is handled by mutating a std::vector<int> of indices rather than physically copying the underlying data. 
This makes subsetting inexpensive and thread‑safe. <br>
-	Simple structs: The inputs to the AFT and Cox model fits are stored in plain C++ structs (trial_data, aftparams_pure, coxparams_purecpp), avoiding the need to construct or destruct Rcpp objects inside threads. <br>

## Parallel bootstrap
- run_replicate(): A new function encapsulating all per replicate operations. It orders the data by treatement, samples with replacement within treatement groups using precomputed bootstrap indices, calling the COX and AFT fits. Returns only the values needed in
the bootstrap to construct the confidence interval. Function entirely removed all R API calls, making it safe to invoke from multiple threads.

## Removed or deferred features
For the minimum viable product, some low priority features from the lower level helper functions were disabled using compiler flags:
-	Baseline hazard estimation and associated return objects.
-	Firth correction for the Cox model.
-	Residuals and rep‑grouping logic.
-	Profile likelihood confidence intervals.
These can be reintroduced once the core multi‑threaded pipeline is stable (and if needed).

## Why make the changes?
The original trtswitch implementation made heavy use of Rcpp and R data structures. While convenient, those structures impose two major constraints:
1.	Thread safety – The R interpreter is not thread‑safe; any calls into R (printing to the console, generating random numbers, or constructing R objects) must occur on the main thread.
2.	Running the bootstrap in parallel would cause sporadic errors, stack imbalance warnings, or crashes if Rcpp functions were invoked in worker threads. By eliminating Rcpp dependencies from the hot loop we can safely useOpenMP to run replicates concurrently.
3.	Performance – Repeatedly subsetting and copying DataFrame objects inside the bootstrap loop introduced significant overhead. Packing covariates into contiguous row‑major matrices and using a single row‑index vector avoids these copies and improves cache locality.
The refactor also clarifies the data flow: shared inputs such as covariate names, distribution choice, and alpha are separated from replicate‑specific scratch buffers. This design makes it easier to reason about memory ownership and
to extend the pipeline (e.g. adding more covariates) without introducing hidden dependencies.

## Impact
The pure C++ refactor yields the following benefits:
- Parallel speedup: On a dataset with 10,000 observations and 500 bootstrap replicates, the multithreaded version reduces wall-clock runtime by roughly 4x  compared to the original serial Rcpp implementation. Performance gains from both multi threading
and more efficent memory accessing w the views and flat matrix.
- Numerical Parity: Unit tests confirm that the point estimates for psi and HR match the CRAN implementation within 1e-6 tolerance, while the bootstrap confidence interval were identitcal in a larger tolerance.
- Better maintainability: With Rcpp dependancies removed from the inner loop, the code seperates clearly into C++ numerical code and a R wrapper. Future optimizations become much easier with this!


## Benchmarking
### E2E runtime with baseline comparison
(Based off bench-tsesimpmulti.R)
n_boot = 500
| Threads | Median time (s) | Speedup vs 1T | Speedup vs baseline (CRAN) |
|--------:|----------------:|--------------:|----------------------------:|
| baseline (CRAN) | 10.7s     | 1.35x (base)       | 1.00x |
| 1       | 14.5s             | 1.00(base)         | 0.74x |
| 2       | 8.25s             | 1.81x              | 1.26x |
| 4       | 4.36s             | 3.43x              | 2.42x |
| 8       | 2.86s             | 5.25x              | 3.52x |
| 12      | 2.78s             | 5.33x              | 3.77x |

The multithreaded implementation performs additional one-time setup (e.g., packing covariates into a row-major matrix and initializing per-replicate scratch buffers). <br>
For small n_boot, this fixed overhead can reduce apparent speedup or even make runs slightly slower. <br>
For larger n_boot, runtime is dominated by the replicate loop, and parallel speedups become more representative. <br>

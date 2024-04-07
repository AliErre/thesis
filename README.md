# wfda: Reconstruction Package

## Description

The `wfda` package provides tools for the reconstruction of partially observed curves. The reconstruction methods are implemented in C++ for efficiency, with wrapper functions provided for easy usage in R.

The package includes a `ReconstructionBase` class, from which three other classes are derived: `Kraus`, `KLAl`, `KLNoAl`, and `Extrapolation`. The user selects the derived class by providing a string identifier. The possible values for this identifier are:

- "Kraus": if `alpha` is `NULL` then it is set through gcv using the complete curves observations.
- "KLAl": if `K` is `NULL` then it is set through gcv using the complete curves observations.
- "KLNoAl": if `K` is `NULL` then it is set through gcv using the complete curves observations.
- "Extrapolation"

The vector t_points must be provided to every method except "Kraus".



## Installation

You can install the package directly from GitHub using the `devtools` package:

```R
# Install devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the wfda package from GitHub
devtools::install_github("username/wfda")
```
## Usage
```R
library(Rcpp)
library(wfda)
load("curves.Rdata")
load("t_points.Rdata")

# expose ReconstructionBase class
reconstruction <- Module("reconstruction", PACKAGE = "wfda")
ReconstructionBase <- reconstruction$ReconstructionBase

# instantiate objects with string identifier
klal <- new(ReconstructionBase, "KLAl", matrix_data)
klnoal <- new(ReconstructionBase, "KLNoAl", matrix_data)
kraus <- new(ReconstructionBase, "Kraus", matrix_data)
extrapolationo <- new(ReconstructionBase, "Extrapolation", matrix_data)

# call reconstruct method. Arguments: alpha, all, t_points, K, maxBins, nRegGrid
klal <- klal$reconstruct(NULL, FALSE, t.points, NULL, NULL, NULL)
klnoal <- klnoal$reconstruct(NULL, FALSE, t.points, NULL, NULL, NULL)
kraus <- kraus$reconstruct(NULL, FALSE, NULL, NULL, NULL, NULL)
extrapolationo <- extrapolationo$reconstruct(NULL, FALSE, t.points, NULL, NULL, NULL)

# reconstruct_wrapper example. Missing arguments default to NULL.
kl <- reconstruct_wrapper("KLAl", Y = curves, t.points = t.points)
```

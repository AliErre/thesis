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
### In case the repo is public
This is valid only for when the repo will be made public.

You can install the package directly from GitHub using the `devtools` package:

```R
# Install devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install the wfda package from GitHub
devtools::install_github("AliErre/thesis/wfda")
```

### In case you've been given access as a contributor: set up Rcpp
Since the `wfda` package uses Rcpp extensively, you need to ensure you have Rcpp installed and a compatible compiler set up. Here's how to prepare:

- **Install Rtools**:
  Rtools provides the necessary tools, including the C/C++ compiler, to build R packages with C++ code. You can download and install Rtools from the official website: [Rtools Download Page](https://cran.r-project.org/bin/windows/Rtools/).

- **Set Up PATH Environment Variable**:
  After installing Rtools, you need to add its bin directory to your PATH environment variable. This allows your system to locate the Rtools compiler when building packages. Here's how to do it on Windows:
  - Search for "Environment Variables" in your computer's search bar and open the "Edit the system environment variables" option.
  - Click on the "Environment Variables" button.
  - In the "System variables" section, find the "Path" variable and click "Edit".
  - Add the path to the Rtools bin directory (e.g., `C:\Rtools\bin`) to the list of paths. Make sure to separate it from other paths with a semicolon (`;`).

### Building and Installing the Package

3. Once you've cloned the repository and set up Rcpp, navigate to the `wfda` directory in your terminal.

After you clone the repository and have a compiler (which Rtools provides), you can do the following:
```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
install.packages("Rcpp")

setwd("path_to_your_local_repo")
Rcpp::compileAttributes()
devtools::build()
# devtools::check()
install.packages("../wfda_1.0.tar.gz")
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

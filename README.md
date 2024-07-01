# wfda: Reconstruction Package

## Description
The `wfda` package revolves around the reconstruction and estimation of complex data. Specifically, it can be split in two areas: one for the reconstruction of partially observed curves over a domain (see `ReconstructionBase` class) and the other one for tuning parameters for the estimation of a spatially adaptive estimator in a Function on Function (FoF) regression problem. For this latter task, the evolutionary algorithm for the adaptive smoothing spline estimator (Centofanti et al. 2023) is implemented (see `Genetic` class). 

The reconstruction methods and the genetic algorithm are implemented in C++ for efficiency (a wrapper function is provided for easy usage in R). This work builds upon the R code in [Teresa Bortolotti's git hub repo](https://github.com/tbortolotti/WFDA.git).

The package includes a `ReconstructionBase` class, from which three other classes are derived: `Kraus`, `KLAl/KLNoAl`, and `Extrapolation`. Each class specializes a `reconstructCurve` method. The user selects the derived class by providing a string identifier to a factory function, together with a matrix containing the observed curves. The possible values for the string identifier are:

- "Kraus": if `alpha` is `NULL` then it is set through gcv using the complete curves observations.
- "KLAl": if `K` is `NULL` then it is set through gcv using the complete curves observations.
- "KLNoAl": if `K` is `NULL` then it is set through gcv using the complete curves observations.
- "Extrapolation"

The vector t_points must be provided non `NULL` to every method except "Kraus" which does not require it.
Since variadic templating is not possible within `Rcpp` the `reconstructCurve` method has the same number of arguments for every method; moreover the only way to have default values for arguments is to create an R wrapper function. The latter is located in the R directory and wraps `reconstructCurve` so that the calls look more intuitive and every method call only needs to make explicit only the arguments it requires.

### Install from source
Download wfda_1.0.tar.gz, and all the .RData.
On R, you can do the following (navigate to the working directory where the source is located):
```R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
install.packages("Rcpp")
library(Rcpp)

setwd("path_to_wfda_1.0.tar.gz")
Rcpp::compileAttributes()
devtools::build()
# devtools::check()
install.packages("./wfda_1.0.tar.gz")
```
## Usage
```R
library(wfda)
load("curves.Rdata")
load("t_points.Rdata")

# expose ReconstructionBase class
reconstruction <- Module("reconstruction", PACKAGE = "wfda")
ReconstructionBase <- reconstruction$ReconstructionBase

# instantiate objects with new() and string identifier
klal <- new(ReconstructionBase, "KLAl", matrix_data)
klnoal <- new(ReconstructionBase, "KLNoAl", matrix_data)
kraus <- new(ReconstructionBase, "Kraus", matrix_data)
extrapolationo <- new(ReconstructionBase, "Extrapolation", matrix_data)

# call reconstruct method. Arguments: alpha, all, t_points, K, maxBins, nRegGrid
klal <- klal$reconstruct(NULL, FALSE, t.points, NULL, NULL, NULL)
klnoal <- klnoal$reconstruct(NULL, FALSE, t.points, NULL, NULL, NULL)
kraus <- kraus$reconstruct(NULL, FALSE, NULL, NULL, NULL, NULL)
extrapolationo <- extrapolationo$reconstruct(NULL, FALSE, t.points, NULL, NULL, NULL)

#access results
reconst_kraus <- kraus$Y_reconst
kraus$alpha
weights_kraus <- kraus$W_reconst

reconstruct_klal <- klal$Y_reconst_list
weights_klal <- klal$W_reconst_list

reconstruct_klnoal <- klnoal$Y_reconst_list
weights_klnoal <- klnoal$W_reconst_list

reconst_e <- extrapolationo$Y_reconst


# reconstruct_wrapper example. Missing arguments default to NULL.
kr <- reconstruct_wrapper("Kraus", curves)
kl <- reconstruct_wrapper("KLAl", Y = curves, t.points = t.points)
ex <- reconstruct_wrapper("Extrapolation", Y = curves, t.points = t.points)

## load genetic module
genetic_module <- Module("Genetic", PACKAGE = "wfda")
genetic <- genetic_module$Genetic
object_genetic <- new(genetic, xlist.full , blist.default, curves.full, curves.fd, t.points, event.id) #mettere xlist
object_genetic$multistart()

#get results with getters
P.full  <- list()
v.full  <- list()
blist.full  <- list()
for(k in 1:5){
  P.full[[k]]= object_genetic$get_P(k)
  v.full[[k]]= object_genetic$get_v(k)
  blist.full[[k]]=object_genetic$get_blist(k)
}

```

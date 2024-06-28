\name{reconstruction}
\alias{reconstruction}
\title{
    Rcpp module: reconstruction
}
\description{
    Rcpp module provides four methods, identified through a string, to reconstruct partially observed functional data.
}
\details{The module exposes the ReconstructionBase class. Through a factory function which is called with a string identifier and a numeric matrix, an object of a class derived from ReconstructionBase is instantiated.}
\arguments{
\itemize{
\item{\code{id}: string identifier. Possible values are "Kraus", "KLAl", "KLNoAl" and "Extrapolation".}
\item{\code{Y}: matrix containing the curves. Each column is a curve. NA values indicate the curve has not been observed.}}
}
\source{
    https://github.com/AliErre/thesis.git
}
\author{Alice Ragonesi}
\references{
    %% Possibly secondary sources and usages
    %% Perhaps references to the C++ code that the module exposes
}
\examples{
    library(wfda)
    matrix_data <- matrix(c(
        -2.267744, -2.468645, -1.6856891, -2.429462, -0.7875974,
        -2.619638, -2.225056, -1.3639621, -2.129250, -1.1095009,
        -2.260310, -2.000934, -1.3973601, -2.357720, -1.0714788,
        -2.284882, -2.208072, -1.5417680, -2.477875, -0.8794883,
        -2.512043, -2.160219, -1.6673675, -2.414988, -0.8976942,
        -2.749826, -2.284415, -1.7724204, -2.402217, -1.1590291,
        -2.910132, -2.205017, -1.6372608, -2.613263, -1.2320787,
        -2.614979, -1.956221, -1.2267533, -2.296165, NA,
        -2.248436, -1.700328, -1.0281668, -1.953007, NA,
        -1.753496, -1.750807, -1.0005058, NA, NA,
        -1.580548, -1.459076, -0.3417352, NA, NA,
        -1.736311, -1.656483, -0.5288931, NA, NA,
        -1.483448, -1.604342, -0.5389337, NA, NA,
        -1.693828, -2.046650, -0.7654954, NA, NA,
        -2.267157, -2.127552, -0.7805399, NA, NA
    ), nrow = 15, byrow = TRUE)
    t.points <- c(0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 
		2.25, 2.50, 2.75, 3.00, 3.25, 3.50)
    reconstruction <- Rcpp::Module("reconstruction", PACKAGE = "wfda")
    ReconstructionBase <- reconstruction$ReconstructionBase
    object <- new(ReconstructionBase, "KLAl", matrix_data)  
}
\keyword{datasets}

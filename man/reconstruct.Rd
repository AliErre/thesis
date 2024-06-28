\name{reconstruct}
\alias{reconstruct}
\title{Reconstruction function}
\description{
    reconstruct() method for the ReconstructionBase class in wfda package.
}
\details{
    The method is specialized based on the derived class the ReconstructionBase object belongs to. Arguments needed for the reconstruction function change with the type of method.
}
\arguments{
    \itemize{
        \item{\code{alpha}: used in "Kraus". If \code{NULL} then set through gcv.}
        \item{\code{all}: indicates whether all data needs to be reconstructed. If \code{FALSE} then only incomplete ones are reconstructed.}
        \item{\code{t_points}: used in "KLAl", "KLNoAl" and "Extrapolation". Must be provided for these methods.}
        \item{\code{K}: integer. Used in "KLAl" and "KLNoAl". If \code{NULL} then set through gcv.}
        \item{\code{maxBins}: integer. Used in "KLAl" and "KLNoAl". If \code{NULL} then set to 1000}
        \item{\code{nRegGrid}: integer. Used in "KLAl" and "KLNoAl". Can stay NULL.}
    }
}
\value{
    All methods return all the reconstructed curves.
    If "Kraus", "KLAl" or "KLNoAl" then it returns also the weights computed in the reconstruction.
}
\examples{
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
    t.points <- c(0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 
		2.50, 2.75, 3.00, 3.25, 3.50)
    reconstruction <- Rcpp::Module("reconstruction", PACKAGE = "wfda")
    ReconstructionBase <- reconstruction$ReconstructionBase
    object <- new(ReconstructionBase, "KLAl", matrix_data)
    reconstruct <- object$reconstruct( alpha = NULL, all = FALSE, t_points = t.points, K = NULL, 
                                       maxBins = NULL, nRegGrid = NULL) 
}
\source{
    https://github.com/AliErre/thesis.git
}

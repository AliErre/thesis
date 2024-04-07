\name{reconstruct_wrapper}
\alias{reconstruct_wrapper}
\title{Wrapper function for reconstruct method exposed with Rcpp Modules}
\description{This function allows for default values of arguments. This is useful since not all reconstruction methods use the same arguments and it may be more readable to omit the unnecessary ones.}
\details{The named arguments are: \code{id}, \code{Y}, \code{alpha}, \code{t.points}, \code{all}, \code{K}, \code{maxBins}, \code{nRegGrid}. For further information on their use in the different methods run \code{help(reconstruct)}. Last six arguments can be omitted and have a default.}
\arguments{
\itemize{
\item{\code{id} (mandatory): string identifier for the reconstruction method to be employed. Possible values "Kraus", "KLAl", "KLNoAl", "Extrapolation".}
\item{\code{Y} (mandatory): matrix data for the curves. Each column is a curve, complete or partially observed.}
\item{\code{alpha} (optional): defaults to \code{NULL}.}
\item{\code{all} (optional): defaults to \code{FALSE}.}
\item{\code{t.points} (optional): defaults to \code{NULL}.}
\item{\code{K} (optional): defaults to \code{NULL}.}
\item{\code{maxBins} (optional): defaults to \code{NULL}.}
\item{\code{nRegGrid} (optional): defaults to \code{NULL}.}
}
}
\author{
Alice Ragonesi}
\source{
    https://github.com/AliErre/thesis.git
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
reconstruct_wrapper("Kraus", matrix_data)
}


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
reconstruct_wrapper("Kraus", matrix_data)
}


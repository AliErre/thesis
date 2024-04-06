\name{reconstruct}
\alias{reconstruct}
\title{Reconstruction function}
\description{
	reconstruct() method for the ReconstructionBase class in wfda package.
}
\details{
	The method is specialized based on the derived class the ReconstructionBase object belongs to. Arguments needed for the reconstruction function change with the type of method.
}
\usage{
	\code{object$reconstruct(alpha = NULL, all = FALSE, t_points = t.points, K = NULL, maxBins = NULL, nRegGrid =  NULL)}
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
    reconstruct <- object$reconstruct( alpha = NULL, all = FALSE, t_points = t.points, K = NULL, 
                                       maxBins = NULL, nRegGrid = NULL) }
\source{
	https://github.com/AliErre/thesis.git
}
\name{wfda-package}
\alias{wfda-package}
\alias{wfda}
\docType{package}
\title{
    wfda: Reconstruct Partially Observed Curves
}
\description{
    This package provides a class to reconstruct partially observed curves through the class method reconstruct. Four methods are provided: "Kraus", "KLAl", "KLNoAl" and "Extrapolation". These are the string identifiers that need to be used to pick the reconstruction method. Use the string identifier when instantiating an object of the exposed ReconstructionBase class. Specifications of the reconstruct method are invoked based on the string identifier. An alternative to mask all this code is provided by the function \code{reconstruct_wrapper}, a wrapper thanks to which it is possible to have default values for arguments. When not using the wrapper, all arguments must be provided to the \code{reconstruct} method exposed by the Rcpp Modules.
}
\details{
    When not using the wrapper: first load the "reconstruction" module into the R session. Then use new() method to create an object of the exposed class.
}
\examples{
reconstruction <- Module("reconstruction", PACKAGE = "wfda")
}
\author{
    Alice Ragonesi.

    Maintainer: Alice Ragonesi <alice.ragonesi@mail.polimi.it>
}
\references{
    %%This optional section can contain literature or other references for background information.
}
\keyword{ package }
\seealso{
    %%Optional links to other man pages
}

\keyword{datasets}


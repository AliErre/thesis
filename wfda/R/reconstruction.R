#' Wrapper function for reconstruction method
#' @param my_wfda An instance of ReconstructionBase
#' @param alpha Numeric, default is NULL
#' @param all Logical, default is FALSE
#' @param t.points Numeric, default is NULL
#' @param K Numeric, default is NULL
#' @param maxBins Numeric, default is NULL
#' @param nRegGrid Numeric, default is NULL
#' @return Result of reconstruct method
reconstruct_wrapper <- function(
    id,
    Y,
    alpha = NULL,
    all = FALSE,
    t.points = NULL,
    K = NULL,
    maxBins = NULL,
    nRegGrid = NULL
) {
  reconstruction <- Module("reconstruction", PACKAGE = "wfda")
  ReconstructionBase <- reconstruction$ReconstructionBase
  
  # Create instance of ReconstructionBase
  my_wfda <- new(ReconstructionBase, id, Y)
  
  # Call reconstruct method
  my_wfda$reconstruct(alpha, all, t.points, K, maxBins, nRegGrid)
}

# Attach wrapper function to your package
assign("reconstruct_wrapper", reconstruct_wrapper, envir = asNamespace("wfda"))
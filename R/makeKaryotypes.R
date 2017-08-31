#' Make a matrix of karyotypes
#'
#' This function will create a matrix of karyotypes that is compatible with the CINsim pipeline. Currently supported species are mouse (20 chromosomes) and human (23 chromsomes).
#'
#' @param numCell The desired number of cells in the matrix.
#' @param species The desired species, either mouse or human (20 or 23 chromosomes).
#' @param copies The number of copies for all chromosomes in the karyotype matrix.
#' @return A matrix of karyotypes according to the input parameters.
#' @author Bjorn Bakker
#' @export

makeKaryotypes <- function(numCell = 100, species = "mouse", copies = 2) {

  if(species == "mouse") {
    numChr <- 20
  } else if(species == "human") {
    numChr <- 23
  } else {
    numChr <- species
  }

  karyotype <- matrix(rep(copies, times = numCell*numChr), ncol = numChr)
  colnames(karyotype) <- c(1:(numChr - 1), "X")
  row.names(karyotype) <- paste0("cell_", 1:nrow(karyotype))
  return(karyotype)

}

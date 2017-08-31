#' Quantify levels of aneuploidy
#'
#' This function will quantify the degree of aneuploidy in a sample (calculated as the absolute deviation from euploid). Higher scores indicate more aneuploidy.
#'
#' @param karyoMat A karyotype matrix with cells in rows and chromosomes in columns.
#' @param numberOfCells The number of cells in the karyotype matrix.
#' @param euploidRef The reference euploid copy number. Defaults to 2 if none is provided.
#' @return A vector of aneuploidy scores per chromosome.
#' @author Bjorn Bakker

qAneuploidy <- function(karyoMat, numberOfCells = NULL, euploidRef = NULL) {

  # check user input and create alternate variables if necessary
  if(is.null(numberOfCells)) {
    numberOfCells <- nrow(karyoMat)
  }
  if(is.null(euploidRef)) {
    euploidRef <- 2
  }

  # quantify and return level of aneuploidy
  if(numberOfCells > 1) {
    aneuploidyScore <- colMeans(abs(karyoMat - euploidRef))
  } else {
    aneuploidyScore <- mean(abs(karyoMat - euploidRef))
  }
  return(aneuploidyScore)

}

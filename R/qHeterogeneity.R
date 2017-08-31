#' Quantify levels of heterogeneity
#'
#' This function will caculate the degree of heterogeneity in a sample. Higher scores indicate greater heterogeneity.
#'
#' @param karyoMat A karyotype matrix with cells in rows and chromosomes in columns.
#' @return A vector of heterogeneity scores per chromosome.
#' @author Bjorn Bakker

qHeterogeneity <- function(karyoMat) {

  # set number of cells
  numberOfCells <- nrow(karyoMat)

  # calculate absolute frequency of copy number states per chromosome
  if(numberOfCells > 1) {
    tabs <- apply(karyoMat, 2, function(x) {sort(table(x), decreasing = TRUE)})
  } else {
    tabs <- sapply(karyoMat, function(x) {sort(table(x), decreasing = TRUE)})
  }

  # calculate and return heterogeneity score
  if(is.list(tabs)) {
    heterogeneityScore <- unlist(lapply(tabs, function(x) {sum(x * 0:(length(x)-1))}))/numberOfCells
  } else if(is.null(dim(tabs))) {
    heterogeneityScore <- sapply(tabs, function(x) {sum(x * 0:(length(x)-1))})/numberOfCells
  } else {
    heterogeneityScore <- apply(tabs, 2, function(x) {sum(x * 0:(length(x)-1))})/numberOfCells
  }
  return(heterogeneityScore)

}

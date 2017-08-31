#' Quantify levels of heterogeneity (parallel variant)
#'
#' This function will caculate the degree of heterogeneity in a sample. Higher scores indicate greater heterogeneity. This is the same function as qHeterogeneity, but uses parallel computing to speed up calculations.
#'
#' @param karyoMat A karyotype matrix with cells in rows and chromosomes in columns.
#' @param numberOfCells The number of cells in the karyotype matrix.
#' @param cl A defined cluster that will be used to compute heterogeneity in a parallel fashion.
#' @return A vector of heterogeneity scores per chromosome.
#' @import parallel
#' @import doParallel
#' @author Bjorn Bakker

qHeterogeneityPar <- function(karyoMat, numberOfCells = NULL, cl) {

    # verify and check user input
    if(is.null(numberOfCells)) {
        numberOfCells <- nrow(karyoMat)
    }

    # calculate absolute frequency of copy number states per chromosome
    if(numberOfCells > 1) {
        tabs <- parApply(cl, karyoMat, 2, function(x) {sort(table(x), decreasing = TRUE)})
    } else {
        tabs <- parSapply(cl, karyoMat, function(x) {sort(table(x), decreasing = TRUE)})
    }

    if(is.list(tabs)) {
      heterogeneityScore <- unlist(parLapply(cl, tabs, function(x) {sum(x * 0:(length(x)-1))}))/numberOfCells
    } else if(is.null(dim(tabs))) {
      heterogeneityScore <- parSapply(cl, tabs, function(x) {sum(x * 0:(length(x)-1))})/numberOfCells
    } else {
      heterogeneityScore <- parApply(cl, tabs, 2, function(x) {sum(x * 0:(length(x)-1))})/numberOfCells
    }
    return(heterogeneityScore)

}

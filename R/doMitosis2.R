#' Perform a round of mitosis with possible mis-segregations
#'
#' This function will perform mitosis on a given karyotype, apply chromosome mis-segregation events if the criteria are met, and returns a matrix with two new karytoypes.
#'
#' @param karyotype A karyotype provided as a vector on which to apply the mitosis.
#' @param pMisseg The probability of mis-segregation per chromosome copy number.
#' @return A 2-row matrix with karyotypes from daughter cells.
#' @author Aaron Taudt
#' @export


doMitosis2 <- function(karyotypes, pMissegs) {

    # recalculate weighted probabilities
    pMissegsMatrix <- array(pMissegs, dim = dim(karyotypes), dimnames = dimnames(karyotypes))
    pChromWeighted <- 1 - (1 - pMissegsMatrix)^(karyotypes)

    # determine missegregation event per chromosome
    gainLoss <- array(runif(length(pChromWeighted)) < pChromWeighted, dim = dim(karyotypes), dimnames = dimnames(karyotypes))

    # apply missegregations and setup daughter mat
    daughterMat <- rbind(karyotypes + gainLoss, karyotypes - gainLoss)

    # return ouput
    return(daughterMat)

}

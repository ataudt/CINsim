#' Calculate adjusted probability of mis-segregation
#'
#' This function calculates the probability of mis-segregation, such that cells with higher fitness have a lower probability of mis-segregating.
#'
#' @param karyotype A single karyotype, supplied as a vector.
#' @param pMisseg The probability of mis-segregation per chromosome copy.
#' @param pMissegCoef The probability co-efficient, calculated from the euploid reference karyotype.
#' @param numberOfChromosomes The number of chromosomes in the provided karyotype.
#' @return A probability of chromosome mis-segregation.
#' @author Bjorn Bakker

calcMisseg2 <- function(karyotypes, pMisseg, pMissegCoef, numberOfChromosomes=NULL, probDf=NULL) {

    # check user input
    if(is.null(numberOfChromosomes)) {
        numberOfChromosomes <- ncol(karyotypes)
    }
    missegFactor <- pMissegCoef / calcSurProb2(karyotypes, numberOfChromosomes, probDf)
    pMissegAdjusted <- pMisseg * (missegFactor)^2
    return(pMissegAdjusted)

}

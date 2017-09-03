#' Calculate probability of division
#'
#' This function will return a probability of division, based on a probability fitness matrix, such that cells with higher fitness have a greater chance of dividing. Whenever the calculated probability is less than 0.7, a value of 0 is returned instead.
#'
#' @param karyotype A single karyotype as a vector.
#' @param numberOfChromosomes The number of chromosomes in the vector.
#' @param probDf A probability fitness matrix.
#' @return A probability of cell division.
#' @author Bjorn Bakker

calcDivProb2 <- function(karyotypes, numberOfChromosomes = NULL, probDf) {

    # check user input
    if(is.null(numberOfChromosomes)) {
        numberOfChromosomes <- ncol(karyotypes)
    }

    probs <- array(NA, dim=dim(karyotypes), dimnames=dimnames(karyotypes))
    for (i in colnames(karyotypes)) {
        probs[,i] <- probDf[cbind(as.character(karyotypes[,i]), as.character(i))]
    }
    divProb <- rowSums(probs)/sum(apply(probDf, 2, max))
    divProb[divProb < 0.7] <- 0

    return(divProb)
}

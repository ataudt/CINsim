#' Calculate probability of division
#'
#' This function will return a probability of division, based on a probability fitness matrix, such that cells with higher fitness have a greater chance of dividing. Whenever the calculated probability is less than 0.7, a value of 0 is returned instead.
#'
#' @param karyotype A single karyotype as a vector.
#' @param numberOfChromosomes The number of chromosomes in the vector.
#' @param probDf A probability fitness matrix.
#' @return A probability of cell division.
#' @author Bjorn Bakker

calcDivProb <- function(karyotype, numberOfChromosomes = NULL, probDf) {

    # check user input
    if(is.null(numberOfChromosomes)) {
        numberOfChromosomes <- length(karyotype)
    }

    probs <- vector("numeric", numberOfChromosomes)
    for(i in 1:numberOfChromosomes) {
        probs[i] <- probDf[as.numeric(karyotype[i]), i]
    }
    divProb <- sum(probs)/sum(apply(probDf, 2, max))
    if(divProb < 0.7) {
        divProb <- 0
        return(divProb)
    } else {
        return(divProb)
    }

}

#' Calculate probability of survival in a fitness-dependent manner
#'
#' This function will calculate a probability of survival based on karyotypic fitness.
#'
#' @param karyotype A karyotype provided as a vector.
#' @param numberOfChromosomes The number of chromosomes in the karyotype.
#' @param probDf A fitness probability matrix.
#' @return A probability of survival.
#' @author Bjorn Bakker

calcSurProb <- function(karyotype, numberOfChromosomes = NULL, probDf = NULL) {

  # set-up collection vector
  if(is.null(numberOfChromosomes)) {
    numberOfChromosomes <- length(karyotype)
  }
  probs <- vector("numeric", numberOfChromosomes)

  # loop through karyotype to get the copy numbers
  for(i in 1:numberOfChromosomes) {
    probs[i] <- probDf[as.numeric(karyotype[i]), i]
  }

  # calculate and return survival probability
  surProb <-  sum(probs)/sum(apply(probDf, 2, max))
  return(surProb)

}

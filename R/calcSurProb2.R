#' Calculate probability of survival in a fitness-dependent manner
#'
#' This function will calculate a probability of survival based on karyotypic fitness.
#'
#' @param karyotype A karyotype provided as a vector.
#' @param numberOfChromosomes The number of chromosomes in the karyotype.
#' @param probDf A fitness probability matrix.
#' @return A probability of survival.
#' @author Bjorn Bakker
calcSurProb2 <- function(karyotypes, numberOfChromosomes = NULL, probDf = NULL) {

  # set-up collection vector
  if(is.null(numberOfChromosomes)) {
    numberOfChromosomes <- ncol(karyotypes)
  }
  probs <- array(NA, dim=dim(karyotypes), dimnames=dimnames(karyotypes))
  for (i in colnames(karyotypes)) {
      probs[,i] <- probDf[cbind(as.character(karyotypes[,i]), as.character(i))]
  }
  # calculate and return survival probability
  surProb <- rowSums(probs) / sum(apply(probDf, 2, max))

  return(surProb)

}

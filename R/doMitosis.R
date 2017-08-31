#' Perform a round of mitosis with possible mis-segregations
#'
#' This function will perform mitosis on a given karyotype, apply chromosome mis-segregation events if the criteria are met, and returns a matrix with two new karytoypes.
#'
#' @param karyotype A karyotype provided as a vector on which to apply the mitosis.
#' @param pMisseg The probability of mis-segregation per chromosome copy number.
#' @param cellID The cell ID of the provided karyotype to ensure faithful clone tracking.
#' @return A 2-row matrix with karyotypes from daughter cells.
#' @author Bjorn Bakker
#' @export

doMitosis <- function(karyotype, pMisseg, cellID) {

  # check user input
  numberOfChromosomes <- length(karyotype)

  # setup daughter mat and duplicate the karyotypes, collect clone ID
  daughterMat <- rbind(karyotype, karyotype)
  rownames(daughterMat) <- rep(cellID, times = 2)

  # recalculate weighted probabilities
  pChromWeighted <- vector("numeric", numberOfChromosomes)
  for(i in 1:numberOfChromosomes) {
    pChromWeighted[i] <- 1 - (1 - pMisseg)^(as.numeric(daughterMat[1, i]))
    #pChromWeighted[i] <- as.numeric(daughterMat[1, i])*pMisseg
  }

  # determine missegregation event per chromosome
  gainLoss <- vector("numeric", numberOfChromosomes)
  for(i in 1:numberOfChromosomes) {
    misseg <- runif(n = 1)
    if(misseg < pChromWeighted[i]) {
      gainLoss[i] <- sample(c(-1, 1), size = 1)
    } else {
      gainLoss[i] <- 0
    }
  }

  # apply missegregations
  daughterMat[1, ] <- as.numeric(daughterMat[1, ]) + gainLoss
  daughterMat[2, ] <- as.numeric(daughterMat[2, ]) - gainLoss

  # return ouput
  return(daughterMat)

}

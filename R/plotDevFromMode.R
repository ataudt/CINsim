#' Plot deviation from mode per chromosome.
#'
#' This function will plot the frequency at which cells deviate from the modal copy number per chromosome.
#'
#' @param karayoSim A karyoSim object that contains a karyotype matrix.
#' @return A plot of chromosomal deviation from mode.
#' @author Bjorn Bakker
#' @export

plotDevFromMode <- function(karyoSim) {

  # get table
  devTable <- devFromMode(karyoSim$karyotypeMatrix)

  # plot all chromosomes deviance
  yRangeVector <- seq(from = range(devTable)[1], to = range(devTable)[2], length.out = nrow(devTable))
  plot(x = 0:(nrow(devTable)-1), y = yRangeVector, type = "n",
       xlab = "Generation", ylab = "Deviance from mode",
       main = "Deviance from mode per chromosome")
  lineCol <- rainbow(n = ncol(devTable))
  for(i in 1:ncol(devTable)) {
    lines(x = 0:(nrow(devTable)-1), y = devTable[, i],
          type = "b", col = lineCol[i], lwd = 1, pch = i, lty = i)
  }
  legend("topleft", pch = 1:ncol(devTable), lty = 1:ncol(devTable), cex = 0.6,
         col = lineCol, legend = as.character(colnames(devTable)), ncol = 2)

}

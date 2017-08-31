#' Calculate the fraction of cells that deviate from the modal copy number
#'
#' This function will determine the fraction of cells that deviate from the model copy number of that chromosome.
#'
#' @param karyoSim A karyoSim object.
#' @param plot Whether to plot the results or return a matrix.
#' @param plotTitle A title for the plot.
#' @return A named vector or plot showing the deviation from modal chromosome copy number.
#' @author Bjorn Bakker
#' @export

devFromMode <- function(karyoSim, plot = FALSE, plotTitle = NULL) {

    # get karyotype matrix
    kldf <- karyoSim$karyotypeMatrix
    # no built-in function to get the mode in R, hence this function
    # get modal copy number per chromosome
    modalCopy <- apply(kldf, 2, function(x) {getMode(x)})
    # compare modal copy number with given copy number
    devianceKldf <- sweep(kldf, 2, modalCopy, '!=')
    # conver
    devianceKldf <- apply(devianceKldf, 2, function(x) {sum(x)/length(x)})
    # plot the fraction of cells deviating from mode
    if(is.null(plotTitle)) {
      plotTitle <- "Fraction of cells deviating from modal copy number"
    }
    if(plot) {
      barplot(devianceKldf, col = "black", xlab = "Chromosome", ylab = "Fraction of cells",
              main = plotTitle, cex.names = 0.8, ylim = c(0, 1))
    } else {
      return(devianceKldf)
    }

}

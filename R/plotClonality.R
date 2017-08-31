#' Plot clonal frequencies
#'
#' This function will allow you to plot the frequencies of the clones over time.
#'
#' @param karyoSim A karyoSim object that contains a summary matrix for clonality frequency per generations.
#' @param plot A logical whether to plot the results or return a table.
#' @param file File name (without extension) to which to print the plot.
#' @return A table or plot of the clonal frequencies.
#' @author Bjorn Bakker
#' @export

plotClonality <- function(karyoSim, plot = TRUE, file = NULL) {

  if(!is.null(file)) {
      pdf(file = paste0(file, ".pdf"), width = 16.69, height = 8.27)
  }

  clonality <- karyoSim$clonality
  row.names(clonality) <- 1:nrow(clonality)
  barColors <- rainbow(n = nrow(clonality))
  set.seed(41)
  barColors <- barColors[sample(1:length(barColors), replace = FALSE, size = length(barColors))]
  finalGen <- names(sort(clonality[, ncol(clonality)], decreasing = TRUE))
  finalGenCol <- barColors[as.numeric(finalGen)]
  finalGen <- head(finalGen, 5)
  finalGenCol <- head(finalGenCol, n = 5)
  if(plot) {
  barplot(clonality, col = barColors, names.arg = 0:(ncol(clonality)-1), border = NA,
          xlim = c(0, ncol(clonality) + 2.5),
          cex.names = 0.8, main = "Relative abundance of clones", space = 0,
          xlab = "Generation", ylab = "Frequency", legend.text = TRUE,
          args.legend = list(x = ncol(clonality) + 2.5, bty = "o",
                             legend = finalGen,
                             fill = finalGenCol))
  } else {
      return(clonality)
  }

  if(!is.null(file)) {
      dev.off()
  }

}

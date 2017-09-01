#' Simulation metrics plot
#'
#' This function will plot various simulation metrics into a single overview. Metrics include the cummulative cell number, fraction of cells surviving per round of cell division, karyotype measures (aneuploidy and heterogeneity), and clonal frequencies when the simulation starting population was greater than 1.
#'
#' @param karyoSim A karyoSim object that contains all the required information and metrics.
#' @param file A file name (without extension) to which the data shall be printed.
#' @return A plot showing the dynamics in various metrics over time.
#' @author Bjorn Bakker
#' @export

simPlots <- function(karyoSim, file = NULL) {

    if(class(karyoSim) != "karyoSim") {
        message("An object of class karyoSim is required")
        stop
    }

    if(!is.null(file)) {
        pdf(file = paste0(file, ".pdf"), width = 16.69, height = 8.27)
    }

    par(mfrow = c(2, 3), mar = c(5, 5, 2, 5))
    par(mar = c(5, 5, 2, 5))

    # survival curves and number of cells
    plotNumSur <- karyoSim$simulationDataFrame[, 1:2]
    plotNumSur$Generation <- 0:(nrow(plotNumSur)-1)
    plotNumSur$numberOfCells <- plotNumSur$numberOfCells
    firstYlab <- "Fraction surviving"
    if(!is.null(karyoSim$simulationDataFrame$meanSurvivalProb)) {
        meanSurvivalProb <- karyoSim$simulationDataFrame$meanSurvivalProb
        plotNumSur$meanSurProb <- meanSurvivalProb
        firstYlab <- "Fraction Surviving | mean pSurvival"
    }
    with(plotNumSur, plot(Generation, fractionSurviving, type = "l", lwd = 2,
                          col = "red3", ylab = firstYlab,
                          ylim = c(0, 1)))
    if(!is.null(plotNumSur$meanSurProb)) {
        par(new = T)
        with(plotNumSur, plot(Generation, meanSurProb, type = "l", lwd = 2,
                              col = "darkgreen", ylab = NA, ylim = c(0, 1)))
    }
    par(new = T)
    with(plotNumSur, plot(Generation, log = "y", numberOfCells, type = "l", lwd = 2,
                          col = "blue3", ylab = NA, xlab = NA, axes = F))
    # create logarithmic ticks for the right y-axis
    ticks <- seq(0, round(max(log(plotNumSur$numberOfCells, 10)), digits = 0), by=1)
    labels <- sapply(ticks, function(x) as.expression(bquote(10^ .(x))))
    axis(4, at=c(10^ticks), labels=labels)
    mtext(side = 4, line = 3, "Number of cells", cex = 0.8)
    # conditional legend
    if(!is.null(plotNumSur$meanSurProb)) {
        legend("bottomright", legend = c("Fraction surviving", "Mean pSurvival", "Number of cells"),
               col = c("red3", "darkgreen", "blue3"), lty = c(1, 1, 1), lwd = c(1, 1, 1), pch = c(NA, NA, NA))
    } else {
        legend("bottomright", legend = c("Fraction surviving", "Number of cells"),
               col = c("red3", "blue3"), lty = c(1, 1), lwd = c(1, 1), pch = c(NA, NA))
    }
    title(main = "Cell survival and number of cells")

    # copy number frequency
    cnFreq(karyoSim, plot = TRUE)

    # plot clonality
    if(nrow(karyoSim$clonality) > 1) {
        plotClonality(karyoSim)
    }

    # genome-wide karyotype measures over time
    aneuploidy <- karyoSim$simulationDataFrame$aneuploidy
    heterogeneity <- karyoSim$simulationDataFrame$heterogeneity
    yRangeVector <- seq(from = range(c(aneuploidy, heterogeneity))[1], to = range(c(aneuploidy, heterogeneity))[2], length.out = length(aneuploidy))
    plot(x = 0:(length(aneuploidy)-1), y = yRangeVector, type = "n",
         xlab = "Generation", ylab = "Score",
         main = "Genome-wide karyotype measures over time")
    lines(x = 0:(length(aneuploidy)-1), y = aneuploidy,
          type = "l", col = "blue", lwd = 2, lty = 1)
    lines(x = 0:(length(heterogeneity)-1), y = heterogeneity,
          type = "l", col = "red", lwd = 2, lty = 1)
    legend("topleft", pch = c(NA, NA), lty = c(1, 1), col = c("blue", "red"),
           legend = c("Aneuploidy", "Heterogeneity"))

    # aneuploidy per chromosome over time
    aneuScores <- karyoSim$aneuploidyPerChromosome
    yRangeVector <- seq(from = range(aneuScores)[1], to = range(aneuScores)[2], length.out = nrow(aneuScores))
    plot(x = 0:(nrow(aneuScores)-1), y = yRangeVector, type = "n",
         xlab = "Generation", ylab = "Aneuploidy score",
         main = "Aneuploidy score per chromosome over time")
    lineCol <- rainbow(n = ncol(aneuScores))
    for(i in 1:ncol(aneuScores)) {
        lines(x = 0:(nrow(aneuScores)-1), y = aneuScores[, i],
              type = "b", col = lineCol[i], lwd = 1, pch = i, lty = i)
    }
    legend("topleft", pch = 1:ncol(aneuScores), lty = 1:ncol(aneuScores), cex = 0.6,
           col = lineCol, legend = as.character(colnames(aneuScores)), ncol = 2)

    # heterogeneity per chromosome over time
    hetScores <- karyoSim$heterogeneityPerChromosome
    yRangeVector <- seq(from = range(hetScores)[1], to = range(hetScores)[2], length.out = nrow(hetScores))
    plot(x = 0:(nrow(hetScores)-1), y = yRangeVector, type = "n",
         xlab = "Generation", ylab = "Heterogeneity score",
         main = "Heterogeneity score per chromosome over time")
    lineCol <- rainbow(n = ncol(hetScores))
    for(i in 1:ncol(hetScores)) {
        lines(x = 0:(nrow(hetScores)-1), y = hetScores[, i],
              type = "b", col = lineCol[i], lwd = 1, pch = i, lty = i)
    }
    legend("topleft", pch = 1:ncol(hetScores), lty = 1:ncol(hetScores), cex = 0.6,
           col = lineCol, legend = as.character(colnames(hetScores)), ncol = 2)

    if(!is.null(file)) {
        dev.off()
    }
    par(mfrow = c(1, 1))
}

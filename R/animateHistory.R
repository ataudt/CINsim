#' Animate karyotype history
#'
#' This function will return a looping plot of the karyotype matrices over time.
#'
#' @param karyoSimHistory A list containing the karyotype matrices at every generation
#' @param subsetSize The maximum number of cells per generation to display
#' @param generations A subset of generations to plot, rather than the entire history. The 0th generation will always be plotted.
#' @param A file name with either a .gif, .mp4, .swf, or .html extension to which the plot will be saved.
#' @param The interval between the generation in seconds.
#' @return An animated ggplot.
#' @author Bjorn Bakker
#' @import reshape2
#' @import ggplot2
#' @import gganimate
#' @export

animateHistory <- function(karyoSimHistory, subsetSize = 1000, generations = NULL, file = NULL, interval = 1) {

    # check user input
    if(class(karyoSimHistory) != "karyoSimHistory") {
        message("An object of class karyoSimHistory must be provided")
        stop
    }

    # subselect generations if desired
    if(!is.null(generations)) {
        generationNames <- c("gen_0", paste0("gen_", generations))
        karyoSimHistory <- karyoSimHistory[generationNames]
    }

    # generate a karyotype data frame composite with generation indexes
    karyoMatDf <- data.frame()
    chromosomes <- colnames(karyoSimHistory[[1]])

    # compile each generation
    for(generationName in generationNames) {

        karyotypeTemp <- as.data.frame(karyoSimHistory[[generationName]])
        row.names(karyotypeTemp) <- NULL

        # subset for a particular size
        if(!is.null(subsetSize)) {
            if(nrow(karyotypeTemp) > subsetSize) {
                subset <- sample(1:nrow(karyotypeTemp), size = subsetSize, replace = FALSE)
                karyotypeTemp <- karyotypeTemp[subset, ]
            }
        }

        if(nrow(karyotypeTemp) < 2) {
            ord <- 1
        } else {
            ord <- hclust(dist(karyotypeTemp[, 1:length(chromosomes)], method = "euclidean"), method = "ward.D")$order
        }
        karyotypeTemp$generation <- rep(generationName, times = nrow(karyotypeTemp))
        karyotypeTemp$cellID <- paste0("cell_", 1:nrow(karyotypeTemp))
        karyotypeTemp$cellID <- factor(karyotypeTemp$cellID, levels = karyotypeTemp$cellID[ord])
        karyoMatDf <- rbind(karyoMatDf, karyotypeTemp)
    }

    # reshape the data frame for compatibility with ggplot2
    kldfPlot <- melt(data = karyoMatDf, measure.vars = chromosomes, variable.name = "chromosome")
    colnames(kldfPlot) <- c("generation", "cellID", "chromosome", "copyNumber")

    # factor generation, cellID and copy number
    kldfPlot$generation <- factor(kldfPlot$generation, levels = names(karyoSimHistory))
    kldfPlot$copyNumber <- factor(kldfPlot$copyNumber)

    # set copy number colors
    copyNumberColors <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy", "lemonchiffon", "dodgerblue", "chartreuse4", "lightcoral", "aquamarine2")
    names(copyNumberColors) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

    # plot
    gPlot <- ggplot(kldfPlot, aes(x = chromosome, y = cellID, frame = generation)) +
        geom_tile(aes(fill = copyNumber)) +
        scale_fill_manual(values = copyNumberColors) +
        theme_classic() +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())

    if(is.null(file)) {
        gganimate(gPlot, interval = interval)
    } else {
        gganimate(gPlot, interval = interval, filename = file)
    }

}

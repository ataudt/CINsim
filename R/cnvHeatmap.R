#' Plot a copy number heatmap.
#'
#' This function will plot a copy number heatmap based on a karyotype matrix.
#'
#' @param karyoSim A karyoSim object (final output of CINsim).
#' @param subsetSize The number of cells to subset and plot from the karyotype matrix.
#' @param clones A character vector of the clone identifiers to plot.
#' @param file A file name (without extension) to which the plot will be printed.
#' @return A heatmap with copy numbers shown in colours.
#' @import ggplot2
#' @import reshape2
#' @author Bjorn Bakker
#' @export

cnvHeatmap <- function(karyoSim = NULL, subsetSize = 1000, clones = NULL, file = NULL) {

    # check user input
    if(class(karyoSim) != "karyoSim") {
        message("An object of class karyoSim is required")
        stop
    }

    if(is.null(file)) {
        plotToFile <- FALSE
    } else {
        plotToFile <- TRUE
    }

    # subselect karyotype data frame and get chromosome names
    kldfPlot <- as.data.frame(karyoSim$karyotypeMatrix)
    chromosomes <- colnames(kldfPlot)

    # check subset size
    if(!is.null(subsetSize)) {
        if(nrow(kldfPlot) < subsetSize) {
            subsetSize <- NULL
        }
    }

    # subset a particular clone
    if(!is.null(clones)) {
        clones <- as.character(clones)
        cloneIDs <- sapply(row.names(kldfPlot), function(x) {unlist(strsplit(x, split = "_"))[2]})
        names(cloneIDs) <- NULL
        if(all(clones %in% cloneIDs)) {
            kldfPlotClones <- data.frame()
            for(clone in clones) {
                kldfPlotClones <- rbind(kldfPlotClones, kldfPlot[clone == cloneIDs, ])
            }
            if(nrow(kldfPlotClones) == 0) {
                message("Clone not present - plotting all clones")
                remove(kldfPlotClones)
            } else {
                kldfPlot <- kldfPlotClones
                remove(kldfPlotClones)
                numCells <- paste0(nrow(kldfPlot), "/", nrow(karyoSim$karyotypeMatrix))
            }
        } else {
            message("Some or all clone IDs not found in list - plotting all clones")
        }
    }

    # subset for a particular size
    if(is.null(subsetSize)) {
        numCells <- as.character(nrow(kldfPlot))
    } else {
        subset <- sample(1:nrow(kldfPlot), size = subsetSize, replace = FALSE)
        kldfPlot <- kldfPlot[subset, ]
        numCells <- paste0(length(subset), "/", nrow(karyoSim$karyotypeMatrix))
    }

    # add cell IDs to data frame, and determine euclidean distance/order by cell ID
    ord <- hclust(dist(kldfPlot[, 1:length(chromosomes)], method = "euclidean"), method = "ward.D")$order
    kldfPlot$cellID <- paste0("cell_", 1:nrow(kldfPlot))

    # reshape df to allow for plotting using ggplot2
    library(reshape2)
    kldfPlot <- melt(data = kldfPlot, measure.vars = chromosomes, variable.name = "chromosome")
    colnames(kldfPlot) <- c("cellID", "chromosome", "copyNumber")

    # factor cell IDs and copy numbers
    kldfPlot$cellID <- factor(kldfPlot$cellID, levels = kldfPlot$cellID[ord])
    kldfPlot$copyNumber <- factor(kldfPlot$copyNumber)

    # set copy number colors
    copyNumberColors <- c("gray90", "darkorchid3", "springgreen2", "red3", "gold2", "navy", "lemonchiffon", "dodgerblue", "chartreuse4", "lightcoral", "aquamarine2")
    names(copyNumberColors) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

    # set plot title
    plotTitle <- paste0("pMisseg: ", karyoSim$simulationParameters["pMisseg"],
                        ", generations: ", nrow(karyoSim$simulationDataFrame)-1,
                        ", number of cells: ", numCells)
    if(!is.null(clones)) {
        plotTitle <- paste0(plotTitle, ", clones: ", clones)
    }

    # plot
    gPlot <- ggplot(kldfPlot, aes(x = chromosome, y = cellID)) +
        geom_tile(aes(fill = copyNumber)) +
        scale_fill_manual(values = copyNumberColors) +
        theme_classic() +
        theme(axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
        labs(title = plotTitle)

    if(plotToFile) {
        ggsave(plot = gPlot, filename = paste0(file, ".pdf"), width = 10, height = 8)
    } else {
        return(gPlot)
    }

}

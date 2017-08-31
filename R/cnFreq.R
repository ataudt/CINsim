#' Calculate copy number frequency
#'
#' This function will determine the frequencies of copy numbers per chromosomes, and return either a table or a barplot.
#'
#' @param karyoSim A karyoSim object (final output of the CINsim function).
#' @param plot A logical, whether to plot the result or not.
#' @return A table or plot with the frequencies of the chromosome copy number states.
#' @author Bjorn Bakker
#' @export

cnFreq <- function(karyoSim, plot = TRUE) {

    # check user input
    if(class(karyoSim) != "karyoSim") {
        message("An object of class karyoSim is required")
        stop
    }

    # get karyotype data frame
    kldf <- as.data.frame(karyoSim$karyotypeMatrix)

    # define chromosomes and set up data frame
    chromosomes <- colnames(kldf)
    countDf <- data.frame()

    # collect the frequencies and output a data frame
    for(chromosome in chromosomes) {
        countVec <- c()
        for(i in 1:8) {
            countVec <- c(countVec, sum(i == kldf[, chromosome, drop = FALSE]))
        }
        countDf <- rbind(countDf, countVec)
    }
    countDf <- t(countDf)

    # tranlate to relative frequencies
    countDf <- apply(countDf, 2, function(x) { x/sum(x) } )

    # redefine column names
    rownames(countDf) <- as.character(1:8)

    # return plot, else the data frame
    if(plot) {
        copyNumberColors <- c("darkorchid3", "springgreen2",
                              "red3", "gold2", "navy", "lemonchiffon",
                              "dodgerblue", "chartreuse4")
        barplot(countDf, col = copyNumberColors, names.arg = chromosomes,
                xlim = c(0, ncol(countDf) + 6), cex.names = 0.8,
                main = "Relative copy number state frequency",
                xlab = "Chromosome", ylab = "Frequency", legend.text = TRUE,
                args.legend = list(x = ncol(countDf) + 7, y = max(colSums(countDf)), bty = "n"))
    } else {
        return(countDf)
    }

}

#' Check viability based on karyotype
#'
#' This examines a karyotype and then determines whehter the cell survives.
#'
#' @param karyotype A single karyotype, supplied as a vector.
#' @param minMonosomyLethal The minimum number of monosomies that will cause the cell to die.
#' @param minNumEuploidChr The minimum number of euploid chromosomes that must remain before the cell dies.
#' @param numberOfChromosomes The number of chromosomes in the karyotype.
#' @param probDf A probability fitness matrix.
#' @return A logical, whether the cell survives or not.
#' @author Bjorn Bakker

checkViability2 <- function(karyotype, minMonosomyLethal, minNumEuploidChr, numberOfChromosomes = NULL, probDf = NULL) {

    # check user input
    if(is.null(numberOfChromosomes)) {
        numberOfChromosomes <- length(karyotype)
    }
    cellSurvival <- NULL

    chrCounts <- table(as.numeric(karyotype))
    monosomies <- chrCounts["1"]
    euploidies <- chrCounts["2"]
    if(!is.na(chrCounts["0"])) {
        return(FALSE)
    }
    if(!is.na(chrCounts["9"])) {
        return(FALSE)
    }
    if(!is.na(monosomies)) {
        if(monosomies >= minMonosomyLethal) {
            return(FALSE)
        }
    }
    if(!is.na(euploidies)) {
        if(euploidies < minNumEuploidChr) {
            return(FALSE)
        }
    }
    if(!is.null(probDf)) {
        survivalThreshold <- runif(n = 1)
        pSurvival <- calcSurProb(karyotype, numberOfChromosomes, probDf)
        if(pSurvival < survivalThreshold) {
            return(FALSE)
        }
    }

    return(TRUE)

}

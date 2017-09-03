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

checkViability2 <- function(karyotypes, minMonosomyLethal, minNumEuploidChr, numberOfChromosomes = NULL, probDf = NULL) {

    # check user input
    if(is.null(numberOfChromosomes)) {
        numberOfChromosomes <- ncol(karyotypes)
    }
    cellSurvival <- rep(TRUE, nrow(karyotypes))
    names(cellSurvival) <- rownames(karyotypes)

    is0or9 <- apply(karyotypes, 2, function(x) { x == 0 | x == 9 })
    num0or9 <- rowSums(is0or9)

    isMonosomy <- apply(karyotypes, 2, function(x) { x == 1 })
    numMonosomy <- rowSums(isMonosomy)

    isEuploid <- apply(karyotypes, 2, function(x) { x == 2 })
    numEuploid <- rowSums(isEuploid)

    cellSurvival <- (num0or9 == 0) & (numMonosomy < minMonosomyLethal) & (numEuploid >= minNumEuploidChr)

    if(!is.null(probDf)) {
        survivalThreshold <- runif(n = nrow(karyotypes))
        pSurvival <- calcSurProb2(karyotypes, numberOfChromosomes, probDf)

        survivesProbdf <- pSurvival >= survivalThreshold
        cellSurvival <- cellSurvival & survivesProbdf
    }


    return(cellSurvival)

}

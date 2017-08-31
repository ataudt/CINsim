#' Simulate chromosomal instability.
#'
#' This the the main function for the CINsim package, that will allow for simulation of cell divisions and chromosomal instability.
#'
#' @param karyotypes A matrix with cells in rows and chromosomes in colums. Rows and colums must be labelled appropriately.
#' @param euploidRef The euploid chromosome copy number (used for calculating levels of aneuploidy).
#' @param g The maximum number of generations to be simulated.
#' @param pMisseg The mis-segregation probability per chromosome copy.
#' @param pMissegG A vector of generation numbers during which pMisseg is set at the desired values. Will otherwise default to 0.
#' @param pDivision The base probability of cell division.
#' @param fitDivision A logical whether pDivision should be karyotype (i.e. fitness) dependent.
#' @param probDf A probability fitness matrix (max 8 rows for the max 8 copy number states).
#' @param karyotypeSelection Whether karyotype selection should take place.
#' @param minMonosomyLethal The minimum number of monosomies that will cause the cell to die.
#' @param minNumEuploidChr The minimum number of euploid chromosomes that must remain before the cell dies.
#' @param downSample The maximum size of the simulated population before down-sampling to 12.5% occurs.
#' @param maxNumCells The maximum number of (theoretical) cells to simulate before the simulation is terminated.
#' @param simTitle A custom title linked to the final output when saving to file.
#' @param saveToFile A logical whether to the final simulation ouput to file.
#' @param numFreeCPU The number of CPUs to keep free for use when the simulation is running.
#' @return A karyoSim object containing all relevant information of the simulation.
#' @author Bjorn Bakker
#' @import parallel
#' @import foreach
#' @import doParallel
#' @export

Cinsim <- function(karyotypes = NULL,
                        euploidRef = NULL,
                        g = 12,
                        pMisseg = 0.01,
                        pMissegG = NULL,
                        fitMisseg = FALSE,
                        pDivision = 0.8,
                        fitDivision = FALSE,
                        probDf = NULL,
                        karyotypeSelection = FALSE,
                        minMonosomyLethal = 3,
                        minNumEuploidChr = 0,
                        downSample = 2.5e+04,
                        maxNumCells = 2e+09,
                        simTitle = "CINsim",
                        saveToFile = FALSE,
                        numFreeCPU = 1) {

  # start timed message
  message("==> Starting simulation <==")
  time0 <- proc.time()

  # check user input for karyotypes
  if (is.null(karyotypes)) {
    # generate a mouse karyotype if none is provided
    message("No karyotype provided - using default female mouse karyotype")
    karyoMat <- makeKaryotypes(n = 1)
    numberOfChromosomes <- ncol(karyoMat)
  } else {
    karyoMat <- karyotypes
    numberOfChromosomes <- ncol(karyoMat)
  }

  # start timed message
  ptm <- startTimedMessage("Initializing simulation ...")

  # check user input (karyotypes, euploidRef) and set initial size of cell population
  if(is.null(euploidRef)) {
    euploidRef <- 2
  }
  chromosomes <- colnames(karyoMat)
  numberOfChromosomes <- length(chromosomes)
  sizeInitial <- nrow(karyoMat[])
  indexSampled <- 0

  # check status of pdivision, and whether to apply asynchronous cell divisions
  if(pDivision == 1) {
    asynch <- FALSE
  } else if(pDivision ==0) {
    message("pDivision must be greater than 0, applying default of 0.8")
    pDivision <- 0.8
    asynch <- TRUE
  } else {
    asynch <- TRUE
  }

  # set-up collection vectors and matrices
  fractionSurviving <- c(1, vector("numeric", g))
  numberOfCells <- c(sizeInitial, vector("numeric", g))
  numberOfDaughters <- c(0, vector("numeric", g))
  numberOfSurvivors <- c(sizeInitial, vector("numeric", g))

  aneuploidyMat <- matrix(NA, nrow = g + 1, ncol = numberOfChromosomes)
  aneuploidyMat[1, ] <- qAneuploidy(karyoMat, numberOfCells[1], euploidRef)
  heterogeneityMat <- matrix(NA, nrow = g + 1, ncol = numberOfChromosomes)
  heterogeneityMat[1, ] <- qHeterogeneity(karyoMat)

  kmGenomeWide <- matrix(NA, nrow = g + 1, ncol = 2)
  kmGenomeWide[1, ] <- c(mean(aneuploidyMat[1, ]), mean(heterogeneityMat[1, ]))

  # collect mean survival probability if probDf is defined
  if(!is.null(probDf)) {
    meanSurvivalProb <- vector("numeric", g + 1)
    if(sizeInitial > 1) {
      meanSurvivalProb[1] <- mean(apply(karyoMat, 1, function(x) {calcSurProb(x, numberOfChromosomes, probDf)}))
    } else {
      meanSurvivalProb[1] <- calcSurProb(karyoMat, numberOfChromosomes, probDf)
    }
  } else {
    meanSurvivalProb <- "undefined"
  }

  # set-up clonality tracking
  clonalityMat <- matrix(NA, nrow = g + 1, ncol = sizeInitial)
  clonalityMat[1, ] <- calcClonality(karyoMat, sizeInitial)

  # if pMisseg is fit-dependent, set pMisseg co-efficient
  if(fitMisseg) {
    euploidRefKaryotype <- rep(euploidRef, times = numberOfChromosomes)
    pMissegCoef <- calcSurProb(euploidRefKaryotype, numberOfChromosomes, probDf)
  } else {
    pMissegCoef <- 1
  }

  # set multi clore processing for improved speed at high cell numbers
  numCores <- detectCores() - numFreeCPU
  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  clusterExport(cl, varlist = c("calcDivProb", "numberOfChromosomes", "probDf", "pDivision", "checkViability", "calcSurProb",
                                "pMisseg", "pMissegCoef", "minMonosomyLethal", "minNumEuploidChr"), envir = environment())

  stopTimedMessage(ptm)

  # run simulation for g number of generations
  for(j in 1:g) {

    # start timed message
    ptm <- startTimedMessage("Processing generation ", j, " out of ", g, " ...")

    # check whether pMisseg is still applicable
    if(!is.null(pMissegG)) {
      if(j %in% pMissegG) {
        pMisseg <- 0
      }
    }

    # set number of cells at start of generation
    numCells <- nrow(karyoMat)

    # asynchronous cell division pipeline
    if(asynch) {

      # set pDivision as a function of karyotypic fitness
      if(fitDivision) {
        if(numCells > 1) {
          dividers <- parApply(cl, karyoMat, 1, function(x) {return(runif(n = 1) < calcDivProb(x, numberOfChromosomes, probDf))})
        } else {
          dividers <- runif(n = 1) < calcDivProb(karyoMat, numberOfChromosomes, probDf)
        }
      } else {
        if(numCells > 1) {
          dividers <- parApply(cl, karyoMat, 1, function(x) {return(runif(n = 1) < pDivision)})
        } else {
          dividers <- runif(n = 1) < pDivision
        }
      }
      nonDividers <- !(dividers)
      karyoMatDividers <- karyoMat[dividers, , drop = FALSE]
      karyoMatNonDividers <- karyoMat[nonDividers, , drop = FALSE]

      # apply mitoses if more than 1 dividing cell is present
      numDividers <- sum(dividers)
      if(numDividers > 0) {

        # determine pMisseg when karyotype dependent
        if(fitMisseg) {
          pMissegs <- parApply(cl, karyoMatDividers, 1, function(x) {calcMisseg(x, pMisseg, pMissegCoef, numberOfChromosomes)})
        } else {
          pMissegs <- rep(pMisseg, times = numDividers)
        }

        cellIDs <- rownames(karyoMat)
        karyoMatDividers <- foreach(cell = 1:numDividers,
                                    .combine = rbind,
                                    .export = c('doMitosis')) %dopar%
          doMitosis(karyoMatDividers[cell, , drop = FALSE], pMissegs[cell], cellIDs[cell])

        karyoMat <- rbind(karyoMatNonDividers, karyoMatDividers)

      }

    } else {

      # determine pMisseg when karyotype dependent
      if(fitMisseg) {
        pMissegs <- parApply(cl, karyoMat, 1, function(x) {calcMisseg(x, pMisseg, pMissegCoef)})
      } else {
        pMissegs <- rep(pMisseg, times = numCells)
      }

      cellIDs <- rownames(karyoMat)
      karyoMat <- foreach(cell = 1:numCells,
                          .combine = rbind,
                          .export = c('doMitosis', 'karyoMat', 'pMissegs', 'cellIDs')) %dopar%
        doMitosis(karyoMat[cell, , drop = FALSE], pMissegs, cellIDs[cell])

    }

    # determine number of daughters after mitosis
    numDaughters <- nrow(karyoMat)

    # remove unviable cells
    if(karyotypeSelection) {
      viableCells <- unlist(parApply(cl, karyoMat, 1, function(x) {checkViability(x, minMonosomyLethal, minNumEuploidChr, numberOfChromosomes, probDf)}))
    } else {
      viableCells <- unlist(parApply(cl, karyoMat, 1, function(x) {checkViability(x, minMonosomyLethal, minNumEuploidChr, numberOfChromosomes)}))
    }

    # count the number of viable cells
    numSurvivingCells <- sum(viableCells)

    # break loop if no surviving cells exist
    if(numSurvivingCells == 0) {

      message("No more surviving cells - exiting simulation")
      j <- j - 1
      break

    } else {

      # sub-select the viable cells
      karyoMat <- karyoMat[viableCells, , drop = FALSE]

      # down sample if the cell count exceeds the predetermined threshold
      if(!is.null(downSample)) {

        if(numSurvivingCells > downSample) {
          numSampled <- round(0.125 * numSurvivingCells)
          sampleSubset <- sample(1:numSurvivingCells, size = numSampled)
          karyoMat <- karyoMat[sampleSubset, , drop = FALSE]
          indexSampled <- indexSampled + 1
          trueCellCount <- numSampled * (1 / 0.125)^(indexSampled)
        } else {
          if(indexSampled > 0) {
            trueCellCount <- numSurvivingCells * (1 / 0.125)^(indexSampled)
          } else {
            trueCellCount <- numSurvivingCells
          }
        }
      } else {
        trueCellCount <- numSurvivingCells
      }

      # collect and store variables
      fractionSurviving[j+1] <- numSurvivingCells/numDaughters
      numberOfCells[j+1] <- trueCellCount
      numberOfDaughters[j+1] <- numDaughters
      numberOfSurvivors[j+1] <- numSurvivingCells
      aneuploidyMat[j+1, ] <- qAneuploidy(karyoMat, numSurvivingCells, euploidRef)
      heterogeneityMat[j+1, ] <- qHeterogeneityPar(karyoMat, nrow(karyoMat), cl)
      kmGenomeWide[j+1, ] <- c(mean(aneuploidyMat[j+1, ]), mean(heterogeneityMat[j+1, ]))
      clonalityMat[j+1, ] <- calcClonality(karyoMat, sizeInitial)
      if(!is.null(probDf)) {
        meanSurvivalProb[j+1] <- mean(parApply(cl, karyoMat, 1, function(x) {calcSurProb(x, numberOfChromosomes, probDf)}))
      }

      # set-up clonality tracking if multiple start cells are used
      if(sizeInitial > 1) {
        clonalityMat[j+1, ] <- calcClonality(karyoMat, sizeInitial)
      } else {
        clonalityMat[j+1, ] <- 1
      }

    }

    stopTimedMessage(ptm)

    # break if maximum number of cells is reached
    if(trueCellCount > maxNumCells) {
      message("True cell count greater than threshold: ", trueCellCount, " - exiting simulation")
      break
    }

  }

  # setting up timed message for parameter compiline
  ptm <- startTimedMessage("Compiling all simulation parameters ...")

  # compile generation specific parameters
  generationRowNames <- paste0("gen_", 0:j)

  genDataFrame <- data.frame(numberOfCells = numberOfCells[1:(j+1)],
                             fractionSurviving = fractionSurviving[1:(j+1)],
                             numberOfDaughters = numberOfDaughters[1:(j+1)],
                             numberOfSurvivors = numberOfSurvivors[1:(j+1)])
  colnames(kmGenomeWide) <- c("aneuploidy", "heterogeneity")
  kmGenomeWide <- kmGenomeWide[1:(j+1), ]
  genDataFrame <- cbind(genDataFrame, kmGenomeWide)
  if(!(any(meanSurvivalProb == "undefined"))) {
    meanSurvivalProb <- meanSurvivalProb[1:(j+1)]
    genDataFrame <- cbind(genDataFrame, meanSurvivalProb)
  }
  row.names(genDataFrame) <- generationRowNames

  # compile karyotype measures
  aneuploidyMat <- aneuploidyMat[1:(j+1), ]
  rownames(aneuploidyMat) <- generationRowNames
  colnames(aneuploidyMat) <- chromosomes
  heterogeneityMat <- heterogeneityMat[1:(j+1), ]
  rownames(heterogeneityMat) <- generationRowNames
  colnames(heterogeneityMat) <- chromosomes

  # compule clonality metrics if applicable
  clonalityMat <- t(clonalityMat[1:(j+1), ])
  colnames(clonalityMat) <- generationRowNames
  rownames(clonalityMat) <- paste0("clone_", 1:sizeInitial)

  # compile the simulation parameters
  simulationTime <- proc.time() - time0
  simParameters <- c(sizeInitial,
                     g,
                     j,
                     pMisseg,
                     fitMisseg,
                     asynch,
                     pDivision,
                     fitDivision,
                     maxNumCells,
                     minNumEuploidChr,
                     minMonosomyLethal,
                     tail(fractionSurviving, n = 1),
                     round(simulationTime[3], 2))

  names(simParameters) <- c("sizeInitialPopulation",
                            "maxGenerations",
                            "simulatedGenerations",
                            "pMisseg",
                            "pMissegFitnessDependent",
                            "asynchronousDivisions",
                            "pDivision",
                            "pDivisionFitnessDependent",
                            "maxNumberOfCells",
                            "minNumberOfEuploidChromosomes",
                            "minNumberOfLethalMonosomies",
                            "finalFractionSurviving",
                            "simulationTime (s)")

  # compile final output list
  karyoSim <- list(karyoMat,
                   genDataFrame,
                   simParameters,
                   aneuploidyMat,
                   heterogeneityMat,
                   clonalityMat)

  names(karyoSim) <- c("karyotypeMatrix",
                       "simulationDataFrame",
                       "simulationParameters",
                       "aneuploidyPerChromosome",
                       "heterogeneityPerChromosome",
                       "clonality")

  # append probDf if provided
  if(!is.null(probDf)) {
    karyoSim[["probabilityDataFrame"]] <- probDf
  }

  if(saveToFile) {
    save(karyoSim, file = paste(simTitle, pMisseg, paste0(j, "g"), "final_output.RData", sep = "_"))
  }

  stopTimedMessage(ptm)

  time1 <- proc.time() - time0
  message("==| Simulation complete - final time: ", round(time1[3], 2), "s |==")

  stopImplicitCluster()
  return(karyoSim)

}

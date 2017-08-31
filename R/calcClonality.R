#' Calculate frequency of clones in a population
#'
#' This function determines the frequency of clones in a karyotype matrix, using the matrix row names as identifiers for the parental clones
#'
#' @param karyoMat A karyotype matrix, with cells in rows and chromosomes in columns.
#' @param sizeInitial The initial size of the total cell population at the start of the simulation.
#' @return A table containting the frequencies of clones.
#' @author Bjorn Bakker


calcClonality <- function(karyoMat, sizeInitial = NULL) {

  # check user input
  if(sizeInitial == 1) {
    return(1)
  } else {

      # get clone IDs and their frequencies
      cloneIDs <- sapply(rownames(karyoMat), function(x) {as.numeric(unlist(strsplit(x, split = "_"))[2])})
      cloneTable <- table(cloneIDs)
      cloneIDs <- unique(names(cloneTable))
      cloneIDsNum <- as.numeric(cloneIDs)

      # calculate frequency of all clones: actual and theoretical
      missingCloneIDs <- c()
      missingClonesNum <- 0
      for(i in 1:sizeInitial) {
        if(i %in% cloneIDsNum) {
          next
        } else {
          missingCloneIDs <- c(missingCloneIDs, i)
          missingClonesNum <- missingClonesNum + 1
        }
      }

      # fill in the missing or extinct clones
      missingClones <- rep(0, times = missingClonesNum)
      names(missingClones) <- missingCloneIDs
      cloneTable <- c(cloneTable, missingClones)

      # determine relative frequencies
      cloneTable <- cloneTable / sum(cloneTable)
      cloneTable <- cloneTable[as.character(1:sizeInitial)]

  }

  # return table
  return(cloneTable)

}

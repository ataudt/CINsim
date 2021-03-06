% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkViability.R
\name{checkViability}
\alias{checkViability}
\title{Check viability based on karyotype}
\usage{
checkViability(karyotype, minMonosomyLethal, minNumEuploidChr,
  numberOfChromosomes = NULL, probDf = NULL)
}
\arguments{
\item{karyotype}{A single karyotype, supplied as a vector.}

\item{minMonosomyLethal}{The minimum number of monosomies that will cause the cell to die.}

\item{minNumEuploidChr}{The minimum number of euploid chromosomes that must remain before the cell dies.}

\item{numberOfChromosomes}{The number of chromosomes in the karyotype.}

\item{probDf}{A probability fitness matrix.}
}
\value{
A logical, whether the cell survives or not.
}
\description{
This examines a karyotype and then determines whehter the cell survives.
}
\author{
Bjorn Bakker
}

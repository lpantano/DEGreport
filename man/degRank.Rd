\name{degRank}
\alias{degRank}
\title{Get rank data frame with best score on the top}
\usage{
degRan(g1,g2,counts,fc,popsize)
}
\arguments{
  \item{g1}{list of samples in group 1}

  \item{g2}{list of samples in group 2}

  \item{counts}{count matrix for each gene and each sample
  that is deregulated}

  \item{fc}{list of FC of deregulated genes. Should be same
  length than counts \code{row.names}}

  \item{popsize}{number of combinations to generate}
}
\value{
data frame with the output of \link{degBIcmd} for each gene
}
\description{
Get rank data frame with best score on the top
}


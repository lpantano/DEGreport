\name{degVar}
\alias{degVar}
\title{Distribution of pvalues by SD range}
\usage{
degVar(pvalues,counts)
}
\arguments{
  \item{pvalues}{pvalues of DEG analysis}

  \item{counts}{matrix with counts for each samples and
  each gene. row number should be the same length than
  pvalues vector.}
}
\value{
ggplot2 object
}
\description{
Distribution of pvalues by SD range
}


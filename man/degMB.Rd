\name{degMB}
\alias{degMB}
\title{Distribution of expression of DE genes compared to the background}
\usage{
degMB(tags,g1,g2,counts,pop=400)
}
\arguments{
  \item{tags}{list of genes that are DE}

  \item{g1}{list of samples in group 1}

  \item{g2}{list of samples in group 2}

  \item{counts}{matrix with counts for each samples and
  each gene. Should be same length than pvalues vector.}

  \item{pop}{number of random samples taken for background
  comparison}
}
\value{
ggplot2 object
}
\description{
Distribution of expression of DE genes compared to the
background
}
\examples{
data(DEGreportSet)
detag<-row.names(DEGreportSet$deg[1:10,])
degMB(detag,DEGreportSet$g1,DEGreportSet$g2,DEGreportSet$counts)
}


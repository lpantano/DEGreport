\name{figurebyvar}
\alias{figurebyvar}
\title{Wrap figure from \code{degVB} into a Nozzle object}
\usage{
figurebyvar(tags, g1, g2, counts, out, pop = 400)
}
\arguments{
  \item{tags}{genes of DEG analysis}

  \item{g1}{group 1}

  \item{g2}{group 2}

  \item{counts}{matrix with counts for each samples and
  each gene. row number should be the same length than
  pvalues vector.}

  \item{out}{path to save the figure}

  \item{pop}{random genes for background}
}
\value{
Nozzle object
}
\description{
Wrap figure from \code{degVB} into a Nozzle object
}


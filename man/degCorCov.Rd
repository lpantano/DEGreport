% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariate.R
\name{degCorCov}
\alias{degCorCov}
\title{Calculate the correlation relationshipt among all covariates
in the metadata table}
\usage{
degCorCov(metadata, fdr = 0.05, use_pval = FALSE, ...)
}
\arguments{
\item{metadata}{data.frame with samples metadata.}

\item{fdr}{numeric value to use as cutoff to determine
the minimum fdr to consider significant correlations
between pcs and covariates.}

\item{use_pval}{boolean to indicate to use p-value instead of FDR to
hide non-significant correlation.}

\item{...}{Parameters to pass to \code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}}.}
}
\value{
: list:
a) cor, data.frame with pair-wise correlations, pvalues, FDR
b) corMat, data.frame with correlation matrix
c) fdrMat, data.frame with FDR matrix
b) plot, Heatmap plot of correlation matrix
}
\description{
This function will calculate the correlation among
all columns in the metadata
}
\examples{
data(humanGender)
library(DESeq2)
idx <- c(1:10, 75:85)
dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
  colData(humanGender)[idx,], design=~group)
cor <- degCorCov(colData(dse))
}
\author{
: Lorena Pantano, Kenneth Daily and Thanneer Malai Perumal
}

#' Deprecated functions in package DEGreport
#' 
#' These functions are provided for compatibility with older versions
#' of DEGreport only and will be defunct at the next release.
#' 
#' @rdname DEGreport-deprecated
#' 
#' @details
#' The following functions are deprecated and will be made defunct; 
#' use the replacement indicated below:
#' * degRank, degPR, degBIcmd, degBI, degFC, degComb, degNcomb: DESeq2::lcfShrink. 
#' This function was trying to avoid big FoldChange
#' in variable genes. There are other methods nowadays like lcfShrink function.


#' DEGreport
#' 
#' @importFrom S4Vectors SimpleList DataFrame
#' @importMethodsFrom BiocGenerics plotMA
#' @importFrom "grDevices" "dev.off" "jpeg" "pdf" "rgb" "colorRampPalette"
#' @importFrom "methods" "slotNames"
#' @importFrom "stats" "as.dist" "as.hclust" "cor" "cor.test"
#'             "cutree" "density" "median.default" "quantile" "sd"
#'             "time" "var" "na.omit" "aov" "prcomp" "median"
#'             "cmdscale" "dist" "hclust" "p.adjust" "xtabs"
#'             "as.dendrogram"
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus calcICL
#' @import ggrepel
#' @import ggplot2
#' @importFrom ggplot2 layer ggproto
#' @importFrom ggdendro dendro_data
#' @import utils
#' @import Nozzle.R1
#' @import edgeR
#' @import cluster
#' @import logging
#' @importFrom Biobase rowMax rowMin
#' @importFrom DESeq2 plotCounts rlog results resultsNames lfcShrink
#' @importFrom SummarizedExperiment SummarizedExperiment
#'             colData assay assays rowData
#' @importFrom psych corr.test
#' @importFrom grid textGrob
#' @importFrom broom tidy
#' @importFrom cowplot draw_plot
#' @importFrom cowplot ggdraw
#' @importFrom cowplot plot_grid
#' @importFrom dplyr arrange
#' @importFrom dplyr select select_if
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate mutate_all mutate_if
#' @importFrom dplyr filter rename
#' @importFrom dplyr group_by ungroup
#' @importFrom dplyr "%>%"
#' @importFrom dplyr n
#' @importFrom dplyr distinct
#' @importFrom dplyr left_join inner_join right_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom tidyr spread gather unite
#' @importFrom tibble column_to_rownames remove_rownames rownames_to_column
#' @importFrom tibble as_tibble
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape melt melt.data.frame
#' @importFrom knitr kable
#' @importFrom methods new slot
#' @importFrom rlang sym
#' @importFrom magrittr set_colnames set_rownames
#' @importFrom stringr str_split
#' @import lasso2
"_PACKAGE"

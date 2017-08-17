#' Deprecated functions in package DEGreport
#' 
#' These functions are provided for compatibility with older versions
#' of DEGreport only, and will be defunct at the next release.
#' 
#' @rdname DEGreport-deprecated
#' 
#' @details
#' The following functions are deprecated and will be made defunct; 
#' use the replacement indicated below:
#' * degRank, degPR, degBIcmd, degBI, degFC, degComb, degNcomb: DESeq2::lcfShrink. 
#' This function was trying to avoid big FoldChange
#' in variable genes. There are other methods nowadays like lcfShrink function.
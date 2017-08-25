#' Method to get the default table to use.
#' 
#' @param object [DEGSet]
#' @author Lorena Pantano
#' @rdname degDefault
#' @export
setMethod("degDefault", signature("DEGSet"), function(object) {
    slot(object, name = "default")
})

 
#' Method to get all table stored for an specific comparison
#'
#' @inheritParams degDefault
#' @param value Character to specify which table to use.
#' @author Lorena Pantano
#' @rdname degTable
#' @export
setMethod("degTable", signature("DEGSet"), function(object, value=NULL) {
    if (is.null(value))
        return(object[[degDefault(object)]])
    if (value == "raw")
        return(object[["raw"]])
    if (value == "shrunk")
        return(object[["shrunk"]])
    stop(value, "not found.")
})

#' MA-plot from base means and log fold changes
#'
#' MA-plot
#'
#' @author Victor Barrera
#' @author Rory Kirchner
#' @author Lorena Pantano
#' 
#' @param object [DESeqResults] or a data frame of DESeqResults.
#' @param title *Optional*. Plot title.
#' @param label_points Optionally label these particular points.
#' @param label_column Match label_points to this column in the results.
#' @param limit Absolute maximum to plot on the log2FoldChange.
#' @param diff Minimum difference between logFoldChange before and
#'   after shrinking.
#' @param ... Optional parameters to pass.
#' 
#' @docType methods
#' @rdname plotMA
#' @name plotMA
#' @aliases plotMA plotMA,DEGSet-method
#' 
#' @return MA-plot [ggplot].
#' @examples 
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(betaSD=1)
#' dds <- DESeq(dds)
#' res <- degComps(dds, contrast = list("condition_B_vs_A"))
#' plotMA(res[["condition_B_vs_A"]])
#' 
#' @export
setMethod("plotMA", signature(object = "DEGSet"), function(object, title = NULL,
                                                           label_points = NULL,
                                                           label_column = "symbol",
                                                           limit = NULL,
                                                           diff = 5, ...){
    .plotMA(object, title, label_points, label_column, limit, diff)
})
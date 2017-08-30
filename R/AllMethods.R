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
#' @param tidy Return data.frame, tibble or original class.
#' @param top Limit number of rows to return. Default: All.
#' @param ... Other parameters to pass for other methods.
#' @author Lorena Pantano
#' @rdname deg
#' @references  
#' * Testing it `top` is whole number or not comes from:
#'   https://stackoverflow.com/a/3477158
#' @export
setMethod("deg", signature("DEGSet"),
          function(object, value=NULL, tidy = NULL, top = NULL, ...) {
              df <- NULL

              if (is.null(value))
                  df <- object[[degDefault(object)]]
              else if (value == "raw")
                  df <- object[["raw"]]
              else if (value == "shrunk")
                  df <- (object[["shrunk"]])
              stopifnot(value %in% c(NULL,"raw", "shrunk"))
              
              if (is.null(top))
                  top <- nrow(df)
              stopifnot(top%%1 == 0)
              
              if (is.null(tidy))
                  return(df %>% .[order(.[["padj"]]),] %>%
                             .[1:top,])
              if (tidy == "data.frame")
                  return(as.data.frame(df) %>% .[order(.[["padj"]]),] %>%
                             .[1:top,])
              if (tidy == "tibble")
                  return(as.data.frame(df) %>%
                             rownames_to_column("gene") %>%
                             .[order(.[["padj"]]),] %>%
                             .[1:top,] %>%
                             as_tibble)
              stop("Not supported format, ", tidy)
          })

#' Method to get the significant genes
#' 
#' Function to get the features that are significant
#' according to some thresholds from a [DEGSet] class.
#' 
#' @author Lorena Pantano
#' @inheritParams degDefault
#' @param padj Cutoff for the FDR column.
#' @param fc Cutoff for the log2FC column.
#' @param ... Passed to [deg]. Default: value = NULL.
#' @rdname significants
#' @export
setMethod("significants", signature("DEGSet"),
          function(object, padj = 0.05, fc = 0, ...){
              df <-  as.data.frame(deg(object, ...))
              filterOut <- abs(df[["log2FoldChange"]]) > fc & df[["padj"]] < padj
              df %>%
                  rownames_to_column("gene") %>%
                  subset(., filterOut) %>%
                  .[order(abs(.[["log2FoldChange"]]), decreasing = TRUE),] %>%
                  .[["gene"]]
              
          })

#' MA-plot from base means and log fold changes
#'
#' MA-plot addaptation to show the shrinking effect.
#'
#' @author Victor Barrera
#' @author Rory Kirchner
#' @author Lorena Pantano
#' 
#' @param object [DEGSet] class.
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
#' @export
setMethod("plotMA", signature(object = "DEGSet"), function(object, title = NULL,
                                                           label_points = NULL,
                                                           label_column = "symbol",
                                                           limit = NULL,
                                                           diff = 5, ...){
    .plotMA(object, title, label_points, label_column, limit, diff)
})
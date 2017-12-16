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
              stopifnot(value %in% c(NULL, names(object)))
              
              if (is.null(value))
                  df <- object[[degDefault(object)]]
              else df <- object[[value]]
              
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
#' according to some thresholds from a [DEGSet],
#' [DESeq2::DESeqResults] and [edgeR::topTags].
#' 
#' @author Lorena Pantano
#' @inheritParams degDefault
#' @param padj Cutoff for the FDR column.
#' @param fc Cutoff for the log2FC column.
#' @param direction Whether to take down/up/ignore. Valid arguments are
#'   down, up and NULL.
#' @param ... Passed to [deg]. Default: value = NULL.
#'   Value can be 'raw', 'shrunken'.
#' @rdname significants
#' @export
setMethod("significants", signature("DEGSet"),
          function(object, padj = 0.05, fc = 0, direction = NULL, ...){
              fc <- abs(fc)
              df <-  as.data.frame(deg(object, ...))
              if (is.null(direction))
                  filterOut <- abs(df[["log2FoldChange"]]) > fc & df[["padj"]] < padj
              else if (direction == "up")
                  filterOut <- df[["log2FoldChange"]] > fc & df[["padj"]] < padj
              else if (direction == "down")
                  filterOut <- df[["log2FoldChange"]] < (-1L * fc) & df[["padj"]] < padj
              else
                  stop("Value ", direction, " is not valid: NULL, down, up.")
              
              df %>%
                  rownames_to_column("gene") %>%
                  subset(., filterOut) %>%
                  .[order(abs(.[["log2FoldChange"]]), decreasing = TRUE),] %>%
                  .[["gene"]]
          })

#' @rdname significants
#' @export
setMethod("significants", signature("DESeqResults"),
          function(object, padj = 0.05, fc = 0, direction = NULL, ...){
              fc <- abs(fc)
              df <-  as.data.frame(object)
              if (is.null(direction))
                  filterOut <- abs(df[["log2FoldChange"]]) > fc & df[["padj"]] < padj
              else if (direction == "up")
                  filterOut <- df[["log2FoldChange"]] > fc & df[["padj"]] < padj
              else if (direction == "down")
                  filterOut <- df[["log2FoldChange"]] < (-1L * fc) & df[["padj"]] < padj
              else
                  stop("Value ", direction, " is not valid: NULL, down, up.")
              
              df %>%
                  rownames_to_column("gene") %>%
                  subset(., filterOut) %>%
                  .[order(abs(.[["log2FoldChange"]]), decreasing = TRUE),] %>%
                  .[["gene"]]
          })

#' @rdname significants
#' @export
setMethod("significants", signature("TopTags"),
          function(object, padj = 0.05, fc = 0, direction = NULL, ...){
              fc <- abs(fc)
              df <-  as.data.frame(object)
              if (is.null(direction))
                filterOut <- abs(df[["logFC"]]) > fc & df[["FDR"]] < padj
              else if (direction == "up")
                  filterOut <- df[["logFC"]] > fc & df[["FDR"]] < padj
              else if (direction == "down")
                  filterOut <- df[["logFC"]] < (-1L * fc) & df[["FDR"]] < padj
              else
                  stop("Value ", direction, " is not valid: NULL, down, up.")
              
              df %>%
                  rownames_to_column("gene") %>%
                  subset(., filterOut) %>%
                  .[order(abs(.[["logFC"]]), decreasing = TRUE),] %>%
                  .[["gene"]]
          })

.supported <- function(x){
    return(class(x) %in% c("DEGSet", "DESeqResults"))
}

#' @rdname significants
#' @export
setMethod("significants", signature("list"),
          function(object, padj = 0.05, fc = 0, direction = NULL, ...){
              object <- object[sapply(object, .supported)]
              if (length(object) == 0){
                  message("Only DEGSet and DESeqResults objects are used.")
                  stop("No compatible objects remained.")
              }
              df <- lapply(object, deg, tidy = "tibble") %>%
                  bind_rows() %>% as.data.frame()
              df[["fdr"]] <- p.adjust(df[["pvalue"]], method = "fdr")
              fc <- abs(fc)
              if (is.null(direction))
                  filterOut <- abs(df[["log2FoldChange"]]) > fc & df[["fdr"]] < padj
              else if (direction == "up")
                  filterOut <- df[["log2FoldChange"]] > fc & df[["fdr"]] < padj
              else if (direction == "down")
                  filterOut <- df[["log2FoldChange"]] < (-1L * fc) & df[["fdr"]] < padj
              else
                  stop("Value ", direction, " is not valid, different from NULL, down, up.")
              
              df %>%
                  subset(., filterOut) %>%
                  .[order(abs(.[["log2FoldChange"]]), decreasing = TRUE),] %>%
                  .[["gene"]] %>% 
                  unique
              
          })

#' @rdname DEGSet
#' @export
setMethod("DEGSetFromEdgeR", signature("TopTags"),
          function(object, default = "shrunken", extras = NULL){
    name <- paste("contrast_", object[["comparison"]], collapse = "-")
    df <- as.data.frame(object) %>% set_colnames(c("log2FoldChange",
                                                "baseMean",
                                                "pvalue",
                                                "padj")) %>%
        DataFrame
    news <- lapply(extras, function(e){
        if (class(e) == "TopTags")
            return(as.data.frame(e) %>% set_colnames(c("log2FoldChange",
                                                       "baseMean",
                                                       "pvalue",
                                                       "padj")) %>%
                       DataFrame)
        stop(class(e), " is not a TopTags object.")
    })
    names(news) <- names(extras)
    l <- c(list(raw = df), news)
    dge <- new("DEGSet", l,
               default = default)
    attr(dge, "comparison") <- name
    dge
})

#' @rdname DEGSet
#' @export
setMethod("DEGSetFromDESeq2", signature("DESeqResults"),
          function(object, default = "shrunken", extras = NULL){
    name <- slot(slot(object, "elementMetadata"), "listData")[[2L]][2L]
    name <- strsplit(name, ":")[[1L]][2L] %>%
        gsub(" ", "", .)
    df <- DataFrame(object)
    news <- lapply(extras, function(e){
        if (class(e) == "DESeqResults")
            return(DataFrame(e))
        stop(class(e), " is not a DESeqResults object.")
    })
    names(news) <- names(extras)
    l <- c(list(raw = df), news)
    dge <- new("DEGSet", l,
               default = default)
    attr(dge, "comparison") <- name
    dge
})

# setMethod("DEGSetFromList", signature("list"),
#           function(object, default = "shrunken"){
#               dgeList <- lapply(object, function(o){
#                   if (class(o) == "DESeqResults")
#                       return(DEGSetFromDESeq2())
#                   else if (class(o) == "TopTags")
#                       return(DEGSetFromDESeq2())
#                   else
#                       message(class(o), " not supported, skipping.")
#               })
#               if (!is.null(names(list)))
#                 names(dgeList) <- names(list)
#               dgeList
# })


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
#' @param raw Whether to plot just the unshrunken log2FC.
#' @param correlation Whether to plot the correlation of the two logFCs.
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
                                                           diff = 5,
                                                           raw = FALSE,
                                                           correlation = FALSE,
                                                           ...){
    .plotMA(object, title, label_points, label_column,
            limit, diff, raw, correlation)
})
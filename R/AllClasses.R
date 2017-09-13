#' DEGSet
#' 
#' S4 class to store data from differentially expression analysis.
#' It should be compatible with different package and stores the information
#' in a way the methods will work with all of them.
#' 
#' @details
#' For now supporting only [DESeq2::results()] output.
#' Use constructor [degComps()] to create the object.
#' 
#' The list will contain one element for each comparison done.
#' Each element has the following structure:
#' 
#' * DEG table
#' * Optional table with shrunk Fold Change when it has been done.
#' 
#' To access the raw table use `deg(dgs, "raw")`, to access the 
#' shrunken table use `deg(dgs, "shrunken")` or just `deg(dgs)`. 
#' 
#' @rdname DEGSet
#' @name DEGSet
#' @param resList List with results as elements containing log2FoldChange,
#'   pvalues and padj as column. Rownames should be feature names. Elements
#'   should have names.
#' @param default The name of the element to use by default.
#' @param extras List of extra tables related to the same comparison.
#' @param object Different objects to be transformed to DEGSet.
#' @param ... Optional parameters of the generic.
#' @aliases DEGSet-class DEGSet
#' @author Lorena Pantano
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(betaSD = 1)
#' colData(dds)[["treatment"]] <- sample(colData(dds)[["condition"]], 12)
#' design(dds) <-  ~ condition + treatment
#' dds <- DESeq(dds)
#' res <- degComps(dds, combs = c("condition"))
#' deg(res[[1]])
#' deg(res[[1]], tidy = "tibble")
#' @export
DEGSet <- setClass("DEGSet",
                          contains = "list",
                          slots = c(default = "character"))

setValidity("DEGSet", function(object) {
    stopifnot(!is.null(names(object)))
    stopifnot(degDefault(object) %in% names(object))
    if (sum(c("raw", "shrunken") %in% names(object)) < 2)
        message("Some functions won't work without
                 'raw' and 'shrunken' elements in the object.")
    TRUE
})

#' @rdname DEGSet
#' @export
DEGSet <- function(resList, default){
    stopifnot(!is.null(names(resList)))
    stopifnot(default %in% names(resList))
    
    new("DEGSet", resList,
        default = default)
}
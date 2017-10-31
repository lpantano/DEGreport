.get_counts <- function(c){
    if (class(c) == "bcbioRNASeq")
        return(counts(c, "rlog"))
    if (class(c) == "DESeqDataSet")
        return(counts(c, normalized = TRUE))
    if (class(c) == "SummarizedExperiment")
        return(assay(c))
    if (class(c) == "data.frame")
        return(c)
    if (class(c) == "matrix")
        return(as.data.frame(c))
    stop("class ", class(c), " no supported.")
}

.get_meta <- function(c){
    if (class(c) %in% c("bcbioRNASeq", "DESeqDataSet", "SummarizedExperiment"))
        return(data.frame(colData(c), stringsAsFactors = FALSE))
    return(NULL)
}

#' Plot gene signature for each group and signature
#' 
#' Given a list of genes beloging to a different classes, like
#' markers, plot for each group, the expression values for all the samples.
#' 
#' @param counts expression data. It accepts bcbioRNASeq, DESeqDataSet and
#'   SummarizedExperiment. As well, data.frame or matrix is supported, but
#'   it requires metadata in that case.
#' @param signature data.frame with two columns: a) genes that match
#'   row.names of counts, b) label to classify the gene inside a group.
#'   Normally, cell tissue name.
#' @param group character in metadata used to split data into different
#'   groups.
#' @param metadata data frame with sample information. Rownames
#'   should match \code{ma} column names
#'   row number should be the same length than p-values vector.
#' @return ggplot plot.
#' @examples
#' data(humanGender)
#' data(geneInfo)
#' library(SummarizedExperiment)
#' degSignature(humanGender, geneInfo, group = "group")
#' @export
degSignature <- function(counts, signature, group = NULL, metadata = NULL){
    c <- .get_counts(counts)
    meta <- .get_meta(counts)
    if (is.null(meta))
        meta <- metadata
    stopifnot(group %in% colnames(meta))
    
    names(signature) <- c("gene", "signature")
    common <- intersect(row.names(c), signature[["gene"]])
    c[common, ] %>%  melt() %>%  data.frame(., stringsAsFactors = FALSE) %>% 
        set_colnames(c("gene", "sample", "expression")) %>% 
        mutate_if(is.factor, as.character) %>% 
        left_join(meta %>% rownames_to_column("sample") %>% 
                      mutate_if(is.factor, as.character),
                  by = "sample") %>%
        left_join(signature %>% 
                      mutate_if(is.factor, as.character), by = "gene") %>% 
        group_by(group, signature, sample) %>% 
        summarise(median = median(expression)) %>% 
        ungroup %>% 
        ggplot(aes_string(x = group, y = "median", color = "signature")) +
        geom_boxplot() +
        geom_jitter() +
        facet_wrap(~signature)
}

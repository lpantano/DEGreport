#' Plot top genes allowing more variables to color and shape points
#'
#' @param dds [DESeq2::DESeqDataSet] object or SummarizedExperiment
#'   or Matrix or data.frame.
#' @param res [DESeq2::DESeqResults] object.
#' @param n Integer number of genes to plot.
#' @param genes Character of gene names matching rownames of count data.
#' @param xs Character, colname in colData that will be used as X-axes.
#' @param group Character, colname in colData to color points and add different
#'   lines for each level.
#' @param batch Character, colname in colData to shape points, normally used by
#'   batch effect visualization.
#' @param ann Columns in rowData (if available) used to print gene names.
#' @param metadata Metadata in case dds is a matrix.
#' @param slot Name of the slot to use to get count data.
#' @param log2 Whether to apply or not log2 transformation.
#' @param xsLab Character, alternative label for x-axis (default: same as xs).
#' @param groupLab Character, alternative label for group (default: same as group).
#' @param batchLab Character, alternative label for batch (default: same as batch).
#' @return ggplot showing the expresison of the genes
#' @examples 
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dse <- DESeq(dse)
#' degPlot(dse, genes = rownames(dse)[1:10], xs = "group")
#' @export
degPlot = function(dds, xs, res=NULL, n=9, genes=NULL,
                   group=NULL, batch=NULL,
                   metadata = NULL,
                   ann=c("external_gene_name", "symbol"),
                   slot=1L,
                   log2 = TRUE,
                   xsLab=xs, groupLab=group, batchLab=batch){
    if (class(dds) %in% c("data.frame", "matrix"))
        dds = SummarizedExperiment(assays = SimpleList(counts = as.matrix(dds)),
                                   colData = metadata)
    
    if ( !("assays" %in% slotNames(dds)) )
        stop("dds object doesn't have assays slot")
    
    if (!is.null(res))
        res <- res[order(res$padj),] %>% .[!is.na(res$padj),]
    
    if (is.null(genes))
        genes =  row.names(res)[1L:n]
    
    anno <- as.data.frame(rowData(dds))
    
    metadata = data.frame(colData(dds))
    if (class(dds) == "DESeqDataSet")
        counts <- counts(dds, normalized = TRUE)
    else counts <- assays(dds)[[slot]]
    
    stopifnot(class(counts) == "data.frame" | class(counts) == "matrix")
    
    if (log2 & max(counts) < 500L)
        warning("Data seems to be already in log2. Please use log2 = FALSE.")
    if (log2)
        counts <- log2(counts + 0.2)
    
    newgenes <- genes
    if (ncol(anno) > 0) {
        name <- intersect(names(anno), ann)
        if (length(name) != 2L)
            message("No genes were mapped to rowData. check ann parameter values.")
        if (length(name) == 2L)
            newgenes <- anno[match(genes, anno[, name[1L]]), name[2L], drop = FALSE]
    }
    
    dd = melt(as.data.frame(counts[genes, , drop = FALSE]) %>%
                  mutate(gene = newgenes))
    colnames(dd) = c("gene", "sample", "count")
    
    dd$xs = as.factor(metadata[as.character(dd$sample), xs])
    
    if (!is.null(group)) {
        dd[, groupLab] = as.factor(metadata[as.character(dd$sample), group])
    }else {
        groupLab = "fake"
        dd[, groupLab] = "fake"
    }
    
    if (!is.null(batch)) {
        dd[, batchLab] = as.factor(metadata[as.character(dd$sample), batch])
        
        p = ggplot(dd, aes_string(x = "xs", y = "count", color = groupLab,
                                  shape = batchLab))
    }else{
        p = ggplot(dd, aes_string(x = "xs", y = "count", color = groupLab))
    }
    
    p = p +
        # geom_violin(alpha=0.3) +
        stat_smooth(fill = "grey80", method = 'loess') +
        geom_point(size = 1, alpha = 0.7,
                   position = position_jitterdodge(dodge.width = 0.9)) +
        facet_wrap(~gene, scales = "free_y") +
        xlab(xsLab)
    if (length(unique(dd[, groupLab])) == 1L) {
        p = p +
            scale_color_manual(guide = FALSE, values = "black") +
            scale_fill_manual(guide = FALSE, values = "black")
    }else{
        p = p +
            scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")
    }
    p = p + theme_bw() +
        theme(strip.background = element_rect(fill = "white"),
              strip.text = element_text(colour = "black"))
    
    suppressWarnings(p)
}

#' Plot selected genes on a wide format
#'
#' @param counts [DESeq2::DESeqDataSet] object or expression matrix
#' @param genes character genes to plot.
#' @param group character, colname in colData to color points and add different
#' lines for each level
#' @param metadata data.frame, information for each sample. Not needed if
#' [DESeq2::DESeqDataSet] given as counts.
#' @param batch character, colname in colData to shape points, normally used by
#' batch effect visualization
#' @return ggplot showing the expresison of the genes on the x
#' axis
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dse <- DESeq(dse)
#' degPlotWide(dse, rownames(dse)[1:10], group = "group")
#' @export
degPlotWide <- function(counts, genes, group, metadata=NULL, batch=NULL){
    if (is.null(metadata))
        metadata = data.frame(colData(counts))
    metadata = data.frame(metadata)
    if (class(counts) == "DESeqDataSet") {
        dd = bind_rows(lapply(genes,function(gene){
            plotCounts(counts, gene,
                       intgroup = group, returnData = TRUE) %>%
                mutate(count = log2(count + 1)) %>%
                mutate(gene = gene, sample = row.names(metadata))}))
    }else if (class(counts) == "matrix") {
        dd = melt(counts[genes, ])
        colnames(dd) = c("gene", "sample", "count")
        dd$group = as.factor(metadata[as.character(dd$sample), group])
    }else{
        stop("No supported for class", class(counts))
    }
    if (is.null(group)) {
        dd$treatment = "one_group"
    }else{
        dd$treatment = dd[,"group"]
    }
    p = ggplot(dd, aes_string(x = "gene", y = "count", color = "treatment"))
    if (!is.null(batch)) {
        dd$batch = as.factor(metadata[dd$sample, batch])
        p = ggplot(dd, aes_string(x = "gene", y = "count",
                                  color = "treatment", shape = "batch"))
    }
    
    p = p +
        geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
        xlab("Genes") +
        ylab("Normalized Counts") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45L, hjust = 1L))
    p
}

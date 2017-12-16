
fmt <- function(){
    function(x) format(x,nsmall = 1,scientific = FALSE)
}
#' Create volcano plot from log2FC and adjusted pvalues data frame
#'
#' @param stats data.frame with two columns: logFC and Adjusted.Pvalue
#' @param side plot UP, DOWN or BOTH de-regulated points
#' @param title title for the figure
#' @param pval.cutoff cutoff for the adjusted pvalue. Default 0.05
#' @param lfc.cutoff cutoff for the log2FC. Default 1
#' @param shade.alpha transparency value. Default 0.25
#' @param shade.colour background color. Default orange.
#' @param point.colour colours for points. Default gray
#' @param point.alpha transparency for points. Default 0.75
#' @param point.outline.colour Default darkgray
#' @param line.colour Defaul gray
#' @param plot_text data.frame with three columns: logFC, Pvalue, Gene name
#' @return The function will plot volcano plot together
#' with density of the fold change and p-values on the top and the
#' right side of the volcano plot.
#' @details 
#' This function was mainly developed by @jnhutchinson.
#' @author Lorena Pantano, John Hutchinson
#' @examples 
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(betaSD = 1)
#' dds <- DESeq(dds)
#' stats <- results(dds)[,c("log2FoldChange", "padj")]
#' stats[["name"]] <- row.names(stats)
#' degVolcano(stats, plot_text = stats[1:10,])
#' @export
degVolcano <- function(stats, side="both", title="Volcano Plot with Marginal Distributions",
                                 pval.cutoff=0.05, lfc.cutoff=1, shade.colour="orange",
                                 shade.alpha=0.25, point.colour="gray", point.alpha=0.75,
                                 point.outline.colour="darkgray", line.colour="gray",
                                 plot_text=NULL) {
    if (class(stats) == "DESeqResults")
        stats <- stats[, c("log2FoldChange", "padj")]
    if (class(stats) == "DEGSet"){
        stats <- deg(stats, tidy = "data.frame")
        stats <- stats[, c("log2FoldChange", "padj")]
    }
    stats <- as.data.frame(stats)
    if (!any(side %in% c("both","down","up")) | length(side)>1)
        stop("side parameter should be: both, up or down.")
    if (ncol(stats)<2)
        stop("Need a data.frame with two columns: logFC and Adjusted.Pvalue")
    if (sum( rowSums(is.na(stats)) ) > 0)
        stats = stats[rowSums(is.na(stats)) == 0,]
    if (any(stats[,2]>1) | any(stats[,2]<0))
        stop("pvalues needs to be >0 and <1")
    names(stats) = c("logFC","adj.P.Val")
    stats[,2] = -log10(stats[,2] + 1e-10)
    
    # get range of log fold change and p-value values to setup plot borders
    range.lfc <- c(floor(min(stats$logFC)), ceiling(max(stats$logFC)))
    range.pval <- c(floor(min(stats$adj.P.Val)), ceiling(max(stats$adj.P.Val)))
    pval.cutoff <- -log10(pval.cutoff)

    #make scatter volcano plot
    scat.poly.up <- with(stats, data.frame(x=as.numeric(c(lfc.cutoff,  lfc.cutoff, max(range.lfc),max(range.lfc))), 
                                           y=as.numeric(c(pval.cutoff, max(range.pval), max(range.pval),pval.cutoff))))
    scat.poly.down <- with(stats, data.frame(x=as.numeric(c(-lfc.cutoff,  -lfc.cutoff, min(range.lfc),min(range.lfc))), 
                                             y=as.numeric(c(pval.cutoff, max(range.pval), max(range.pval),pval.cutoff))))
    
    scatter <- ggplot(stats, aes_string(x="logFC", y="adj.P.Val")) +
        geom_point(alpha=point.alpha, pch=21, fill=point.colour, color=point.outline.colour) +
        xlab("log2 fold change") + ylab("-log10(adjusted p-value)") +
        theme_bw()+
        theme(legend.position="none") +
        theme(plot.margin=unit(c(3,-5.5,4,3), "mm") )+
        scale_x_continuous(limits = range.lfc, breaks = range.lfc[1]:range.lfc[2], expand = c(.05,.05))+
        scale_y_continuous(labels=fmt(), limits = range.pval)+ labs(list(title="Volcano plot"))
    if (side=="both" | side=="up")
        scatter = scatter + geom_polygon(data=scat.poly.up, aes_string(x="x",y="y"), fill=shade.colour, alpha=shade.alpha)
    if (side=="both" | side=="down")
        scatter = scatter + geom_polygon(data=scat.poly.down, aes_string(x="x",y="y"), fill=shade.colour, alpha=shade.alpha)

    if (!is.null(plot_text)){
        plot_text <- as.data.frame(plot_text)
        names(plot_text) = c("logFC", "adj.P.Val", "name")
        plot_text[,2] <- -log10(plot_text[,2] + 1e-10)
        scatter <- scatter + 
            geom_text_repel(data=plot_text, aes_string(x="logFC", y="adj.P.Val", label="name"), size=3)
    }
    scatter
}

# plot group of genes according time and group
.plot_cluster  = function(norm_sign, g_in_c, xs ,groups, title) {
    ma = as.data.frame(norm_sign)[g_in_c,]
    ma_long = suppressMessages(melt(cbind(gene=row.names(ma), ma), variable_name = "sample"))
    ma_long$x = xs[ma_long$sample]
    ma_long$group = groups[ma_long$sample]
    ggplot(ma_long, aes(x=x, y=value, fill=group, color=group)) + 
        geom_violin(alpha=0.3) + 
        stat_smooth(aes(x=x, y=value, group=group, color=group),method = "lm",formula = y~poly(x,2)) +
        ggtitle(paste("Group:", title, "(", length(g_in_c), " genes )")) +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Scale from 1 to 0 the expression genes
.scale <- function(e){
    (e - min(e))/(max(e) - min(e))
}

# use diana package to detect clusters
.make_clusters <- function(counts_group){
    m = (1-cor(t(counts_group), method = "kendall"))
    d = as.dist(m^2)
    c = diana(d, diss = TRUE, stand = FALSE)
    
    cutree(as.hclust(c), h = c$dc)
}

#' Make groups of genes using expression profile
#' @aliases degPatterns
#' @param ma  log2 normalized count matrix
#' @param metadata  data frame with sample information. Rownames
#' should match \code{ma} column names
#' row number should be the same length than pvalues vector.
#' @param minc integer minimum number of genes in a group that
#' will be return
#' @param summarize column name in metadata that group replicates.
#' For instance, a merge between summarize and time parameters:
#' control_point0 ... etc
#' @param time column name in metadata that will be used as
#' variable that changes, normally time
#' @param col column name in metadata to separate
#' samples. Normally control/mutant
#' @return ggplot2 object
degPatterns = function(ma, metadata, minc=15, summarize="group", time="time", col="condition"){
    
    if (!col %in% names(metadata))
        metadata[,col] = paste0(metadata[,summarize], metadata[,time])
    stopifnot(class(metadata) == "data.frame")
    stopifnot(class(ma) == "matrix")
    stopifnot(summarize %in% names(metadata))
    stopifnot(time %in% names(metadata))
    if (nrow(ma)>3000)
        message("Large number of genes given. Please,",
                "make sure is not an error. Normally",
                "Only DE genes are useful for this function.")
    
    counts_group = t(sapply(rownames(ma), function(g){
        sapply(unique(metadata[,summarize]), function(i){
            idx = which(metadata[,summarize] == i)
            mean(ma[g,idx], na.rm=TRUE)
        })
    }))
    
    groups = .make_clusters(counts_group)
    
    norm_sign = t(apply(counts_group, 1, .scale))
    metadata_groups = metadata %>% distinct(metadata[,summarize])
    rownames(metadata_groups) = metadata_groups[,summarize]

    to_plot = names(table(groups))[table(groups) > minc]
    plots = lapply(to_plot, function(x){
        .plot_cluster(norm_sign, as.character(names(groups[groups==x])), 
                     metadata_groups[,time], metadata_groups[,col], x)
    })
    
    do.call(grid.arrange, plots)
    data.frame(genes=names(groups),clutser=groups)
}

# output plots as a way to avoid running the function again
# add plotGenes to documentation. Option to plot by names

.logger = function(toout,msg=""){
    logdebug(paste("\n\nchecking data:" , msg, "\n\n"))
    if (getLogger()[['level']]!=20)
        print(toout)
}

# It takes a table with gene and metadata information and convert
# it to be compatible with degPlotCluster
.process <- function(table, time, color){
    stopifnot(c("genes", "sample", time, "cluster", "expression")  %in% names(table))

        if (!is.null(color))
        stopifnot(color  %in% names(table))
    if (is.null(color)){
        grouped <- group_by(table, genes,
                            !!sym(time), cluster)
    }else{
        grouped <- group_by(table, genes,
                            !!sym(time), cluster, !!sym(color))
    }
    
    grouped %>% 
        summarise(expression=mean(expression)) %>% 
        group_by(!!sym("genes")) %>% 
        mutate(value = scale(expression)) %>% 
        ungroup()
}

#' Plot clusters from degPattern function output
#' 
#' This function helps to format the cluster plots from [degPatterns()].
#' It allows to control the layers and it returns a ggplot object that
#' can accept more ggplot functions to allow customization.
#' 
#' @param table `normalized` element from [degPatterns()] output.
#'   It can be a data.frame with the following columns in there:
#'   `genes, sample, expression, cluster, xaxis_column, color_column`.
#' @param time column name to use in the x-axis.
#' @param process whether to process the table if it is not
#'   ready for plotting.
#' @param color column name to use to color and divide the samples.
#' @param min_genes minimum number of genes to be added to the plot.
#' @param points Add points to the plot.
#' @param boxes Add boxplot to the plot.
#' @param smooth Add regression line to the plot.
#' @param lines Add gene lines to the plot.
#' @param facet Split figures based on cluster ID.
#' @param cluster_column column name if cluster is in a column
#'   with a different name. Usefull, to plot cluster with different
#'   cutoffs used when grouping genes from the clustering step.
#' @param prefix_title text to add before the cluster ID in the figure title.
#' @return [ggplot2] object.
#' @examples
#' data(humanGender)
#' library(SummarizedExperiment)
#' library(ggplot2)
#' ma <- assays(humanGender)[[1]][1:100,]
#' des <- colData(humanGender)
#' des[["other"]] <- sample(c("a", "b"), 85, replace = TRUE)
#' res <- degPatterns(ma, des, time="group", col = "other", plot = FALSE)
#' degPlotCluster(res$normalized, "group", "other")
#' degPlotCluster(res$normalized, "group", "other", lines = FALSE)
#' 
#' library(dplyr)
#' library(tidyr)
#' library(tibble)
#' table <- rownames_to_column(as.data.frame(ma), "genes") %>%
#'     gather("sample", "expression", -genes) %>%
#'     right_join(distinct(res$df[,c("genes", "cluster")]),
#'                by = "genes") %>%
#'     left_join(rownames_to_column(as.data.frame(des), "sample"),
#'               by = "sample") %>% 
#'               as.data.frame()
#' degPlotCluster(table, "group", "other", process = TRUE)
#' @export
degPlotCluster <- function(table, time, color = NULL,
                           min_genes = 10,
                           process = FALSE,
                           points = TRUE,
                           boxes = TRUE,
                           smooth = TRUE,
                           lines = TRUE,
                           facet = TRUE,
                           cluster_column = "cluster",
                           prefix_title = "Group: "){
    stopifnot(class(table)[1] == "data.frame")

    if (cluster_column  %in% colnames(table)){
        table[["cluster"]] = table[[cluster_column]]
    }
    if (process){
        table <- .process(table, time, color)
    }

    if ("cluster"  %in% colnames(table)){
        counts <- table(distinct(table, genes, cluster)[["cluster"]])
        counts <- counts[counts>=min_genes]
        if (length(counts)==0)
            stop("No clusters with min_genes > ", min_genes)
        table <- inner_join(table,
                            data.frame(cluster = as.integer(names(counts)),
                                       title = paste(prefix_title,
                                                      names(counts),
                                                      "- genes:" ,
                                                      counts),
                                       stringsAsFactors = FALSE),
                            by = "cluster")
    }

    if (is.null(color)){
        color = "dummy"
        table[[color]] = ""
        lines = FALSE
    }
    table[["line_group"]] = paste(table[["genes"]],
                                  table[[color]])
    
    
    splan <- length(unique(table[[time]])) - 1L
    table$title <- table$title %>% as.factor()
    old <- table$title %>% levels() %>% stringi::stri_extract_first(., regex="\\d+") %>% as.integer()
    index <- sort(old, index.return = TRUE)[[2]]
    table$title <- factor(table$title, levels = levels(table$title)[index])
    p <- ggplot(table, aes_string(x = time, y = "value",
                                  fill = color, color = color))
    
    if (boxes)
        p <- p + geom_boxplot(alpha = 0,
                              outlier.size = 0,
                              outlier.shape = NA, )
    if (points)
        p <- p + 
        geom_point(alpha = 0.4, size = 1,
                   position = position_jitterdodge(dodge.width = 0.9))
    if (smooth)
        p <- p + 
        stat_smooth(aes_string(x = time, y = "value",
                               group = color, color = color),
                    se = FALSE,
                    method = "lm", formula = y~poly(x, splan))
    if (lines)
        p <- p + geom_line(aes_string(group = "line_group"), alpha = 0.1)
    if (facet)
        p <- p + facet_wrap(~title)
    
    p <- p + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ylab("Z-score of gene abundance") +
        xlab("")
    p + theme_bw() 
    
}

.plot_benchmarking_curve <- function(benchmarking){
    if (is.null(benchmarking))
        return(NULL)
    df <- benchmarking[["genes"]]
    nc <- apply(df[2:ncol(df)], 2, function(c){
        length(unique(c))
    })
    ng <- apply(df[2:ncol(df)], 2, function(c){
        sum(!is.na(c))
    })
    data.frame(cutoff = names(nc),
               cluster = rev(nc),
               genes = rev(ng),
               pct_variance = unlist(benchmarking[["pcts"]])) %>% 
        ggplot(aes_string("cutoff", "pct_variance", group = 1)) +
        geom_line() +
        geom_text(aes_string(label="nc"), nudge_y = 1) +
        geom_text(aes_string(label="ng"), nudge_y = -1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust=1))
}

.plot_benchmarking <- function(normalized, benchmarking,
                               time, color){

    if (is.null(benchmarking))
        return(NULL)
    benchmarking[["pcts"]] <- benchmarking[["pcts"]][lapply(benchmarking[["pcts"]],length)>0] 
    p <- lapply(names(benchmarking[["pcts"]]), function(serie){
        table <- normalized %>% group_by(!!sym(color),
                                         !!sym(time),
                                         !!sym(serie)) %>% 
            summarise(value = median(value))
        n = length(unique(table[[serie]]))
        p <- ggplot(table, aes_string(x = time, y = "value",
                                 color = color, group = serie)) +
            geom_line() +
            guides(color="none") +
            ggtitle(paste0(serie, ": ",
                          round(benchmarking[["pcts"]][[serie]]),
                          "% clusters:",
                          n)) +
            theme_minimal() +
            theme(title = element_text(size=6),
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
        
        if (color == "colored"){
            p <- p + scale_color_manual(values = "black")
        }
        p
    }) %>% plot_grid(plotlist = .)
    p
}

# plot group of genes according time and group
.plot_cluster  <- function(norm_sign, g_in_c, xs , groups, title, fixy=NULL) {
    ma <- as.data.frame(norm_sign)[g_in_c,]
    ma_long <- suppressMessages(melt(cbind(genes = row.names(ma), ma),
                                     variable_name = "sample"))
    ma_long$x <- xs[ma_long$sample]
    ma_long$group <- groups[ma_long$sample]
    plotting_data <- ma_long
    p <- suppressWarnings(degPlotCluster(ma_long, "x", "group", facet = FALSE))
    p <- p +
        ggtitle(paste("Group:", title, "(", length(g_in_c), " genes )"))
    if (!is.null(fixy))
        p <- p + ylim(fixy[1], fixy[2])
    if (length(unique(groups)) == 1) {
        p <- p + scale_color_brewer(guide = "none", palette = "Set2") +
            scale_fill_brewer(guide = "none", palette = "Set2")
    }else{
        p <- p + scale_color_brewer(palette = "Set2") +
            scale_fill_brewer(palette = "Set2")
    }
    p
}

# Scale from 1 to 0 the expression genes
.scale <- function(e){
    scale(e)
}

.pct_var <- function(counts, clusters){
    clusters_sd <- lapply(unique(clusters), function(c){
        counts[names(clusters[clusters == c]),, drop = FALSE] %>% 
            apply(., 2, sd) %>% 
            median(.)
    }) %>% unlist() %>%  median
    (1 - clusters_sd / median(apply(counts, 2, sd))) * 100
}

.reduce <- function(groups, counts_group){
    lapply(unique(groups), function(g){
        ma <- counts_group[names(groups)[groups==g],]
        nokeep <- apply(ma, 2L, function(x){
            out <- boxplot(x, plot = FALSE)$out
            rownames(ma)[which(x %in% out)]
        }) %>% as.vector() %>% unlist() %>% unique()
        keep <- setdiff(rownames(ma), nokeep)
        group <- rep(g, length(keep))
        names(group) <- keep
        group
    }) %>%  unlist()
}

.select_pattern <- function(df){
    if(!("bycor" %in% colnames(df)))
        selected = df[df[["fdr"]] < 0.05, "gene"]
    selected = df[df[["cor"]] > 0.7, "gene"]
    selected = selected[!is.na(selected)]
    groups = rep(1, length(selected))
    names(groups) = selected
    groups
}

# find correlation to a pattern
.find_pattern <- function(counts_group, reference){
    cor_ma = apply(counts_group, 1, function(x){
        if(length(x) > 2){
            c = cor.test(x, reference, method = "kendall")
            return(data.frame(cor = c$estimate, pval=c$p.value))
        }else{
            c = cor(x, reference, method = "kendall")
            return(data.frame(cor = c[1], pval=0, bycor = 1))
        }
    }) %>%  bind_rows()
    cor_ma[["fdr"]] = p.adjust(cor_ma[["pval"]], "fdr") 
    cor_ma[["gene"]] = rownames(counts_group)
    cor_ma
}

# use diana package to detect clusters
.make_clusters <- function(counts_group){
    m <- (1 - cor(t(counts_group), method = "kendall"))
    d <- as.dist(m^2)
    c <- diana(d, diss = TRUE, stand = FALSE)
    c
}

.make_concensus_cluster <- function(counts_group, ...){
    ConsensusClusterPlus(t(counts_group),
                         reps = 500, maxK = 10,
                         distance = "pearson", finalLinkage = "ward.D2",
                         ...)
}

.select_concensus_genes <- function(c){
    icl <- calcICL(c)
    consensus <- icl[["itemConsensus"]] %>%
        group_by(!!!sym("k")) %>% 
        summarise(score = quantile(!!!sym("itemConsensus"),
                                   0.75,
                                   na.rm = TRUE)) %>% 
        arrange(desc(score), k) %>% .[["k"]] %>% .[1]
    
    df <- icl[["itemConsensus"]] %>%
        filter(k == consensus, itemConsensus > 0.9) %>% 
        .[,c("item", "cluster")]
    
    genes <- df[["cluster"]]
    names(genes) <- df[["item"]]
    
    genes
}

.select_genes <- function(c, counts_group, minc=15,
                          reduce=FALSE, cutoff=0.30, nClusters = NULL){
    if (is.null(nClusters)){
        h <- c$dc
        select <- cutree(as.hclust(c), h = h)
        select <- select[select %in% names(table(select))[table(select) > minc]]}
    else {
        select <- cutree(as.hclust(c), k = nClusters)
        select <- select[select %in% names(table(select))[table(select) > minc]]}
    
    if (reduce)
        select <- .reduce(select, counts_group)
        select <- select[select %in% names(table(select))[table(select) > minc]]
        message("Working with ", length(select), " genes after filtering: minc > ",minc)
    return(select)
}

.benckmark_cutoff <- function(tree, counts, minc = 15){
    series <- unique(round(tree$height, digits = 3))
    series <- sort(series[2:length(series)])
    list_clusters <- lapply(series, function(s) {
        select <- cutree(as.hclust(tree), h = s)
        select <- select[select %in% names(table(select))[table(select) > minc]]
    })
    list_pct <- lapply(list_clusters, function(c) {
        .pct_var(counts, c)
    })
    names(list_clusters) <- paste0("cutoff", series)
    names(list_pct) <- paste0("cutoff", series)
    df_clusters <- lapply(names(list_clusters), function(s){
        c = list_clusters[[s]]
        if (length(c) == 0) {return(data.frame())}
        data.frame(genes = make.names(names(c)),
                   cluster = c,
                   cutoff = s,
                   stringsAsFactors = F)
    }) %>% bind_rows() %>% 
        distinct() %>% 
        spread(cutoff, cluster)
    list(genes=df_clusters, pcts=list_pct)
}

# summarize matrix into groups and scale if needed
.summarize_scale <- function(ma, group, scale = TRUE){
    counts_group = t(sapply(rownames(ma), function(g){
        sapply(levels(group), function(i){
            idx = which(group == i)
            mean(ma[g, idx], na.rm = TRUE)
        })
    }))
    colnames(counts_group) = levels(group)
    .logger(head(counts_group), "summarize_scale::counts_group")
    .logger(head(group), "summarize_scale::group")
    if (scale) {
        norm_sign <- t(apply(counts_group, 1, .scale))
    }else{
        norm_sign <- counts_group
    }
    colnames(norm_sign) = colnames(counts_group)
    norm_sign
}

.median_per_cluster <- function(ma, clusters){
    # matrix, and df
    .logger(table(clusters[["cluster"]]),
            ".median_per_cluster::table=clusters")
    .logger(head(clusters), ".median_per_cluster::clusters")
    .t = do.call(rbind, lapply(unique(clusters$cluster), function(nc){
        .g = as.character(clusters$genes[clusters$cluster == nc])
        .e = apply(ma[.g, ], 2, median.default, na.rm = TRUE)
        .e
    }))
    rownames(.t) = unique(clusters$cluster)
    colnames(.t) = colnames(ma)
    .t
}

.filter <- function(df, co = 4){
    .sum <- table(df[["cluster"]])
    .pass <- names(.sum)[.sum > co]
    df[df[["cluster"]] %in% .pass,]
}

.get_features <- function(genes, norm, mapping){
    if (is.null(mapping))
        return(data.frame(touse = intersect(genes, rownames(norm)), pair = "."))
    df <- mapping[match(genes, mapping[, 1]),]
    names(df) = c("pair","touse")
    df
}

.plot_base <- function(genes, norm, raw, metadata, time, col, nc, name){
    .logger(length(intersect(genes, rownames(norm))), "num genes")
    .logger(head(norm), "Plot expression")
    .logger(head(metadata), "Plot expression")
    p <- .plot_cluster(norm, genes,
                       metadata[,time],
                       metadata[,col], nc) +
        ggtitle( paste("Pattern of the Genes in", name) )
    print(p)
    .logger(head(melt(as.data.frame(raw[genes,]))), "Plot expression")
    p <- ggplot(melt(as.data.frame(raw[genes,])),
               aes_string(x = "value", color = "variable")) + geom_density() +
        ggtitle( paste("Expression of the Genes in", name) )
    print(p)

}

.integrate <- function(nc1, .exp_base, .clus_base, .norm_base, .metadata_group,
                       .maraw_base, .ma, .norm, .metadata, time, col, summarize,
                       col_transformed, name, mapping=NULL){
    .genes_nc1 <- as.character(.clus_base$genes[.clus_base$cluster == nc1])
    keep_info <-  .get_features(.genes_nc1, .norm, mapping)
    keep <- as.character(unique(keep_info$touse))
    if (length(keep) < 5)
        return(NULL)
    clusters <- degPatterns(as.matrix(.ma[keep,]),
                           .metadata,
                           minc = 5, summarize = summarize,
                           time = time, col = col)
    if (length(clusters$pass) == 0)
        return(NULL)
    .exp <- .median_per_cluster(.norm,
                               clusters$df %>%
                               .[clusters$df$cluster %in% clusters$pass,])
    df <- lapply(rownames(.exp), function(nc2){
        .genes_nc2 = as.character(clusters$df$genes[clusters$df$cluster==nc2])
        .est = cor.test(.exp_base[nc1,], .exp[nc2,colnames(.exp_base)])
        data.frame(nc1=nc1,
                   nc2=paste0(nc1, ".", nc2),
                   cor=.est$estimate,
                   pval=.est$p.value,
                   gene=keep_info$touse[keep_info$touse %in% .genes_nc2],
                   pair=keep_info$pair[keep_info$touse %in% .genes_nc2])
    })
    .logger(head(df[[1]]), ".integrate::final_df")
    do.call(rbind, df)
}

.group_metadata <- function(metadata, time, col, summarize){
    if (is.null(col)) {
        col_transformed = "condition"
        metadata[,col_transformed] = rep("one_group", nrow(metadata))
    }
    if (!summarize %in% names(metadata))
        metadata[,summarize] = paste0(metadata[,col_transformed], metadata[,time])

    metadata_groups = metadata %>% dplyr::distinct(!!sym(summarize), .keep_all = TRUE)
    rownames(metadata_groups) = metadata_groups[,summarize]
    return(metadata_groups)
}


#' Integrate data comming from degPattern into one data object
#'
#' The simplest case is if you want to convine the pattern profile
#' for gene expression data and proteomic data. It will use the first element
#' as the base for the integration. Then, it will loop through clusters
#' and run [degPatterns] in the second data set to detect patterns that match
#' this one.
#' @param matrix_list list expression data for each element
#' @param cluster_list list df item from degPattern output
#' @param metadata_list list data.frames from each element
#' with design experiment. Normally \code{colData} output
#' @param summarize character column to use to group samples
#' @param time character column to use as x-axes in figures
#' @param col character column to color samples in figures
#' @param scale boolean scale by row expression matrix
#' @param mapping data.frame mapping table in case elements use
#' different ID in the row.names of expression matrix. For instance,
#' when integrating miRNA/mRNA.
#' @return A data.frame with information on what genes are in each cluster in
#' all data set, and the correlation value for each pair cluster comparison.
degMerge <- function(matrix_list, cluster_list, metadata_list,
                     summarize="group", time="time", col="condition",
                     scale=TRUE, mapping=NULL){
    # basicConfig(level='FINEST')
    stopifnot(length(matrix_list) > 1 | length(metadata_list) > 1)
    stopifnot(length(matrix_list) == length(metadata_list))
    matrixnorm_list = list()
    cluster_expression = list()
    metadata_groups_list = list()
    matrixabs_list = list()
    for (name in names(matrix_list)) {
        .logger(name, msg = "Name list")
        metadata = metadata_list[[name]]
        matrixnorm_list[[name]] = .summarize_scale( as.matrix(matrix_list[[name]]),
                          group = metadata[, summarize],
                          scale = scale)
        .logger(head(matrix_list[[name]]), msg = "cluster Matrix")
        matrixabs_list[[name]] = .summarize_scale(as.matrix(matrix_list[[name]]),
                                                    group = metadata[,summarize],
                                                    scale = FALSE)
        .logger(head(matrixnorm_list[[name]]), msg = "Matrix norm DF")

        if (is.null(col))
            col_transformed = "condition"

        metadata_groups_list[[name]] = .group_metadata(metadata, time, col, summarize)
        .logger(metadata_groups_list[[name]], "Metadata")
        if (name %in% names(cluster_list)) {
            cluster = cluster_list[[name]]
            cluster_list[[name]] = .filter(cluster[as.character(cluster$genes) %in% rownames(matrix_list[[name]]), ])
            cluster_expression[[name]] = .median_per_cluster(
                matrixnorm_list[[name]],
                cluster_list[[name]])
        }
        .logger(head(cluster_expression[[name]]), "Cluster summarized")
    }

    base = names(matrix_list)[1]
    .norm_base = matrixnorm_list[[base]]
    .exp_base = cluster_expression[[base]]
    .clus_base = cluster_list[[base]]
    .metadata_group = metadata_groups_list[[base]]
    .maraw_base = matrixabs_list[[base]]
    cat("\n\n## Integration of :", paste(names(matrix_list)),"{.tabset}\n\n")
    df = lapply(rownames(.exp_base), function(nc1){
        .genes_nc1 = as.character(.clus_base$genes[.clus_base$cluster == nc1])
        cat("\n\n### Cluster number ", nc1, "\n\n")
        .plot_base(.genes_nc1, .norm_base, .maraw_base, .metadata_group,
                   time, col_transformed, nc1, paste(base, nc1))
        df = lapply(names(matrix_list), function(name) {
            if (name == base) return(NULL)
            .norm = matrixnorm_list[[name]]
            .ma =  matrix_list[[name]]
            .metadata = metadata_list[[name]]
            map = NULL
            if (name %in% names(mapping))
                map = mapping[[name]]
            .logger(head(map), "Mapping")
            .integrate(nc1, .exp_base, .clus_base, .norm_base, .metadata_group,
                       .maraw_base, .ma, .norm, .metadata,
                       time, col, summarize, col_transformed, name, map)
        })
        do.call(rbind, df)
    })
    df

}

.table_w_fc <- function(dds, contrast){
    if (!contrast[[1]][1] %in% names(colData(dds))) {
        stop("column not in dds object:", contrast[[1]][1])
    }
    fc_df <- do.call(rbind, lapply(contrast, function(cntr){
        tb <- results(dds, contrast = c(cntr[1], cntr[2], cntr[3]),
                      tidy = "TRUE")
        tb  %>% .[, c("row", "log2FoldChange")] %>%
            mutate(comp = paste0(cntr[2], "vs", cntr[3]))
    }))
    fc_df <- fc_df %>% tidyr::spread(comp, log2FoldChange)
    rownames(fc_df) <- fc_df$row
    fc_df[,2:ncol(fc_df)]
}

.run_cluster_profiler <- function(out_df, FDR, FC, org, minc=30){
    .res = as.data.frame(out_df)
    .idx = .res$padj < FDR & .res$absMaxLog2FC > FC
    .idx[is.na(.idx)] = FALSE
    cat("\n\n### GO ontology of DE genes (log2FC >" , FC, " and FDR <", FDR, "): ", sum(.idx),"\n\n")
    ego <- enrichGO(gene = row.names(out_df[.idx,]), keytype = "ENSEMBL",
                    OrgDb = org, ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

    if ("result" %in%  slotNames(ego)) {
        print(knitr::kable(simplify(ego)@result[,1:7]))
        cat("\n\n")
        return(ego)
    }
    return(NULL)
}

.convertIDs <- function(ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
    stopifnot( inherits( db, "AnnotationDb" ) )
    ifMultiple <- match.arg( ifMultiple )
    if (sum(ids %in% keys(db, from))==0)
        return(ids)
    suppressMessages( selRes <- AnnotationDbi::select(
        db, keys = ids, keytype = from, columns = c(from,to) ) )
    if ( ifMultiple == "putNA" ) {
        duplicatedIds <- selRes[ duplicated( selRes[,from] ), from ]
        selRes <- selRes[ !(selRes[,from] %in% duplicatedIds), ]
    }
    return( selRes[ match( ids, selRes[,from] ), to ] )
}

.save_file <- function(dat, fn, basedir=".", csv=TRUE){
    if (is.null(basedir))
        return(invisible(NULL))
    tab <- cbind(id=data.frame(id=row.names(dat)), as.data.frame(dat))
    sep="\t"
    quote=FALSE
    if (csv){
        sep=","
        quote=TRUE
    }
    write.table(tab, file.path(basedir, fn), quote=quote, sep=sep, row.names=F)
}


.pca_loadings = function(counts, ntop=500) {
    pca <- prcomp((t(counts)))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    names(percentVar) = colnames(pca$x)
    pca$percentVar = percentVar
    pca
}

#' Complete report from DESeq2 analysis
#'
#' @param res  output from [DESeq2::results()] function.
#' @param dds  [DESeq2::DESeqDataSet()] object.
#' @param rlogMat matrix from [DESeq2::rlog()] function.
#' @param name string to identify results
#' @param org an organism annotation object, like org.Mm.eg.db.
#' NULL if you want to skip this step.
#' @param FDR int cutoff for false discovery rate.
#' @param do_go boolean if GO enrichment is done.
#' @param FC int cutoff for log2 fold change.
#' @param group string column name in colData(dds) that
#' separates samples in meaninful groups.
#' @param xs string column name in colData(dss)
#' that will be used as X axes in plots (i.e time)
#' @param path_results character path where files are stored.
#' NULL if you don't want to save any file.
#' @param contrast list with character vector indicating the
#' fold change values from different comparisons to add to the output table.
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dse <- DESeq(dse)
#' res <- degResults(dds = dse, name = "test", org = NULL,
#'   do_go = FALSE, group = "group", xs = "group", path_results = NULL)
#' @export
degResults <- function(res=NULL, dds, rlogMat=NULL, name,
                                org=NULL, FDR=0.05, do_go=FALSE,
                                FC=0.1, group="condition", xs="time",
                                path_results =".",
                                contrast=NULL){
    if (is.null(rlogMat)){
        message("Doing rlog...")
        rlogMat <- assay(rlog(dds))
    }
    if (is.null(res)){
        message("Getting result...")
        res <- results(dds)
    }
    cat(paste("## Comparison: ", name, "{.tabset} \n\n"))
    metadata = as.data.frame(colData(dds))
    metadata = metadata[,!grepl("replaceable", names(metadata)),
                        drop=FALSE]
    out_df = as.data.frame(res)
    out_df = out_df[!is.na(out_df$padj),]
    out_df = out_df[order(out_df$padj),]
    if (! is.null(org)){
        out_df$symbol = .convertIDs(rownames(out_df),
                                    "ENSEMBL", "SYMBOL", org, "useFirst")
        out_df$description = .convertIDs(rownames(out_df),
                                         "ENSEMBL", "GENENAME", org, "useFirst")
    }

    cat("\n",paste(capture.output(summary(res)[1:8]), collapse = "<br>"),"\n")

    if (!is.null(contrast)){
        fc_df <- .table_w_fc(dds, contrast)
        cols_names <- names(fc_df)
        out_df <- cbind(out_df, fc_df[rownames(out_df),])
        out_df$absMaxLog2FC <- rowMax(abs(as.matrix(out_df[, cols_names])))
    }else{
        out_df$absMaxLog2FC <- abs(out_df$log2FoldChange)
    }

    fn = paste(name, "_de.csv", sep="")
    fn_log = paste(name, "_log2_counts.csv", sep="")
    .save_file(out_df, fn, path_results, csv=TRUE)
    .save_file(rlogMat, fn_log, path_results)
    cat("\n\nDifferential expression file at: ", fn)
    cat("\n\nNormalized counts matrix file at: ", fn_log)

    cat("\n\n### MA plot plot\n\n")
    DESeq2::plotMA(res)

    cat("\n\n### Volcano plot\n\n")
    stats = as.data.frame(res[,c(2,6)])
    degVolcano(stats, title=name, lfc.cutoff=1.5)

    cat("\n\n### QC for DE genes\n")
    show= !is.na(out_df$pvalue)
    p = degQC(rlogMat[row.names(out_df)[show],],
              metadata[,xs],
              pvalue = out_df[["pvalue"]][show])
    print(p)

    sign = row.names(out_df)[out_df$padj<FDR & !is.na(out_df$padj) & out_df$absMaxLog2FC > FC]
    cat("\n\n### Most significants, FDR<", FDR,
        " and log2FC > ", FC, ": ", length(sign),"\n")

    if (length(sign) < 2){
        cat("Too few genes to plot.")
    }else{
        Heatmap(rlogMat[sign, ], show_row_names = FALSE,
                top_annotation = HeatmapAnnotation(df = metadata),
                clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2",
                clustering_distance_columns = "correlation"
        )
        print(degMDS(rlogMat[sign,],condition = metadata[,xs])+
                  theme_minimal())
    }
    cat("\n")

    cat("\n\n### Plots top 9 most significants\n")
    degPlot(dds, out_df, xs = xs, group = group)
    cat("\n")

    cat("\n\n### Top DE table\n\n")
    print(kable(head(out_df, 20)))
    cat("\n\n")
    goterm = ""
    if (do_go){
        df <- .run_cluster_profiler(out_df, FDR, FC, org)
        if (!is.null(df)){
            .save_file(df@result, paste0(name, "_goenrich.csv"), path_results)
            goterm = df
        }


    }
    return(list(sign=sign, table=out_df, go_res=goterm))
}


#' smart PCA from count matrix data
#'
#' nice plot using ggplot2 from prcomp function
#'
#' @param counts matrix with count data
#' @param metadata dara.frame with sample information
#' @param pc1 character PC to plot on x-axis
#' @param pc2 character PC to plot on y-axis
#' @param condition character column in metadata to use to color samples
#' @param name character if given, column in metadata to print label
#' @param shape character if given, column in metadata to shape points
#' @param data Whether return PCA data or just plot the PCA.
#' @author Lorena Pantano, Rory Kirchner, Michael Steinbaugh
#' @return if `results <-` used, the function return the output
#'   of [prcomp()].
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#' colData(humanGender)[idx,], design=~group)
#' degPCA(log2(counts(dse)+0.5), colData(dse),
#'   condition="group", name="group", shape="group")
#' @export
degPCA <- function(counts, metadata = NULL, condition=NULL,
                   pc1 = "PC1", pc2 = "PC2",
                   name = NULL, shape = NULL,
                   data = FALSE){
    pc <- .pca_loadings(counts)
    # error if pc1/2 are not in columns of pc
    idx1 <- which(names(pc[["percentVar"]]) == pc1)
    idx2 <- which(names(pc[["percentVar"]]) == pc2)
    comps <- data.frame(pc[["x"]])
    comps[["Name"]] <- rownames(comps)
    
    if (!is.null(metadata))
        comps <- bind_cols(comps,
                           as.data.frame(metadata)[as.character(comps$Name), ,
                                                   drop=FALSE])
    # [Feature] check metadata has name, shape, condition
    
    if (!is.null(metadata))
        stopifnot(all.equal(colnames(counts), rownames(metadata)))
    
    p <- ggplot(comps, aes_string(pc1, pc2,
                                  color = condition,
                                  shape = shape)) +
        geom_point(size = 3)
    
    if (!is.null(condition)){
        if (is.factor(comps[[condition]]))
            p <- p + scale_color_brewer(palette = "Set2")
    }
    
    if (!is.null(name))
        p <- p + geom_text(aes_string(label = name),
                           nudge_x = 1, nudge_y = 1)

    p <- p +
        xlab(paste0(pc1, ": ",
                    round(pc[["percentVar"]][idx1] * 100),
                    "% variance")) +
        ylab(paste0(pc2, ": ",
                    round(pc[["percentVar"]][idx2] * 100),
                    "% variance"))
    
    if(data)
        return(list(pca = pc, plot = p))
    p
}

#' Plot MDS from normalized count data
#'
#' Uses cmdscale to get multidimensional scaling of data matrix,
#' and plot the samples with ggplot2.
#' @param counts matrix samples in columns, features in rows
#' @param condition vector define groups of samples in counts.
#' It has to be same order than the count matrix for columns.
#' @param k integer number of dimensions to get
#' @param d type of distance to use, c("euclidian", "cor").
#' @param xi number of component to plot in x-axis
#' @param yi number of component to plot in y-axis
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' degMDS(counts(dse), condition = colData(dse)[["group"]])
#' @export
degMDS = function(counts, condition=NULL, k=2, d="euclidian", xi=1, yi=2) {
    if (d == "euclidian"){
        distances = dist(t(counts))
    }else if (d == "cor"){
        distances = as.dist((1-cor(counts))^2)
    } else {
        stop(d, " is not implemented.")
    }
    fit = cmdscale(distances, eig = TRUE, k = k)
    eigs = data.frame(variance_explained = fit$eig / sum(fit$eig))
    xnames = paste0("PC", 1:k, " ", round(eigs[1:k, 1] * 100, digits = 2), "%")
    df = as.data.frame(fit$points[, c(xi,yi)])
    names(df) = c("one", "two")
    df[["label"]] = rownames(df)
    if(!is.null(condition)) {
        df[["condition"]] = condition
        p = ggplot(df, aes_string("one", "two",
                                  label = "label", color = "condition")) +
            geom_text(aes_string("one", "two",
                                 label = "label"), size = 3) +
            labs(list(x = xnames[xi], y = xnames[yi])) +
            scale_x_continuous(expand = c(0.3, 0.3))
    } else {
        p = ggplot(df, aes_string("one", "two")) +
            geom_text(aes_string("one", "two",
                                 label = "label"), size = 3) +
            labs(list(x = xnames[xi], y = xnames[yi])) +
            scale_x_continuous(expand = c(0.3, 0.3))

    }

    p
}

.remove_low_difference <- function(ma, groupDifference, each_step){
    if (!each_step){
        keep <- rowMax(ma) - rowMin(ma) > groupDifference
    } else {
        keep <- sapply(1:(ncol(ma) -1), function(i){
            abs(ma[,i] - ma[,i + 1]) > groupDifference
        }) %>% apply(., 1, function(x) all(x))
    }
    message("Filtering ", sum(!keep), " after groupDifference applied.")
    if (sum(keep) < 10)
        stop("After applying groupDifference: ", groupDifference,
             " number of features in matrix is less than 10.",
             " Reduce the value and try again.")
    ma[keep,]
}

#' Make groups of genes using expression profile.
#' 
#'  
#' Note that this function doesn't calculate significant
#' difference between groups, so the
#' matrix used as input should be already filtered to contain only
#' genes that are significantly different or the most interesting genes
#' to study.
#'
#' @aliases degPatterns
#' @param ma  log2 normalized count matrix
#' @param metadata  data frame with sample information. Rownames
#'   should match \code{ma} column names
#'   row number should be the same length than p-values vector.
#' @param minc integer minimum number of genes in a group that
#' will be return
#' @param summarize character column name in metadata that will be used to group
#'   replicates. If the column doesn't exist it'll merge the `time` and
#'   the `col` columns, if `col` doesn't exist it'll use `time` only.
#'   For instance, a merge between summarize and time parameters:
#'   control_point0 ... etc
#' @param time character column name in metadata that will be used as
#'   variable that changes, normally a time variable.
#' @param col character column name in metadata to separate
#'   samples. Normally control/mutant
#' @param consensusCluster Indicates whether using [ConsensusClusterPlus]
#'   or [cluster::diana()]
#' @param reduce boolean remove genes that are outliers of the cluster
#'   distribution. `boxplot` function is used to flag a gene in any
#'   group defined by `time` and `col` as outlier and it is removed
#'   from the cluster. Not used if `consensusCluster` is TRUE.
#' @param cutoff This is deprecated.
#' @param scale boolean scale the \code{ma} values by row
#' @param pattern numeric vector to be used to find patterns like this
#'   from the count matrix. As well, it can be a character indicating the
#'   genes inside the count matrix to be used as reference.
#' @param groupDifference Minimum abundance difference between the
#'   maximum value and minimum value for each feature. Please, 
#'   provide the value in the same range than the `ma` value
#'    ( if `ma` is in log2, `groupDifference` should be inside that range).
#' @param eachStep Whether apply `groupDifference` at each stem over
#'   `time` variable. **This only work properly for one group
#'    with multiple time points**.
#' @param plot boolean plot the clusters found
#' @param fixy vector integers used as ylim in plot
#' @param nClusters an integer scalar or vector with the desired number of groups
#' @param skipDendrogram a boolean to run or not dendextend. Temporary fix to memory
#'                  issue in linux.
#' @details 
#' It can work with one or more groups with 2 or
#' more several time points. 
#' Before calculating the genes similarity among samples,
#' all samples inside the same time point (`time` parameter) and
#' group (`col` parameter) are collapsed together, and the `mean`
#' value is the representation of the group for the gene abundance.
#' Then, all pair-wise gene expression is calculated using
#' `cor.test`  R function using kendall as the statistical
#' method. A distance matrix is created from those values.
#' After that, [cluster::diana()] is used for the 
#' clustering of gene-gene distance matrix and cut the tree using
#' the divisive coefficient of the clustering, giving as well by diana.
#' Alternatively, if `consensusCluster` is on, it would use
#' [ConsensusClusterPlus] to cut the tree in stable clusters.
#' Finally, for each group of genes, only the ones that have genes
#' higher than `minc` parameter will be added to the figure.
#' The y-axis in the figure is the results of applying `scale()` 
#' R function, what is similar to creating a 
#' `Z-score` where values are centered to the `mean`  and
#' scaled to the `standard desviation` by each gene.
#' 
#' The different patterns can be merged
#' to get similar ones into only one pattern. The expression
#' correlation of the patterns will be used to decide whether
#' some need to be merged or not.
#' @return list wiht two items:
#' *  `df` is a data.frame
#' with two columns. The first one with genes, the second
#' with the clusters they belong. 
#' * `pass` is a vector of the clusters that pass the `minc` cutoff.
#' * `plot` ggplot figure.
#' * `hr` clustering of the genes in hclust format.
#' * `profile` normalized count data used in the plot.
#' * `raw` data.frame with gene values summarized by biological replicates and
#'   with metadata information attached.
#' * `summarise` data.frame with clusters values summarized by group and
#'   with the metadata information attached.
#' * `normalized` data.frame with the clusters values 
#'   as used in the plot.
#' * `benchmarking` plot showing the different patterns at different
#'   values for clustering cuttree function.
#' * `benchmarking_curve` plot showing how the numbers of clusters and genes
#'   changed at different values for clustering cuttree function.
#' @examples
#' data(humanGender)
#' library(SummarizedExperiment)
#' library(ggplot2)
#' ma <- assays(humanGender)[[1]][1:100,]
#' des <- colData(humanGender)
#' des[["other"]] <- sample(c("a", "b"), 85, replace = TRUE)
#' res <- degPatterns(ma, des, time="group", col = "other")

#' # Use the data yourself for custom figures
#'  ggplot(res[["normalized"]],
#'         aes(group, value, color = other, fill = other)) +
#'   geom_boxplot() +
#'    geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
#'    # change the method to make it smoother
#'    geom_smooth(aes(group=other), method = "lm")
#' @export
degPatterns = function(ma, metadata, minc=15, summarize="merge",
                       time="time", col=NULL,
                       consensusCluster = FALSE,
                       reduce=FALSE, cutoff=0.70,
                       scale=TRUE,
                       pattern = NULL,
                       groupDifference = NULL,
                       eachStep = FALSE,
                       plot=TRUE, fixy=NULL, nClusters = NULL,
                       skipDendrogram=TRUE){
    benchmarking <- NULL
    metadata <- as.data.frame(metadata)
    ma = ma[, row.names(metadata)]
    rownames(ma) = make.names(rownames(ma))
    if (is.null(col)){
        col = "colored"
        metadata[,col] = rep("one_group", nrow(metadata))
    }

    if (!summarize %in% names(metadata))
        metadata[,summarize] = as.factor(paste0(metadata[,col],
                                                metadata[,time]))

    # ensure there are no missing levels in summarize
    metadata[,summarize] = droplevels(metadata[,summarize])

    stopifnot(class(metadata)[1] == "data.frame")
    stopifnot(class(ma)[1] == "matrix" | class(ma)[1] == "data.frame")
    stopifnot(summarize %in% names(metadata))
    stopifnot(time %in% names(metadata))

    ma <- as.matrix(ma)
    if (!is.null(fixy))
        stopifnot(length(fixy) == 2)

    if (nrow(ma)>3000 & is.null(pattern))
        message("A large number of genes was given-- please, ",
                "make sure this is not an error. Normally, ",
                "only DE genes will be useful for this function.")
    message("Working with ", nrow(ma), " genes.")
    
    counts_group <- .summarize_scale(ma,
                                     metadata[[summarize]],
                                     FALSE)
    if (!is.null(groupDifference))
        counts_group <- .remove_low_difference(counts_group,
                                               groupDifference,
                                               eachStep)
    if (scale)
        norm_sign <- t(apply(counts_group, 1, .scale))
    else norm_sign <- counts_group
    
    colnames(norm_sign) <- colnames(counts_group)
    metadata_groups <- metadata %>%
        dplyr::distinct(!!sym(summarize), .keep_all = TRUE)
    rownames(metadata_groups) = metadata_groups[,summarize]
    norm_sign <- norm_sign[, row.names(metadata_groups), drop = FALSE]
    if (nrow(ma) == 1){
        p = .plot_cluster(norm_sign,
                      as.character(rownames(norm_sign)),
                      metadata_groups[,time],
                      metadata_groups[,col], rownames(norm_sign), fixy)
        cluster_genes <- c()
        groups <- as.factor("1")
        names(groups) <- c(row.names(ma))

    }else if (!consensusCluster & is.null(pattern)){
        cluster_genes = .make_clusters(counts_group)
        groups <- .select_genes(cluster_genes, norm_sign, minc,
                               reduce = reduce,
                               cutoff = cutoff, nClusters = nClusters)
        benchmarking <- .benckmark_cutoff(cluster_genes, norm_sign, minc)
    }else if (consensusCluster & is.null(pattern)){
        cluster_genes <- .make_concensus_cluster(counts_group)
        groups <- .select_concensus_genes(cluster_genes)
    }else if(is.character(pattern)){
        stopifnot(pattern %in% rownames(counts_group))
        if (length(pattern) > 1)
            reference <- rowMedians(counts_group[pattern, ])
        else reference <- counts_group[pattern, ]
        cluster_genes <- .find_pattern(counts_group, reference)
        groups <- .select_pattern(cluster_genes)
    }else if(is.numeric(pattern)){
        stopifnot(length(pattern) == ncol(counts_group))
        cluster_genes <- .find_pattern(counts_group, pattern)
        groups <- .select_pattern(cluster_genes)
    }
    temp <- names(groups)
    
    dend_plot <- NA
    if (length(unique(groups)) > 0 & is.null(nClusters) & !skipDendrogram){
        dend <- cluster_genes 
        h = dend$dc
        clust <- cutree(as.hclust(dend), h = h)
        clust.cutree <- dendextend:::cutree(dend, h = h, order_clusters_as_data = FALSE)
        dend <- as.dendrogram(dend, h = h)
        idx <- order(names(clust.cutree))
        clust.cutree <- clust.cutree[idx]
        df.merge <- merge(clust,clust.cutree,by='row.names')
        df.merge.sorted <- df.merge[order(df.merge$y),]
        lbls <- unique(df.merge.sorted$x)
        dend_plot <- dendextend::color_branches(dend, h = h, groupLabels = TRUE,
                                                warn = FALSE) %>% 
            dendextend::set("labels", "") %>% suppressWarnings()
        
        groups <- match(groups, lbls)
        
        if (plot)
            plot(dend_plot, xlab="", ylab="", main="", sub="", axes=FALSE, cex = 2)
    }
    if (length(unique(groups)) > 0 & is.numeric(nClusters) & !skipDendrogram){
        dend <- cluster_genes 
        clust <- cutree(as.hclust(dend), k = nClusters)
        clust.cutree <- dendextend:::cutree(dend, k = nClusters, order_clusters_as_data = FALSE)
        dend <- as.dendrogram(dend, k = nClusters)
        idx <- order(names(clust.cutree))
        clust.cutree <- clust.cutree[idx]
        df.merge <- merge(clust,clust.cutree,by='row.names')
        df.merge.sorted <- df.merge[order(df.merge$y),]
        lbls <- unique(df.merge.sorted$x)
        dend_plot <- dendextend::color_branches(dend, k = nClusters,
                                                groupLabels = TRUE, 
                                                warn = FALSE ) %>%
                                                dendextend::set("labels", "") %>% 
                                                suppressWarnings()
        groups <- match(groups, lbls)
        
        if (plot)
            plot(dend_plot, xlab="", ylab="", main="", sub="", axes=FALSE, cex = 2)
    }
    
    names(groups) <- temp
    df <- data.frame(genes = names(groups), 
                    cluster = groups, stringsAsFactors = FALSE)

    raw <- counts_group %>% as.data.frame %>% 
        rownames_to_column("genes") %>%
        gather(!!sym(summarize), "value", -genes) %>%
        inner_join(metadata_groups %>%
                       mutate_if(is.factor, as.character)) %>%
        inner_join(df, by = "genes")
    
    summarise <- raw %>%
        group_by(!!sym(summarize), !!sym("cluster"),
                 !!sym(time), !!sym(col)) %>%
        summarise(abundance = median(value),
                  n_genes = n()) %>% 
        ungroup()
    
    normalized <- norm_sign %>% as.data.frame() %>% 
        rownames_to_column("genes") %>%
        gather(!!sym(summarize), "value", -genes) %>%
        inner_join(metadata_groups %>%
                       mutate_if(is.factor, as.character)) %>%
        inner_join(df, by = "genes") 

    if (!is.null(benchmarking))
        normalized <- normalized %>% left_join(benchmarking[["genes"]], by = "genes")
    normalized[[time]] = factor(normalized[[time]],
                                levels = levels(metadata[[time]]))
    
    plot_benchmarking <- .plot_benchmarking(normalized, benchmarking, time, col)
    plot_benchmarking_curve <- .plot_benchmarking_curve(benchmarking)
    plotting_data <- list(norm = normalized, time = time, col = col, min_genes = minc) 
    if (length(unique(groups)) > 0){
        p <- degPlotCluster(normalized, time, col, min_genes = minc)
        if (!is.null(fixy))
            p <- p + ylim(fixy[1], fixy[2])
        if (plot)
            print(p)
    }
    
    invisible(list(df = df,
         pass = unique(groups),
         plot = p,
         dend = dend_plot,
         hr = cluster_genes,
         normalized = normalized,
         summarise = summarise,
         raw = raw,
         counts = ma[raw[["genes"]],],
         benchmarking = plot_benchmarking,
         benchmarking_curve = plot_benchmarking_curve))
}

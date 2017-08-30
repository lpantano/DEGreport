# output plots as a way to avoid running the function again
# add plotGenes to documentation. Option to plot by names

.logger = function(toout,msg=""){
    logdebug(paste("\n\nchecking data:" , msg, "\n\n"))
    if (getLogger()[['level']]!=20)
        print(toout)
}


# plot group of genes according time and group
.plot_cluster  = function(norm_sign, g_in_c, xs , groups, title, fixy=NULL) {
    ma <- as.data.frame(norm_sign)[g_in_c,]
    ma_long <- suppressMessages(melt(cbind(gene = row.names(ma), ma),
                                     variable_name = "sample"))
    ma_long$x <- xs[ma_long$sample]
    ma_long$group <- groups[ma_long$sample]
    splan <- max(c(round(length(unique(ma_long$x))/3*2,0), 1))
    p <- suppressWarnings(
        ggplot(ma_long, aes_string(x = "x", y = "value",
                                   fill = "group", color = "group")) +
        geom_boxplot(alpha = 0.3, outlier.size = 0, outlier.shape = NA) +
        geom_point(alpha = 0.4, width = 0.2, size = 1,
                    position = position_jitterdodge(dodge.width = 0.9)) +
        stat_smooth(aes_string(x = "x", y = "value",
                               group = "group", color = "group"),
                    method = "lm", formula = y~poly(x, splan)) +
        ggtitle(paste("Group:", title, "(", length(g_in_c), " genes )")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ylab("scaled expression") + xlab("") )
    if (!is.null(fixy))
        p <- p + ylim(fixy[1], fixy[2])
    if (length(unique(groups)) == 1) {
        p <- p + scale_color_brewer(guide = FALSE, palette = "Set1") +
                scale_fill_brewer(guide = FALSE, palette = "Set1")
    }else{
        p <- p + scale_color_brewer(palette = "Set1") +
            scale_fill_brewer(palette = "Set1")
    }
    p
}

# Scale from 1 to 0 the expression genes
.scale <- function(e){
    .max <- max(e)
    .min <- min(e)
    scale(e)
}

.reduce <- function(group, counts_group, cutoff=0.70){
    ngroup <- unique(group)
    cor <- lapply(ngroup, function(nc1){
        sapply(ngroup, function(nc2){
            g1 = colMeans(counts_group[names(group[group == nc1]),])
            g2 = colMeans(counts_group[names(group[group == nc2]),])
            (1 - cor.test(g1, g2)$estimate)^2
        })
    })
    cor <- do.call(rbind, cor)
    colnames(cor) <- ngroup
    rownames(cor) <- ngroup
    h <- hclust(as.dist(cor), method = "ward.D2")
    c <- cutree(h, h = (1 - cutoff)^2)
    new <- c[as.character(group)]
    names(new) <- names(group)
    new
}


# use diana package to detect clusters
.make_clusters <- function(counts_group){
    m <- (1 - cor(t(counts_group), method = "kendall"))
    d <- as.dist(m^2)
    c <- diana(d, diss = TRUE, stand = FALSE)
    c
}

.select_genes <- function(c, counts_group, minc=15,
                          reduce=FALSE, cutoff=0.30){
    select <- cutree(as.hclust(c), h = c$dc)
    select <- select[select %in% names(table(select))[table(select) > minc]]
    message("Working with ", length(select), " genes after filtering: minc > ",minc)
    if (reduce & length(unique(select) > 1) & ncol(counts_group) > 2)
        select <- .reduce(select, counts_group, cutoff)
    return(select)
}

# summarize matrix into groups and scale if needed
.summarize_scale <- function(ma, group, scale=TRUE){
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
        norm_sign = t(apply(counts_group, 1, .scale))
    }else{
        norm_sign = counts_group
    }
    colnames(norm_sign) = colnames(counts_group)
    norm_sign
}


.median_per_cluster <- function(ma, clusters){
    # matrix, and df
    .logger(table(clusters$cluster),
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

.filter <- function(df){
    .sum <- table(df$cluster)
    .pass <- names(.sum)[.sum > 4]
    df[df$cluster %in% .pass,]
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
    keep = as.character(unique(keep_info$touse))
    if (length(keep) < 5)
        return(NULL)
    clusters <- degPatterns(as.matrix(.ma[keep,]),
                           .metadata,
                           minc = 5, summarize = summarize,
                           time = time, col = col)
    # print(clusters$pass)
    if (length(clusters$pass) == 0)
        return(NULL)
    .exp <- .median_per_cluster(.norm,
                               clusters$df %>%
                               .[clusters$df$cluster %in% clusters$pass,])
    # print(.exp)
    # for (nc2 in rownames(.exp)){
    df = lapply(rownames(.exp), function(nc2){
        #print(nc2)
        # print(.exp_base[nc1,])
        # print(.exp[nc2,colnames(.exp_base)])
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

    metadata_groups = metadata %>% dplyr::distinct_(summarize, .keep_all = TRUE)
    rownames(metadata_groups) = metadata_groups[,summarize]
    return(metadata_groups)
}


#' Integrate data comming from degPattern into one data object
#'
#' The simplest case is if you want to convine the pattern profile
#' for gene expression data and proteomic data. It will use the first element
#' as the base for the integration. Then, it will loop through clusters
#' and run \link{degPatterns} in the second data set to detect patterns that match
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
#' @param res  output from \link[DESeq2]{results} function.
#' @param dds  \link[DESeq2]{DESeqDataSet} object.
#' @param rlogMat matrix from \link[DESeq2]{rlog} function.
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
    p = degQC(out_df$pvalue[show],
              rlogMat[row.names(out_df)[show],],
              metadata[,xs])
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
#' @author Lorena Pantano, Rory Kirchner, Michael Steinbaugh
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#' colData(humanGender)[idx,], design=~group)
#' degPCA(log2(counts(dse)+0.5), colData(dse),
#'   condition="group", name="group", shape="group")
#' @export
degPCA <- function(counts, metadata, condition="condition",
                   pc1="PC1", pc2="PC2",
                   name=NULL, shape=NULL){
    pc = .pca_loadings(counts)
    idx1 = which(names(pc[["percentVar"]]) == pc1)
    idx2 = which(names(pc[["percentVar"]]) == pc2)
    comps = data.frame(pc[["x"]])
    comps$Name = rownames(comps)
    comps = bind_cols(comps,
                      as.data.frame(metadata)[as.character(comps$Name), ,
                                                     drop=FALSE])
    # [Feature] check metadata has name, shape, condition
    # [Feature] check counts has same samples than metadata
    p <- ggplot(comps, aes_string(pc1, pc2, color=condition))
    if (!is.null(shape))
        p <- ggplot(comps, aes_string(pc1, pc2, color=condition, shape=shape))
    if (!is.null(name))
        p <- p + geom_text(aes_string(label=name), nudge_x = 1, nudge_y = 1)
    else
        p <- p + geom_point(size=3)
    p +
        scale_color_brewer(palette = "Set1") +
        xlab(paste0(pc1, ": ",
                    round(pc[["percentVar"]][idx1] * 100),
                    "% variance")) +
        ylab(paste0(pc2, ": ",
                    round(pc[["percentVar"]][idx2] * 100),
                    "% variance")) +
        theme_minimal()
}

#' Plot MDS from normalized count data
#'
#' Uses cmdscale to get multidimensional scaling of data matrix,
#' and plot the samples with ggplot2.
#' @param counts matrix samples in columns, features in rows
#' @param condition vector define groups of samples in counts.
#' It has to be same order than the count matrix for columns.
#' @param k integer number of dimensions to get
#' @param d type of distance to use
#' @param xi number of component to plot in x-axis
#' @param yi number of component to plot in y-axis
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' degMDS(counts(dse), condition=colData(dse)[["group"]])
#' @export
degMDS = function(counts, condition=NULL, k=2, d="euclidian", xi=1, yi=2) {
    if (d=="euclidian"){
        distances = dist(t(counts))
    }else if (d=="cor"){
        distances = as.dist((1-cor(counts))^2)
    }
    fit = cmdscale(distances, eig = TRUE, k = k)
    eigs = data.frame(variance_explained=fit$eig / sum(fit$eig))
    xnames = paste0("PC", 1:k, " ", round(eigs[1:k,1]*100,digits=2), "%")
    df = as.data.frame(fit$points[,c(xi,yi)])
    names(df) = c("one", "two")
    df$label = rownames(df)
    if(!is.null(condition)) {
        df$condition = condition
        p = ggplot(df, aes_string("one", "two",
                                  label = "label", color = "condition")) +
            geom_text(aes_string("one", "two",
                                 label = "label"), size = 3) +
            labs(list(x = xnames[xi], y = xnames[yi])) +
            scale_x_continuous(expand = c(0.3, 0.3))
    }
    else {
        p = ggplot(df, aes_string("one", "two")) +
            geom_text(aes_string("one", "two",
                                 label = "label"), size = 3) +
            labs(list(x = xnames[xi], y = xnames[yi])) +
            scale_x_continuous(expand = c(0.3, 0.3))

    }

    p
}


#' Make groups of genes using expression profile
#'
#' @aliases degPatterns
#' @param ma  log2 normalized count matrix
#' @param metadata  data frame with sample information. Rownames
#' should match \code{ma} column names
#' row number should be the same length than p-values vector.
#' @param minc integer minimum number of genes in a group that
#' will be return
#' @param summarize character column name in metadata that will be used to gorup
#' replicates.
#' For instance, a merge between summarize and time parameters:
#' control_point0 ... etc
#' @param time character column name in metadata that will be used as
#' variable that changes, normally a time variable.
#' @param col character column name in metadata to separate
#' samples. Normally control/mutant
#' @param reduce boolean reduce number of clusters using
#' correlation values between them.
#' @param cutoff integer threshold for correlation
#' expression to merge clusters (0 - 1)
#' @param scale boolean scale the \code{ma} values by row
#' @param plot boolean plot the clusters found
#' @param fixy vector integers used as ylim in plot
#' @details It would be used \link[cluster]{diana} function
#' to detect a value to cut the expression based clustering
#' at certain height. It can work with one or more groups with 2 or
#' more several time points. The different patterns can be merged
#' to get similar ones into only one pattern. The expression
#' correlation of the patterns will be used to decide whether
#' some need to be merged or not.
#' @return list wiht two items. \code{df} is a data.frame
#' with two columns. The first one with genes, the second
#' with the clusters they belong. \code{pass_to_plot} is a vector
#' of the clusters that pass the \code{minc} cutoff.
#' @examples
#' data(humanGender)
#' library(SummarizedExperiment)
#' ma <- assays(humanGender)[[1]][1:100,]
#' des <- colData(humanGender)
#' res <- degPatterns(ma, des, time="group")
#' @export
degPatterns = function(ma, metadata, minc=15, summarize="merge",
                       time="time", col=NULL,
                       reduce=FALSE, cutoff=0.70,
                       scale=TRUE, plot=TRUE, fixy=NULL){
    metadata <- as.data.frame(metadata)
    ma = ma[, row.names(metadata)]
    if (is.null(col)){
        col = "colored"
        metadata[,col] = rep("one_group", nrow(metadata))
    }
    if (!summarize %in% names(metadata))
        metadata[,summarize] = as.factor(paste0(metadata[,col], metadata[,time]))
    stopifnot(class(metadata) == "data.frame")
    stopifnot(class(ma) == "matrix")
    stopifnot(summarize %in% names(metadata))
    stopifnot(time %in% names(metadata))

    if (!is.null(fixy))
        stopifnot(length(fixy) == 2)

    if (nrow(ma)>3000)
        message("Large number of genes given. Please,",
                "make sure is not an error. Normally",
                "Only DE genes are useful for this function.")
    message("Working with ", nrow(ma), " genes.")
    counts_group = t(sapply(rownames(ma), function(g){
        sapply(levels(metadata[,summarize]), function(i){
            idx = which(metadata[,summarize] == i)
            mean(ma[g, idx], na.rm=TRUE)
        })
    }))
    # colnames(counts_group) = unique(metadata[,summarize])

    cluster_genes = .make_clusters(counts_group)
    groups = .select_genes(cluster_genes, counts_group, minc,
                           reduce = reduce,
                           cutoff = cutoff)

    if (scale)
        norm_sign = t(apply(counts_group, 1, .scale))
    else norm_sign = counts_group

    colnames(norm_sign) = colnames(counts_group)
    metadata_groups = metadata %>%
        dplyr::distinct_(summarize, .keep_all = TRUE)
    rownames(metadata_groups) = metadata_groups[,summarize]
    norm_sign = norm_sign[, row.names(metadata_groups)]
    to_plot = unique(groups)
    plots = lapply(to_plot, function(x){
        .plot_cluster(norm_sign,
                      as.character(names(groups[groups == x])),
                      metadata_groups[,time],
                      metadata_groups[,col], x, fixy)
    })
    nc = 3
    all <- NULL
    if (length(plots) < 3)
        nc = length(plots)
    if (length(plots) > 0){
        all <- plot_grid(plotlist = plots, ncol = nc)
        print(all)
    }

    invisible(list(df = data.frame(genes = names(groups), cluster = groups),
         pass = to_plot, plot = all, hr = as.hclust(cluster_genes),
         profile = counts_group))
}

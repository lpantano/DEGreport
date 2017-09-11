# Add function for PCA of discarded genes in the independing filtering by deseq2

#' Filter genes by group
#' 
#' This function will keep only rows that have a minimum counts of
#' 1 at least in a \code{min} number of samples (default 80%).
#' 
#' @param counts Matrix with expression data, columns are samples
#'   and rows are genes or other feature.
#' @param metadata Data.frame with information about
#'   each column in counts matrix. Rownames should match
#' \code{colnames(counts)}.
#' @param group Character column in metadata used to
#'   group samples and applied the cutoff.
#' @param min Numeric value indicating the minimum
#'   number of samples in each group that should have
#'   more than 0 in count matrix.
#' @param minreads Integer minimum number of reads to consider 
#'   a feature expressed.
#' @return count \code{matrix} after filtering genes (features)
#'   with not enough expression in any group.
#' @examples
#' data(humanGender)
#' library(SummarizedExperiment)
#' idx <- c(1:10, 75:85)
#' c <- degFilter(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], "group", min=1)
#' @export
degFilter <- function(counts, metadata, group, min=0.8, minreads=0){
    .unique_group <- as.character(unique(metadata[,group]))
    .keep = sapply(.unique_group, function(g){
        .samples = rownames(metadata)[metadata[, group] == g]
        .n = length(.samples)
        rowSums(counts[, .samples] > minreads) >= .n * min
    })
    counts[rowSums(.keep) > 0, ]
}
    

#' Plot main figures showing p-values distribution and mean-variance correlation
#' 
#' This function joins the output of \code{\link[DEGreport]{degMean}}, 
#' \code{\link[DEGreport]{degVar}} and \code{\link[DEGreport]{degMV}} in a
#' single plot. See these functions for further information.
#' 
#' @param counts Matrix with counts for each samples and each gene.
#' @param groups Character vector with group name for each sample in the
#'   same order than counts column names.
#' @param object [DEGSet] oobject.
#' @param pvalue pvalues of DEG analysis.
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' degQC(counts(dds, normalized=TRUE), colData(dds)[["group"]],
#'   pvalue = res[["pvalue"]])
#' @export
degQC <- function(counts, groups, object=NULL, pvalue=NULL){
    if (is.null(pvalue) & is.null(object))
        stop("You need to provide DEGset object or pvalue.")
    if (!is.null(object) & class(object) != "DEGSet")
        stop("Object should be a DEGSet class.")
    stopifnot(class(counts) %in% c("matrix", "data.frame"))
    
    if (class(object) == "DEGSet"){
        df <- deg(object, tidy = "tibble")
        if (length(intersect(rownames(counts), df[["gene"]])) < nrow(df))
            stop("Not all features in DEGSet are in counts table.")
        counts <- counts[df[["gene"]],]
        pvalue <- df[["pvalue"]]
        padj <- df[["padj"]]
    }else{
        padj <- p.adjust(pvalue, method = "fdr")
    }

    counts <- counts[!is.na(padj),]
    pvalue <- pvalue[!is.na(padj)]
    padj <- padj[!is.na(padj)]
    pmean <- degMean(pvalue, counts) +
        xlab("pvalues along expression quantiles") +
        guides(fill = FALSE)
    pvar <- degVar(pvalue, counts) +
        xlab("pvalues along variance quantiles") +
        guides("") +
        theme(legend.position="top") +
        theme(legend.text=element_text(size=3),
              legend.title=element_blank(),
              legend.key.size = unit(0.5,"line"))
    pmv <- degMV(groups, padj, counts)
    suppressWarnings(ggdraw() + 
        draw_plot(pvar, 0, 0.4, 0.6, 0.6 ) +
        draw_plot(pmv, 0.6, 0, 0.4, 0.7) +
        draw_plot(pmean, 0, 0, 0.6, 0.4))
}

#' Distribution of gene ratios used to calculate Size Factors.
#' 
#' This function check the median ratio normalization used by
#' DESeq2 and similarly by edgeR to visualy check whether
#' the median is the best size factor to represent depth.
#' 
#' @aliases degCheckFactors
#' @param counts Matrix with counts for each samples and each gene.
#'   row number should be the same length than pvalues vector.
#' @param each Plot each sample separately.
#' @return ggplot2 object
#' @details This function will plot the gene ratios for each sample. To calculate
#' the ratios, it follows the simliar logic than DESeq2/edgeR uses, where the expression
#' of each gene is divided by the mean expression of that gene. The distribution
#' of the ratios should approximate to a normal shape and the factors should be similar
#' to the median of distributions. If some samples show different distribution,
#' the factor may be bias due to some biological or technical factor.
#' @references 
#' * Code to calculate size factors comes from
#'   [DESeq2::estimateSizeFactorsForMatrix()].
#' @examples
#' data(humanGender)
#' library(SummarizedExperiment)
#' degCheckFactors(assays(humanGender)[[1]][, 1:10])
#' @export
degCheckFactors <-
    function(counts, each = FALSE)
    {
        counts <- as.data.frame(counts)
        geoMeanNZ <- function(x) {
            if (all(x == 0)) { 0 } else { exp( sum(log(x[x > 0])) / length(x) ) }
        }
        geoMeans <- apply(counts, 1, geoMeanNZ)
        loggeomeans <- log(geoMeans)
        df <- lapply(colnames(counts), function(s) {
            cnts <- counts[[s]]
            r <- (log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]
            data.frame(ratios = r,
                       sample = s)
        }) %>% bind_rows()
        p <- ggplot(df, aes(ratios, group = sample)) +
            geom_density() +
            theme_bw() +
            xlim(-4, 4)
        if (each)
            p <- p + facet_wrap(~sample)
        p
    }


#' Distribution of pvalues by expression range
#' 
#' This function plot the p-values distribution colored by
#' the quantiles of the average count data.
#' 
#' @param pvalues pvalues of DEG analysis.
#' @param counts  Matrix with counts for each samples and each gene.
#'   row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' degMean(res[, 4], counts(dds))
#' @export
degMean <-
    function(pvalues, counts)
{
    meanv  <-  apply(counts, 1, mean)
    q <- quantile(meanv, seq(.1, 1, .1))
    q <- q[!duplicated(q)]
    meanvfac <- cut(meanv,
                breaks=c(0, q),
                labels=names(q),
                right=TRUE)
    pvalfac  <-  cut(pvalues,
               breaks=c(-1, seq(.1, 1, .1)),
               labels=seq(0.1, 1, 0.1), right=TRUE)
    d  <-  data.frame(pvalues=factor(pvalfac),
                meanv=factor(meanvfac))
    suppressWarnings(ggplot(d, aes(pvalues, fill=meanv))+
        geom_bar()+
        theme_bw() +
        scale_fill_brewer("mean\nquantiles", palette="RdYlBu")) +
        labs(list(x="p-values", y="# genes"))
}

#' Distribution of pvalues by standard desviation range
#' 
#' This function pot the p-valyes distribution colored by
#' the quantiles of the standard desviation of count data.
#' 
#' @param pvalues pvalues of DEG analysis
#' @param counts Matrix with counts for each samples and each gene.
#'   row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' degVar(res[, 4], counts(dds))
#' @export
degVar <-
    function(pvalues, counts)
{
    sdv <- apply(counts, 1, sd)
    q <- quantile(sdv, seq(.1, 1, .1))
    q <- q[!duplicated(q)]
    sdvfac <- cut(sdv,
                breaks=c(0, q),
                labels=names(q),
                right=TRUE)
    pvalfac <- cut(pvalues,
               breaks=c(-1, seq(.1, 1, .1)),
               labels=seq(0.1, 1, 0.1), right=TRUE)
    d <- data.frame(pvalues=factor(pvalfac), 
                sdv=factor(sdvfac))
    suppressWarnings(ggplot(d, aes(pvalues, fill=sdv))+
        geom_bar()+
        theme_bw() +
        scale_fill_brewer("variance\nquantiles", palette="RdYlBu") +
        labs(list(x="p-values", y="# genes")))
}

#' Correlation of the standard desviation and the mean of the abundance of a
#' set of genes.
#' 
#' @param group Character vector with group name for each sample in the
#'   same order than counts column names.
#' @param pvalues pvalues of DEG analysis.
#' @param counts Matrix with counts for each samples and each gene.
#' @param sign Defining the cutoff to label significant features.
#'   row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' degMV(colData(dds)[["group"]],
#'       res[, 4],
#'       counts(dds, normalized = TRUE))
#' @export
degMV <-
    function(group, pvalues, counts, sign=0.01)
{
    var_ma <- sapply(unique(as.character(group)), function(g){
        apply(counts[, group==g], 1, sd, na.rm=TRUE)
    })    
    sdv <- apply(var_ma, 1, max)
    mean_ma <- sapply(unique(as.character(group)), function(g){
        apply(counts[, group==g], 1, median, na.rm=TRUE)
    })    
    meanv <- apply(mean_ma, 1, min)
    pv <- cut(pvalues, breaks=c(-1, sign, 1.1),
          labels=c("Sign", "NoSig"))
    # pv[is.na(pvalues)] <- "NoSig"
    d <- data.frame(pvalues=pv, max_sd=log2(sdv), 
                min_median=log2(meanv))
    red <- d[d[["pvalues"]] == "Sign" & !is.na(d[["pvalues"]]),]
    black <- rgb(0.9, 0.9, 0.9, 0.6)
    suppressWarnings(ggplot(d, aes_string("min_median", "max_sd"))+
                         geom_point(data = d, color = black) +
                         geom_point(data = red, color = "red") +
                         scale_color_manual(values=c("red", black))+
                         theme_bw()+
                         stat_quantile(aes(min_median, max_sd), colour="blue",
                                       quantiles = c(0.025, 0.975), 
                                       linetype=2, formula=y ~ x)) +
        theme(legend.position="top")
    }

#' Distribution of expression of DE genes compared to the background
#'
#' @aliases degMB
#' @param tags  List of genes that are DE.
#' @param group Character vector with group name for each sample in the
#'   same order than counts column names.
#' @param counts  Matrix with counts for each samples and each gene
#'   Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' degMB(row.names(res)[1:20], colData(dds)[["group"]],
#'   counts(dds, normalized = TRUE))
#' @export
degMB <-
    function(tags, group, counts, pop=400)
{
    delen <- length(tags)
    g <- ""

    # if counts is tiny, sample with replacement
    if(nrow(counts) < pop) {
        warning("The number of genes < samples requested, sampling with replacement.")
        replace = TRUE
    }
    else {
        replace = FALSE
    }
    rand <- sample(row.names(counts), pop, replace=replace)
    
    sign_group <- lapply(unique(as.character(group)), function(g){
        m <- apply(counts[tags, group==g, drop=FALSE], 1, mean, na.rm=TRUE)
        data.frame(mean = m, group = g,
                   type = paste("significants -", g),
                   stringsAsFactors = FALSE)
    })    
    
    rand_group <- lapply(unique(as.character(group)), function(g){
        m <- apply(counts[rand, group==g, drop=FALSE], 1, mean, na.rm=TRUE)
        data.frame(mean = m, group = g,
                   type = paste("background -", g),
                   stringsAsFactors = FALSE)
    }) 

    res <- bind_rows(c(sign_group, rand_group))

    suppressWarnings(ggplot(res, aes_string("type", "mean", 
        fill="group", colour="group"))+
        geom_violin(alpha=0.2)+
        scale_y_log10()+
        theme_bw())
}

#' Distribution of the standard desviation of
#'     DE genes compared to the background
#' @aliases degVB
#' @param tags List of genes that are DE.
#' @param group Character vector with group name for each sample in the
#'   same order than counts column names.
#' @param counts  matrix with counts for each samples and each gene.
#'     Should be same length than pvalues vector.
#' @param pop Number of random samples taken for background comparison.
#' @return ggplot2 object
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' degVB(row.names(res)[1:20], colData(dds)[["group"]],
#'   counts(dds, normalized = TRUE))
#' @export
degVB <-
    function(tags, group, counts, pop=400)
{
    delen <- length(tags)
    g <- ""
    # if counts is tiny, sample with replacement
    if(nrow(counts) < pop) {
        warning("The number of genes < samples requested, sampling with replacement.")
        replace = TRUE
    }
    else {
        replace = FALSE
    }
    rand <- sample(row.names(counts), pop)
    sign_group <- lapply(unique(as.character(group)), function(g){
        m <- apply(counts[tags, group==g, drop=FALSE], 1, sd, na.rm=TRUE)
        data.frame(sd = m, group = g,
                   type = paste("significants -", g),
                   stringsAsFactors = FALSE)
    })    
    
    rand_group <- lapply(unique(as.character(group)), function(g){
        m <- apply(counts[rand, group==g, drop=FALSE], 1, sd, na.rm=TRUE)
        data.frame(sd = m, group = g,
                   type = paste("background -", g),
                   stringsAsFactors = FALSE)
    }) 
    
    res <- bind_rows(c(sign_group, rand_group))
    
    suppressWarnings(ggplot(res, aes_string("type", "sd", 
                                            fill="group", colour="group"))+
                         geom_violin(alpha=0.2)+
                         scale_y_log10()+
                         theme_bw())
    }

degNcomb <- function() {
    .Deprecated("DESeq2::lcfShrink")
}

degComb <- function()
{
    .Deprecated("DESeq2::lcfShrink")
}

degFC <- function()
{
    .Deprecated("DESeq2::lcfShrink")
}

degBI <- function()
{
    .Deprecated("DESeq2::lcfShrink")
}

degBIcmd <- function()
{
    .Deprecated("DESeq2::lcfShrink")
}

degRank <- function()
{
    .Deprecated("DESeq2::lcfShrink")
}


degPR <- function()
{
    .Deprecated("DESeq2::lcfShrink")
}

#' Create a deg object that can be used to plot expression values
#'     at shiny server:runGist(9930881)
#' @aliases degObj
#' @param counts Output from get_rank function.
#' @param design Colour used for each gene.
#' @param outfile File that will contain the object.
#' @return R object to be load into vizExp.
#' @examples 
#' data(humanGender)
#' library(SummarizedExperiment)
#' degObj(assays(humanGender)[[1]], colData(humanGender), NULL)
#' @export
degObj <-
    function(counts, design, outfile)
{
    deg <- NULL
    deg <- list(counts, design)
    if (!is.null(outfile))
        save(deg, file=outfile)
    if (is.null(outfile))
        message("please, give an output file name.")
}


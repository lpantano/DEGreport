# add function for PCA of discarded genes in the independing filtering by deseq2

#' Filter genes by group
#' 
#' This function will keep only rows that have a minimum counts of
#' 1 at least in a \code{min} number of samples (default 80%).
#' 
#' @param counts matrix with expression data, columns are samples
#' and rows are genes or other feature
#' @param metadata data.frame with information about
#' each column in counts matrix. Rownames should match
#' \code{colnames(counts)}.
#' @param group character column in metadata used to
#' group samples and applied the cutoff
#' @param min numeric value indicating the minimum
#' number of samples in each group that should have
#' more than 0 in count matrix.
#' @param minreads integer minimum number of reads to consider 
#' a feature expressed.
#' @return count \code{matrix} after filtering genes (features)
#' with not enough expression in any group.
#' @examples
#' data(humanSexDEedgeR)
#' idx <- c(1:10, 75:85)
#' c <- degFilter(humanSexDEedgeR$counts[1:1000, idx],
#' humanSexDEedgeR$samples[idx,], "group", min=1)
degFilter <- function(counts, metadata, group, min=0.8, minreads=0){
    .unique_group <- as.character(unique(metadata[,group]))
    .keep = sapply(.unique_group, function(g){
        .samples = rownames(metadata)[metadata[,group] == g]
        .n = length(.samples)
        rowSums(counts[, .samples] > minreads) >= .n * min
    })
    counts[rowSums(.keep) > 0,]
}
    

#' Plot main figures showing p-values distribution and mean-variance correlation
#' 
#' This function joins the output of \code{\link[DEGreport]{degMean}}, 
#' \code{\link[DEGreport]{degVar}} and \code{\link[DEGreport]{degMV}} in a
#' single plot. See these functions for further information.
#' 
#' @param pvalue  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene.
#' @param groups character vector with group name for each sample in the
#' same order than counts column names.
#' @return ggplot2 object
#' @examples
#' library(DESeq2)
#' data(humanSexDEedgeR)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(humanSexDEedgeR$counts[1:1000, idx], 
#' humanSexDEedgeR$samples[idx,], design=~group)
#' dse <- DESeq(dse)
#' res <- results(dse)
#' degQC(res$pvalue, counts(dse, normalized=TRUE),colData(dse)$group)
degQC <- function(pvalue, counts, groups){
    pmean <- degMean(pvalue, counts) + guides(fill=FALSE)
    pvar <- degVar(pvalue, counts)+ theme(legend.position="top")
    pmv <- degMV(groups, pvalue, counts)
    suppressWarnings(ggdraw() + 
        draw_plot(pmean, 0, 0.5, 0.6, 0.4 ) +
        draw_plot(pmv, 0.6, 0, 0.4, 0.7) +
        draw_plot(pvar, 0, 0, 0.6, 0.6))
}

#' Distribution of gene ratios used to calculate Size Factors.
#' @aliases degCheckFactors
#' @param counts  matrix with counts for each samples and each gene.
#' row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @details This function will plot the gene ratios for each sample. To calculate
#' the ratios, it follows the simliar logic than DESeq2/edgeR uses, where the expression
#' of each gene is divided by the mean expression of that gene. The distribution
#' of the ratios should approximate to a normal shape and the factors should be similar
#' to the median of distributions. If some samples show different distribution,
#' the factor may be bias due to some biological or technical factor.
#' @examples
#' data(DEGreportSet)
#' degCheckFactors(DEGreportSet$counts[,1:10])
degCheckFactors <-
    function(counts)
    {
        meanv  <-  rowMeans(counts)
        ratios <- sweep(counts, 1, meanv, "/")
        df <- suppressWarnings(reshape::melt.data.frame(as.data.frame(ratios)))
        suppressWarnings(ggplot(df, aes(value))+
                             geom_histogram(binwidth = 0.3)+
                             theme_bw() +
                             facet_wrap(~variable) +
                             xlim(-4,4))
    }


#' Distribution of pvalues by expression range
#' @aliases degMean
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene.
#' row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' degMean(DEGreportSet$deg[,4],DEGreportSet$counts)
degMean <-
    function(pvalues,counts)
{
    meanv  <-  apply(counts,1,mean)
    q <- quantile(meanv,seq(.1,1,.1))
    q <- q[!duplicated(q)]
    meanvfac <- cut(meanv,
                breaks=c(0,q),
                labels=names(q),
                right=TRUE)
    pvalfac  <-  cut(pvalues,
               breaks=c(-1,seq(.1,1,.1)),
               labels=seq(0.1,1,0.1),right=TRUE)
    d  <-  data.frame(pvalues=factor(pvalfac),
                meanv=factor(meanvfac))
    suppressWarnings(ggplot(d,aes(pvalues,fill=meanv))+
        geom_bar()+
        theme_bw() +
        scale_fill_brewer("mean\nquantiles",palette="RdYlBu")) +
        labs(list(x="p-values",y="# genes"))
}

#' Distribution of pvalues by standard desviation range
#' @aliases degVar
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene.
#' row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' degVar(DEGreportSet$deg[,4],DEGreportSet$counts)
degVar <-
    function(pvalues,counts)
{
    sdv <- apply(counts,1,sd)
    q <- quantile(sdv,seq(.1,1,.1))
    q <- q[!duplicated(q)]
    sdvfac <- cut(sdv,
                breaks=c(0,q),
                labels=names(q),
                right=TRUE)
    pvalfac <- cut(pvalues,
               breaks=c(-1,seq(.1,1,.1)),
               labels=seq(0.1,1,0.1),right=TRUE)
    d <- data.frame(pvalues=factor(pvalfac),
                sdv=factor(sdvfac))
    suppressWarnings(ggplot(d,aes(pvalues,fill=sdv))+
        geom_bar()+
        theme_bw() +
        scale_fill_brewer("variance\nquantiles",palette="RdYlBu") +
        labs(list(x="p-values",y="# genes")))
}

#' Correlation of the standard desviation and the mean of the abundance of a
#' set of genes.
#' @aliases degMV
#' @param group character vector with group name for each sample in the
#' same order than counts column names.
#' @param pvalues  pvalues of DEG analysis
#' @param counts  matrix with counts for each samples and each gene.
#' @param sign defining the cutoff to label significant features.
#' row number should be the same length than pvalues vector.
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' degMV(c(rep("M", length(DEGreportSet$g1)), rep("F", length(DEGreportSet$g2))),
#'       DEGreportSet$deg[,4],
#'       DEGreportSet$counts)
degMV <-
    function(group, pvalues, counts, sign=0.01)
{
    var_ma <- sapply(unique(as.character(group)), function(g){
        apply(counts[,group==g], 1, sd, na.rm=TRUE)
    })    
    sdv <- apply(var_ma, 1, max)
    mean_ma <- sapply(unique(as.character(group)), function(g){
        apply(counts[,group==g], 1, mean, na.rm=TRUE)
    })    
    meanv <- apply(mean_ma, 1, min)
    pv <- cut(pvalues,breaks=c(-1,sign,1.1),
          labels=c("Sign","NoSig"))
    d <- data.frame(pvalues=pv,sdv=log2(sdv),
                meanv=log2(meanv))
    suppressWarnings(ggplot(d,aes(meanv,sdv,
    colour=pvalues))+
    geom_point()+
    scale_color_manual(values=c("red",rgb(0.9,0.9,0.9,0.6)))+
    theme_bw()+
    stat_quantile(aes(meanv,sdv),colour="blue",
        quantiles = c(0.025,0.975),
        linetype=2,formula=y ~ x)) +
        theme(legend.position="top")
}

#' Distribution of expression of DE genes compared to the background
#'
#' @aliases degMB
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene
#' Should be same length than pvalues vector
#' @param pop number of random samples taken for background comparison
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' detag <- row.names(DEGreportSet$deg[1:10,])
#' degMB(detag,DEGreportSet$g1,DEGreportSet$g2,DEGreportSet$counts)
degMB <-
    function(tags,g1,g2,counts,pop=400)
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
    rand <- sample(row.names(counts),pop, replace=replace)
    g1var <- apply(counts[tags,g1,drop=FALSE],1,mean)    
    g2var <- apply(counts[tags,g2,drop=FALSE],1,mean)    
  
    rand.s1 <- apply(counts[rand,g1,drop=FALSE],1,mean)
    rand.s2 <- apply(counts[rand,g2,drop=FALSE],1,mean)
    res <- rbind(data.frame(g="g1",mean=g1var),
             data.frame(g="g2",mean=g2var))
    res <- rbind(res,data.frame(g="r1",mean=rand.s1),
             data.frame(g="r2",mean=rand.s2))

    suppressWarnings(ggplot(res,aes(g,mean,
        fill=g,colour=g))+
        geom_violin(alpha=0.2)+
        scale_y_log10()+
        theme_bw())
}

#' Distribution of the standard desviation of
#'     DE genes compared to the background
#' @aliases degVB
#' @param tags  list of genes that are DE
#' @param g1 list of samples in group 1
#' @param g2 list of samples in group 2
#' @param counts  matrix with counts for each samples and each gene.
#'     Should be same length than pvalues vector.
#' @param pop number of random samples taken for background comparison
#' @return ggplot2 object
#' @examples
#' data(DEGreportSet)
#' detag <- row.names(DEGreportSet$deg[1:10,])
#' degVB(detag,DEGreportSet$g1,DEGreportSet$g2,DEGreportSet$counts)
degVB <-
    function(tags,g1,g2,counts,pop=400)
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
    rand <- sample(row.names(counts),pop)
    g1var <- apply(counts[tags,g1,drop=FALSE],1,sd)
    g2var <- apply(counts[tags,g2,drop=FALSE],1,sd)

    rand.s1 <- apply(counts[rand,g1,drop=FALSE],1,sd)
    rand.s2 <- apply(counts[rand,g2,drop=FALSE],1,sd)
    res <- rbind(data.frame(g="g1",var=g1var),
             data.frame(g="g2",var=g2var))
    res <- rbind(res,data.frame(g="r1",var=rand.s1),
             data.frame(g="r2",var=rand.s2))

    suppressWarnings(ggplot(res,aes(g,var,fill=g,colour=g))+
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

degBIcmd <- function(x,iter=1000)
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
#' @param counts output from get_rank function
#' @param design colour used for each gene
#' @param outfile file that will contain the object
#' @return R object to be load into vizExp
#' @examples 
#' data(DEGreportSet)
#' de = data.frame(row.names=colnames(DEGreportSet$counts),
#' sex = c(rep("M", length(DEGreportSet$g1)),
#'         rep("F", length(DEGreportSet$g2))))
#' degObj(DEGreportSet$counts, de, NULL)
degObj <-
    function(counts,design,outfile)
{
    deg <- NULL
    deg <- list(counts, design)
    if (!is.null(outfile))
        save(deg, file=outfile)
    if (is.null(outfile))
        message("please, give an output file name.")
}


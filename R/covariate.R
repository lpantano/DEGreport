# find inter class correlation between factor and continuous covariates
# inspired from http://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables
.getfactorcontassociationstatistics <- function(factorcontnames,
                                                covariates,
                                                na.action='remove',
                                                alpha = 0.05){
    if (na.action == "remove")
        covariates = na.omit(covariates[,factorcontnames])

    covariates[,2] <- as.numeric(covariates[,2])
    data <- list(x = covariates[, 1],
                y = covariates[, 2])
    stats <- cor.test(scale(data[[1]]),
                     scale(data[[2]]),
                     method = "kendall", exact = FALSE)
    cor <- stats$estimate
    pval <- stats$p.value

    return(c(estimate = cor, pval = pval))
}

# function to run principal component analysis
.runpca <- function(genesbysamples, scale_data_for_pca = TRUE,
                    min_pve_pct_pc = 1.0) {

    # estimate variance in data by pc:
    pca.res <- prcomp(t(genesbysamples), center=TRUE,
                      scale.=scale_data_for_pca, retx=TRUE)

    # examine how much variance is explained by pcs,
    # and consider those with pve >= (min_pve_pct_pc %):
    pc.var <- pca.res$sdev^2
    pve <- 100 * (pc.var / sum(pc.var))
    npca <- max(1, length(which(pve >= min_pve_pct_pc)))

    samplepcvals <- pca.res$x[, 1:npca, drop=FALSE]

    list(samplepcvals=samplepcvals, pve=pve)
}

# function to calculate correlation and plot
.calccompletecorandplot <- function(compare_data, covar_data,
                                   correlationtype, title,
                                   weights = NULL,
                                   exclude_vars_from_fdr=NULL,
                                   max_fdr = 0.1) {
    # get factor and continuous covariates
    character_vars <- lapply(covar_data, class) == "character"
    covar_data[, character_vars] <- apply(covar_data[, character_vars,
                                                     drop = FALSE],
                                          1L, as.factor)
    
    factorcovariates <- colnames(covar_data)[sapply(covar_data, is.factor)]
    contcovariates <- colnames(covar_data)[sapply(covar_data, is.numeric)]
    
    all_covariates <- cbind(covar_data[, contcovariates, drop = FALSE],
                            covar_data[, factorcovariates, drop = FALSE] %>%
                                mutate_all(as.numeric))
    
    cov_cor <- corr.test(compare_data,
                         all_covariates,
                         use = 'pairwise.complete.obs',
                         method = correlationtype,
                         adjust = "none")
    all_cor_vals <- cov_cor$r
    all_cor_p <- cov_cor$p
    
    rownames(all_cor_vals) <- colnames(compare_data)
    colnames(all_cor_vals) <- colnames(all_covariates)
    rownames(all_cor_p) <- colnames(compare_data)
    colnames(all_cor_p) <- colnames(all_covariates)
    
    effects.significantcovars <- all_cor_vals
    effects.significantcovars[all_cor_p > max_fdr] <- 0
    effects.significantcovars <- colSums(abs(effects.significantcovars) * replicate(dim(effects.significantcovars)[2L], weights / sum(weights)))
    effects.significantcovars <- effects.significantcovars[order(abs(effects.significantcovars), decreasing = TRUE)]

    cor_mat <- melt(all_cor_p, varnames = c("compare", "covar"))
    colnames(cor_mat)[colnames(cor_mat) == "value"] <- "pvalue"

    cor_mat$compare <- factor(cor_mat$compare, levels = rownames(all_cor_p))
    cor_mat$covar <- factor(cor_mat$covar, levels = colnames(all_cor_p))

    cor_mat$r <- melt(all_cor_vals)$value
    cor_mat$fdr <- p.adjust(cor_mat$pvalue, method = "fdr")
    return(list(mat = cor_mat,
                effects.significantcovars = effects.significantcovars))
}

#' Find correlation between pcs and covariates
#'
#' This function will calculate the pcs using prcomp function,
#' and correlate categorical and numerical variables from
#' metadata.
#'
#' @author: Lorena Pantano, Kenneth Daily and Thanneer Malai Perumal
#'
#' @param counts normalized counts matrix
#' @param metadata data.frame with samples metadata.
#' @param fdr numeric value to use as cutoff to determine
#' the minimum fdr to consider significant correlations
#' between pcs and covariates.
#' @param scale boolean to determine wether counts matrix should be
#' scaled for pca. default FALSE.
#' @param min_pc_pct numeric value that will be used as cutoff
#' to select only pcs that explain more variability than this.
#' @param correlation character determining the method for the
#' correlation between pcs and covariates.
#'
#' @return: list:
#' a)significantcovars,
#' b)plot,
#' c)cor_matrix,
#' d)effects.significantcovars: that is pcs pct * absolute correlation between covariate and pcs,
#' e)pcs_matrix: pcs loading for each sample
#' @examples
#' data(humanSexDEedgeR)
#' library(DESeq2)
#' idx <- c(1:5, 75:80)
#' dse <- DESeqDataSetFromMatrix(humanSexDEedgeR$counts[1:1000, idx],
#' humanSexDEedgeR$samples[idx,], design=~group)
#' res <- degCovariates(log2(counts(dse)+0.5),
#' colData(dse))
#' res$plot
degCovariates <- function(counts, metadata,
                                      fdr = 0.1,
                                      scale = FALSE,
                                      min_pc_pct = 5.0,
                                      correlation = "spearman") {
    title <- paste(ifelse(scale, "s", "un-s"), "caled ",
                  " data in pca; pve >= ",
                  min_pc_pct, "%; ", correlation,
                  " correlations ", sep = "")
    message(paste("\nrunning pca and calculating correlations for:\n",
                  title, sep = ""))

    metadata <- as.data.frame(metadata)
    pcares <- .runpca(genesbysamples = counts,
                      scale_data_for_pca = scale,
                      min_pve_pct_pc = min_pc_pct)

    samplepcvals <- pcares$samplepcvals
    pve <- pcares$pve
    npca <- ncol(samplepcvals)
    colnames(samplepcvals) = paste(colnames(samplepcvals), " (",
                                   sprintf("%.2f", pve[1L:npca]), "%)",
                                   sep = "")

    # find covariates without any missing data
    samplesbyfullcovariates <- metadata[, which(apply(metadata, 2L, function(dat) all(!is.na(dat))))]

    exclude_vars_from_fdr <- setdiff(colnames(metadata),
                                    colnames(samplesbyfullcovariates))

    corrres <- .calccompletecorandplot(samplepcvals,
                                      samplesbyfullcovariates,
                                      correlation,
                                      title,
                                      weights = pve[1:dim(samplepcvals)[2L]],
                                      exclude_vars_from_fdr)

    significantcovars <- corrres$mat[corrres$mat$fdr < fdr,"covar"]
    ma = corrres$mat
    ma$r[ma$fdr > fdr] <- NA
    p = ggplot(ma, aes_(fill = ~r,x = ~covar,y = ~compare)) +
        geom_tile() +
        theme_minimal() +
        ggtitle(title) +
        scale_fill_gradient2(low = "darkblue", high = "darkorange",
                             guide = "colorbar", na.value = "grey90",
                             limits = c(-1,1)) +
        theme(axis.text.x = element_text(angle = 90,
                                         hjust = 1,
                                         vjust = 0.5))
    print(p)
    invisible(list(significantcovars = significantcovars,
                plot = p,
                cor_matrix = corrres$mat,
                effects.significantcovars = corrres$effects.significantcovars,
                pcs_matrix = samplepcvals))
}

degClean <- function(ma){
    lapply(ma, function(x) {
        if (length(unique(x)) < length(x) * 0.20 & is.numeric(x))
            x <- as.factor(x)
        if ((class(x) %in% c("factor", "character"))) {
            .f = as.factor(x)
            if (length(levels(.f)) < 2)
                return(NULL)
            if (length(levels(.f)) > length(.f) * 0.80)
                return(NULL)
            return(.f)
        }
        if (sd(x, na.rm = TRUE) == 0)
            return(NULL)
        as.numeric(x)
    })  %>% Filter(Negate(is.null), .) %>% bind_cols
}

#' Calculate the correlation relationshipt among all covariates
#' in the metadata table
#'
#' This function will calculate the correlation among
#' all columns in the metadata
#'
#' @author: Lorena Pantano, Kenneth Daily and Thanneer Malai Perumal
#'
#' @param metadata data.frame with samples metadata.
#' @param fdr numeric value to use as cutoff to determine
#' the minimum fdr to consider significant correlations
#' between pcs and covariates.
#'
#' @return: list:
#' a) cor, data.frame with pair-wise correlations, pvalues, FDR
#' b) corMat, data.frame with correlation matrix
#' c) fdrMat, data.frame with FDR matrix
#' b) plot, Heatmap plot of correlation matrix
#' @examples
#' data(humanSexDEedgeR)
#' library(DESeq2)
#' idx <- c(1:5, 75:80)
#' dse <- DESeqDataSetFromMatrix(humanSexDEedgeR$counts[1:1000, idx],
#' humanSexDEedgeR$samples[idx,], design=~group)
#' cor <- degCorCov(colData(dse))
degCorCov <- function(metadata, fdr=0.05){
    clean <- degClean(metadata) %>%
        mutate_all(as.numeric)
    cor = .calccompletecorandplot(clean,
                                      clean,
                                      "kendall",
                                      "",
                                      weights = 1)
    cor
    corMat <- cor$mat[, c("r", "compare", "covar")] %>%
        spread(compare, r) %>% remove_rownames() %>%
        column_to_rownames("covar")
    fdrMat <- cor$mat[, c("fdr", "compare", "covar")] %>%
        spread(compare, fdr) %>% remove_rownames() %>%
        column_to_rownames("covar")
    corMat[fdrMat > fdr] <- 0
    # corMat[fdrMat > 0.05] <- NA
    cols <- colorRampPalette(c("steelblue", "white", "orange"))(10)
    if (sum(!is.na(corMat)) > 2) {
        p <- Heatmap(corMat, col = cols)
        print(p)

    }
    invisible(list(cor = cor$mat, corMat = corMat, fdrMat = fdrMat, plot = p))

}
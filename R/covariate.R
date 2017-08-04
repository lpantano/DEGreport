# find inter class correlation between factor and continuous covariates
# inspired from http://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables
.getfactorcontassociationstatistics <- function(factorcontnames,covariates,
                                               na.action='remove',
                                               alpha = 0.05){
    if (na.action == "remove")
        covariates = na.omit(covariates[,factorcontnames])

    covariates[,2] <- as.numeric(covariates[,2])
    stats = ICC(covariates[,factorcontnames], alpha = alpha)
    pval = summary(aov(covariates[,factorcontnames[1]]~covariates[,factorcontnames[2]]))[[1]][["Pr(>F)"]][1]


    return(c(estimate = stats$results['Single_raters_absolute','ICC'],
             pval = pval))
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
    npca <- max(1,length(which(pve >= min_pve_pct_pc)))

    samplepcvals <- pca.res$x[, 1:npca, drop=FALSE]

    list(samplepcvals=samplepcvals, pve=pve)
}

# function to calculate correlation and plot
.calccompletecorandplot <- function(compare_data, covar_data,
                                   correlationtype, title,
                                   weights = NULL,
                                   exclude_vars_from_fdr=NULL,
                                   max_fdr = 0.1) {

    # require(plyr)

    # get factor and continuous covariates
    character_vars <- lapply(covar_data, class) == "character"
    covar_data[, character_vars] <- apply(covar_data[, character_vars, drop=FALSE], 1, as.factor)
    factorcovariates <- colnames(covar_data)[sapply(covar_data,is.factor)]
    contcovariates <- colnames(covar_data)[sapply(covar_data,is.numeric)]



    # calculate correlation between compare_data and factor covariates
    if (length(factorcovariates) > 0){
        comb <- expand.grid(colnames(compare_data),factorcovariates)
        factcont_cor <- apply(comb,1,
                              .getfactorcontassociationstatistics,
                              cbind(compare_data,covar_data[rownames(compare_data),factorcovariates, drop=FALSE]),
                              alpha=max_fdr)
        rownames(factcont_cor) = c("Estimate", "Pval")
        factcont_cor_vals <- matrix(factcont_cor['Estimate',],
                                    nrow = length(colnames(compare_data)),
                                    ncol = length(factorcovariates))
        factcont_cor_p <- matrix(factcont_cor['Pval',],
                                 nrow = length(colnames(compare_data)),
                                 ncol = length(factorcovariates))

        rownames(factcont_cor_vals) <- colnames(compare_data)
        colnames(factcont_cor_vals) <- factorcovariates

        rownames(factcont_cor_p) <- colnames(compare_data)
        colnames(factcont_cor_p) <- factorcovariates
    } else {
        factcont_cor_vals <- NULL
        factcont_cor_p <- NULL
    }

    # calculate correlation between compare_data and factor covariates
    if (length(contcovariates) > 0){
    cont_cor <- corr.test(compare_data,
                              covar_data[,contcovariates, drop=FALSE],
                              use='pairwise.complete.obs',
                              method=correlationtype,
                              adjust="none")
        cont_cor_vals <- cont_cor$r
        cont_cor_p <- cont_cor$p

        rownames(cont_cor_vals) <- colnames(compare_data)
        colnames(cont_cor_vals) <- contcovariates

        rownames(cont_cor_p) <- colnames(compare_data)
        colnames(cont_cor_p) <- contcovariates
    } else {
        cont_cor_vals <- NULL
        cont_cor_p <- NULL
    }

    all_cor_vals = cbind(factcont_cor_vals,cont_cor_vals)
    all_cor_p = cbind(factcont_cor_p,cont_cor_p)

    effects.significantcovars = all_cor_vals
    effects.significantcovars[all_cor_p>max_fdr] = 0
    effects.significantcovars = colSums(abs(effects.significantcovars)*replicate(dim(effects.significantcovars)[2],weights/sum(weights)))
    effects.significantcovars = effects.significantcovars[order(abs(effects.significantcovars),decreasing=TRUE)]

    cor_mat = melt(all_cor_p, varnames=c("compare", "covar"))
    colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"

    cor_mat$compare = factor(cor_mat$compare, levels=rownames(all_cor_p))
    cor_mat$covar = factor(cor_mat$covar, levels=colnames(all_cor_p))

    cor_mat$r = melt(all_cor_vals)$value
    cor_mat$fdr = p.adjust(cor_mat$pvalue, method="fdr")
    return(list(mat=cor_mat,
                effects.significantcovars = effects.significantcovars))
}

#' Find correlation between pcs and covariates
#'
#' This function will calculate the pcs using prcomp function,
#' and correlate categorical and numerical variables from
#' metadata.
#'
#' @author: Kenneth Daily and Thanneer Malai Perumal
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
    title = paste(ifelse(scale, "s", "un-s"), "caled ",
                  " data in pca; pve >= ",
                  min_pc_pct, "%; ", correlation,
                  " correlations ", sep="")
    message(paste("\nrunning pca and calculating correlations for:\n",
                  title, sep=""))

    metadata = as.data.frame(metadata)
    pcares <- .runpca(genesbysamples=counts,
                      scale_data_for_pca=scale,
                      min_pve_pct_pc=min_pc_pct)

    samplepcvals <- pcares$samplepcvals
    pve <- pcares$pve

    npca <- ncol(samplepcvals)

    colnames(samplepcvals) = paste(colnames(samplepcvals), " (",
                                   sprintf("%.2f", pve[1:npca]), "%)",
                                   sep="")

    # find covariates without any missing data
    samplesbyfullcovariates = metadata[, which(apply(metadata, 2,function(dat) all(!is.na(dat))))]

    exclude_vars_from_fdr = setdiff(colnames(metadata),
                                    colnames(samplesbyfullcovariates))

    corrres = .calccompletecorandplot(samplepcvals,
                                      samplesbyfullcovariates,
                                      correlation,
                                      title,
                                      weights = pve[1:dim(samplepcvals)[2]],
                                      exclude_vars_from_fdr)

    significantcovars = corrres$mat[corrres$mat$fdr<fdr,"covar"]
    ma = corrres$mat
    ma$r[ma$fdr > fdr] <- NA
    p = ggplot(ma, aes_(fill= ~r,x= ~covar,y= ~compare)) + 
        geom_tile() +
        theme_minimal() +
        ggtitle(title) +
        scale_fill_gradient2(low="darkblue", high="darkorange",
                             guide="colorbar", na.value = "grey90",
                             limits=c(-1,1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)
    invisible(list(significantcovars=significantcovars,
                plot=p,
                cor_matrix=corrres$mat,
                effects.significantcovars = corrres$effects.significantcovars,
                pcs_matrix=samplepcvals))
}


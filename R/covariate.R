# Find Inter Class Correlation between factor and continuous covariates
# Inspired from http://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables
.getFactorContAssociationStatistics <- function(factorContNames,COVARIATES,
                                               na.action='remove',
                                               alpha = 0.05){
    if (na.action == "remove")
        COVARIATES = na.omit(COVARIATES[,factorContNames])

    COVARIATES[,2] <- as.numeric(COVARIATES[,2])
    stats = ICC(COVARIATES[,factorContNames], alpha = alpha)
    Pval = summary(aov(COVARIATES[,factorContNames[1]]~COVARIATES[,factorContNames[2]]))[[1]][["Pr(>F)"]][1]


    return(c(Estimate = stats$results['Single_raters_absolute','ICC'],
             Pval = Pval))
}

# Function to run principal component analysis
.runPCA <- function(genesBySamples, SCALE_DATA_FOR_PCA = TRUE, MIN_PVE_PCT_PC = 1.0) {

    # estimate variance in data by PC:
    pca.res <- prcomp(t(genesBySamples), center=TRUE, scale.=SCALE_DATA_FOR_PCA, retx=TRUE)

    # examine how much variance is explained by PCs, and consider those with PVE >= (MIN_PVE_PCT_PC %):
    pc.var <- pca.res$sdev^2
    pve <- 100 * (pc.var / sum(pc.var))
    npca <- max(1,length(which(pve >= MIN_PVE_PCT_PC)))

    samplePCvals <- pca.res$x[, 1:npca, drop=FALSE]

    list(samplePCvals=samplePCvals, pve=pve)
}

# Function to calculate correlation and plot
.calcCompleteCorAndPlot <- function(COMPARE_data, COVAR_data,
                                   correlationType, title,
                                   WEIGHTS = NULL,
                                   EXCLUDE_VARS_FROM_FDR=NULL,
                                   MAX_FDR = 0.1) {

    # require(plyr)

    # Get factor and continuous covariates
    character_vars <- lapply(COVAR_data, class) == "character"
    COVAR_data[, character_vars] <- apply(COVAR_data[, character_vars, drop=FALSE], 1, as.factor)
    FactorCovariates <- colnames(COVAR_data)[sapply(COVAR_data,is.factor)]
    ContCovariates <- colnames(COVAR_data)[sapply(COVAR_data,is.numeric)]



    # Calculate correlation between compare_data and factor covariates
    if (length(FactorCovariates) > 0){
        comb <- expand.grid(colnames(COMPARE_data),FactorCovariates)
        factCont_cor <- apply(comb,1,
                              .getFactorContAssociationStatistics,
                              cbind(COMPARE_data,COVAR_data[rownames(COMPARE_data),FactorCovariates, drop=FALSE]),
                              alpha=MAX_FDR)
        rownames(factCont_cor) = c("Estimate", "Pval")
        factCont_cor_vals <- matrix(factCont_cor['Estimate',],
                                    nrow = length(colnames(COMPARE_data)),
                                    ncol = length(FactorCovariates))
        factCont_cor_p <- matrix(factCont_cor['Pval',],
                                 nrow = length(colnames(COMPARE_data)),
                                 ncol = length(FactorCovariates))

        rownames(factCont_cor_vals) <- colnames(COMPARE_data)
        colnames(factCont_cor_vals) <- FactorCovariates

        rownames(factCont_cor_p) <- colnames(COMPARE_data)
        colnames(factCont_cor_p) <- FactorCovariates
    } else {
        factCont_cor_vals <- NULL
        factCont_cor_p <- NULL
    }

    # Calculate correlation between compare_data and factor covariates
    if (length(ContCovariates) > 0){
    cont_cor <- corr.test(COMPARE_data,
                              COVAR_data[,ContCovariates, drop=FALSE],
                              use='pairwise.complete.obs',
                              method=correlationType,
                              adjust="none")
        cont_cor_vals <- cont_cor$r
        cont_cor_p <- cont_cor$p

        rownames(cont_cor_vals) <- colnames(COMPARE_data)
        colnames(cont_cor_vals) <- ContCovariates

        rownames(cont_cor_p) <- colnames(COMPARE_data)
        colnames(cont_cor_p) <- ContCovariates
    } else {
        cont_cor_vals <- NULL
        cont_cor_p <- NULL
    }

    all_cor_vals = cbind(factCont_cor_vals,cont_cor_vals)
    all_cor_p = cbind(factCont_cor_p,cont_cor_p)

    Effects.significantCovars = all_cor_vals
    Effects.significantCovars[all_cor_p>MAX_FDR] = 0
    Effects.significantCovars = colSums(abs(Effects.significantCovars)*replicate(dim(Effects.significantCovars)[2],WEIGHTS/sum(WEIGHTS)))
    Effects.significantCovars = Effects.significantCovars[order(abs(Effects.significantCovars),decreasing=T)]

    cor_mat = melt(all_cor_p, varnames=c("COMPARE", "COVAR"))
    colnames(cor_mat)[colnames(cor_mat) == "value"] = "pvalue"

    cor_mat$COMPARE = factor(cor_mat$COMPARE, levels=rownames(all_cor_p))
    cor_mat$COVAR = factor(cor_mat$COVAR, levels=colnames(all_cor_p))

    cor_mat$r = melt(all_cor_vals)$value
    cor_mat$fdr = p.adjust(cor_mat$pvalue, method="fdr")
    return(list(mat=cor_mat,
                Effects.significantCovars = Effects.significantCovars))
}

#' Find correlation between PCs and covariates
#'
#' This function will calculate the PCs using prcomp function,
#' and correlate categorical and numerical variables from
#' metadata.
#'
#' @author: Kenneth Daily and Thanneer Malai Perumal
#'
#' @param counts normalized counts matrix
#' @param metadata data.frame with samples metadata.
#' @param fdr numeric value to use as cutoff to determine
#' the minimum FDR to consider significant correlations
#' between PCs and covariates.
#' @param scale boolean to determine wether counts matrix should be
#' scaled for PCA. Default FALSE.
#' @param min_pc_pct numeric value that will be used as cutoff
#' to select only PCs that explain more variability than this.
#' @param correlation character determining the method for the
#' correlation between PCs and covariates.
#'
#' @return: list:
#' a)significantCovars,
#' b)plot,
#' c)cor_matrix,
#' d)Effects.significantCovars: that is PCs pct * absolute correlation between covariate and PCs,
#' e)PCs_matrix: PCs loading for each sample
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
    title = paste(ifelse(scale, "S", "Un-s"), "caled ",
                  " data in PCA; PVE >= ",
                  min_pc_pct, "%; ", correlation,
                  " correlations ", sep="")
    message(paste("\nRunning PCA and calculating correlations for:\n",
                  title, sep=""))

    metadata = as.data.frame(metadata)
    pcaRes <- .runPCA(genesBySamples=counts,
                     SCALE_DATA_FOR_PCA=scale,
                     MIN_PVE_PCT_PC=min_pc_pct)

    samplePCvals <- pcaRes$samplePCvals
    pve <- pcaRes$pve

    npca <- ncol(samplePCvals)

    colnames(samplePCvals) = paste(colnames(samplePCvals), " (",
                                   sprintf("%.2f", pve[1:npca]), "%)",
                                   sep="")

    # Find covariates without any missing data
    samplesByFullCovariates = metadata[, which(apply(metadata, 2,function(dat) all(!is.na(dat))))]

    EXCLUDE_VARS_FROM_FDR = setdiff(colnames(metadata),
                                    colnames(samplesByFullCovariates))

    corrRes = .calcCompleteCorAndPlot(samplePCvals,
                                     samplesByFullCovariates,
                                     correlation,
                                     title,
                                     WEIGHTS = pve[1:dim(samplePCvals)[2]],
                                     EXCLUDE_VARS_FROM_FDR)

    significantCovars = corrRes$mat[corrRes$mat$fdr<fdr,"COVAR"]
    ma = corrRes$mat
    ma$r[ma$fdr > fdr] <- NA
    p = ggplot(ma, aes_(fill= ~r,x= ~COVAR,y= ~COMPARE)) + geom_tile() +
        theme_minimal() +
        ggtitle(title) +
        scale_fill_gradient2(low="darkblue", high="darkorange",
                             guide="colorbar", na.value = "grey90",
                             limits=c(-1,1)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return(list(significantCovars=significantCovars,
                plot=p,
                cor_matrix=corrRes$mat,
                Effects.significantCovars = corrRes$Effects.significantCovars,
                PCs_matrix=samplePCvals))
}


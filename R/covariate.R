# find inter class correlation between factor and continuous covariates
# inspired from http://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables
.getfactorcontassociationstatistics <- function(factorcontnames,
                                                covariates,
                                                na.action='remove',
                                                alpha = 0.05){
    if (na.action == "remove")
        covariates = na.omit(covariates[, factorcontnames])

    covariates[,2] <- as.numeric(covariates[,2])
    data <- list(x = covariates[, 1L],
                y = covariates[, 2L])
    stats <- cor.test(scale(data[[1L]]),
                     scale(data[[2L]]),
                     method = "kendall", exact = FALSE)
    cor <- stats[["estimate"]]
    pval <- stats[["p.value"]]

    return(c(estimate = cor, pval = pval))
}

# function to run principal component analysis
.runpca <- function(genesbysamples, scale_data_for_pca = TRUE,
                    min_pve_pct_pc = 1.0) {

    # estimate variance in data by pc:
    pca.res <- prcomp(t(genesbysamples), center = TRUE,
                      scale. = scale_data_for_pca, retx = TRUE)

    # examine how much variance is explained by pcs,
    # and consider those with pve >= (min_pve_pct_pc %):
    pc.var <- pca.res$sdev^2L
    pve <- 100L * (pc.var / sum(pc.var))
    npca <- max(1L, length(which(pve >= min_pve_pct_pc)))

    samplepcvals <- pca.res$x[, 1L:npca, drop = FALSE]

    list(samplepcvals = samplepcvals, pve = pve)
}

# function to calculate correlation and plot
.calccompletecorandplot <- function(compare_data, covar_data,
                                    correlationtype, title,
                                    weights = NULL,
                                    exclude_vars_from_fdr=NULL,
                                    max_fdr = 0.1) {
    # get factor and continuous covariates
    character_vars <- lapply(covar_data, class) == "character"
    if (sum(character_vars) > 0 )
        covar_data[, character_vars] <- apply(covar_data[, character_vars,
                                                         drop = FALSE],
                                              1L, as.factor)

    factorcovariates <- select_if(covar_data, is.factor) %>% colnames
    contcovariates <- select_if(covar_data, is.numeric) %>% colnames
    
    all_covariates <- cbind(covar_data[, contcovariates, drop = FALSE],
                            covar_data[, factorcovariates, drop = FALSE] %>%
                                mutate_all(as.numeric))
    
    cov_cor <- corr.test(compare_data,
                         all_covariates,
                         use = 'pairwise.complete.obs',
                         method = correlationtype,
                         adjust = "none")
    all_cor_vals <- cov_cor[["r"]]
    all_cor_p <- cov_cor[["p"]]
    
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

    cor_mat[["compare"]] <- factor(cor_mat[["compare"]],
                                   levels = rownames(all_cor_p))
    cor_mat[["covar"]] <- factor(cor_mat[["covar"]],
                                 levels = colnames(all_cor_p))

    cor_mat[["r"]] <- melt(all_cor_vals)[["value"]]
    cor_mat[["fdr"]] <- p.adjust(cor_mat[["pvalue"]], method = "fdr")
    return(list(mat = cor_mat,
                effects.significantcovars = effects.significantcovars))
}

.generate_scatter_plot <- function(metadata, corMat){
    if (sum(corMat[["fdr"]] < 0.1, na.rm = T) == 0) return(NULL)
    plist <- apply(corMat[corMat[["fdr"]] < 0.1, ], 1, function(row){
        xs <- strsplit(row[1], " ")[[1]][1]
        ys <- row[2]
        ggplot(metadata, aes_string(x = xs, y = ys)) +
            geom_point() +
            ggtitle(paste(xs, ys))
    })
    names <- corMat[corMat[["fdr"]] < 0.1, ] %>%
        mutate(name = paste(compare, covar, sep = ":")) %>% 
        .[["name"]]
    names(plist) <- names
    plist
}

#' Find correlation between pcs and covariates
#'
#' This function will calculate the pcs using prcomp function,
#' and correlate categorical and numerical variables from
#' metadata. The size of the dots indicates the importance of the
#' metadata, for instance, when the range of the values is pretty
#' small (from 0.001 to 0.002 in ribosimal content),
#' the correlation results is not important. See details to know
#' how this is calculated.
#'
#' @author: Lorena Pantano, Kenneth Daily and Thanneer Malai Perumal
#'
#' @param counts normalized counts matrix
#' @param metadata data.frame with samples metadata.
#' @param fdr numeric value to use as cutoff to determine
#'   the minimum fdr to consider significant correlations
#'   between pcs and covariates.
#' @param scale boolean to determine wether counts matrix should be
#'   scaled for pca. default FALSE.
#' @param minPC numeric value that will be used as cutoff
#'   to select only pcs that explain more variability than this.
#' @param correlation character determining the method for the
#'   correlation between pcs and covariates.
#' @param addCovDen boolean. Whether to add the covariates
#'   dendograme to the plot to see covariates relationship.
#'   It will show [degCorCov()] dendograme on top of the columns of
#'   the heatmap.
#' @param legacy boolean. Whether to plot the legacy version.
#' @param plot Whether to plot or not the correlation matrix.
#' @details This method is adapeted from Daily et al 2017 article.
#'   Principal components from PCA analysis are correlated with 
#'   covariates metadata. Factors are transformed to numeric variables.
#'   Correlation is measured by cor.test function with Kendall method
#'   by default.
#'   
#'   The size of the dot, or effect size, indicates the importance of 
#'   the covariate based on the range of the values. Covariates
#'   where the range is very small (like a % of mapped reads that
#'   varies between 0.001 to 0.002) will have a very small dot size.
#'   To get to this value, each covariate is normalized using this 
#'   equation: v/max(v). Then the standard deviation is calculated,
#'   and the minimum and maximum values are set to
#'   1 and 3 respectively.
#'   
#' @references 
#' Daily, K. et al.  Molecular, phenotypic, and sample-associated data to describe pluripotent stem cell lines and derivatives. Sci Data 4, 170030 (2017).
#'  
#' @return: list:
#' a) significantCovars, covariates with FDR below the cutoff.
#' b) plot, heatmap of the correlation found. * means pvalue < 0.05.
#'    Only variables with FDR value lower than the cutoff are colored.
#' c) corMatrix, correlation, p-value, FDR values
#'    for each covariate and PCA pais
#' d) effectsSignificantcovars: that is PCs % * absolute
#'    correlation between covariate and PCs,
#' e) pcsMatrix: PCs loading for each sample
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' res <- degCovariates(log2(counts(dse)+0.5), colData(dse))
#' res <- degCovariates(log2(counts(dse)+0.5),
#'   colData(dse), legacy = TRUE)
#' res$plot
#' res$scatterPlot[[1]]
#' @export
degCovariates <- function(counts, metadata,
                          fdr = 0.1,
                          scale = FALSE,
                          minPC = 5.0,
                          correlation = "kendall",
                          addCovDen = TRUE,
                          legacy = FALSE,
                          plot = TRUE) {
    title <- paste(ifelse(scale, "s", "un-s"), "caled ",
                   " data in pca;\npve >= ",
                   minPC, "%;\n", correlation,
                   " cor ", sep = "")
    message(paste("\nrunning pca and calculating correlations for:\n",
                  title, sep = ""))
    metadata <- as.data.frame(metadata)
    covar_class <- sapply(metadata[1,], class)
    
    metadata <- degClean(metadata) %>%
        mutate_all(as.numeric) %>% 
        as.data.frame() %>% 
        set_rownames(row.names(metadata))

    stopifnot(identical(colnames(counts), rownames(metadata)))
    
    pcares <- .runpca(genesbysamples = counts,
                      scale_data_for_pca = scale,
                      min_pve_pct_pc = minPC)
    
    samplepcvals <- pcares[["samplepcvals"]]
    pve <- pcares[["pve"]]
    npca <- ncol(samplepcvals)
    original_names <- colnames(samplepcvals)
    colnames(samplepcvals) <- paste(colnames(samplepcvals), " (",
                                    sprintf("%.2f", pve[1L:npca]), "%)",
                                    sep = "")
    
    # find covariates without any missing data
    samplesbyfullcovariates <- metadata[, which(apply(metadata, 2L,
                                                      function(dat) all(!is.na(dat)))), drop = FALSE]
    covar_class <- covar_class[colnames(samplesbyfullcovariates)]
    
    exclude_vars_from_fdr <- setdiff(colnames(metadata),
                                     colnames(samplesbyfullcovariates))
    
    covar_factors <- samplesbyfullcovariates[,names(covar_class)[covar_class != "numeric"], drop = FALSE]
    covar_numeric <- samplesbyfullcovariates[,names(covar_class)[covar_class == "numeric"]]
    samplesbyfullcovariates = cbind(covar_factors, covar_numeric)
    
    corrRes <- .calccompletecorandplot(samplepcvals,
                                       samplesbyfullcovariates,
                                       correlation,
                                       title,
                                       weights =
                                           pve[1L:dim(samplepcvals)[2L]],
                                       exclude_vars_from_fdr)
    
    significantcovars <- corrRes[["mat"]][corrRes[["mat"]][["fdr"]] < fdr,"covar"]
    ma <- corrRes[["mat"]]
    ma[["r"]][ma[["fdr"]] > fdr] <- NA
    corMa <- ma[, c("r", "compare", "covar")] %>%
        spread(!!sym("compare"), !!sym("r")) %>%
        remove_rownames() %>%
        column_to_rownames("covar")
    
    cor_meta <- .calccompletecorandplot(samplesbyfullcovariates,
                                        samplesbyfullcovariates,
                                        correlation,
                                        "",
                                        weights = 1L)
    
    corMeta <- cor_meta[["mat"]][, c("r", "compare", "covar")] %>%
        spread(!!sym("compare"), !!sym("r"), fill = 0) %>%
        remove_rownames() %>%
        column_to_rownames("covar") 
    
    hc <-  hclust(as.dist((1-corMeta)^2),
                method = "ward.D")
    ma[["covar"]] = as.character(ma[["covar"]])
    ma_sd <- bind_rows(
        data.frame(covar = colnames(covar_numeric),
                   effect_size = apply(covar_numeric,
                                       2,
                                       function(v) {
                                           sd(v/max(v))
                                       }),
                   stringsAsFactors = FALSE),
        data.frame(covar = colnames(covar_factors),
                   effect_size = apply(covar_factors,
                                       2,
                                       sd),
                   stringsAsFactors = FALSE),
    )
    ma <- left_join(ma,
                    ma_sd,
              by = "covar")
    ma[["effect_size"]][ma[["effect_size"]] < 0.05] <- 0.05
    ma[["effect_size"]][ma[["effect_size"]] > 3] <- 3

    if (legacy){
        if (addCovDen){
            p <- Heatmap(t(corMa),
                         name = "cor",
                         row_title = title,
                         cluster_rows = FALSE,
                         cluster_columns = hclust(as.dist((1-corMeta)^2),
                                                  method = "ward.D"),
                         col = colorRamp2(c(1, 0, -1),
                                          c("darkorange", "white", "darkblue")))
        }else{
            p <- ggplot(ma, aes_(fill = ~r,x = ~covar,y = ~compare)) +
                geom_tile() +
                theme_minimal() +
                ggtitle(title) +
                scale_fill_gradient2(low = "darkblue", high = "darkorange",
                                     guide = "colorbar", na.value = "grey90",
                                     limits = c(-1L, 1L)) +
                theme(axis.text.x = element_text(angle = 90L,
                                                 hjust = 1L,
                                                 vjust = 0.5))
            if (sum(ma[["r"]] > 0, na.rm = TRUE))
                p <- p + geom_text(data = ma[ma[["pvalue"]] < 0.05, ], aes(label="*"))
        }
        
    }else{
        
        dhc <- as.dendrogram(hc)
        # Rectangular lines
        ddata <- dendro_data(dhc, type = "rectangle")
        dn = ggplot(ddata[["segments"]]) + 
            geom_segment(aes(x = x, y = -y, xend = xend, yend = -yend)) + 
            scale_y_reverse(expand = c(0, 0)) +
            scale_x_continuous(expand=c(0.01, 0.01)) +
            xlab("") +
            ylab("") +
            theme_minimal() +
            theme(axis.text.x = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
        
        ma[["covar"]] <- factor(ma[["covar"]], levels = hc$labels[hc$order])
        tile = ggplot(ma, aes_(x = ~covar,y = ~compare,
                               size = ~effect_size, color = ~r)) +
            geom_point() +
            theme_minimal() +
            
            scale_size_continuous(range = c(0.01, 3)) + 
            scale_color_gradient2(low = "darkblue", high = "darkorange",
                                  guide = "colorbar", na.value = "grey90",
                                  limits = c(-1L, 1L)) +
            xlab("") +
            theme(axis.text.x = element_text(angle = 90L,
                                             hjust = 1L,
                                             vjust = 0.5),
                  legend.position = "bottom") +
            scale_x_discrete(expand=c(0.01,0.01))
        
        if (addCovDen){
            p <- plot_grid(
                dn, tile, align = "v",
                nrow = 2, rel_heights = c(1,4)
            ) + ggtitle(title)
        }else{
            p <- tile  + ggtitle(title)
        }
    }

    samplepcvals <- as.data.frame(samplepcvals) %>% 
        set_colnames(original_names) %>% 
        rownames_to_column("samples")
    samplepcvals <- bind_cols(samplepcvals,
                              metadata)
    scatterPlot <- .generate_scatter_plot(samplepcvals, corrRes[["mat"]])
    
    if (plot) print(p)
    
    invisible(list(significantCovars = significantcovars,
                plot = p,
                corMatrix = ma,
                pcsMatrix = samplepcvals,
                scatterPlot = scatterPlot,
                effectsSignificantCovars =
                    corrRes[["effects.significantcovars"]]
                ))
}

degClean <- function(ma){
    ma <- ma[,colSums(is.na(ma))<nrow(ma)]
    lapply(ma, function(x) {
        if (length(unique(x)) < length(x) * 0.20 & is.numeric(x))
            x <- as.factor(x)
        if ((class(x) %in% c("factor", "character"))) {
            .f = as.factor(x)
            if (length(levels(.f)) < 2L)
                return(NULL)
            if (length(levels(.f)) > length(.f) * 0.80)
                return(NULL)
            return(.f)
        }
        if (sd(x, na.rm = TRUE) == 0)
            return(NULL)
        as.numeric(x)
    })  %>% Filter(Negate(is.null), .) %>% bind_cols %>%
        as.data.frame() %>%
        set_rownames(row.names(ma))
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
#'   the minimum fdr to consider significant correlations
#'   between pcs and covariates.
#' @param ... Parameters to pass to [ComplexHeatmap::Heatmap()].
#'
#' @return: list:
#' a) cor, data.frame with pair-wise correlations, pvalues, FDR
#' b) corMat, data.frame with correlation matrix
#' c) fdrMat, data.frame with FDR matrix
#' b) plot, Heatmap plot of correlation matrix
#' @examples
#' data(humanGender)
#' library(DESeq2)
#' idx <- c(1:10, 75:85)
#' dse <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
#'   colData(humanGender)[idx,], design=~group)
#' cor <- degCorCov(colData(dse))
#' @export
degCorCov <- function(metadata, fdr=0.05, ...){
    clean <- degClean(metadata) %>%
        mutate_all(as.numeric)
    cor <- .calccompletecorandplot(clean,
                                      clean,
                                      "kendall",
                                      "",
                                      weights = 1L)
    
    corMat <- cor[["mat"]][, c("r", "compare", "covar")] %>%
        spread(!!sym("compare"), !!sym("r")) %>% remove_rownames() %>%
        column_to_rownames("covar")
    fdrMat <- cor[["mat"]][, c("fdr", "compare", "covar")] %>%
        spread(compare, fdr) %>% remove_rownames() %>%
        column_to_rownames("covar")
    corMat[fdrMat > fdr] <- 0
    # corMat[fdrMat > 0.05] <- NA
    cols <- colorRampPalette(c("steelblue", "white", "orange"))(10)
    if (sum(!is.na(corMat)) > 2) {
        p <- Heatmap(corMat, col = cols, name = "cor-value",
                     row_dend_reorder = FALSE,
                     column_dend_reorder = FALSE, ...)
        print(p)
    }
    invisible(list(cor = cor[["mat"]], corMat = corMat,
                   fdrMat = fdrMat, plot = p))

}
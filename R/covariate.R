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
# cat vs continous
# ICCest(Chick, weight, data = ChickWeight, CI.type = "S")
# broom::tidy(aov(weight~Chick, data = ChickWeight))
# data("Arthritis") # cat vs cat
# tab <- xtabs(~Improved + Treatment, data = Arthritis)
# summary(assocstats(tab))
# 
# assocstats(UCBAdmissions)
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
    #effects.significantcovars[all_cor_p > max_fdr] <- 0
    effects.significantcovars <- colSums(abs(effects.significantcovars) * replicate(dim(effects.significantcovars)[2L], weights / sum(weights)))
    effects.significantcovars <- effects.significantcovars[order(abs(effects.significantcovars), decreasing = TRUE)]

    cor_mat <- melt(all_cor_p, varnames = c("compare", "covar"))
    colnames(cor_mat)[colnames(cor_mat) == "value"] <- "pvalue"

    cor_mat[["compare"]] <- factor(cor_mat[["compare"]],
                                   levels = rownames(all_cor_p))
    cor_mat[["covar"]] <- factor(cor_mat[["covar"]],
                                 levels = colnames(all_cor_p))

    # cor_mat[["r"]] <- melt(all_cor_vals)[["value"]]
    cor_mat[["r"]] <- as.vector(all_cor_vals)
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

.numeric_effect_size <- function(covar_numeric, smart = TRUE){
    .smart <- function(v) {
        if (min(v) >= 0 & max(v)<=1)
            return(max(v)-min(v))
        if (min(v) >= 0 & max(v)<=100)
            return(max(v)-min(v))
        1 - min(v/max(v))
    }
    apply(covar_numeric,
          2, function(v){
              if (smart)
                  return(.smart(v))
              max(v)-min(v)
          }
    )
}

.effect_size <- function(ma, covar_numeric, covar_factors,
                         smart = TRUE){
    ma_sd <- data.frame(covar = colnames(covar_numeric),
                   effect_size = .numeric_effect_size(covar_numeric,
                                                      smart),
                   stringsAsFactors = FALSE)
    if (ncol(covar_factors) > 0){
        ma_sd <- bind_rows(ma_sd,
                           data.frame(covar = colnames(covar_factors),
                   effect_size = 1,
                   stringsAsFactors = FALSE)
        )
    }

    ma <- left_join(ma, ma_sd, by = "covar")
    ma[["effect_size"]][ma[["effect_size"]] < 0.01] <- 0.01
    ma[["effect_size"]][ma[["effect_size"]] > 1] <- 1
    ma[["type_variable"]] <- "categorical"
    ma[["type_variable"]][ma[["covar"]]  %in% colnames(covar_numeric)] <- "numeric"
    ma
}

.model <- function(data, method = "lm"){
    if (method == "lm"){
        pc_sig <- lm(PC~., data=data) %>% 
            broom::tidy()
    }
    # else if (method == "lasso"){
    #     pc_sig <- lm.lasso <- l1ce(PC ~ 0 + ., data=data, sweep.out = NULL) %>% 
    #         summary() %>% 
    #         .[["coefficients"]] %>% 
    #         broom::tidy()
    # }
    names(pc_sig) <- c("term", "estimate", "std.error", "statistic", "p.value")
    return(pc_sig)
}

# reduce covariates to significant ones that predict PCs
.reduce_covariates <- function(corMatrix, pcsMatrix, method = "lm"){
    pcs <- colnames(pcsMatrix)[grepl("PC[0-9]+", colnames(pcsMatrix))]
    significants <- lapply(pcs, function(pc){
        pc_var <- corMatrix %>% 
            filter(fdr < 0.01, grepl(pc, compare)) %>% 
            arrange(abs(r)) %>% 
            .[["covar"]] %>% 
            as.character()
        if (length(pc_var) == 0)
            return(NULL)
        data <- pcsMatrix[,c(pc, pc_var)]
        colnames(data)[1] = "PC"
        data[,2:ncol(data)] <- apply(data[,2:ncol(data), drop = FALSE], 2, scale)
        pc_sig <- .model(data, method)
        pc_sig[["PC"]] = pc
        pc_sig %>% filter(p.value < 0.05)
    }) %>% bind_rows()

    if (nrow(significants) == 0)
        return(data.frame(estimate=0, p.value=0, PC="1", term="noterm"))
    significants[,c("estimate", "p.value", "PC", "term")]
}

#' Find correlation between pcs and covariates
#'
#' This function will calculate the pcs using prcomp function,
#' and correlate categorical and numerical variables from
#' metadata. The size of the dots indicates the importance of the
#' metadata, for instance, when the range of the values is pretty
#' small (from 0.001 to 0.002 in ribosimal content),
#' the correlation results is not important. If black stroke lines
#' are shown, the correlation analysis has a FDR < 0.05 for that
#' variable and PC.
#' Only significant variables according the linear model are colored.
#' See details to know
#' how this is calculated. 
#' 
#' @author: Lorena Pantano, Victor Barrera, 
#'          Kenneth Daily and Thanneer Malai Perumal
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
#' @param smart boolean. Whether to avoid normalization of the
#'   numeric covariates when calculating importance. This is not
#'   used if `legacy = TRUE`. See @details for more information.
#' @param method character. Whether to use `lm` to
#'   calculate the significance of the variable during reduction
#'   step. See @details for more information.
#' @param plot Whether to plot or not the correlation matrix.
#' @details This method is adapeted from Daily et al 2017 article.
#'   Principal components from PCA analysis are correlated with 
#'   covariates metadata. Factors are transformed to numeric variables.
#'   Correlation is measured by `cor.test` function with Kendall method
#'   by default.
#'   
#'   The size of the dot, or importance, indicates the importance of 
#'   the covariate based on the range of the values. Covariates
#'   where the range is very small (like a % of mapped reads that
#'   varies between 0.001 to 0.002) will have a very small size (0.1*max_size).
#'   The maximum value is set to 5 units.
#'   To get to importance, each covariate is normalized using this 
#'   equation: `1 - min(v/max(v))`,
#'   and the minimum and maximum values are set to
#'   0.01 and 1 respectively. For instance, 0.5 would mean there is at least
#'   50% of difference between the minimum value and the maximum value.
#'   Categorical variables are plot using the maximum size always, since
#'   it is not possible to estimate the variability. By default, it 
#'   won't do `v/max(v)` if the values are already between 0-1 or
#'   0-100 (already normalized values as rates and percentages).
#'   If you want to ignore the importance, use `legacy = TRUE`.
#'   
#'   Finally, a linear model is used to calculate the significance
#'   of the covariates effect on the PCs. For that, this function
#'   uses `lm` to regress the data and uses the p-value calculated by
#'   each variable in the model to define significance (pvalue < 0.05).
#'   Variables with a black stroke
#'   are significant after this step. Variables with grey stroke 
#'   are significant at the first pass considering p.value < 0.05
#'   for the correlation analysis.
#' @references 
#' Daily, K. et al.  Molecular, phenotypic, and sample-associated data to describe pluripotent stem cell lines and derivatives. Sci Data 4, 170030 (2017).
#'  
#' @return: list:
#' 
#' * plot, heatmap showing the signifcance of the variables. 
#' * corMatrix, correlation, p-value, FDR values
#'   for each covariate and PCA pais
#' * pcsMatrix: PCs loading for each sample
#' * scatterPlot: plot for each significant covariate and
#'   the PC values.
#' * significants: contains the significant covariates
#'   using a linear model to predict the coefficient
#'   of covariates that have some color in the plot.
#'   All the significant covariates from the liner model analysis
#'   are returned.
#' 
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
                          smart = TRUE,
                          method = "lm",
                          plot = TRUE) {
    title <- paste(ifelse(scale, "s", "un-s"), "caled ",
                   " data in pca;\npve >= ",
                   minPC, "%;\n", correlation,
                   " cor ", sep = "")
    message(paste("\nrunning pca and calculating correlations for:\n",
                  title, sep = ""))
    metadata <- as.data.frame(metadata)
    
    metadata <- degClean(metadata) 
    
    covar_class <- sapply(metadata[1,], class)

    metadata <- metadata %>%
        mutate_all(as.numeric) %>% 
        as.data.frame() %>% 
        set_rownames(row.names(metadata))
    
    stopifnot(identical(colnames(counts), rownames(metadata)))
    
    pcares <- .runpca(genesbysamples = counts,
                      scale_data_for_pca = scale,
                      min_pve_pct_pc = minPC)
    
    samplepcvals <- pcares[["samplepcvals"]]
    pve <- pcares[["pve"]]
    original_names <- colnames(samplepcvals)
    pc_pct <- data.frame(
        pc = colnames(samplepcvals),
        pct = paste(" (",
                    sprintf("%.2f", pve[1L:ncol(samplepcvals)]), "%)",
                    sep = ""),
        stringsAsFactors = FALSE
        )

        # find covariates without any missing data
    samplesbyfullcovariates <- metadata[, which(apply(metadata, 2L,
                                                      function(dat) all(!is.na(dat)))), drop = FALSE]
    covar_class <- covar_class[colnames(samplesbyfullcovariates)]
    
    exclude_vars_from_fdr <- setdiff(colnames(metadata),
                                     colnames(samplesbyfullcovariates))
    
    covar_factors <- samplesbyfullcovariates[,names(covar_class)[covar_class != "numeric"], drop = FALSE]
    covar_numeric <- samplesbyfullcovariates[,names(covar_class)[covar_class == "numeric"], drop = FALSE]
    
    samplesbyfullcovariates = cbind(covar_factors, covar_numeric)
    
    corrRes <- .calccompletecorandplot(samplepcvals,
                                       samplesbyfullcovariates,
                                       correlation,
                                       title,
                                       weights =
                                           pve[1L:dim(samplepcvals)[2L]],
                                       exclude_vars_from_fdr)
    
    ma <- corrRes[["mat"]]
    ma[["r"]][ma[["fdr"]] > fdr] <- NA
    ma[["fdr"]][ma[["fdr"]] > fdr] <- NA
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
    ma[["compare"]] = as.character(ma[["compare"]])

    ma <- .effect_size(ma, covar_numeric, covar_factors)
    
    samplepcvals <- as.data.frame(samplepcvals) %>% 
        set_colnames(original_names) %>% 
        rownames_to_column("samples")
    samplepcvals <- bind_cols(samplepcvals,
                              metadata)
    significants <- .reduce_covariates(ma, samplepcvals, method)
    scatterPlot <- .generate_scatter_plot(samplepcvals, corrRes[["mat"]])
    
    ma_plot <- left_join(ma,
                         significants,
              by = c("compare" = "PC", 
                     "covar" = "term")) %>% 
        left_join(pc_pct, by = c("compare" = "pc")) %>% 
        unite(col = "compare", !!sym("compare"), !!sym("pct"),
              sep = " ")
    ma_plot[["r"]][is.na(ma_plot[["p.value"]])] = NA
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
        
        
        ma_plot[["covar"]] <- factor(ma_plot[["covar"]],
                                levels = hc$labels[hc$order])

        tile <- ggplot(ma_plot, aes_(x = ~covar,y = ~compare,
                             size = ~effect_size,
                             color = ~r,
                             fill = ~r,
                             shape = ~type_variable)) +
            geom_point() +
            geom_point(data = filter(ma_plot, !is.na(fdr)),
                       stroke = 1, color="black") +
            geom_point(data = filter(ma_plot, pvalue<0.05),
                       stroke = 1, color="grey60") +
            theme_minimal() +
            scale_shape_manual(values = c(22, 21)) + 
            scale_size_continuous(name = "importance",
                                  limits = c(0.01, 1),
                                  range = c(0.01, 5)) + 
            scale_color_gradient2(low = "darkblue", high = "darkorange",
                                  guide = "colorbar", na.value = "grey90",
                                  limits = c(-1L, 1L)) +
            scale_fill_gradient2(low = "darkblue", high = "darkorange",
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

    
    if (plot) print(p)

    invisible(list(
                plot = p,
                corMatrix = ma,
                pcsMatrix = samplepcvals,
                scatterPlot = scatterPlot,
                significants = significants
                ))
}

degClean <- function(ma){
    ma <- ma[,colSums(is.na(ma))<nrow(ma)]
    lapply(ma, function(x) {
        if (length(unique(x)) < length(x) * 0.20 & is.numeric(x))
            x <- as.factor(x)
        if ((class(x)[1] %in% c("factor", "character"))) {
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
#' @param use_pval boolean to indicate to use p-value instead of FDR to
#'   hide non-significant correlation.
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
degCorCov <- function(metadata, fdr=0.05, use_pval = FALSE, ...){
    clean <- degClean(metadata) %>%
        mutate_all(as.numeric)
    cor <- .calccompletecorandplot(clean,
                                   clean,
                                   "kendall",
                                   "",
                                   max_fdr = fdr,
                                   weights = 1L)
    #browser()
    if (use_pval){
        cor[["mat"]]['fdr'] <- cor[["mat"]]['pvalue']
    }
    corMat <- cor[["mat"]][, c("r", "compare", "covar")] %>%
        spread(!!sym("compare"), !!sym("r")) %>% remove_rownames() %>%
        column_to_rownames("covar")
    fdrMat <- cor[["mat"]][, c("fdr", "compare", "covar")] %>%
        spread(compare, fdr) %>% remove_rownames() %>%
        column_to_rownames("covar")

    corMat[fdrMat > fdr] <- 0
    # corMat[fdrMat > 0.05] <- NA
    if (sum(!is.na(corMat)) > 2) {
        p <- Heatmap(corMat, name = "cor-value",
                     col = colorRamp2(c(-1, 0, 1),
                                      c("steelblue", "white", "orange")),
                     row_dend_reorder = FALSE,
                     column_dend_reorder = FALSE, ...)
        print(p)
    }
    invisible(list(cor = cor[["mat"]], corMat = corMat,
                   fdrMat = fdrMat, plot = p))

}
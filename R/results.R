.guessResults <- function(object, what, alpha){
    coef <- match(what[[1]], resultsNames(object))
    if (is.na(coef) & length(what) == 1L)
        stop("coef ", what, " not found in resultsNames(dds).")
    if (length(what) == 1)
        res <- results(object, name = what, alpha = alpha)
    else res <- results(object, contrast = what, alpha = alpha)
    res[order(res[["padj"]]),]
}

.guessShrunken <- function(object, what, unShrunken){
    coef <- match(what[[1]], resultsNames(object))
    if (is.na(coef) & length(what) == 1L)
        stop("coef ", what, " not found in resultsNames().")
    if (length(what) == 1)
        res <- lfcShrink(
            dds = object,
            coef = coef,
            res = unShrunken)
    else res <- lfcShrink(object, contrast = what, res = unShrunken)
    res[order(res[["padj"]]),]
}

.createComb <- function(dds, combs){
    lapply(combs, function(comb){
        colData(dds)[[comb]] %>%
            unique %>%
            combn(., 2, simplify = FALSE) %>%
            lapply(., function(x) c(comb, x))
    }) %>% unlist(., recursive = FALSE)
}

.normalizeNames <- function(combs, dds) {
    lapply(combs, function(x) {
        if (grepl("^[0-9]+$", x))
            return(colnames(colData(dds))[[as.numeric(x)]])
        else return(x)
    }) %>% unlist(., recursive = FALSE)
}

.guessComb <- function(dds, combs, contrast, pairs = FALSE){
    all_combs <- list()
    if (!is.null(combs)){
        if (length(combs) > 1)
            combs <- .normalizeNames(combs, dds)
        if (sum(combs %in% names(colData(dds))) > 0 & !pairs)
            all_combs <- lapply(combs, function(comb)
                resultsNames(dds) %>%
                .[grepl(comb, resultsNames(dds))]) %>%
                unlist(., recursive = FALSE)
        if (sum(combs %in% names(colData(dds))) > 0 & pairs)
            all_combs <- .createComb(dds, combs)
        if (sum(combs %in% resultsNames(dds)) > 1)
            all_combs <- c(all_combs, combs)
    }

    if (!is.null(contrast))
        all_combs <- c(all_combs, contrast)
    all_combs <- unique(all_combs)
    contrast_string <- lapply(all_combs,
                              function(x) {
                                  if (length(x) == 3)
                                      return(paste0(x[1L], "_",
                                                    paste(x[2L:3L],
                                                          collapse = "_vs_")))
                                  paste(x, collapse = ":")
                              })
    names(all_combs) <- contrast_string
    all_combs
}

#' Automatize the use of `results()` for multiple comparisons
#'
#' This function will extract the output of [DEseq2::results()]
#' and [DESeq2::lcfSrink()] for multiple comparison using:
#'
#' * coefficients
#' * contrast
#' * Multiple columns in `colData` that match coefficients
#' * Multiple columns in `colData` to create all possible
#' contrasts
#' @param dds [DESeq2::DESeqDataSet] obcject.
#' @param combs Optional vector indicating the coefficients or columns
#'   fom `colData(dds)` to create group comparisons.
#' @param contrast Optional vector to specify contrast. See [DESeq2::results()].
#' @param alpha Numeric value used in independent filtering in [DESeq2::results()].
#' @param pairs Boolean to indicate whether create all comparisons or only
#'   use the coefficient already created from [DESeq2::resultsNames()].
#'   
#' @author Lorena Pantano
#'
#' @return [DEGSet] with unSrunken and Srunken results.
#' @examples
#' library(DESeq2)
#' dds <- makeExampleDESeqDataSet(betaSD=1)
#' colData(dds)[["treatment"]] <- sample(colData(dds)[["condition"]], 12)
#'   design(dds) <-  ~ condition + treatment
#' dds <- DESeq(dds)
#' res <- degComps(dds, combs = c("condition", 2),
#'               contrast = list("treatment_B_vs_A", c("condition", "A", "B")))
#' @export
degComps <- function(dds, combs = NULL, contrast = NULL,
                          alpha = 0.05, pairs = FALSE) {
    stopifnot(class(dds) == "DESeqDataSet")

    all_combs <- .guessComb(dds,
                            combs = combs, contrast = contrast,
                            pairs = pairs)
    message("Doing ", length(all_combs), " element(s).")

    default <- "shrunk"
    attributes(default) <- list(default = "shrunken")
    message("Doing results() for each element.")
    resUnshrunken <- lapply(all_combs, function(x) .guessResults(dds, x, alpha))
    message("Doing lcfSrink() for each element.")
    resShrunken <- lapply(names(all_combs), function(x)
        .guessShrunken(dds,
                       all_combs[[x]],
                       resUnshrunken[[x]]))
    
    names(resShrunken) <- names(all_combs)
    rdsList <- lapply(names(all_combs), function(x)
        new("DEGSet", list(raw = resUnshrunken[[x]],
                                 shrunk = resShrunken[[x]]),
            default = default)
        )
    names(rdsList) <- names(all_combs)
    rdsList
}

# plotMA for DEGSet object
.plotMA <- function(results, 
                    title = NULL,
                    label_points = NULL,
                    label_column = "symbol",
                    limit = NULL,
                    diff = 5) {
    
    res_unshrunken <- degTable(results, "raw")
    res_shrunken <- degTable(results)
    
    res_all <- rownames_to_column(as.data.frame(res_unshrunken),
                                  var = "ID") %>%
        inner_join(rownames_to_column(as.data.frame(res_shrunken),
                                      var = "ID"), 
                   by = "ID", suffix = c("_unshrunken","_shrunken"))
    
    res_all[["sign"]] <- res_all[["padj_shrunken"]] < 0.05 * 1
    res_all <- res_all %>% filter(!is.na(.[["padj_shrunken"]]))
    
    toplot <- (abs(res_all[["log2FoldChange_shrunken"]] - res_all[["log2FoldChange_unshrunken"]])) >= diff
    
    if (!is.null(limit)){
            res_all[["log2FoldChange_shrunken"]][res_all[["log2FoldChange_shrunken"]] < -1 * limit] <-  -1 * limit
            res_all[["log2FoldChange_shrunken"]][res_all[["log2FoldChange_shrunken"]] > 1 * limit] <-  1 * limit
            res_all[["log2FoldChange_unshrunken"]][res_all[["log2FoldChange_unshrunken"]] < -1 * limit] <- -1 * limit
            res_all[["log2FoldChange_unshrunken"]][res_all[["log2FoldChange_unshrunken"]] > 1 * limit] <- 1 * limit
    }
    
    res_all_subset <- res_all[toplot,]
    
    p <- ggplot(res_all, aes_string("baseMean_shrunken",
                                    "log2FoldChange_shrunken")) +
        geom_point(size = 0.8, color = "black") +
        geom_point(data = res_all[res_all[["sign"]],], 
                   aes_string("baseMean_shrunken",
                              "log2FoldChange_shrunken"),
                   color = "red", size = 0.9) + 
        scale_x_log10(
            breaks = trans_breaks("log10", function(x) 10L ^ x),
            labels = trans_format("log10", math_format(10L ^ .x))) + # nolint
        annotation_logticks(sides = "b") +
        xlab("mean expression across all samples") +
        ylab(expression(log[2]*" fold change")) + # nolint
        scale_color_manual(values = c("black","red", "green")) +
        guides(color = FALSE) + 
        geom_point(data = res_all_subset,
                   aes_string("baseMean_unshrunken",
                              "log2FoldChange_unshrunken"),
                   color = "blue") +
        geom_segment(data = res_all_subset,
                     aes_string(x = "baseMean_unshrunken",
                                xend = "baseMean_shrunken",
                                y = "log2FoldChange_unshrunken", 
                                yend = "log2FoldChange_shrunken"),
                     color = "blue",
                     arrow = arrow(length = unit(0.05, "inches"),
                                   type = "closed"))
    
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    if (!is.null(label_points)) {
        labels <- res_all[res_all[[label_column]] %in% label_points, ]
        p <- p + geom_text(data = labels,
                           aes_string("baseMean_shrunken",
                                      "log2FoldChange_shrunken",
                                      label = "label_column"), size = 3L)
    }
    p
}

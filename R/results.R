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
    object <- object[rownames(unShrunken),]
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
#' This function will extract the output of [DESeq2::results()]
#' and [DESeq2::lfcShrink()] for multiple comparison using:
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

    default <- "shrunken"
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
                           shrunken = resShrunken[[x]]),
            default = default)
        )
    names(rdsList) <- names(all_combs)
    rdsList
}

.plot_shrunken <- function(res_all, res_all_subset){
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
    p
}

.plot_raw <- function(res_all){
    ggplot(res_all, aes_string("baseMean_unshrunken",
                               "log2FoldChange_unshrunken")) +
        geom_point(size = 0.8, color = "black") +
        geom_point(data = res_all[res_all[["sign"]],], 
                   aes_string("baseMean_unshrunken",
                              "log2FoldChange_unshrunken"),
                   color = "red", size = 0.9) + 
        scale_x_log10(
            breaks = trans_breaks("log10", function(x) 10L ^ x),
            labels = trans_format("log10", math_format(10L ^ .x))) + # nolint
        annotation_logticks(sides = "b") +
        xlab("mean expression across all samples") +
        ylab(expression(log[2]*" fold change")) + # nolint
        scale_color_manual(values = c("black","red", "green")) +
        guides(color = FALSE)
}

.plot_correlation <- function(res_all){
    ggplot(res_all, aes_string("log2FoldChange_unshrunken",
                               "log2FoldChange_shrunken")) +
        geom_point(size = 0.8, color = "black") +
        geom_point(data = res_all[res_all[["sign"]],], 
                   aes_string("log2FoldChange_unshrunken",
                              "log2FoldChange_shrunken"),
                   color = "red", size = 0.9) + 
        xlab("fold change - unshrunken") +
        ylab(expression(log[2]*" fold change - shrunken")) + # nolint
        scale_color_manual(values = c("black","red", "green")) +
        guides(color = FALSE)
}

.merge_results <- function(results, raw){
    res_unshrunken <- deg(results, "raw")
    res_shrunken <- deg(results, "shrunken")
    base_mean <- "baseMean_shrunken"
    if (raw)
        base_mean <- "baseMean_unshrunken"
    log2fc <- "log2FoldChange_shrunken"
    if (raw)
        log2fc <- "log2FoldChange_unshrunken"
    res_all <- rownames_to_column(as.data.frame(res_unshrunken),
                                  var = "ID") %>%
        inner_join(rownames_to_column(as.data.frame(res_shrunken),
                                      var = "ID"), 
                   by = "ID", suffix = c("_unshrunken","_shrunken"))
    
    res_all[["sign"]] <- res_all[["padj_shrunken"]] < 0.05 * 1
    res_all <- res_all %>% filter(!is.na(.[["padj_shrunken"]]))
    res_all
}

# plotMA for DEGSet object
.plotMA <- function(results, 
                    title = NULL,
                    label_points = NULL,
                    label_column = "symbol",
                    limit = NULL,
                    diff = 5,
                    raw = FALSE,
                    correlation = FALSE) {
    
    if (raw & correlation)
        stop("Use one or another. Incompatible parameters.")
    
    res_all <- .merge_results(results, raw)    
    toplot <- (abs(res_all[["log2FoldChange_shrunken"]] - res_all[["log2FoldChange_unshrunken"]])) >= diff
    
    if (!is.null(limit)){
            res_all[["log2FoldChange_shrunken"]][res_all[["log2FoldChange_shrunken"]] < -1 * limit] <-  -1 * limit
            res_all[["log2FoldChange_shrunken"]][res_all[["log2FoldChange_shrunken"]] > 1 * limit] <-  1 * limit
            res_all[["log2FoldChange_unshrunken"]][res_all[["log2FoldChange_unshrunken"]] < -1 * limit] <- -1 * limit
            res_all[["log2FoldChange_unshrunken"]][res_all[["log2FoldChange_unshrunken"]] > 1 * limit] <- 1 * limit
    }
    
    res_all_subset <- res_all[toplot,]
    
    if (correlation)
        p <- .plot_correlation(res_all)
    else if (raw)
        p <- .plot_raw(res_all)
    else p <- .plot_shrunken(res_all, res_all_subset)
    
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    
    if (!is.null(label_points)) {
        labels <- res_all[res_all[[label_column]] %in% label_points, ]
        p <- p + geom_text(data = labels,
                           aes_string(base_mean,
                                      log2fc,
                                      label = "label_column"), size = 3L)
    }
    p
}


# helper ====
.summary <- function(object, contrast, alpha){
    if (class(object) == "DESeqDataSet"){
        if (is.null(contrast)) {
            contrast <- resultsNames(object)[[2L]]
        }
        return(summary(.guessResults(object, contrast, alpha)))
    }
    if (class(object) == "DESeqResults")
        return(summary(object))
    if (class(object) == "DEGSet")
        return(summary(deg(object)))
    stop("No class supported.")
}
#' Print Summary Statistics of Alpha Level Cutoffs
#'
#' @rdname degSummary
#' @name degSummary
#' 
#' @param object Can be [DEGSet] or [DESeqDataSet] or [DESeqResults].
#' @param alpha Numeric vector of desired alpha cutoffs.
#' @param contrast Character vector to use with [results()] function.
#' @param caption Character vector to add as caption to the table.
#' @param kable Whether return a [knitr::kable()] output. Default is data.frame.
#' @return [data.frame] or [knitr::kable()].
#' 
#' @author Lorena Pantano
#' @references 
#' * original idea of multiple alpha values and code syntax
#'   from  Michael Steinbaugh.
#' @examples
#' library(DESeq2)
#' data(humanGender)
#' idx <- c(1:5, 75:80)
#' counts <- assays(humanGender)[[1]]
#' dse <- DESeqDataSetFromMatrix(counts[1:1000, idx],
#'                               colData(humanGender)[idx,],
#'                               design = ~group)
#' dse <- DESeq(dse)
#' res1 <- results(dse)
#' res2 <- degComps(dse, contrast = c("group_Male_vs_Female"))
#' degSummary(dse, contrast = "group_Male_vs_Female")
#' degSummary(res1)
#' degSummary(res1, kable = TRUE)
#' degSummary(res2[[1]])
#' @export
degSummary <- function(
    object,
    alpha = c(0.1, 0.05, 0.01),
    contrast = NULL,
    caption = "",
    kable = FALSE) {

    if (!is.null(contrast)) {
        caption <- contrast
    }
    df <- lapply(seq_along(alpha), function(a) {
        info <- capture.output(
            .summary(object, contrast, alpha[a])
        ) %>%
            # Get the lines of interest from summary
            .[4L:8L]
        parse <- info[1L:5L] %>%
            # Extract the values after the colon in summary
            sapply(function(a) {
                gsub("^.+\\:\\s(.+)\\s$", "\\1", a)
            }) %>%
            # Coerce to character here to remove names
            as.character
        data.frame(alpha = parse)
    }) %>%
        bind_cols %>%
        set_colnames(alpha) %>%
        set_rownames(c("LFC > 0 (up)",
                       "LFC < 0 (down)",
                       "outliers",
                       "low counts",
                       "cutoff")) 
    
    if (kable) 
        df %>%
        kable(caption = paste(caption)) %>%
        show
    else return(df)
}


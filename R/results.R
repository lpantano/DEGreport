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
#' design(dds) <-  ~ condition + treatment
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
        new("DEGSet", SimpleList(raw = resUnshrunken[[x]],
                                 shrunk = resShrunken[[x]]),
            default = default)
        )
    names(rdsList) <- names(all_combs)
    rdsList
}

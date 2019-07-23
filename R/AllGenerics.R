#' @rdname deg
#' @export
setGeneric("deg", function(object, value=NULL, tidy = NULL, top = NULL, ...)
    standardGeneric("deg"))

#' @rdname degDefault
#' @export
setGeneric("degDefault", function(object)
    standardGeneric("degDefault"))

#' @rdname degDefault
#' @export
setGeneric("degCorrect", function(object, fdr)
    standardGeneric("degCorrect"))


#' @rdname significants
#' @export
setGeneric("significants", function(object, padj = 0.05, fc = 0,
                                    direction = NULL, full = FALSE, ...)
    standardGeneric("significants"))

#' @rdname DEGSet
#' @export
setGeneric("as.DEGSet", function(object, ...)
    standardGeneric("as.DEGSet"))

#' @rdname deg
#' @export
setGeneric("deg", function(object, ...)
    standardGeneric("deg"))

#' @rdname degDefault
#' @export
setGeneric("degDefault", function(object)
    standardGeneric("degDefault"))

#' @rdname significants
#' @export
setGeneric("significants", function(object, ...)
    standardGeneric("significants"))

#' @rdname DEGSet
#' @export
setGeneric("DEGSetFromEdgeR", function(object, ...)
    standardGeneric("DEGSetFromEdgeR"))

#' @rdname DEGSet
#' @export
setGeneric("DEGSetFromDESeq2", function(object, ...)
    standardGeneric("DEGSetFromDESeq2"))

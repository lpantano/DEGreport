#' Method to get the default table to use.
#' 
#' @param object [DEGSet]
#' @author Lorena Pantano
#' @rdname degDefault
#' @export
setMethod("degDefault", signature("DEGSet"), function(object) {
    slot(object, name = "default")
})

 
#' Method to get all table stored for an specific comparison
#'
#' @inheritParams degDefault
#' @param value Character to specify which table to use.
#' @author Lorena Pantano
#' @rdname degTable
#' @export
setMethod("degTable", signature("DEGSet"), function(object, value=NULL) {
    if (is.null(value))
        return(object[[degDefault(object)]])
    if (value == "raw")
        return(object[["raw"]])
    if (value == "shrunk")
        return(object[["shrunk"]])
    stop(value, "not found.")
})

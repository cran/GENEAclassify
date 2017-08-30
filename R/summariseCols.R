

#' @title Perform summary of vectors in the parent frame
#'
#' @description Parses a matrix of object/summary instructions,
#' and applies the specified function to the object.
#' The object(s) named in colfun must exist in a parent frame.
#'
#' @param colfun a character matrix with two columns.
#'     The first should be a function name, the second a column in data.
#' @param nfr single positive integer number of frames up to look for named variables (default 3)
#' @return data frame with column names constructed from colfun
#' @keywords internal
#' @examples
#' A <- split(1:10, rep(1, each = 10))
#' dataCols <- matrix(c("A", "mean",
#'         "A", "median",
#'         "A", "var"), ncol = 2, byrow = TRUE)
#' GENEAclassify:::summariseCols(colfun = dataCols)
#' A <- split(1:10, rep(1:2, each = 5))
#' GENEAclassify:::summariseCols(colfun = dataCols)

summariseCols <- function(colfun, nfr = 3L) {

    if (!( is(object = colfun, class2 = "matrix") && is.character(colfun) &&
        ncol(colfun) == 2 && nrow(colfun) > 0 )) {
        stop("dataCols must be a character matrix with two columns") }

    if (!all(sapply(colfun[, 1], exists, envir = parent.frame()))) {
        stop("variables named in colfun do not exist") }

    if (!all(sapply(colfun[, 2], exists))) {
        stop("functions named in colfun do not exist") }

    sumfun <- function(x, nfr) {
        ex <- paste0("sapply(", x[1], ", ", x[2], ")")
        val <- eval(parse(text = ex), envir = parent.frame(n = nfr))
        return(val)
    }

    vals <- apply(X = colfun, MARGIN = 1, FUN = sumfun, nfr = nfr)

    if (!is.matrix(vals)) { vals <- matrix(vals, ncol = nrow(colfun), byrow = FALSE) }

    colnames(vals) <- apply(colfun, MARGIN = 1, paste, collapse = ".")

    vals <- as.data.frame(vals)

    rownames(vals) <- seq_len(nrow(vals))

    return(vals)
}

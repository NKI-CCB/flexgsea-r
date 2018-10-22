#' Score genes with limma with a design matrix.
#'
#' Scores genes by their coefficients in from limma. 
#' The argument \code{y} should be a design matrix, as produced by
#' \code{model.matrix}. 
#' For RNA-seq data, use the
#' \code{EList} object produced by the \code{\link[limma]{voom}} function.
#' Do not call directly, but give as the \option{gene.score.fn}  argument to
#' \code{\link{flexgsea}}.
#'
#' @family gene scoring functions
#' @usage flexgsea_limma
#'
#' @export
flexgsea_limma <- function (x, y, abs=F) {
    stopifnot(is.model.matrix(y))

    fit <- limma::lmFit(x, y)
    fit <- limma::eBayes(fit)
    t_stat <- fit$t

    if (colnames(t_stat)[[1]] == "(Intercept)") {
        t_stat = t_stat[, -1, drop=F]
    }

    if (abs) {
        abs(t_stat)
    } else{
        t_stat
    }
}

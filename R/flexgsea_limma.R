#' Score genes with limma with a design matrix.
#'
#' Scores genes by their coefficients in from limma. 
#' The argument \code{y} should be a design matrix, as produced by
#' \code{model.matrix}. 
#' For RNA-seq data, use the
#' \code{EList} object produced by the \code{\link[limma]{voom}} function.
#' Do not call directly, but give as the \option{gene.score.fn}  argument to
#' \code{\link{flexgsea}}.
#' Old function, use flexgsea_limma_voom or flexgsea_limma_trend instead.
#'
#'
#' @family gene scoring functions
#' @usage flexgsea_limma
#'
#' @export
flexgsea_limma <- list(
    prepare = function(x, y, gene.names, abs) {
        warning('flexgsea_limma has been deprecated, use flexgsea_limma_voom instead')
        if (!(methods::is(x, 'EList'))) {
            stop('x should be an limma EList from limma::voom')
        }
        n.samples <- ncol(x)
        if (!is.null(gene.names)) {
            stopifnot(nrow(x) == length(gene.names))
            rownames(x) <- gene.names
        } else {
            if (is.null(rownames(x))) {
                stop("Gene names should be given in gene.name or as row",
                     "names of x")
            }
            gene.names <- rownames(x)
        }
        list(x = x, y = y, gene.names = gene.names, n.samples = n.samples, abs = abs)
    },
    score = function (x, y, abs=F) {
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
)


prep_limma <- function(x, y, gene.names, abs) {
    prep <- list(x=x, y=y, gene.names=gene.names, abs=abs)
    if (!is.matrix(x)) {
        stop('x should be a matrix')
    }
    if (!is.matrix(y)) {
        stop('y should be a design matrix')
    }
    prep$n.samples <- nrow(x)
    if (!is.null(gene.names)) {
        stopifnot(ncol(x) == length(gene.names))
        colnames(prep$x) <- gene.names
    } else {
        if (is.null(colnames(x))) {
            stop("Gene names should be given in gene.name or as col",
                 "names of x")
        }
        prep$gene.names <- colnames(x)
    }
    prep
}

flexgsea_limma_voom <- list(
    prepare = function(x, y, gene.names, abs) {
        prep <- prep_limma(x, y, gene.names, abs)
        dge = edgeR::calcNormFactors(edgeR::DGEList(t(prep$x)))
        prep$x = limma::voom(dge, prep$y)
        prep
    },
    score = function (x, y, abs) {
        fit <- limma::lmFit(x, y)
        fit <- limma::eBayes(fit)
        t_stat <- fit$t
        if (colnames(t_stat)[[1]] == "(Intercept)") {
            t_stat = t_stat[, -1, drop=F]
        }
        rownames(t_stat) = rownames(x)
        colnames(t_stat) = colnames(y)[-1]
        if (abs) {
            abs(t_stat)
        } else {
            t_stat
        }
    }
)

flexgsea_limma_trend <- list(
    prepare = function(x, y, gene.names, abs) {
        prep = prep_limma(x, y, gene.names, abs)
        dge = edgeR::calcNormFactors(edgeR::DGEList(t(prep$x)))
        prep$x = edgeR::cpm(dge, log=TRUE, prior.count=3)
        prep
    },
    score = function (x, y, abs) {
        fit = limma::lmFit(x, y)
        fit = limma::eBayes(fit, trend=TRUE)
        t_stat <- fit$t
        if (colnames(t_stat)[[1]] == "(Intercept)") {
            t_stat = t_stat[, -1, drop=F]
        }
        rownames(t_stat) = rownames(x)
        colnames(t_stat) = colnames(y)[-1]
        if (abs) {
            abs(t_stat)
        } else {
            t_stat
        }
    }
)

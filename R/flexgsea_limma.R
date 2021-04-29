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
flexgsea_limma <- list(
    prepare = function(x, y) {
        warning('flexgsea_limma has been deprecated, use flexgsea_limma_voom instead')
        x
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

flexgsea_limma_voom <- list(
    prepare = function(x, y, abs) {
        dge = edgeR::calcNormFactors(edgeR::DGEList(t(x)))
        str(dge)
        str(y)
        v = limma::voom(dge, y)
        list(
            x = v,
            y = y,
            abs=abs)
    },
    score = function (x, y, abs) {
        fit <- limma::lmFit(x, y)
        fit <- limma::eBayes(fit)
        t_stat <- fit$t[, -1, drop=F]
        if (abs) {
            abs(t_stat)
        } else {
            t_stat
        }
    }
)

flexgsea_limma_trend <- list(
    prepare = function(x, y, abs) {
        dge = edgeR::calcNormFactors(edgeR::DGEList(t(x)))
        list(
            x = edgeR::cpm(dge, log=TRUE, prior.count=3),
            y = y,
            abs=abs)
    },
    score = function (x, y, abs) {
        fit <- limma::lmFit(x, y)
        fit <- limma::eBayes(fit, trend=TRUE)
        t_stat <- fit$t[, -1, drop=F]
        if (abs) {
            abs(t_stat)
        } else {
            t_stat
        }
    }
)

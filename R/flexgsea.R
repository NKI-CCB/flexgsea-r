#' @importFrom stats lm sd p.adjust
#' @importFrom abind abind adrop
#' @importFrom tibble tibble

#' @useDynLib flexgsea
#' @importFrom Rcpp sourceCpp

named_full_list <- function(value, names) {
    structure(rep(list(value), length(names)), names=names)
}
named_empty_list <- function(names) {
    structure(vector('list', length(names)), names=names)
}

#' Flexible Gene Set Enrichment Analysis.
#'
#' \code{flexgsea()} does a gene set enrichment analysis, calculating significance
#' by sample permutation. Functions to score genes, calculate enrichment
#' statistic (ES), or calculate significance can be user defined and several
#' options are supplied in the \pkg{flexgsea} package.
#'
#' Gene sets are filtered. First, only genes which exist in the data set
#' \option{x} are kept. Then, gene sets smaller than \option{gs.size.min} or
#' larger than \option{gs.size.max} are filtered out.
#'
#' Runs in parallel by default if \pkg{foreach} environment is setup and
#' \option{block.size} is smaller than the number of permutations.
#'
#' @section Possible values for \option{return_values}:
#' \describe{
#'   \item{\code{es_null}:}{Null distribution of ES.}
#'   \item{\code{gene_names}:}{Gene names, as supplied to this function.}
#' }
#' Additional return values might be available when using specific gene set
#' enrichment functions.
#'
#' @section User-defined gene scoring function \code{gene.score.fn}:
#' A gene score calculation function should take the following arguments:
#' \describe{
#'   \item{\code{x}:}{The data matrix \code{x}, exactly as given to the
#'      \code{gsea} function.}
#'   \item{\code{y}:}{Response variables to test for gene set enrichment.
#'      The \code{y} given to the \code{gsea} function or a permutation of
#'      \code{y}.
#'      This is a matrix with samples in the rows, and output variables in
#'      the columns.}
#' }
#' It should return a matrix with genes in the rows and responses in the columns. The number of
#' responses is determined by this function, but must be the same for permuted y. The number of
#' responses can be simply one. A simple example is \code{\link{flexgsea_lm}}.
#'
#' @section User-defined gene set enrichment function \code{es.fn}:
#' A list of two functions (\code{prepare} and \code{run}) and two character
#' vectors (\code{extra_stats} and \code{extra}). The code{prepare} function
#' can be used to do calculations that are the same for all gene sets. It takes
#' a single argument \code{gene.score}  and can return anything, which is
#' passed to the \code{run} function. This function can be called one or
#' multiple times on any subset of permutations, so this function  should not
#' modify global state.
#' The \code{run} function should take the following arguments:
#' \describe{
#'   \item{gene.score}{Gene scores of one or more permutations in an array
#'      (genes x response variable x permutation).}
#'   \item{gene.set}{Gene set as an integer vector which indexes the first
#'      dimension of the \code{gene.score} array.}
#'   \item{prep}{Whatever the \code{prepare} function returned for this
#'      \code{gene.score}.}
#'   \item{return_stats}{A character vector of statistics to return. This
#'      function can advertise which stats are available trough
#'      \code{extra_stats} in the list. Should default to \code{c()}.}
#'   \item{return}{A character vector of other extra values to return. This
#'      function can advertise which values are available trough
#'      \code{extra} in the list. Should default to \code{c()}.}
#' }
#' It should return a list with \code{es} and any requested extra statistics
#' and other values. The extra statistics are put into the results table, while
#' the other extra values are added to the list returned by \code{flexgsea}. The
#' \code{es} element should be a matrix (response x permutation).
#' A simple example is \code{\link{flexgsea_mean}}.
#'
#' @section User-defined significance calculation \code{sig.fun}:
#' A significance calculation function should take the following arguments:
#' \describe{
#'   \item{\code{es}:}{Enrichment scores for a single output variable, a
#'     numeric vector with a length equal to the number of gene sets.}
#'   \item{\code{es_null}:}{Enrichment scores from permuted labels, a numeric
#'     array with dimensions number of gene sets by number of permutations.}
#'   \item{\code{verbose}:}{Passed from main \code{flexgsea} function.}
#'   \item{\code{abs}:}{Passed from main \code{flexgsea} function.}
#' }
#' It should return a data frame with a row for every gene set, and a column
#' for every statistic. This data frame is returned by the main \code{flexgsea}
#' function in the \code{table} list after appending gene set names.
#'
#' @seealso Gene scoring functions: \code{\link{flexgsea_s2n}},
#'   \code{\link{flexgsea_lm}}.
#' @seealso Gene set enrichment functions: \code{\link{flexgsea_mean}},
#'   \code{\link{flexgsea_weighted_ks}}, \code{\link{flexgsea_maxmean}}.
#' @seealso Functions for significance calculation:
#'   \code{\link{flexgsea_calc_sig}},\code{\link{flexgsea_calc_sig_simple}}.
#'
#' @param x Gene expression matrix (samples by genes), or EList object
#'   produced by, for example, \code{limma::\link[limma]{voom}}.
#' @param y Classes or other response variables to analyse for gene set enrichment.
#'   Vector with length of the number of features, or sample by variable
#'   matrix.
#' @param gene.sets Gene sets. Either a filename of a gmt file, or gene sets
#'   read by the \code{\link{read_gmt}} function.
#' @param gene.score.fn Function to calculate gene scores. The signal to noise
#'   ratio (\code{\link{flexgsea_s2n}}) is appropriate for comparing two classes.
#'   Correlation (\code{\link{flexgsea_lm}}) can be  used for real valued
#'   variables. Can be user-defined, as documented below.
#' @param es.fn Function to calculate enrichment scores (ES). Default is the
#'   weighted KS statistic by Subramanian et al (2005). Can be user-defined, as
#'   documented below.
#' @param sig.fun Function to calculate significance of results. Using
#'   \code{flexgsea_calc_sig_simple} is recommended for a \code{es.fn} function
#'   other than the default \code{flexgsea_weighted_ks} as the default might not
#'   be appropriate. Can be user-defined, as documented below.
#' @param gene.names Gene identifiers for the genes in the data \code{x} that
#'   match the identifiers in \code{gene.sets}. Defaults to the the
#'   row names of \code{x}.
#' @param nperm Number of permutations to run.
#' @param gs.size.min Minimum number genes in a gene set that are also in
#'   \code{x} for a gene set to be included in the analysis.
#' @param gs.size.max Maximum number genes in a gene set that are also in
#'   \code{x} for a gene set to be included in the analysis.
#' @param verbose Should progress be printed. Progress is never printed when
#'   running in parallel.
#' @param block.size Number of permutations for which gene scoring and
#'   calculation of enrichment statistic is done in one batch. One batch can
#'   use only one thread, so this setting also effects parallel processing.
#'   Lower values use less memory, but might lose performance.
#' @param parallel Should computation be done in parallel.
#' @param abs Should the absolute enrichment score be used. This appropriate
#'   when gene sets have no direction, such as the MsigDB c2.cp gene set
#'   collection.
#' @param return_values Character vector of values to be returned other than
#'   table with statistics. Possible values are documented below, and with
#'   the enrichment function used.
#' @return A list. The \code{table} element is a list with a data frame of
#'   enrichment statistics for each response variable in \option{y}. Other
#'   elements are the values requested in \option{return_values}.
#'
#' @export
flexgsea <- function(x, y, gene.sets, gene.score.fn=flexgsea_s2n,
                  es.fn=flexgsea_weighted_ks, sig.fun=flexgsea_calc_sig,
                  gene.names=NULL, nperm=1000, gs.size.min=10,
                  gs.size.max=300, verbose=TRUE, block.size=100,
                  parallel=NULL, abs=FALSE, return_values=character()) {

    #########################
    # Prepare and check input
    if (is.vector(y) || is.factor(y)) {
        y <- matrix(y, ncol=1)
    }
    stopifnot(is.matrix(y) || is.data.frame(y))

    if (is.matrix(x)) {
        t.x <- F
    } else if (methods::is(x, 'EList')) {
        t.x <- T
    } else {
        stop('x should be a matrix or limma EList')
    }

    if (t.x) {
        n.genes <- nrow(x)
        n.samples <- ncol(x)
    } else {
        n.genes <- ncol(x)
        n.samples <- nrow(x)
    }
    stopifnot(n.samples == nrow(y))

    stopifnot(is.vector(block.size))
    stopifnot(length(block.size) == 1)
    stopifnot(block.size > 0)

    stopifnot(is.vector(nperm))
    stopifnot(length(nperm) == 1)

    stopifnot(is.character(return_values))
    '' %in% return_values # Try this, so it fails fast, not after permutations.
    return_stats <- es.fn$extra_stats

    if (is.null(parallel)) {
        if ((nperm / block.size) > 1 && isNamespaceLoaded('foreach')) {
            parallel = T
        }  else {
            parallel = F
        }
    }
    if (parallel && !requireNamespace('foreach', quietly=T)) {
        stop("foreach package required for parallel computation")
    }

    if (is.character(gene.sets)) {
        if (verbose) {
            message("Reading gene sets from file")
        }
        gene.sets <- read_gmt(gene.sets)
    }
    stopifnot(is.list(gene.sets))
    for (i in 1:length(gene.sets)) {
        stopifnot(is.character(gene.sets[[i]]))
    }

    if (!is.null(gene.names)) {
        if(t.x) {
            stopifnot(nrow(x) == length(gene.names))
            rownames(x) <- gene.names
        } else {
            stopifnot(ncol(x) == length(gene.names))
            colnames(x) <- gene.names
        }
    } else if (t.x) {
        if (is.null(rownames(x))) {
            stop("Gene names should be given in gene.name or as row",
                 "names of x")
        }
        gene.names <- rownames(x)
    } else {
        if (is.null(colnames(x))) {
            stop("Gene names should be given in gene.name or as col",
                 "names of x")
        }
        gene.names <- colnames(x)
    }

    if (verbose) {
        message("Filtering gene sets on size in dataset")
    }
    gene.sets <- filter_gene_sets(gene.sets, gene.names,
        gs.size.min=gs.size.min, gs.size.max=gs.size.max, verbose=verbose)
    n.gene.sets <- length(gene.sets)
    if (n.gene.sets == 0) {
        stop("No valid gene sets after filtering for size.")
    }

    #########################
    # Calculating observed ES
    if (verbose) {
        message("Scoring Genes (Observed)")
    }
    if ('t.x' %in% methods::formalArgs(gene.score.fn)) {
        gene.scores <- gene.score.fn(x, y, t.x=t.x, abs=abs)
    } else {
        gene.scores <- gene.score.fn(x, y, abs=abs)
    }
    if (any(!is.finite(gene.scores))) {
        stop("Gene scoring function returned NA or Inf. Check that all genes in x have a positive variance.")
    }
    stopifnot(!is.null(dim(gene.scores)))
    responses = colnames(gene.scores)
    if (is.null(responses)) {
        responses <- paste0('Response ', seq(ncol(gene.scores)))
    }
    stopifnot(!is.null(responses))
    stopifnot(dim(gene.scores)[1] == n.genes)
    gene.scores <- array(gene.scores, c(dim(gene.scores), 1),
                     dimnames=c(dimnames(gene.scores), list(NULL)))
    if (verbose) {
        message("Calculating ES (Observed)")
    }
    prep <- es.fn$prepare(gene.scores)
    es <- array(NA_real_, c(n.gene.sets, length(responses)))

    extra_stats <- structure(
        rep(list(array(NA_real_, dim(es))), length(return_stats)),
        names=return_stats)
    extra = intersect(es.fn$extra, return_values)
    extra <- named_full_list(
        named_full_list(
            named_empty_list(names=names(gene.sets)),
            names=responses),
        names=es.fn$extra)
    for (gs.i in seq_along(gene.sets)) {
        gs.index <- match(gene.sets[[gs.i]], gene.names)
        s <- es.fn$run(gene.scores, gs.index, prep,
                       return_values=return_values, return_stats=return_stats)
        stopifnot(dim(s$es)[2] == 1)  # one 'permutation'
        es[gs.i, ] <- s$es[, 1]
        for (n in names(extra_stats)) {
            extra_stats[[n]][gs.i, ] <- s[[n]][, 1]
        }
        for (n in names(extra)) {
            for (response.i in seq_along(responses)) {
                extra[[n]][[response.i]][[gs.i]] <- s[[n]][[response.i]][[1]]
            }
        }
    }
    rm(prep, gene.scores)

    ##################
    # Permutation test
    if (verbose) {
        message(paste0("Performing Permutation test"))
    }
    if (parallel) {
        perm.fun <- flexgsea_perm_parallel
    } else {
        perm.fun <- flexgsea_perm_sequential
    }
    es.null <- perm.fun(x, y, gene.sets, gene.names, nperm, block.size,
                        gene.score.fn, es.fn, responses, abs=abs,
                        verbose=verbose)

    ##############
    # Significance
    if (verbose) {
        message(paste0("Calculating Significance"))
    }
    if (parallel) {
        `%dopar%` <- foreach::`%dopar%`
        sig <- foreach::foreach(response.i = seq_along(responses)) %dopar% {
            sig.fun(
                es[, response.i],
                adrop(es.null[, response.i, , drop=F], c(F, T, F)),
                verbose=FALSE, abs=abs)
        }
    } else {
        sig <- vector('list', length(responses))
        for (response.i in seq_along(responses)) {
            sig[[response.i]] <- sig.fun(
                es[, response.i],
                adrop(es.null[, response.i, , drop=F], c(F, T, F)),
                verbose=verbose, abs=abs)
        }
    }

    ################
    # Prepare output
    for (response.i in seq_along(responses)) {
        for (n in names(extra_stats)) {
            sig[[response.i]][[n]] <- extra_stats[[n]][, response.i]
        }
    }
    res_table <- lapply(sig, dplyr::mutate, GeneSet=names(gene.sets))
    names(res_table) <- responses
    res <- list(table=res_table)
    if ('es_null' %in% return_values) {
        res$es_null=es.null
    }
    if ('gene_names' %in% return_values) {
        res$gene_names=gene.names
    }
    for (n in names(extra)) {
        if (!(n %in% names(res)) & (n %in% return_values)) {
            res[[n]] <- extra[[n]]
        }
    }
    res
}

flexgsea_perm_sequential <- function(x, y, gene.sets, gene.names, nperm,
                                  block.size, gene.score.fn, es.fn, responses,
                                  abs=F, verbose=F) {
    n.gene.sets <- length(gene.sets)
    n.genes <- length(gene.names)
    n.samples <- nrow(y)

    es.null <- array(NA_real_,
        dim=c(n.gene.sets, length(responses), nperm),
        dimnames=list(GeneSet=names(gene.sets),
                      Response=responses,
                      perm=NULL)
    )
    block.start <- 1
    while (block.start <= nperm) {
        block.end <- block.start + block.size - 1
        if (block.end > nperm) {
            block.end = nperm
        }
        nperm.block  <- block.end - block.start + 1

        gene.scores.null <- array(0, c(n.genes, length(responses),
                                       nperm.block))
        if (verbose) {
            message(paste0("Scoring Genes (Null) ", block.start, '--',
                           block.end))
        }
        for (perm.i in seq_len(nperm.block)) {
            y.perm <- y[sample.int(nrow(y)), , drop=F]
            y_attrs <- attributes(y)
            for (i in seq_along(y_attrs)) {
                if (names(y_attrs)[[i]] %in% c('dim', 'dimnames', 'names',
                                               'row.names')) {
                    next
                }
                attr(y.perm, names(y_attrs)[[i]]) <- y_attrs[[i]]
            }
            r <- gene.score.fn(x, y.perm, abs=abs)
            gene.scores.null[, , perm.i] <- r
            if (verbose) {
                message(".", appendLF=F)
                utils::flush.console()
            }
        }
        if (any(!is.finite(gene.scores.null))) {
            stop("Gene scoring function returned NA or Inf.")
        }
        if (verbose) {
            message("", appendLF=T)
        }
        if (verbose) {
            message(paste0("Calculating ES (Null) ", block.start, '--',
                           block.end))
        }
        prep <- es.fn$prepare(gene.scores.null)
        for (gs.i in seq_along(gene.sets)) {
            gs.index <- match(gene.sets[[gs.i]], gene.names)
            es.null[gs.i, , seq(block.start, block.end)] <-
                es.fn$run(gene.scores.null, gs.index, prep)$es
        }

        block.start <- block.end + 1
    }
    es.null
}

flexgsea_perm_parallel <- function(x, y, gene.sets, gene.names, nperm,
                                block.size, gene.score.fn, es.fn, responses,
                                abs=F, verbose=F) {
    `%dopar%` <- foreach::`%dopar%`
    n.gene.sets <- length(gene.sets)
    n.genes <- length(gene.names)
    n.samples <- nrow(y)

    n.blocks <- ceiling(nperm / block.size)
    block.i <- 0
    es.null <- foreach::foreach(block.i=seq_len(n.blocks), .combine=abind,
                                .inorder=F, .multicombine=T) %dopar% {
        if (block.i==n.blocks & (nperm %% block.size) > 0) {
            nperm.block <- nperm %% block.size
        } else {
            nperm.block <- block.size
        }

        gene.scores.null <- array(0, c(n.genes, length(responses),
                                       nperm.block))
        for (perm.i in seq_len(nperm.block)) {
            y.perm <- y[sample.int(nrow(y)), , drop=F]
            y_attrs <- attributes(y)
            for (i in seq_along(y_attrs)) {
                if (names(y_attrs)[[i]] %in% c('dim', 'dimnames', 'names',
                                               'row.names')) {
                    next
                }
                attr(y.perm, names(y_attrs)[[i]]) <- y_attrs[[i]]
            }
            gene.scores.null[, , perm.i] <- gene.score.fn(x, y.perm, abs=abs)
        }
        if (any(!is.finite(gene.scores.null))) {
            stop("Gene scoring function returned NA or Inf.")
        }
        es.null.block <- array(NA_real_,
            dim=c(n.gene.sets, length(responses), nperm.block),
            dimnames=list(GeneSet=names(gene.sets),
                          Response=responses,
                          perm=NULL)
        )
        prep <- es.fn$prepare(gene.scores.null)
        for (gs.i in seq_along(gene.sets)) {
            gs.index <- match(gene.sets[[gs.i]], gene.names)
            es.null.block[gs.i, , ] <-
                es.fn$run(gene.scores.null, gs.index, prep)$es
        }
        es.null.block
    }
    es.null
}

adj_fdr_nes <- function (fdr, nes) {
    fdr_adj <- fdr
    min_fdr <- 1.0
    fdr_p <- fdr[nes >= 0]
    for (i in order(nes[nes >= 0], decreasing=F)) {
        if (fdr_p[i] <= min_fdr) {
            min_fdr <- fdr_p[i]
        } else {
            fdr_p[i] <- min_fdr
        }
    }
    fdr_adj[nes >= 0] <- fdr_p

    min_fdr <- 1.0
    fdr_n <- fdr[nes < 0]
    for (i in order(nes[nes < 0], decreasing=T)) {
        if (fdr_n[i] <= min_fdr) {
            min_fdr <- fdr_n[i]
        } else {
            fdr_n[i] <- min_fdr
        }
    }
    fdr_adj[nes < 0] <- fdr_n
    fdr_adj
}

calc_fdr_nes <- function (nes, nes_null, verbose=F, abs=F) {
    nes_count_pos <- sum(nes > 0)
    nes_count_neg <- sum(nes < 0)
    nes_null_count_pos <- sum(nes_null > 0)
    nes_null_count_neg <- sum(nes_null < 0)

    fdr <- numeric(length(nes))

    if (verbose) {
        message("Calulating FDR: ", appendLF=F)
        utils::flush.console()
    }
    for (gs_i in seq_along(nes)) {
        if (verbose) {
            message(".", appendLF=F)
            utils::flush.console()
        }

        if ((nes[gs_i] >= 0) | abs) {
            rel_rank <- (sum(nes > nes[gs_i])+1) / nes_count_pos
            rel_rank_null <- (sum(nes_null > nes[gs_i])+1) /
                (nes_null_count_pos)
        } else {
            rel_rank <- (sum(nes < nes[gs_i])+1) / nes_count_neg
            rel_rank_null <- (sum(nes_null < nes[gs_i])+1) /
                (nes_null_count_neg)
        }
        fdr[gs_i] <- rel_rank_null / rel_rank
    }
    if (verbose) {
        message("Adjusting FDR")
    }
    adj_fdr_nes(fdr, nes)
}


#' Calculate significance of gene set enrichments from permutations.
#'
#' Calculates significance by the rank of ES score in the permuted values.
#' Significance is calculated separately per gene set.
#'
#' @usage flexgsea_calc_sig_simple
#' @family functions for significance calculation
#' @export
flexgsea_calc_sig_simple <- function (es, es.null, verbose=F, abs=F) {
    flexgsea_calc_sig(es, es.null, split.p=F, calc.nes=F, verbose=verbose,
                   abs=abs)
}

#' Calculate significance of gene set enrichments from permutations.
#'
#' Calculates significance like Subramanian et al.. Significance is calculated
#' separately for positive and ES scores. A normalized enrichment score (NES)
#' is calculated. The NES for all gene sets are combined to get more precision
#' with less permutations.
#'
#' @references Subramanian, A. et al. (2005) Gene Set Enrichment Analysis: A
#'   Knowledge-Based Approach for Interpreting Genome-Wide Expression
#'   Profiles. \emph{PNAS} \strong{102} (43): 15545-50.
#'   doi:10.1073/pnas.0506580102.
#' @usage flexgsea_calc_sig
#' @family functions for significance calculation
#' @export
#' @export
flexgsea_calc_sig <- function (es, es_null, split.p=T, calc.nes=T, verbose=F,
                            abs=F) {
    stopifnot(is.numeric(es))
    stopifnot(is.vector(es))
    stopifnot(is.numeric(es_null))
    n.gene.sets <- length(es)
    stopifnot(is.matrix(es_null))
    stopifnot(dim(es_null)[1] == n.gene.sets)
    n.perm <- dim(es_null)[2]

    res <- tibble(es = es)
    if (n.perm == 0) {
        return (res)
    }
    if (split.p) {
        res$p <- sapply(seq_len(n.gene.sets), function (gs.i) {
            if ((es[gs.i] >= 0.0) | abs) {
                n <- es_null[gs.i, es_null[gs.i, ] >= 0.0]
                (sum(es[gs.i] <= n) + 1) / (length(n) + 1)
            } else {
                n <- es_null[gs.i, es_null[gs.i, ] <= 0.0]
                (sum(es[gs.i] >= n) + 1) / (length(n) + 1)
            }
        })
    } else {
        res$p.high <- sapply(seq_len(n.gene.sets), function (gs.i) {
            (sum(es[gs.i] <= es_null[gs.i, ]) + 1) / (n.perm + 1)
        })
        if (abs) {
            res$p <- res$p.high
        } else {
            res$p.low <- sapply(seq_len(n.gene.sets), function (gs.i) {
                (sum(es[gs.i] >= es_null[gs.i, ]) + 1) / (n.perm + 1)
            })
            res$p <- pmin(pmin(res$p.low, res$p.high) * 2, 1.0)
        }
    }
    if (calc.nes) {
        mean_es_null_pos <- apply(es_null, 1, function (e) { mean(e[e>=0]) })
        mean_es_null_pos_nar = is.na(mean_es_null_pos) & es >= 0.0
        mean_es_null_pos[mean_es_null_pos_nar] <- es[mean_es_null_pos_nar]

        if (abs) {
            res$nes <- es / mean_es_null_pos
            nes_null <- apply(es_null, 2, function (esn) {
                esn / mean_es_null_pos
            })
        } else {
            mean_es_null_neg <- apply(es_null, 1, function (e) {
                mean(e[e<0])
            })
            mean_es_null_neg_nar = is.na(mean_es_null_neg) & es < 0.0
            mean_es_null_neg[mean_es_null_neg_nar] <- es[mean_es_null_neg_nar]
            res$nes <- es /
                ifelse(es >= 0, mean_es_null_pos, -mean_es_null_neg)
            nes_null <- apply(es_null, 2, function (esn) {
                esn / ifelse(esn >= 0, mean_es_null_pos, -mean_es_null_neg)
            })
        }
        res$fdr <- calc_fdr_nes(res$nes, nes_null, verbose=verbose, abs=abs)
    } else {
        res$fdr=p.adjust(res$p, 'BH')
    }
    res$fwer=p.adjust(res$p, 'bonferroni')
    res
}

is.model.matrix <- function(x) {
    is.matrix(x) && !is.null(attr(x, 'assign', T))
}

#' Score genes in a linear regression model.
#'
#' Scores genes by their coefficients in a linear regression
#' model of \option{y} onto \option{x}. When there is one response variable
#' this is equivalent to Pearson correlation.
#' Do not call directly, but give as the \option{gene.score.fn}  argument to
#' \code{\link{flexgsea}}.
#'
#' @family gene scoring functions
#' @usage flexgsea_lm
#'
#' @export
flexgsea_lm <- function (x, y, abs=F) {
    if (!is.model.matrix(y)) {
        if (is.null(ncol(y))) {
            y = as.matrix(y, ncol=1)
        }
        if (is.null(colnames(y))) {
            if (ncol(y) > 1) {
                colnames(y) <- paste0("Response", 1:ncol(y))
            } else {
                colnames(y) <- "Response"
            }
        }
        if (!is.data.frame(y)) {
            y = as.data.frame(y, stringsAsFactors=T)
        }
        y = model.matrix(~ ., y)
    }
    n.genes <- ncol(x)
    fit <- lm.fit(y, x)
    coeff <- t(fit$coefficients)
    if (colnames(coeff)[[1]] == "(Intercept)") {
        coeff = coeff[, -1, drop=F]
    }
    rownames(coeff) <- colnames(x)
    t <- apply(coeff, 2, '/', apply(x, 2, sd))
    if (abs) {
        t <- abs(t)
    }
    t
}

#' Score genes by there signal to noise ratio.
#'
#' Scores genes by the signal to noise ration (s2n) a. Requires \option{y}
#' to have two classes.
#' Do not call directly, but give as the \option{gene.score.fn}  argument to
#' \code{\link{flexgsea}}.
#'
#' @family gene scoring functions
#' @usage flexgsea_s2n
#'
#' @export
flexgsea_s2n <- function (x, y, abs=F) {
    if (!is.matrix(y)) {
        y = matrix(y, dimnames=list(names(y), 'Response'))
    }
    n.response <- ncol(y)
    n.genes <- ncol(x)
    if (is.logical(y)) {
        classes = c(TRUE, FALSE)
        per_response_classes = F
    } else if (is.factor(y)) {
        classes = levels(y)
        per_response_classes = F
    } else {
        per_response_classes = T
    }
    y_bin = apply(y, 2, function (phenotype) {
        if (per_response_classes) {
            classes <- sort(unique(phenotype))
        }
        stopifnot(length(classes) == 2)
        phenotype == classes[1]
    })

    coef = s2n_C(x, y_bin)
    # When the variance is zero in both classes the s2n is infinite, so we normalize
    # to a value larger than all others.
    coef[is.infinite(coef)] <- max(coef[is.finite(coef)]) * 1.1
    rownames(coef) = colnames(x)
    colnames(coef) = colnames(y)
    if (abs) {
        coef <- abs(coef)
    }
    coef
}

flexgsea_maxmean_ <- function(gene.score, gene.set, prep,
                           return_values=character(), return_stats=F) {
    total.n.genes <- dim(gene.score)[1]
    n.response <- dim(gene.score)[2]
    n.perm <- dim(gene.score)[3]

    res <- list()
    gs.o <- gene.score[gene.set, , ,drop=F]
    abs.gs.o <- abs(gs.o)
    pos <- apply((gs.o + abs.gs.o) / 2, 3, colMeans)
    neg <- apply((-gs.o + abs.gs.o) / 2, 3, colMeans)
    res$es <- pmax(pos, neg)
    res$es[neg > pos] = -1 * res$es[neg > pos]
    rownames(res$es) <- colnames(gene.score)

    res
}

#' The maxmean statistic for calculating gene set enrichement scores.
#'
#' @usage flexgsea_maxmean
#' @family gene set enrichment functions
#'
#' @export
flexgsea_maxmean <- list(
    run = flexgsea_maxmean_,
    prepare = function(gene.score) { list() },
    extra_stats = character()
)

flexgsea_mean_ <- function(gene.score, gene.set, prep,
                        return_stats=F, return_values=character()) {
    total.n.genes <- dim(gene.score)[1]
    n.response <- dim(gene.score)[2]
    n.perm <- dim(gene.score)[3]

    res <- list()
    gs.o <- gene.score[gene.set, , ,drop=F]
    res$es <- apply(gs.o, 3, colMeans)
    rownames(res$es) <- colnames(gene.score)

    res
}

#' The mean statistic for calculating gene set enrichement scores.
#'
#' @usage flexgsea_mean
#' @family gene set enrichment functions
#'
#' @export
flexgsea_mean <- list(
    run = flexgsea_mean_,
    prepare = function(gene.score) { list() },
    extra_stats = character()
)

flexgsea_weighted_ks_ <- function(gene.score, gene.set, prep, p=1.0,
                               return_stats=c(), return_values=character()) {
    total.n.genes <- dim(gene.score)[1]
    n.response <- dim(gene.score)[2]
    n.perm <- dim(gene.score)[3]

    res <- list()
    res$es <- matrix(0.0, n.response, n.perm)

    ret_extra = F
    if ('max_es_at' %in% return_stats) {
        res$max_es_at <- matrix(NA, n.response, n.perm)
        ret_extra = T
    }
    if ('le_prop' %in% return_stats) {
        res$le_prop <- matrix(NA, n.response, n.perm)
        ret_extra = T
    }
    if ('leading_edge' %in% return_values) {
        res$leading_edge <- list()
        ret_extra = T
    }
    if ('running_es_pos' %in% return_values) {
        res$running_es_pos <- list()
        ret_extra = T
    }
    if ('running_es_neg' %in% return_values) {
        res$running_es_neg <- list()
        ret_extra = T
    }
    if ('running_es_at' %in% return_values) {
        res$running_es_at <- list()
        ret_extra = T
    }
    do_le = 'leading_edge' %in% return_values || 'le_prop' %in% return_stats
    if (!ret_extra) {
        for (i in seq(n.response)) {
            for (j in seq_len(n.perm)) {
                w <- abs(gene.score[gene.set, i, j])**p
                w <- w / sum(w)
                r <- prep$gene.rank[gene.set, i, j]
                o = order(r)
                wo = w[o]
                d_hit <- cumsum(wo)
                d_miss = ((r[o] - seq_along(r)) /
                          (total.n.genes - length(gene.set)))

                running_es_pos = d_hit-d_miss
                running_es_neg = running_es_pos-wo

                es_neg <- min(running_es_neg)
                es_pos <- max(running_es_pos)
                do_pos = es_pos > -es_neg
                if (do_pos) {
                    res$es[i, j] <- es_pos
                } else {
                    res$es[i, j] <- es_neg
                }
            }
        }
    } else {
      for (i in seq(n.response)) {
          if ('leading_edge' %in% return_values) {
              res$leading_edge[[i]] <- list()
          }
          if ('running_es_pos' %in% return_values) {
              res$running_es_pos[[i]] <- list()
          }
          if ('running_es_neg' %in% return_values) {
              res$running_es_neg[[i]] <- list()
          }
          if ('running_es_at' %in% return_values) {
              res$running_es_at[[i]] <- list()
          }
          for (j in seq_len(n.perm)) {
              w <- abs(gene.score[gene.set, i, j])**p
              w <- w / sum(w)
              r <- prep$gene.rank[gene.set, i, j]
              o = order(r)
              wo = w[o]
              d_hit <- cumsum(wo)
              d_miss = ((r[o] - seq_along(r)) /
                        (total.n.genes - length(gene.set)))

              running_es_pos = d_hit-d_miss
              running_es_neg = running_es_pos-wo

              es_neg <- min(running_es_neg)
              es_pos <- max(running_es_pos)
              do_pos = es_pos > -es_neg
              if (do_pos) {
                  res$es[i, j] <- es_pos
              } else {
                  res$es[i, j] <- es_neg
              }

              if ('running_es_pos' %in% return_values) {
                  res$running_es_pos[[i]][[j]] <- running_es_pos
              }
              if ('running_es_neg' %in% return_values) {
                  res$running_es_neg[[i]][[j]] <- running_es_neg
              }
              if ('running_es_at' %in% return_values) {
                  res$running_es_at[[i]][[j]] = gene.set[o]
              }
              if (do_pos) {
                  wmes <- which.max(running_es_pos)
                  if (do_le) {
                      le_idx <- gene.set[o][1:wmes]
                  }
              } else {
                  wmes <- which.min(running_es_neg)
                  if (do_le) {
                      le_idx <- gene.set[o][wmes:length(gene.set)]
                  }
              }
              if ('max_es_at' %in% return_stats) {
                  res$max_es_at[i, j] <- r[o][wmes]
              }
              if ('le_prop' %in% return_stats) {
                  res$le_prop[i, j] <- length(le_idx) / length(gene.set)
              }
              if ('leading_edge' %in% return_values) {
                  res$leading_edge[[i]][[j]] <- le_idx
              }
          }
      }
    }
    res
}

#' The weighted KS statistic for calculating gene set enrichement scores.
#'
#' The weighted KS-like statistic as described by Subramanian et al (2005).
#'
#' @references Subramanian, A. et al. (2005) Gene Set Enrichment Analysis: A
#'   Knowledge-Based Approach for Interpreting Genome-Wide Expression
#'   Profiles. \emph{PNAS} \strong{102} (43): 15545-50.
#'   doi:10.1073/pnas.0506580102.
#' @usage flexgsea_weighted_ks
#' @family gene set enrichment functions
#'
#'
#' @export
flexgsea_weighted_ks <- list(
    run = flexgsea_weighted_ks_,
    prepare = function(gene.score) {
        total.n.genes <- dim(gene.score)[1]
        gene.rank <- apply(gene.score, c(2, 3), function (s) {
            total.n.genes - rank(s, 'keep', 'first') + 1
        })
        list(gene.rank = gene.rank)
    },
    extra_stats = c('max_es_at', 'le_prop'),
    extra = c('running_es_pos', 'running_es_neg', 'running_es_at',
              'leading_edge')
)

filter_gene_sets <- function(gene.sets, gene.names, gs.size.min=10,
                             gs.size.max=300, verbose=TRUE) {
    gs.len.before <- sapply(gene.sets, length)
    gs.filtered <- lapply(gene.sets, function (gs) {
        gs[gs %in% gene.names]
    })
    gs.len <- sapply(gs.filtered, length)
    n.before <- sum(gs.len.before)
    n.after <- sum(gs.len)
    if (verbose) {
        message(paste0("Filtered ", n.before-n.after, " of ", n.before,
                       " genes out"))
    }
    gs.filtered <- gs.filtered[gs.len >= gs.size.min & gs.len <= gs.size.max]
    if (verbose) {
        message(paste0("Filtered ", length(gene.sets) - length(gs.filtered),
                       " of ", length(gene.sets), " gene sets out"))
    }
    gs.filtered
}

#' Read gene sets in gmt format.
#'
#' @param file Either a path to a file, a connection, or literal data. Passed to
#'    the \pkg{readr} package so also supports some compression and urls.
#' @param progress Display a progress bar? Passed to the \pkg{readr} package.
#' @export
read_gmt <- function(file, progress=interactive()) {
    lines <- readr::read_lines(file, progress=progress)
    lines <- stringr::str_split(lines, '\t', 3)
    pw_names <- sapply(lines, `[[`, 1)
    pw_genes <- stringr::str_split(sapply(lines, `[[`, 3), '\t')
    names(pw_genes) <- pw_names
    pw_genes
}

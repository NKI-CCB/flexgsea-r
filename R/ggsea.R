#' @importFrom stats lm sd p.adjust
#' @importFrom abind abind

named_full_list <- function(value, names) {
    structure(rep(list(value), length(names)), names=names)
}
named_empty_list <- function(names) {
    structure(vector('list', length(names)), names=names)
}

#' Flexible Gene Set Enrichment Analysis.
#'
#' \code{ggsea} does a gene set enrichment analysis, calculating significance
#' by sample permutation. Functions to score genes, calculate enrichment
#' statistic (ES), or calculate significance can be user defined and several
#' options are supplied in the \pkg{ggsea} package.
#'
#' Gene sets are filtered. First, only genes with exist in the data set
#' \option{x} are kept. Then, gene sets smaller than \option{gs.size.min} or
#' larger than \option{gs.size.max} are filtered out.
#'
#' Runs in parallel by default if \pkg{foreach} environment is setup and
#' \option{block.size} is smaller than the number of permutations.
#'
#' Possible values for \option{return_values}:
#' \describe{
#'   \item{\code{es_null}:}{Null distribution of ES.}
#'   \item{\code{gene_names}:}{Gene names, as supplied to this function.}
#' }
#'
#' @param x Gene expression matrix (genes by samples), or EList object
#'   produced by, for example, \code{limma::voom()}.
#' @param y Classes or other variables to analyse for gene set enrichment.
#'   Vector with length of the number of features, or sample by variable
#'   matrix.
#' @param gene.sets Gene sets. Either a filename of gmt file, or gene sets
#'   read by \code{read_gmt}.
#' @param gene.score.fn Function to calculate gene scores. The signal to noise
#'   ratio (\code{ggsea_s2n}) is appropriate for comparing two classes.
#'   Correlation (\code{ggsea_lm}) can be  used for real valued variables.
#' @param es.fn Function to calculate ES.
#' @param sig.fun Function to calculate significance of results. Using
#'   \code{ggsea_calc_sig_simple} is recommended for \code{es.fn} other than
#'   \code{ggsea_weighted_ks} as the default might not be appropriate.
#' @param gene.names Gene identifiers for the genes in \code{x} that match the
#'   identifiers in \code{gene.sets}. Can also be given as the row names of
#'   \code{x}.
#' @param nperm Number of permutations to run.
#' @param gs.size.min Minimum number genes in a gene set that are also in
#'   \code{x} for a gene set to be included in the analysis.
#' @param gs.size.max Maximum number genes in a gene set that are also in
#'   \code{x} for a gene set to be included in the analysis.
#' @param verbose Should progress be printed. Progress is never printed when
#'   running in parallel.
#' @param block.size Number of permutations for which gene scoring and
#'   calculation of enrichment statistic is done in one batch. One batch is
#'   can use only thread, so this setting also effects parallel processing.
#' @param parallel Should computation be done in parallel.
#' @param return_values Character vector of values to be returned other than
#'   table with statistics. Possible values are documented below, and with
#'   the enrichment function used.
#' @return A list:
#'   \item{table}{A list with a data frame of enrichment statistics for each
#'     variable in \option{y}.}
#'   \item{...}{Values requested in \code{return_values}}
#'
#' @export
ggsea <- function(x, y, gene.sets, gene.score.fn=ggsea_s2n,
                  es.fn=ggsea_weighted_ks, sig.fun=ggsea_calc_sig,
                  gene.names=NULL, nperm=1000, gs.size.min=10,
                  gs.size.max=300, verbose=TRUE, block.size=100,
                  parallel=NULL, abs=F, return_values=character()) {

    #########################
    # Prepare and check input
    if (is.vector(y)) {
        y <- matrix(y, ncol=1)
    }
    stopifnot(is.matrix(y))
    n.response <- ncol(y)

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

    if (t.x & is.null(rownames(x))) {
        if(is.null(gene.names) | nrow(x) != length(gene.names)) {
            stop("Gene names should be given in gene.name or as column",
                 "names of x")
        }
        rownames(x) <- gene.names
    } else if (!t.x & is.null(colnames(x))) {
        if(is.null(gene.names) | ncol(x) != length(gene.names)) {
            stop("Gene names should be given in gene.name or as column",
                 "names of x")
        }
        colnames(x) <- gene.names
    } else if (t.x) {
        gene.names <- rownames(x)
    } else {
        gene.names <- colnames(x)
    }
    if (is.null(colnames(y))) {
        colnames(y) <- paste0('Response ', seq(ncol(y)))
    }
    if (verbose) {
        message("Filtering gene sets on size in dataset")
    }
    gene.sets <- filter_gene_sets(gene.sets, gene.names,
        gs.size.min=gs.size.min, gs.size.max=gs.size.max, verbose=verbose)
    n.gene.sets <- length(gene.sets)
    if (n.gene.sets == 0) {
        warning("No valid gene sets after filtering for size.")
    }

    #########################
    # Calculating observed ES
    if (verbose) {
        message("Scoring Genes (Observed)")
    }
    if ('t.x' %in% methods::formalArgs(formals(gene.score.fn))) {
        gene.scores <- gene.score.fn(x, y, t.x=t.x, abs=abs)
    } else {
        gene.scores <- gene.score.fn(x, y, abs=abs)
    }
    stopifnot(!is.null(dim(gene.scores)))
    stopifnot(dim(gene.scores) == c(n.genes, n.response))
    gene.scores <- array(gene.scores, c(dim(gene.scores), 1),
                     dimnames=c(dimnames(gene.scores), list(NULL)))
    if (verbose) {
        message("Calculating ES (Observed)")
    }
    prep <- es.fn$prepare(gene.scores)
    es <- array(NA_real_, c(n.gene.sets, n.response))

    extra_stats <- structure(
        rep(list(array(NA_real_, dim(es))), length(return_stats)),
        names=return_stats)
    extra = intersect(es.fn$extra, return_values)
    extra <- named_full_list(
        named_full_list(
            named_empty_list(names=names(gene.sets)),
            names=colnames(y)),
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
            for (response.i in seq(n.response)) {
                extra[[n]][[response.i]][[gs.i]] <- s[[n]][[response.i]][[1]]
            }
        }
    }
    rm(prep, gene.scores)

    ##################
    # Permutation test
    if (parallel) {
        perm.fun <- ggsea_perm_parallel
    } else {
        perm.fun <- ggsea_perm_sequential
    }
    es.null <- perm.fun(x, y, gene.sets, gene.names, nperm, block.size,
                        gene.score.fn, es.fn, abs=abs, verbose=verbose)

    if (verbose) {
        message(paste0("Calculating Significance"))
    }
    sig <- vector('list', n.response)
    for (response.i in seq(n.response)) {
        sig[[response.i]] <- sig.fun(es[, response.i],
                                     es.null[, response.i, ],
                                     verbose=verbose, abs=abs)
        for (n in names(extra_stats)) {
            sig[[response.i]][[n]] <- extra_stats[[n]][, response.i]
        }
    }

    ################
    # Prepare output
    res_table <- lapply(sig, dplyr::mutate_, GeneSet=~names(gene.sets))
    names(res_table) <- colnames(y)
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

ggsea_perm_sequential <- function(x, y, gene.sets, gene.names, nperm,
                                  block.size, gene.score.fn, es.fn,
                                  abs=F, verbose=F) {
    n.gene.sets <- length(gene.sets)
    n.genes <- length(gene.names)
    n.samples <- nrow(y)
    n.response <- ncol(y)

    es.null <- array(NA_real_,
        dim=c(n.gene.sets, n.response, nperm),
        dimnames=list(GeneSet=names(gene.sets),
                      Response=colnames(y),
                      perm=NULL)
    )
    block.start <- 1
    while (block.start <= nperm) {
        block.end <- block.start + block.size - 1
        if (block.end > nperm) {
            block.end = nperm
        }
        nperm.block  <- block.end - block.start + 1

        gene.scores.null <- array(0, c(n.genes, n.response, nperm.block))
        if (verbose) {
            message(paste0("Scoring Genes (Null) ", block.start, '--',
                           block.end))
        }
        for (perm.i in seq_len(nperm.block)) {
            y.perm <- y[sample.int(nrow(y)),]
            gene.scores.null[, , perm.i] <- gene.score.fn(x, y.perm, abs=abs)
            if (verbose) {
                message(".", appendLF=F)
                utils::flush.console()
            }
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

ggsea_perm_parallel <- function(x, y, gene.sets, gene.names, nperm,
                                block.size, gene.score.fn, es.fn, abs=F,
                                verbose=F) {
    `%dopar%` <- foreach::`%dopar%`
    n.gene.sets <- length(gene.sets)
    n.genes <- length(gene.names)
    n.samples <- nrow(y)
    n.response <- ncol(y)

    n.blocks <- ceiling(nperm / block.size)
    block.i <- 0
    es.null <- foreach::foreach(block.i=seq_len(n.blocks), .combine=abind,
                                .inorder=F, .multicombine=T) %dopar% {
        if (block.i==n.blocks & (nperm %% block.size) > 0) {
            nperm.block <- nperm %% block.size
        } else {
            nperm.block <- block.size
        }

        gene.scores.null <- array(0, c(n.genes, n.response, nperm.block))
        for (perm.i in seq_len(nperm.block)) {
            y.perm <- y[sample.int(nrow(y)),]
            gene.scores.null[, , perm.i] <- gene.score.fn(x, y.perm, abs=abs)
        }
        es.null.block <- array(NA_real_,
            dim=c(n.gene.sets, n.response, nperm.block),
            dimnames=list(GeneSet=names(gene.sets),
                          Response=colnames(y),
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
    nes_count_pos <- sum(nes >= 0)
    nes_count_neg <- sum(nes < 0)
    nes_null_count_pos <- apply(nes_null, 2, function (nes) { sum(nes >= 0) })
    nes_null_count_neg <- apply(nes_null, 2, function (nes) { sum(nes < 0) })

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
            if (nes_count_pos > 0) {
                rel_rank <- sum(nes >= nes[gs_i]) / nes_count_pos
            } else {
                rel_rank <- 1.0
            }
            rel_rank_null <- apply(nes_null >= nes[gs_i], 2, sum) /
                nes_null_count_pos
            rel_rank_null[nes_null_count_pos == 0] = 1.0
        } else {
            if (nes_count_neg > 0) {
                rel_rank <- sum(nes < nes[gs_i]) / nes_count_neg
            } else {
                rel_rank <- 1.0
            }
            rel_rank_null <- apply(nes_null < nes[gs_i], 2, sum) /
                nes_null_count_neg
            rel_rank_null[nes_null_count_neg == 0] = 1.0
        }
        if (rel_rank > 0.0) {
            fdr[gs_i] <- mean(rel_rank_null) / rel_rank
        } else {
            fdr[gs_i] <- 0.0
        }
    }
    if (verbose) {
        message("Adjusting FDR")
    }
    adj_fdr_nes(fdr, nes)
}


#' @export
ggsea_calc_sig_simple <- function (es, es.null, verbose=F, abs=F) {
    ggsea_calc_sig(es, es.null, split.p=F, calc.nes=F, verbose=verbose,
                   abs=abs)
}

#' @export
ggsea_calc_sig <- function (es, es_null, split.p=T, calc.nes=T, verbose=F,
                            abs=F) {
    stopifnot(is.numeric(es))
    stopifnot(is.vector(es))
    stopifnot(is.numeric(es_null))
    n.gene.sets <- length(es)
    stopifnot(dim(es_null)[1] == n.gene.sets)
    n.perm <- dim(es_null)[2]

    res <- dplyr::data_frame_(list(es = ~es))
    if (n.perm == 0) {
        return (res)
    }
    if (split.p) {
        res$p <- sapply(seq_len(n.gene.sets), function (gs.i) {
            if ((es[gs.i] >= 0.0) | abs) {
                n <- es_null[gs.i, es_null[gs.i, ] >= 0.0]
                sum(es[gs.i] <= n) / length(n)
            } else {
                n <- es_null[gs.i, es_null[gs.i, ] <= 0.0]
                sum(es[gs.i] >= n) / length(n)
            }
        })
    } else {
        res$p.high <- sapply(seq_len(n.gene.sets), function (gs.i) {
            sum(es[gs.i] <= es_null[gs.i, ]) / n.perm
        })
        if (abs) {
            res$p <- res$p.high
        } else {
            res$p.low <- sapply(seq_len(n.gene.sets), function (gs.i) {
                sum(es[gs.i] >= es_null[gs.i, ]) / n.perm
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

#' @export
ggsea_lm <- function (x, y, abs=F) {
    if (!is.matrix(y)) {
        y = matrix(y, dimnames=list(names(y), 'Response'))
    }
    n.response <- ncol(y)
    n.genes <- ncol(x)
    coef <- t(lm(x ~ y)$coefficients[-1, , drop=F])
    colnames(coef) <- colnames(y)
    rownames(coef) <- colnames(x)
    coef <- apply(coef, 2, '/', apply(x, 2, sd))
    if (abs) {
        coef <- base::abs(coef)
    }
    coef
}

#' @export
ggsea_s2n <- function (x, y, abs=F) {
    if (!is.matrix(y)) {
        y = matrix(y, dimnames=list(names(y), 'Response'))
    }
    n.response <- ncol(y)
    n.genes <- ncol(x)
    coef <- apply(y, 2, function (phenotype) {
        classes <- unique(phenotype)
        stopifnot(length(classes) == 2)
        x1 <- x[phenotype==classes[1], , drop=F]
        x2 <- x[phenotype==classes[2], , drop=F]
        m1 <- apply(x1, 2, mean)
        m2 <- apply(x2, 2, mean)
        sd1 <- apply(x1, 2, sd)
        sd2 <- apply(x2, 2, sd)
        (m1 - m2) / (sd1 + sd2)
    })

    if (abs) {
        coef <- abs(coef)
    }
    coef
}

ggsea_maxmean_ <- function(gene.score, gene.set, prep,
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

#' @export
ggsea_maxmean <- list(
    run = ggsea_maxmean_,
    prepare = function(gene.score) { list() },
    extra_stats = character()
)

ggsea_mean_ <- function(gene.score, gene.set, prep,
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

#' @export
ggsea_mean <- list(
    run = ggsea_mean_,
    prepare = function(gene.score) { list() },
    extra_stats = character()
)

ggsea_weighted_ks_ <- function(gene.score, gene.set, prep, p=1.0,
                               return_stats=F, return_values=character()) {
    total.n.genes <- dim(gene.score)[1]
    n.response <- dim(gene.score)[2]
    n.perm <- dim(gene.score)[3]

    res <- list()
    res$es <- matrix(0.0, n.response, n.perm)

    ret_extra = F
    if ('max_es_at' %in% return_stats) {
        res$max_es_at <- matrix(0.0, n.response, n.perm)
        ret_extra = T
    }
    if ('le_prop' %in% return_stats) {
        res$le_prop <- matrix(0.0, n.response, n.perm)
        ret_extra = T
    }
    if ('leading_edge' %in% return_values) {
        res$leading_edge <- list()
        ret_extra = T
    }
    if ('running_es' %in% return_values) {
        res$running_es <- list()
        ret_extra = T
    }
    for (i in seq(n.response)) {
        if ('leading_edge' %in% return_values) {
            res$leading_edge[[i]] <- list()
        }
        if ('running_es' %in% return_values) {
            res$running_es[[i]] <- list()
        }
        for (j in seq_len(n.perm)) {
            g.o <- prep$gene.order[, i, j]
            gs.i <- rep(0, total.n.genes)
            gs.i[gene.set] <- 1
            gs <- abs(gene.score[, i, j])**p
            p.hit <- (gs.i * gs) / sum(gs[gene.set])
            p.miss <- (1-gs.i) / (total.n.genes - length(gene.set))
            p.r <- cumsum((p.hit - p.miss)[g.o])
            p.min <- min(p.r)
            p.max <- max(p.r)
            if (p.max > abs(p.min)) {
                res$es[i, j] = p.max
            } else {
                res$es[i, j] = p.min
            }

            if (ret_extra) {
                if ('running_es' %in% return_values) {
                    res$running_es[[i]][[j]] <- p.r
                }
                if (p.max > abs(p.min)) {
                    w.p.max <- which.max(p.r)
                    if ('max_es_at' %in% return_values) {
                        res$max_es_at[i, j] <- w.p.max
                    }
                    if ('leading_edge' %in% return_values) {
                        res$leading_edge[[i]][[j]] <-
                            which(gs.i[g.o][1:w.p.max] > 0)
                    }
                } else {
                    w.p.min <- which.min(p.r)
                    if ('max_es_at' %in% return_stats) {
                        res$max_es_at[i, j] <- w.p.min
                    }
                    if ('leading_edge' %in% return_values) {
                        res$leading_edge[[i]][[j]] <-
                            which(gs.i[g.o][w.p.min:length(gs.i)] > 0)
                    }
                }
                if ('le_prop' %in% return_stats) {
                    res$le_prop[i, j] <- length(res$leading_edge[[i]][[j]]) /
                        length(gene.set)
                }
            }
        }
    }
    res
}
#' @export
ggsea_weighted_ks <- list(
    run = ggsea_weighted_ks_,
    prepare = function(gene.score) {
        gene.order <- apply(gene.score, c(2, 3), order, decreasing=T)
        list(gene.order = gene.order)
    },
    extra_stats = c('max_es_at', 'le_prop'),
    extra = c('running_es', 'leading_edge')
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

#' @export
read_gmt <- function(file, progress=interactive()) {
    lines <- readr::read_lines(file, progress=progress)
    lines <- stringr::str_split(lines, '\t', 3)
    pw_names <- sapply(lines, `[[`, 1)
    pw_genes <- stringr::str_split(sapply(lines, `[[`, 3), '\t')
    names(pw_genes) <- pw_names
    pw_genes
}

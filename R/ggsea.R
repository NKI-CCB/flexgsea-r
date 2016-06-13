#' @importFrom stats lm sd p.adjust

#' @export
ggsea <- function(x, y, gene.sets, gene.score.fn=ggsea_lm, es.fn=ggsea_maxmean,
                  sig.fun=ggsea_calc_sig_simple, gene.names=NULL,
                  nperm=1000, gs.size.min=10, gs.size.max=300,
                  verbose=TRUE, block.size=1000) {

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
        gene.scores <- gene.score.fn(x, y, t.x=t.x)
    } else {
        gene.scores <- gene.score.fn(x, y)
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
    for (gs.i in seq_along(gene.sets)) {
        gs.index <- match(gene.sets[[gs.i]], gene.names)
        s <- es.fn$run(gene.scores, gs.index, prep)
        stopifnot(dim(s)[2] == 1)
        es[gs.i, ] <- s[, 1]
    }
    rm(prep, gene.scores)

    ##################
    # Permutation test

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
            gene.scores.null[, , perm.i] <- gene.score.fn(x, y.perm)
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
                es.fn$run(gene.scores.null, gs.index, prep)
        }

        rm(gene.scores.null)
        rm(prep)
        block.start <- block.end + 1
    }

    if (verbose) {
        message(paste0("Calculating Significance"))
    }
    sig <- rep(list(vector('list', n.gene.sets)), n.response)
    for (response.i in seq(n.response)) {
        sig[[response.i]][[gs.i]] <- sig.fun(es[, response.i],
                                             es.null[, response.i, ])
    }

    ################
    # Prepare output
    bind_and_set_names <- function (s) {
        dplyr::mutate_(dplyr::bind_rows(s), GeneSet=~names(gene.sets))
    }
    res <- lapply(sig, bind_and_set_names)
    names(res) <- colnames(y)
    list(table=res, es_null=es.null)
}

#' @export
ggsea_calc_sig_simple <- function (es, es.null) {
    ggsea_calc_sig(es, es.null, split.p=F, calc.nes=F)
}

#' @export
ggsea_calc_sig <- function (es, es.null, split.p=T, calc.nes=T) {
    stopifnot(is.numeric(es))
    stopifnot(is.vector(es))
    stopifnot(is.numeric(es.null))
    n.gene.sets <- length(es)
    stopifnot(dim(es.null)[1] == n.gene.sets)
    n.perm <- dim(es.null)[2]

    res <- dplyr::data_frame_(list(es = ~es))
    if (n.perm == 0) {
        return (res)
    }
    if (split.p) {
        res$p <- sapply(seq_len(n.gene.sets), function (gs.i) {
            if (es[gs.i] >= 0.0) {
                n <- es.null[gs.i, es.null[gs.i, ] >= 0.0]
                sum(es[gs.i] <= n) / length(n)
            } else {
                n <- es.null[gs.i, es.null[gs.i, ] <= 0.0]
                sum(es[gs.i] >= n) / length(n)
            }
        })
    } else {
        res$p.low <- sapply(seq_len(n.gene.sets), function (gs.i) {
            sum(es[gs.i] >= es.null[gs.i, ]) / n.perm
        })
        res$p.high <- sapply(seq_len(n.gene.sets), function (gs.i) {
            sum(es[gs.i] <= es.null[gs.i, ]) / n.perm
        })
        res$p <- pmin(res$p.low, res$p.high)
    }
    if (calc.nes) {
        pos.mean <- apply(es.null, 1, function (x) { mean(x[x>=0]) })
        neg.mean <- apply(es.null, 1, function (x) { -mean(x[x<0]) })
        pos.i = es >= 0
        neg.i = es < 0
        res$nes = rep(NA_real_, n.gene.sets)
        res$nes[pos.i] = es[pos.i] / pos.mean[pos.i]
        res$nes[neg.i] = es[neg.i] / neg.mean[neg.i]
    }
    res$fdr=p.adjust(res$p, 'BH')
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

ggsea_maxmean_ <- function(gene.score, gene.set, prep) {
    total.n.genes <- dim(gene.score)[1]
    n.response <- dim(gene.score)[2]
    n.perm <- dim(gene.score)[3]

    gs.o <- gene.score[gene.set, , ,drop=F]
    abs.gs.o <- abs(gs.o)
    pos <- apply((gs.o + abs.gs.o) / 2, 3, colMeans)
    neg <- apply((-gs.o + abs.gs.o) / 2, 3, colMeans)
    res <- pmax(pos, neg)
    res[neg > pos] = -1 * res[neg > pos]
    rownames(res) <- colnames(gene.score)
    res
}

#' @export
ggsea_maxmean <- list(
    run = ggsea_maxmean_,
    prepare = function(gene.score) { list() }
)

#' @export
ggsea_mean_ <- function(gene.score, gene.set, prep) {
    total.n.genes <- dim(gene.score)[1]
    n.response <- dim(gene.score)[2]
    n.perm <- dim(gene.score)[3]

    gs.o <- gene.score[gene.set, , ,drop=F]
    res <- apply(gs.o, 3, colMeans)
    rownames(res) <- colnames(gene.score)
    res
}

#' @export
ggsea_mean <- list(
    run = ggsea_mean_,
    prepare = function(gene.score) { list() }
)

ggsea_weighted_ks_ <- function(gene.score, gene.set, prep, p=1.0) {
    total.n.genes <- dim(gene.score)[1]
    n.response <- dim(gene.score)[2]
    n.perm <- dim(gene.score)[3]

    es <- matrix(0.0, n.response, n.perm)
    for (i in seq(n.response)) {
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
                es[i, j] = p.max
            } else {
                es[i, j] = p.min
            }

        }
    }
    es
}
#' @export
ggsea_weighted_ks <- list(
    run = ggsea_weighted_ks_,
    prepare = function(gene.score) {
        gene.order <- apply(gene.score, c(2, 3), order, decreasing=T)
        list(gene.order = gene.order)
    }
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

read_gmt <- function(file, progress=interactive()) {
    lines <- readr::read_lines(file, progress=progress)
    lines <- stringr::str_split(lines, '\t', 3)
    pw_names <- sapply(lines, `[[`, 1)
    pw_genes <- stringr::str_split(sapply(lines, `[[`, 3), '\t')
    names(pw_genes) <- pw_names
    pw_genes
}

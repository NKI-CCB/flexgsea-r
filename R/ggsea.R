#' @importFrom stats lm sd

#' @export
ggsea <- function(x, y, gene.sets, gene.score.fn=ggsea_lm, es.fn=ggsea_maxmean,
                  sig.fun=ggsea_calc_sig_simple, gene.names=NULL,
                  nperm=1000, gs.size.min=10, gs.size.max=300,
                  verbose=TRUE, block.size=1000) {
    if (is.vector(y)) {
        y <- matrix(y, ncol=1)
    }
    stopifnot(is.matrix(y))
    n.response <- ncol(y)
    stopifnot(is.matrix(x))
    n.genes <- ncol(x)
    n.samples <- nrow(x)
    stopifnot(n.samples == nrow(y))

    es.fn.prepare = es.fn$prepare
    es.fn.run = es.fn$run

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

    if (is.null(colnames(x))) {
        if(is.null(gene.names) | ncol(x) != length(gene.names)) {
            stop("Gene names should be given in gene.name or as column",
                 "names of x")
        }
        colnames(x) <- gene.names
    }
    if (is.null(colnames(y))) {
        colnames(y) <- paste0('Response ', seq(ncol(y)))
    }
    gene.sets.f <- filter_gene_sets(gene.sets, colnames(x), 
        gs.size.min=gs.size.min, gs.size.max=gs.size.max, verbose=verbose)
    n.gene.sets <- length(gene.sets.f)
    if (n.gene.sets == 0) {
        warning("No valid gene sets after filtering for size.")
    }

    if (verbose) {
        message("Scoring Genes (Observed)")
    }
    gene.scores <- gene.score.fn(x, y)
    gene.scores <- array(gene.scores, c(dim(gene.scores), 1),
                     dimnames=c(dimnames(gene.scores), list(NULL)))

    # Note: could do the order + stat blocked over permutations
    gene.scores.null <- array(0, c(n.genes, n.response, nperm))
    if (verbose) {
        message("Scoring Genes (Null)")
    }
    for (perm.i in seq_len(nperm)) {
        if (verbose) {
            cat('.')
        }
        y.perm <- y[sample.int(nrow(y)),]
        gene.scores.null[, , perm.i] <- gene.score.fn(x, y.perm)
    }
    if (verbose) {
        cat('\n')
    }
    prep <- es.fn.prepare(gene.scores)
    prep.null <- es.fn.prepare(gene.scores.null)
    sig <- rep(list(vector('list', n.gene.sets)), n.response)
    for (gs.i in seq_along(gene.sets.f)) {
        if (verbose) {
            message(paste0("Calculating ES (Gene Set ", gs.i, ")"))
        }
        gs.index <- match(gene.sets.f[[gs.i]], colnames(x))
        s <- es.fn.run(gene.scores, gs.index, prep)
        s.null <- es.fn.run(gene.scores.null, gs.index, prep.null)
        if (verbose) {
            message(paste0("Calculating Significance (Gene Set ", gs.i, ")"))
        }
        for (response.i in seq(n.response)) {
            sig[[response.i]][[gs.i]] <- sig.fun(s[response.i, 1],
                                                 s.null[response.i, ])
        }
    }
    bind_and_set_names <- function (s) {
        dplyr::mutate_(dplyr::bind_rows(s), GeneSet=~names(gene.sets.f))
    }
    res <- lapply(sig, bind_and_set_names) 
    names(res) <- colnames(y)
    res
}

#' @export
ggsea_calc_sig_simple <- function (es, es.null) {
    stopifnot(is.numeric(es))
    stopifnot(length(es) == 1)
    stopifnot(is.numeric(es.null))
    if (length(es.null) == 0) {
        dplyr::data_frame_(list(es = ~es))
    } else {
        dplyr::data_frame_(list(
            es = ~es,
            p.low = ~1-(sum(es > es.null) / length(es.null)),
            p.high = ~1-(sum(es < es.null) / length(es.null)),
            p=~min(p.low, p.high),
            fdr=~p.adjust(p, 'BH'),
            fwer=~p.adjust(p, 'bonferroni')
        ))
    }
}

#' @export
ggsea_calc_sig_split <- function (es, es.null) {
    stopifnot(is.numeric(es))
    stopifnot(length(es) == 1)
    stopifnot(is.numeric(es.null))

    if (length(es.null) == 0) {
        dplyr::data_frame_(list(es = ~es))
    } else {
        if (es >= 0.0) {
            es.null <- es.null[es.null >= 0.0]
            p.perm <- sum(es <= es.null) / length(es.null)
        } else {
            es.null <- es.null[es.null <= 0.0]
            p.perm <- sum(es >= es.null) / length(es.null)
        }
        dplyr::data_frame_(list(
            es = ~es,
            p = ~p.perm,
            fdr=~p.adjust(p, 'BH'),
            fwer=~p.adjust(p, 'bonferroni')
        ))
    }
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

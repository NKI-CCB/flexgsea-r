ggsea <- function(x, y, gene.sets, gene.score.fn=ggsea_lm, es.fn=ggsea_maxmean,
                  sig.fun=ggsea_calc_sig_simple, gene.names=NULL, nperm=1000) {
    stopifnot(is.matrix(y))
    n.response <- ncol(y)
    stopifnot(is.matrix(x))
    n.genes <- ncol(x)
    n.samples <- nrow(x)
    stopifnot(n.samples == nrow(y))

    stopifnot(is.list(gene.sets))
    n.gene.sets <- length(gene.sets)
    for (i in 1:n.gene.sets) {
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

    gene.scores <- gene.score.fn(x, y)
    gene.scores <- array(gene.scores, c(dim(gene.scores), 1),
                     dimnames=c(dimnames(gene.scores), list(NULL)))
    # Note: could do the order + stat blocked over permutations
    gene.scores.null <- array(0, c(n.genes, n.response, nperm))
    for (perm.i in seq(nperm)) {
        y.perm <- y[sample.int(nrow(y)),]
        gene.scores.null[, , perm.i] <- gene.score.fn(x, y.perm)
    }
    prep <- es.fn$prepare(gene.scores)
    prep.null <- es.fn$prepare(gene.scores.null)
    sig <- rep(list(vector('list', n.gene.sets)), n.response)
    for (gs.i in seq(n.gene.sets)) {
        gs.index <- match(gene.sets[[gs.i]], colnames(x))
        s <- es.fn$run(gene.scores, gs.index, prep)
        s.null <- es.fn$run(gene.scores.null, gs.index, prep.null)
        for (response.i in seq(n.response)) {
            sig[[response.i]][[gs.i]] <- sig.fun(s[response.i, 1],
                                                 s.null[response.i, ])
        }
    }
    res <- lapply(sig, . %>%
        bind_rows %>%
        mutate_(GeneSet=~names(gene.sets)))
    
    names(res) <- colnames(y)
    res
}

ggsea_calc_sig_simple <- function (es, es.null) {
    stopifnot(is.numeric(es))
    stopifnot(length(es) == 1)
    stopifnot(is.numeric(es.null))
    stopifnot(length(es.null) > 1)
    data_frame_(list(
        es = ~es,
        p.low = ~1-(sum(es > es.null) / length(es.null)),
        p.high = ~1-(sum(es < es.null) / length(es.null)),
        p=~min(p.low, p.high),
        fdr=~p.adjust(p, 'BH'),
        fwer=~p.adjust(p, 'bonferroni')
    ))
}

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

ggsea_maxmean <- list(
    run = function(gene.score, gene.set, prep) {
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
    },
    prepare = function(gene.score) { list() }
)

ggsea_mean <- list(
    run = function(gene.score, gene.set, prep) {
        total.n.genes <- dim(gene.score)[1]
        n.response <- dim(gene.score)[2]
        n.perm <- dim(gene.score)[3]

        gs.o <- gene.score[gene.set, , ,drop=F]
        res <- apply(gs.o, 3, colMeans)
        rownames(res) <- colnames(gene.score)
        res
    },
    prepare = function(gene.score) { list() }
)

ggsea_weighted_ks <- list(
    run = function(gene.score, gene.set, prep, p=1.0) {
        total.n.genes <- dim(gene.score)[1]
        n.response <- dim(gene.score)[2]
        n.perm <- dim(gene.score)[3]


        es <- matrix(0.0, n.response, n.perm)
        for (i in seq(n.response)) {
            for (j in seq(n.perm)) {
                gs.o <- prep$gene.order[gene.set, i, j]
                s <- gene.score[gs.o, i, j]
                s <- abs(s)**p

                p.hit <- cumsum(s) / sum(s)
                r <- prep$gene.rank[gs.o, i, j]
                f <- (total.n.genes - length(s))
                p.miss = (r - seq_along(r)) / f

                p.max = max(p.hit - p.miss)
                p.min = min(c(0, p.hit) - c(p.miss, 1))
                if (p.max > abs(p.min)) {
                    es[i, j] = p.max
                } else {
                    es[i, j] = p.min
                }
            }
        }
        es
    },
    prepare = function(gene.score) {
        gene.order <- apply(gene.score, c(2, 3), order)
        gene.rank <- apply(gene.score, c(2, 3), rank, ties.method='first')
        list(gene.order = gene.order, gene.rank = gene.rank)
    }
)

ggsea <- function(x, y, gene.sets, gene.score.fn=ggsea_lm, es.fn=ggsea_maxmean,
                  gene.names=NULL, nperm=1000) {
    n.gene.sets <- length(gene.sets)
    n.response <- ncol(y)
    n.genes <- ncol(x)
    n.samples <- nrow(x)
    stopifnot(n.samples == nrow(y))

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

    es <- matrix(0.0, n.gene.sets, n.response)
    colnames(es) <- colnames(y)
    p.low <- matrix(1.0, n.gene.sets, n.response)
    colnames(p.low) <- colnames(y)
    p.high <- matrix(1.0, n.gene.sets, n.response)
    colnames(p.high) <- colnames(y)

    gene.scores <- gene.score.fn(x, y)
    gene.scores <- array(gene.scores, c(dim(gene.scores), 1),
                     dimnames=c(dimnames(gene.scores), list(NULL)))
    # Note: could do the order + stat blocked over permutations
    gene.scores.null <- array(0, c(n.genes, n.response, nperm))
    for (perm.i in seq(nperm)) {
        y.perm <- y[sample.int(nrow(y)),]
        gene.scores.null[, , perm.i] <- gene.score.fn(x, y.perm)
    }
    prep <- es.fn$prep(gene.scores)
    prep.null <- es.fn$prep(gene.scores.null)
    for (gs.i in seq(n.gene.sets)) {
        gs.index <- match(gene.sets[[gs.i]], colnames(x))
        s <- es.fn$run(gene.scores, gs.index, prep)
        s.null <- es.fn$run(gene.scores.null, gs.index, prep.null)
        es[gs.i, ] <- s[,1]
        p.low[gs.i, ] = rowMeans(s[, 1] < s.null)
        p.high[gs.i, ] = rowMeans(s[, 1] > s.null)
    }
    
    res <- lapply(colnames(y), function(n) {
        data_frame(
            GeneSet=names(gene.sets),
            es=es[, n],
            p.low=p.low[, n],
            p.high=p.high[, n],
            p=pmin(p.low, p.high),
            fdr=NA_real_,
            max.es.et=n.genes)
        })
    names(res) <- colnames(y)
    res
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

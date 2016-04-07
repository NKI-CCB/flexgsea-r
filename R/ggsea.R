ggsea <- function(x, y, gene.sets, ranking.fn=ggsea_lm, es.fn=ggsea_maxmean,
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

    ranking <- ranking.fn(x, y)
    ranking <- array(ranking, c(dim(ranking), 1),
                     dimnames=c(dimnames(ranking), list(NULL)))
    # Note: could do the order + stat blocked over permutations
    ranking.null <- array(0, c(n.genes, n.response, nperm))
    for (perm.i in seq(nperm)) {
        ranking.null[, , perm.i] <- ranking.fn(x, y[sample.int(nrow(y)),])
    }
    p <- double(n.gene.sets)
    for (gs.i in seq(n.gene.sets)) {
        gs.index <- match(gene.sets[[gs.i]], colnames(x))
        s <- es.fn(ranking, gs.index)
        s.null <- es.fn(ranking.null, gs.index)
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
        coef <- abs(coef)
    }
    coef
}

ggsea_maxmean <- function(gene.order, gene.set) {
    total.n.genes <- dim(gene.order)[1]
    n.response <- dim(gene.order)[2]
    n.perm <- dim(gene.order)[3]

    gs.o <- gene.order[gene.set, , ,drop=F]
    abs.gs.o <- abs(gs.o)
    pos <- apply((gs.o + abs.gs.o) / 2, 3, colMeans)
    neg <- apply((-gs.o + abs.gs.o) / 2, 3, colMeans)
    res <- pmax(pos, neg)
    res[neg > pos] = -1 * res[neg > pos]
    rownames(res) <- colnames(gene.order)
    res
}

ggsea_mean <- function(gene.order, gene.set) {
    n.genesets <- length(gene.sets)
    n.response <- ncol(gene.order)
    res <- matrix(0.0, n.genesets, n.response)
    for (gs.i in seq_along(gene.sets)) {
        res[gs.i, ] <- colMeans(gene.order[gene.sets[[gs.i]], ])
    }
    colnames(res) <- colnames(gene.order)
    rownames(res) <- names(gene.sets)
    res
}

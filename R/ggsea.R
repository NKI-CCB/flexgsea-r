ggsea <- function(x, y, gene.sets, ranking=ggsea_lm, es=ggsea_maxmean,
                  gene.names=NULL, nperm=1000) {
    n.genesets <- length(gene.sets)
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

    es <- as_data_frame(matrix(0.0, n.genesets, n.response))
    colnames(es) <- colnames(y)
    fdr <- as_data_frame(matrix(1.0, n.genesets, n.response))
    colnames(fdr) <- colnames(y)
    max.es.at <- as_data_frame(matrix(n.genes, n.genesets, n.response))
    colnames(max.es.at) <- colnames(y)
    
    res <- lapply(colnames(y), function(n) {
        data_frame(
            GeneSet=names(gene.sets),
            es=0.0,
            p=1.0,
            fdr=0.0,
            max.es.et=n.genes)
        })
    names(res) <- colnames(y)
    res
}

ggsea_lm <- function (x, y) {
    n.response <- ncol(y)
    n.genes <- ncol(x)
    matrix(rep(seq(n.genes), n.response), n.genes, n.response)
}

ggsea_maxmean <- function(gene.order, gene.sets) {
    n.genesets <- length(gene.sets)
    n.response <- ncol(gene.order)
    matrix(0.0, n.genesets, n.response)
}

requireNamespace('dplyr', quietly=T)
requireNamespace('stringr', quietly=T)

check_gsea_r <- function () {
    if (!dir.exists('../lib/GSEA-P-R')) {
        skip("GSEA-P-R not available")
    }
}

if (requireNamespace('foreach', quietly=T) &&
        requireNamespace('doMC', quietly=T)) {
    parallel=T
    doMC::registerDoMC()
} else {
    warning('Skipping parallel tests')
}

# Write GCT file for GSEA
write_gct <- function(d, file) {
    cat("#1.2\n", file=file)
    catf <- function(...) cat(..., file=file, append=T, sep='\t')
    n_genes <- ncol(d)
    n_samples <- nrow(d)
    sample_names = paste0("Sample", seq(n_samples))
    catf(n_genes, n_samples)
    catf('\n')
    catf("NAME", "Description", sample_names)
    catf('\n')
    for (i in 1:ncol(d)) {
        catf(colnames(d)[i], 'NA', d[[i]])
        catf('\n')
    }
}

# Write CLS file for GSEA
write_cls <- function(phenotype, file) {
    classes <- unique(phenotype)
    cat(length(phenotype), length(classes), 1, sep=' ', file=file)
    catf <- function(...) cat(..., file=file, append=T)
    catf('\n')
    catf("#", classes, sep=' ')
    catf('\n')
    catf(phenotype, sep=' ')
    catf('\n')
}

# Write GMT file for GSEA
write_gmt <- function(gs, file) {
    cat(file=file)
    catf <- function(...) cat(..., file=file, append=T, sep='\t')
    for (i in seq_along(gs)) {
        catf(names(gs)[i], "NA", gs[[i]])
        catf('\n')
    }
}

# Read results from GSEA
read_gsea_res <- function(res_dir, class) {
    fn <- paste0(res_dir, 'GSEA.analysis.SUMMARY.RESULTS.REPORT.',
                 class, '.txt')
    readr::read_tsv(fn)
}

run_gsea <- function(x, y, gs, ...) {
    if (getOption('skip_expensive_tests', F)) skip("Expensive Test")
    check_gsea_r()
    source('../lib/GSEA-P-R/GSEA.1.0.R')
    # Plotting is broken
    GSEA.HeatMapPlot <<- function(...) {}

    gct_fn <- tempfile('gct')
    write_gct(x, gct_fn)
    cls_fn <- tempfile('cls')
    write_cls(y, cls_fn)
    gmt_fn <- tempfile('gmt')
    write_gmt(gs, gmt_fn)

    out_dir <- paste0(tempfile('GSEA'), '/')
    dir.create(out_dir)
    GSEA(gct_fn, cls_fn, gs.db=gmt_fn, output.directory=out_dir,
         fdr.q.val.threshold=-1, ...)
    dplyr::bind_rows(lapply(unique(y), function (cls) {
        read_gsea_res(out_dir, cls) }))
}

context("replicate GSEA")

set.seed(487)

n_genes = 1000
n_samples_a = 10
n_samples_b = 15
n_samples = n_samples_a + n_samples_b
n_gene_sets = 100

x <- matrix(rnorm(n_genes*n_samples), n_samples)
colnames(x) <- paste0('g', seq(n_genes))
y <- sample(c(rep('a', n_samples_a), rep('b', n_samples_b)))
names(y) <- paste0('s', seq(n_samples))
rownames(x) <- names(y)
gs <- lapply(seq(n_gene_sets), function (gsn) {
    gs_size <- 10 + sample.int(50, 1)
    sample(colnames(x), gs_size)
})
names(gs) <- paste0('gs', seq(n_gene_sets))

x[y == 'a', gs[[1]]] = x[y == 'a', gs[[1]]] + 1
x[y == 'b', gs[[2]]] = x[y == 'b', gs[[2]]] + 1
x[y == 'a', gs[[11]]] = x[y == 'a', gs[[11]]] -1
x[y == 'b', gs[[12]]] = x[y == 'b', gs[[12]]] - 1
x[y == 'a', gs[[21]]] = x[y == 'a', gs[[21]]] + 10
x[y == 'b', gs[[22]]] = x[y == 'b', gs[[22]]] + 10
x[y == 'a', gs[[31]]] = x[y == 'a', gs[[31]]] - 10
x[y == 'b', gs[[32]]] = x[y == 'b', gs[[32]]] - 10

x_df <- data.frame(x)

test_that("Same results of complete GSEA", {
    nperm = 1000

    set.seed(7242)
    res_gsea <- run_gsea(x_df, y, gs,
        gs.size.threshold.min=1, gs.size.threshold.max=100,
        topgs=min(length(gs), 10), nperm=nperm)
    set.seed(7242)
    res <- flexgsea(x, y, gs,
        gene.score.fn=flexgsea_s2n, es.fn=flexgsea_weighted_ks,
        sig.fun=flexgsea_calc_sig, nperm=nperm, verbose=F,
        gs.size.min=1, gs.size.max=100, block.size=100, parallel=T)$table[[1]]
    expect_equal(res_gsea[order(res_gsea$GS), 'ES'][[1]],
                 res[order(res$GeneSet), 'es'][[1]], tolerance=.001)
    expect_equal(res_gsea[order(res_gsea$GS), 'NES'][[1]],
                 res[order(res$GeneSet), 'nes'][[1]], tolerance=0.1)
    expect_equal(res_gsea[order(res_gsea$GS), 'NOM p-val'][[1]],
                 res[order(res$GeneSet), 'p'][[1]], tolerance=0.1)
    fdr_gsea <- res_gsea[order(res_gsea$GS), 'FDR q-val'][[1]]
    fdr_gsea[is.na(fdr_gsea)] <- 0.0
    expect_equal(fdr_gsea, res[order(res$GeneSet), 'fdr'][[1]], tolerance=0.1)
})

requireNamespace('dplyr')

check_gsea_r <- function () {
    if (!dir.exists('../lib/GSEA-P-R')) {
        skip("GSEA-P-R not available")
    }
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

d = dplyr::data_frame(
    x1 = c(-0.91036224, -0.23633198,  0.07675882, -0.14270400),
    x2 = c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x2.big = 10*c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x3 = c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    x4 = -c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    x5 = c(0.1979672, -1.0453921, -0.1528823,  0.0897534),
    phenotype = c('a', 'b', 'b',  'a'))

gs = list(pw1=c('x1', 'x2', 'x3'), pw2=c('x2', 'x2.big'),
          pw3=c('x3', 'x5'))

test_that("Same ES", {
    x <- dplyr::select(d, starts_with('x'))
    y <- d$phenotype
    res_gsea <- run_gsea(x, y, gs,
        gs.size.threshold.min=1, topgs=min(length(gs), 10), nperm=5)
    res <- ggsea(data.matrix(x), y, gs,
        gene.score.fn=ggsea_s2n, es.fn=ggsea_weighted_ks, nperm=5, verbose=T,
        gs.size.min=1)[[1]]
    expect_equal(res_gsea[order(res_gsea$GS), 'ES'][[1]],
                 res[order(res$GeneSet), 'es'][[1]])
})

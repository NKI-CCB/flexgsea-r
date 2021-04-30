requireNamespace('limma')


n_genes <- 100
n_samples <- 20

suppressWarnings(RNGversion('3.4'))
set.seed(798)
counts <- matrix(rpois(n_genes * n_samples, lambda = 1000), nrow = n_samples)
colnames(counts) <- paste0('x', 1:n_genes)
y = rnorm(n_samples)

gs = list(pw1=c('x1', 'x2', 'x3'), pw2=c('x2', 'x5'),
          pw3=c('x10', 'x11', 'x12'))

design <- model.matrix(~ y)

test_res <- function(res, gs, ynames) {
    test_that("flexgsea result is a list", {
        expect_true(is.list(res))
    })
    test_that("flexgsea result table is a list", {
        expect_true(is.list(res[['table']]))
    })
    test_that("flexgsea gives results for all predictors in order", {
        expect_equal(names(res[['table']]), ynames)
    })
    test_that("flexgsea gives results for all pathways, in order", {
        for (r in res[['table']]) {
            expect_equal(r$GeneSet, names(gs))
        }
    })
    test_that("flexgsea result null is an array", {
        expect_true(is.array(res[['es_null']]))
    })
    test_that("flexgsea gives results es_null for all pathways, in order", {
        expect_equal(dimnames(res[['es_null']])[[1]], names(gs))
    })
    test_that("flexgsea gives results es_null for all predictors in order", {
        expect_equal(dimnames(res[['es_null']])[[2]], ynames)
    })
    test_that("flexgsea gives p values between 0 and 1", {
        for (i in seq_along(res[['table']])) {
            expect_true(all(res[['table']][[i]]$p >= 0.0))
            expect_true(all(res[['table']][[i]]$p <= 1.0))
        }
    })
}

funs <- alist(
    flexgsea_limma,
    flexgsea_limma_trend,
    flexgsea_limma_voom)

for (gene.score.fn_ in funs) {
    gene.score.fn <- eval(gene.score.fn_)
    context(deparse(gene.score.fn_))
    if (deparse(gene.score.fn_) == 'flexgsea_limma') {
        v <- limma::voom(t(counts), design)
        expect_warning({
            prep <- gene.score.fn$prepare(v, design, gene.names=NULL, abs=F)
        }, '.*flexgsea_limma.*')
    } else {
        prep <- gene.score.fn$prepare(counts, design, gene.names=NULL, abs=F)
    }

    test_that("Number of samples is reported correctly", {
        expect_equal(prep$n.samples, n_samples)
    })
    scores <- gene.score.fn$score(x=prep$x, y=prep$y, abs=F)

    test_that("Output is matrix", {
        expect_true(is.matrix(scores))
    })
    test_that("Number of rows in output is number of genes", {
        expect_equal(nrow(scores), ncol(counts))
    })
    test_that("Number of columns in output is number of response variables", {
        expect_equal(ncol(scores), ncol(design) - 1)
    })
    test_that("Names of columns in output is names of response variables", {
        expect_equal(colnames(scores), 'y')
    })
    test_that("Names of rows in output is names of genes", {
        expect_equal(rownames(scores), colnames(counts))
    })
    test_that("Output is finite", {
        expect_true(all(is.finite(scores)))
    })
    
    context(paste0(deparse(gene.score.fn_), ' full'))
    
    if (deparse(gene.score.fn_) == 'flexgsea_limma') {
        expect_warning({
            res <- flexgsea(v, design, gs, es.fn=flexgsea_weighted_ks, nperm=100, verbose=F,
                            gs.size.min=1,
                            block.size=11, parallel=F,
                            gene.score.fn=gene.score.fn, return_values=c('es_null'))
        }, '.*flexgsea_limma.*')
    } else {
        res <- flexgsea(counts, design, gs, es.fn=flexgsea_weighted_ks, nperm=100, verbose=T,
                        gs.size.min=1,
                        block.size=11, parallel=F,
                        gene.score.fn=gene.score.fn, return_values=c('es_null'))
    }
    

    test_res(res, gs, 'y')

}

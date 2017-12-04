requireNamespace('dplyr', quietly=T)


d = dplyr::data_frame(
    x1 = c(-0.91036224, -0.23633198,  0.07675882, -0.14270400),
    x2 = c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x2.big = 10*c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x3 = c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    x4 = -c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    x5 = c(0.1979672, -1.0453921, -0.1528823,  0.0897534),
    y1 = c(0.19801108, -1.05384953, -0.15942901,  0.09702427),
    y2 = c(1.793303,  1.589571, -1.412860,  2.192915))
x = data.matrix(dplyr::select(d, dplyr::starts_with('x')))
y = data.matrix(dplyr::select(d, dplyr::starts_with('y')))

gs = list(pw1=c('x1', 'x2', 'x3'), pw2=c('x2', 'x2.big'),
          pw3=c('x3', 'x5'))

to_mm <- function(x) {
    x = as.data.frame(x)
    m = model.matrix(~ ., x)
    m
}

to_mm_interactions <- function(x) {
    x = as.data.frame(x)
    m = model.matrix(~ y1 * y2, x)
    m
}

for (y_format_ in alist(as.matrix, as.data.frame, to_mm, to_mm_interactions)) {
for (esf_ in alist(ggsea_maxmean, ggsea_weighted_ks)) {
for (parallel in c(F,T,NULL)) {
    context(paste("ggsea", deparse(y_format_), deparse(esf_),
                  format(parallel)))
    y_format = eval(y_format_)
    esf = eval(esf_)
    if (parallel && !(requireNamespace('foreach', quietly=T) ||
                      requireNamespace('doMC', quietly=T))) {
        next
    } else {
        doMC::registerDoMC(2)
    }
    res <- ggsea(x, y_format(y), gs, es.fn=esf, nperm=100, verbose=F,
                 gs.size.min=1,
                 block.size=11, parallel=parallel, gene.score.fn=ggsea_lm,
                 return_values=c('es_null'))

    test_that("ggsea result is a list", {
        expect_true(is.list(res))
    })
    test_that("ggsea result table is a list", {
        expect_true(is.list(res[['table']]))
    })
    ynames = colnames(y)
    if (deparse(y_format_) == 'to_mm_interactions') {
        ynames = c(ynames, 'y1:y2')
    }
    test_that("ggsea gives results for all predictors in order", {
        expect_equal(names(res[['table']]), ynames)
    })
    test_that("ggsea gives results for all pathways, in order", {
        for (r in res[['table']]) {
            expect_equal(r$GeneSet, names(gs))
        }
    })
    test_that("ggsea result null is an array", {
        expect_true(is.array(res[['es_null']]))
    })
    test_that("ggsea gives results es_null for all pathways, in order", {
        expect_equal(dimnames(res[['es_null']])[[1]], names(gs))
    })
    test_that("ggsea gives results es_null for all predictors in order", {
        expect_equal(dimnames(res[['es_null']])[[2]], ynames)
    })
    test_that("ggsea gives p values between 0 and 1", {
        for (i in seq_along(res[['table']])) {
            expect_true(all(res[['table']][[i]]$p >= 0.0))
            expect_true(all(res[['table']][[i]]$p <= 1.0))
        }
    })
}}}

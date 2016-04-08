
context("ggsea")

d = data_frame(
    x1 = c(-0.91036224, -0.23633198,  0.07675882, -0.14270400),
    x2 = c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x2.big = 10*c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x3 = c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    x4 = -c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    y1 = c(0.19801108, -1.05384953, -0.15942901,  0.09702427),
    y2 = c(1.793303,  1.589571, -1.412860,  2.192915))
x = data.matrix(select(d, starts_with('x')))
y = data.matrix(select(d, starts_with('y')))

gs = list(pw1=c('x1', 'x2', 'x3'), pw2=c('x2', 'x2.big'))

res <- ggsea(x, y, gs, nperm=10)

test_that("ggsea result is a list", {
    expect_true(is.list(res))
})
test_that("ggsea gives results for all predictors in order", {
    expect_equal(names(res), colnames(y))
})
test_that("ggsea gives results for all pathways, in order", {
    for (r in res) {
        expect_equal(r$GeneSet, names(gs))
    }
})
test_that("ggsea gives p values between 0 and 1", {
    for (i in seq_along(res)) {
        expect_true(all(res[[i]]$p >= 0.0), info=names(res)[i])
        expect_true(all(res[[i]]$p <= 1.0), info=names(res)[i])
    }
})

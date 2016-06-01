requireNamespace('dplyr')

context("s2n")

d = dplyr::data_frame(
    x1 = c(-0.91036224, -0.23633198,  0.07675882, -0.14270400),
    x2 = c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x2.big = 10*c(-0.0005694878, -0.9734813268, -2.0230448058,  0.3552175208),
    x3 = c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    x4 = -c(0.1979443, -1.0453220, -0.1528467,  0.0897873),
    y1 = c('a', 'b', 'b',  'a'),
    y2 = c(1,  1, 0,  0))
x = data.matrix(dplyr::select(d, starts_with('x')))
y = as.matrix(dplyr::select(d, starts_with('y')))

scores <- ggsea_s2n(x, y)
scores.abs <- ggsea_s2n(x, y, abs=T)
scores.single <- ggsea_s2n(x, y[,1])
test_that("Output is matrix", {
    expect_true(is.matrix(scores))
    expect_true(is.matrix(scores.abs))
    expect_true(is.matrix(scores.single))
})
test_that("Number of rows in output is number of genes", {
    expect_equal(nrow(scores), ncol(x))
    expect_equal(nrow(scores.abs), ncol(x))
    expect_equal(nrow(scores.single), ncol(x))
})
test_that("Number of columns in output is number of response variables", {
    expect_equal(ncol(scores), ncol(y))
    expect_equal(ncol(scores.abs), ncol(y))
    expect_equal(ncol(scores.single), 1)
})
test_that("Names of columns in output is names of response variables", {
    expect_equal(colnames(scores), colnames(y))
    expect_equal(colnames(scores.single), "Response")
    expect_equal(colnames(scores.abs), colnames(y))
})
test_that("Names of rows in output is names of genes", {
    expect_equal(rownames(scores), colnames(x))
    expect_equal(rownames(scores.single), colnames(x))
    expect_equal(rownames(scores.abs), colnames(x))
})
test_that("Inflated variables get ranked similary to normal one", {
    expect_equal(scores['x2', ], scores['x2.big', ])
    expect_equal(scores.abs['x2', ], scores.abs['x2.big', ])
    expect_equal(scores.single['x2', ], scores.single['x2.big', ])
})

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

ranking <- ggsea_lm(x, y)
ranking.abs <- ggsea_lm(x, y, abs=T)
ranking.single <- ggsea_lm(x, y[,1])
test_that("Output is matrix", {
    expect_true(is.matrix(ranking))
    expect_true(is.matrix(ranking.abs))
    expect_true(is.matrix(ranking.single))
})
test_that("Number of rows in output is number of genes", {
    expect_equal(nrow(ranking), ncol(x))
    expect_equal(nrow(ranking.abs), ncol(x))
    expect_equal(nrow(ranking.single), ncol(x))
})
test_that("Number of columns in output is number of response variables", {
    expect_equal(ncol(ranking), ncol(y))
    expect_equal(ncol(ranking.abs), ncol(y))
    expect_equal(ncol(ranking.single), 1)
})
test_that("Names of columns in output is names of response variables", {
    expect_equal(colnames(ranking), colnames(y))
    expect_equal(colnames(ranking.single), "Response")
    expect_equal(colnames(ranking.abs), colnames(y))
})
test_that("Names of rows in output is names of genes", {
    expect_equal(rownames(ranking), colnames(x))
    expect_equal(rownames(ranking.single), colnames(x))
    expect_equal(rownames(ranking.abs), colnames(x))
})
test_that("Correlated response and variable get ranked at top", {
    expect_equal(names(which.max(ranking[, 'y1'])), 'x3')
    expect_true(names(which.max(ranking.abs[, 'y1'])) %in% c('x3', 'x4'))
    expect_equal(names(which.max(ranking.single[,1])), 'x3')
})
test_that("Anti-correlated response and variable get ranked at bottom", {
    expect_equal(names(which.min(ranking[, 'y1'])), 'x4')
    expect_equal(names(which.min(ranking.single[,1])), 'x4')
})
test_that("Inflated variables get ranked similary to normal one", {
    expect_equal(ranking['x2', ], ranking['x2.big', ])
    expect_equal(ranking.abs['x2', ], ranking.abs['x2.big', ])
    expect_equal(ranking.single['x2', ], ranking.single['x2.big', ])
})

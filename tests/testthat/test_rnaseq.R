requireNamespace('limma')


n_genes <- 1000
n_samples <- 20

suppressWarnings(RNGversion('3.4'))
set.seed(798)
counts <- matrix(rpois(n_genes * n_samples, lambda = 1000), nrow = n_samples)
y = rnorm(n_samples)

design <- model.matrix(~ y)

context("limma-trend")

prep <- flexgsea_limma_trend$prepare(counts, design, abs=F)
scores <- do.call(flexgsea_limma_trend$score, prep)

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
    expect_equal(colnames(scores), colnames(y))
})
test_that("Names of rows in output is names of genes", {
    expect_equal(rownames(scores), colnames(counts))
})
test_that("Output is finite", {
    expect_true(all(is.finite(scores)))
})

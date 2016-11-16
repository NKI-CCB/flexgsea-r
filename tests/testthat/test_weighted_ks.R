
context("weighted_ks")
n_genes <- 1000
n_response <- 3
n_perm <- 1

set.seed(7242)
gene_scores <- array(rnorm(n_genes*n_response*n_perm),
                     c(n_genes, n_response, n_perm))
gene_names <-paste0('g', seq(n_genes))
dimnames(gene_scores)[[1]] <- gene_names
dimnames(gene_scores)[[2]] <- paste0('response', seq(n_response))

gene_set <- sample(seq(n_genes), 20)


prep <- ggsea_weighted_ks$prepare(gene_scores)
res <- ggsea_weighted_ks$run(gene_scores, gene_set, prep,
                             return_stats=c('le_prop', 'max_es_at'))

test_that("Output is data frame", {
    expect_is(res, "list")
})
test_that("Output has ES", {
    expect_is(res$es, "matrix")
    expect_true(is.numeric(res$es))
})
test_that("Output has le_prop", {
    expect_is(res$le_prop, "matrix")
    expect_true(is.numeric(res$le_prop))
    expect_gt(min(res$le_prop), 0)
    expect_lte(max(res$le_prop), 1)
})
test_that("Output has max_es_at", {
    expect_is(res$max_es_at, "matrix")
    expect_true(is.numeric(res$max_es_at))
    expect_gte(min(res$max_es_at), 1)
    expect_lte(max(res$le_prop), n_genes)
})
test_that("Results did not change", {
    expect_equal_to_reference(res$es, 'test_weighted_ks-es.rds')
    expect_equal_to_reference(res$le_prop, 'test_weighted_ks-le_prop.rds')
    expect_equal_to_reference(res$max_es_at, 'test_weighted_ks-max_es_at.rds')
})


context("weighted_ks")
n_genes <- 1000
n_response <- 3
n_perm <- 1

suppressWarnings(RNGversion('3.4'))
set.seed(7242)
gene_scores <- array(rnorm(n_genes*n_response*n_perm),
                     c(n_genes, n_response, n_perm))
gene_names <-paste0('g', seq(n_genes))
dimnames(gene_scores)[[1]] <- gene_names
dimnames(gene_scores)[[2]] <- paste0('response', seq(n_response))

gene_set <- sample(seq(n_genes), 20)


prep <- flexgsea_weighted_ks$prepare(gene_scores)
res <- flexgsea_weighted_ks$run(gene_scores, gene_set, prep,
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

context("weighted_ks_extra_returns")
res <- flexgsea_weighted_ks$run(
    gene_scores, gene_set, prep,
    return_values=c('running_es_pos', 'running_es_neg', 'running_es_at', 'leading_edge'))


test_that("Output has leading edge", {
    expect_is(res$leading_edge, "list")
    expect_length(res$leading_edge, n_response)
    for (i in range(n_response)) {
        expect_length(res$leading_edge[[i]], n_perm)
        expect_lte(length(res$leading_edge[[i]][[1]]), length(gene_set))
        expect_true(is.integer(res$leading_edge[[i]][[1]]))
    }
})

test_that("Output has running ES at", {
    expect_is(res$running_es_at, "list")
    expect_length(res$running_es_at, n_response)
    for (i in range(n_response)) {
        expect_length(res$running_es_at[[i]], n_perm)
        expect_length(res$running_es_at[[i]][[1]], length(gene_set))
        expect_true(is.integer(res$running_es_at[[i]][[1]]))
    }
})

test_that("Output has positive running ES", {
    expect_is(res$running_es_neg, "list")
    expect_length(res$running_es_neg, n_response)
    for (i in range(n_response)) {
        expect_length(res$running_es_neg[[i]], n_perm)
        expect_length(res$running_es_at[[i]][[1]], length(gene_set))
        expect_true(is.double(res$running_es_neg[[i]][[1]]))
    }
})

test_that("Output has negative running ES", {
    expect_is(res$running_es_pos, "list")
    expect_length(res$running_es_pos, n_response)
    for (i in range(n_response)) {
        expect_length(res$running_es_pos[[i]], n_perm)
        expect_length(res$running_es_at[[i]][[1]], length(gene_set))
        expect_true(is.double(res$running_es_pos[[i]][[1]]))
    }
})

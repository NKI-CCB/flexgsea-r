n_gene_sets <- 10
n_perm = 1000

suppressWarnings(RNGversion('3.4'))
set.seed(7242)
es <- rnorm(n_gene_sets)
es_null <- matrix(rnorm(n_gene_sets*n_perm), n_gene_sets)
es[4] <- min(es_null[4, ]) - 1
es[5] <- max(es_null[5, ]) + 1

sig_funs = list(
    flexgsea_calc_sig=flexgsea_calc_sig,
    flexgsea_calc_sig_simple=flexgsea_calc_sig_simple
)

for (i in seq_along(sig_funs)) {
    n = names(sig_funs)[i]
    context(n)
    sig = sig_funs[[i]](es, es_null)
    test_that("Output is data frame with right columns", {
        expect_is(sig$p, "numeric")
        expect_gte(min(sig$p), 0)
        expect_lte(max(sig$p), 1)

        expect_is(sig$fdr, "numeric")
        expect_gte(min(sig$fdr), 0)
        expect_lte(max(sig$fdr), 1)

        expect_is(sig$fwer, "numeric")
        expect_gte(min(sig$fwer), 0)
        expect_lte(max(sig$fwer), 1)

        if (n == "flexgsea_calc_sig") {
            expect_is(sig$nes, "numeric")
        }
    })

    test_that("Significanct ES are significant", {
        expect_lte(sig$p[4], .1)
        expect_lte(sig$p[5], .1)
        expect_lte(sig$fdr[4], .1)
        expect_lte(sig$fdr[5], .1)
    })
    test_that("Nonsignificant ES are not significant", {
        for (i in c(1, 2, 3, 6, 7, 8, 9, 10)) {
            expect_gt(sig$fdr[i], .1)
        }
    })

    test_that("Permutation p-values are strictly larger than zero", {
        expect_gt(min(sig$p), 0.0)
        expect_gt(min(sig$fdr), 0.0)
        expect_gt(min(sig$fwer), 0.0)
    })
}

ranking_example <- function() {
  statistics <- data.frame(
    gene = c("tie_b", "below", "tie_a", "failed", "missing", "partial", "negative", "bad_p"),
    F21 = c(2, 8, 1, 100, 3, 4, -1, 6),
    p21 = c(.2, .01, .3, .001, .2, .1, 1, .02),
    q21 = c(.3, .03, .4, .008, .3, .2, 1, .04),
    F32 = c(20, 1000, 20, 2000, NA, 30, -2, 40),
    p32 = c(.2, .001, .8, .0001, .3, .04, 1, 1.1),
    q32 = c(.3, .004, .9, .001, .4, NA, 1, 1),
    F43 = c(10, 9, 8, 100, 7, NA, -3, 6),
    p43 = c(.01, .02, .03, .001, .04, NA, 1, -.1),
    q43 = c(.03, .04, .05, .01, .06, NA, 1, .5),
    r2_slide_wide_he = c(.1, .2, .3, .4, .1, .2, -.01, .3),
    r2_location_specific_he = c(.03, .029, .1, .2, .5, .03, -.02, .2),
    r2_spatial_beyond_he = c(.03, .1, .02, .5, .04, .2, -.03, .1))
  status <- data.frame(gene = statistics$gene,
    status = c("complete", "complete", "complete", "failed", "partial", "partial", "complete", "complete"))
  structure(list(statistics = statistics, status = status), class = "morphoanova_fit")
}

test_that("ranking uses the inclusive R-squared gate and F order with gene tie breaks", {
  fit <- ranking_example()
  before <- fit
  ranked <- rank_genes(fit, top_n = 100)
  expect_identical(ranked$location_specific_he$gene, c("partial", "tie_a", "tie_b"))
  expect_identical(ranked$spatial_beyond_he$gene, c("tie_b", "below", "missing"))
  expect_identical(ranked$location_specific_he$rank, 1:3)
  expect_named(ranked$location_specific_he, c("gene", "rank", "F", "p", "q", "r2"))
  expect_identical(ranked$location_specific_he$p, c(.04, .8, .2))
  expect_identical(ranked$location_specific_he$q, c(NA_real_, .9, .3))
  expect_identical(rank_genes(fit, top_n = 2)$location_specific_he$gene, c("partial", "tie_a"))
  expect_identical(fit, before)
})

test_that("the universe retains all analyzed genes before test and ranking filters", {
  fit <- ranking_example()
  ranked <- rank_genes(fit, top_n = 1)
  expect_identical(ranked$universe, fit$statistics$gene)
  expect_true(all(c("failed", "missing", "partial") %in% ranked$universe))
  none <- rank_genes(fit, top_n = 0, min_r2 = .5)
  expect_identical(none$universe, ranked$universe)
  expect_equal(nrow(none$location_specific_he), 0L)
  expect_named(none$location_specific_he, c("gene", "rank", "F", "p", "q", "r2"))
  negative <- rank_genes(fit, min_r2 = -.03)
  expect_equal(tail(negative$location_specific_he$F, 1), -2)
  expect_equal(tail(negative$location_specific_he$r2, 1), -.02)
  fit$status <- fit$status[rev(seq_len(nrow(fit$status))), ]
  expect_identical(rank_genes(fit, top_n = 1), ranked)
})

test_that("named settings apply independently and optional slide-wide ranking uses F21", {
  fit <- ranking_example()
  ranked <- rank_genes(fit,
    top_n = c(spatial_beyond_he = 2, location_specific_he = 1),
    min_r2 = c(spatial_beyond_he = .04, location_specific_he = .03))
  expect_identical(ranked$location_specific_he$gene, "partial")
  expect_identical(ranked$spatial_beyond_he$gene, c("below", "missing"))
  expect_identical(ranked$settings$top_n,
    c(location_specific_he = 1, spatial_beyond_he = 2))
  expect_identical(ranked$settings$min_r2,
    c(location_specific_he = .03, spatial_beyond_he = .04))
  one <- rank_genes(fit, components = "slide_wide_he", top_n = 2, min_r2 = .15)
  expect_identical(one$slide_wide_he$gene, c("below", "bad_p"))
  expect_identical(one$slide_wide_he$F, c(8, 6))
  expect_named(one, c("slide_wide_he", "universe", "settings"))
  expect_identical(one$settings$top_n, c(slide_wide_he = 2))
})

test_that("component setting names cannot silently shift thresholds", {
  fit <- ranking_example()
  expect_error(rank_genes(fit, min_r2 = c(.03, .05)), "scalar or a named vector")
  expect_error(rank_genes(fit, min_r2 = c(location_specific_he = .03)), "names must match")
  expect_error(rank_genes(fit, top_n = c(location_specific_he = 1, location_specific_he = 2)), "names must match")
  expect_error(rank_genes(fit, top_n = 1.5), "whole numbers")
  expect_error(rank_genes(fit, min_r2 = Inf), "finite numbers")
  expect_error(rank_genes(fit, components = "total_he"), "Choose unique components")
})

enrichment_plot_example <- function() {
  table <- data.frame(ID = c("b", "a", "c"), Description = c("Repeated term", "Repeated term", "Third term"),
    GeneRatio = c("2/10", "4/10", "1/10"), BgRatio = "5/20", Count = c(2, 4, 1),
    pvalue = c(.01, .01, .1), p.adjust = c(.02, .02, .2), qvalue = c(.02, .02, .2),
    geneID = c("1/2", "1/2/3/4", "1"))
  result <- methods::new("enrichResult", result = table, pvalueCutoff = 1, qvalueCutoff = 1,
                        pAdjustMethod = "BH")
  zero <- result
  zero@result$p.adjust[1] <- 0
  list(location_specific_he = list(GO = result, KEGG = NULL),
       spatial_beyond_he = list(GO = result, KEGG = zero),
       settings = list(components = c("location_specific_he", "spatial_beyond_he"),
                       databases = c("GO", "KEGG")))
}

test_that("enrichment plots rank within each panel and preserve raw statistics", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("clusterProfiler")
  x <- enrichment_plot_example()
  original <- x
  p <- plot_enrichment(x, top_n = 2)
  expect_s3_class(p, "ggplot")
  expect_identical(x, original)
  expect_identical(p$data$ID, c("a", "b", "a", "b", "b", "a"))
  expect_equal(p$data$gene_ratio, c(.4, .2, .4, .2, .2, .4))
  expect_equal(p$data$p.adjust, c(.02, .02, .02, .02, 0, .02))
  expect_true(all(is.finite(p$data$log_p)))
  expect_equal(p$data$log_p[5], -log10(.Machine$double.xmin))
  expect_equal(length(unique(p$data$term)), 6)
  expect_equal(levels(p$data$panel), c("GO\nLocation-specific H&E effect", "KEGG\nLocation-specific H&E effect",
    "GO\nSpatial variation beyond H&E", "KEGG\nSpatial variation beyond H&E"))
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$layout$layout), 4)
  expect_equal(nrow(built$data[[1]]), 6)
  expect_identical(built$data[[2]]$label, "No enriched terms")
  expect_gt(built$data[[2]]$x, min(p$data$gene_ratio))
  expect_lt(built$data[[2]]$x, max(p$data$gene_ratio))
  expect_identical(p$labels$colour, "-log10(BH-adjusted p-value)")
})

test_that("enrichment plots support selecting a database and empty results", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("clusterProfiler")
  x <- enrichment_plot_example()
  p <- plot_enrichment(x, components = "spatial_beyond_he", databases = "KEGG", top_n = 1)
  expect_identical(p$data$ID, "b")
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 1)
  empty <- plot_enrichment(x, components = "location_specific_he", databases = "KEGG")
  expect_equal(nrow(empty$data), 0)
  expect_equal(ggplot2::ggplot_build(empty)$data[[2]]$label, "No enriched terms")
  expect_s3_class(ggplot2::ggplotGrob(empty), "gtable")
})

test_that("enrichment plot inputs have meaningful selections and ratios", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("clusterProfiler")
  x <- enrichment_plot_example()
  expect_error(plot_enrichment(x, components = "total_he"), "components")
  expect_error(plot_enrichment(x, databases = "Reactome"), "GO and/or KEGG")
  expect_error(plot_enrichment(x, top_n = 0), "positive whole")
  x$location_specific_he$GO@result$GeneRatio[1] <- "2/0"
  expect_error(plot_enrichment(x), "GeneRatio")
})


test_that("BH ties use raw p-values before term IDs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("clusterProfiler")
  x <- enrichment_plot_example()
  x$location_specific_he$GO@result$pvalue <- c(.001, .009, .1)
  p <- plot_enrichment(x, components = "location_specific_he", databases = "GO", top_n = 1)
  expect_identical(p$data$ID, "b")
  expect_equal(p$data$pvalue, .001)
  expect_equal(p$data$p.adjust, .02)
  table <- as.data.frame(x$location_specific_he$GO)
  table$pvalue <- NULL
  x$location_specific_he$GO <- table
  fallback <- plot_enrichment(x, components = "location_specific_he", databases = "GO", top_n = 1)
  expect_identical(fallback$data$ID, "a")
  expect_true(is.na(fallback$data$pvalue))
  mixed <- plot_enrichment(x, databases = "GO", top_n = 1)
  expect_identical(mixed$data$ID, c("a", "a"))
})

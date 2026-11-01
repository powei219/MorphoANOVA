composition_summary_fit <- function() {
  ids <- paste0("spot", 1:6)
  effects <- list(a = 1:6 - 3.5, b = 100 * (6:1 - 3.5),
    c = c(2, 3, 1, 5, 6, 4) - 3.5, constant = rep(0, 6))
  spatial <- lapply(effects, function(x) data.frame(spot_id = ids,
    total_he = x, spatial_beyond_he = -x, expression = x + 1000))
  spatial$b <- spatial$b[c(3, 1, 6, 2, 5, 4), ]
  fit <- structure(list(statistics = data.frame(gene = names(spatial)), spatial = spatial,
    status = data.frame(gene = names(spatial), status = "complete", message = "")),
    class = "morphoanova_fit")
  proportions <- data.frame(spot_id = ids, tumor = seq(.1, .6, by = .1),
    immune = c(.6, .2, .3, .1, .4, .5), constant = .5)
  list(fit = fit, proportions = proportions)
}

test_that("composition summaries aggregate signed gene-level correlations", {
  skip_if_not_installed("ggplot2")
  d <- composition_summary_fit()
  genes <- names(d$fit$spatial)
  components <- c("total_he", "spatial_beyond_he", "expression")
  independent <- list()
  for (gene in genes) for (component in components) for (cell_type in c("tumor", "immune", "constant")) {
    spatial <- d$fit$spatial[[gene]]
    x <- spatial[[component]][match(d$proportions$spot_id, spatial$spot_id)]
    y <- d$proportions[[cell_type]]
    independent[[length(independent) + 1L]] <- data.frame(gene = gene, component = component,
      cell_type = cell_type, n = 6L, rho = suppressWarnings(stats::cor(x, y, method = "spearman")))
  }
  independent <- do.call(rbind, independent)
  for (method in c("mean", "median")) {
    p <- plot_composition(d$fit, d$proportions, components = components, summary = method)
    original <- attr(p, "correlations")
    expect_equal(original$rho, independent$rho, tolerance = 1e-12)
    expect_identical(original$gene, independent$gene)
    for (i in seq_len(nrow(p$data))) {
      row <- p$data[i, ]
      rho <- independent$rho[independent$component == row$component & independent$cell_type == row$cell_type]
      rho <- rho[is.finite(rho)]
      expected <- if (!length(rho)) NA_real_ else if (method == "mean") mean(rho) else stats::median(rho)
      expect_equal(row$rho, expected, tolerance = 1e-12)
      expect_equal(row$n_genes, length(rho))
    }
    expect_identical(levels(p$data$component), components)
  }
  # A dominant negatively associated gene reverses correlation of pooled means.
  p <- plot_composition(d$fit, d$proportions, genes = c("a", "b", "c"),
    components = "total_he", summary = "mean")
  pooled <- rowMeans(vapply(c("a", "b", "c"), function(gene) {
    spatial <- d$fit$spatial[[gene]]
    spatial$total_he[match(d$proportions$spot_id, spatial$spot_id)]
  }, numeric(6)))
  expect_true(p$data$rho[p$data$cell_type == "tumor"] > 0)
  expect_true(stats::cor(pooled, d$proportions$tumor, method = "spearman") < 0)
})

test_that("composition summaries preserve ID matching and finite-gene counts", {
  skip_if_not_installed("ggplot2")
  d <- composition_summary_fit()
  d$fit$spatial$c$total_he[1] <- NA_real_
  d$proportions$immune[6] <- NA_real_
  baseline <- plot_composition(d$fit, d$proportions, summary = "median")
  reordered <- plot_composition(d$fit, d$proportions[6:1, ], summary = "median")
  expect_identical(reordered$data, baseline$data)
  expect_identical(attr(reordered, "correlations"), attr(baseline, "correlations"))
  matrix <- as.matrix(d$proportions[-1])
  rownames(matrix) <- d$proportions$spot_id
  expect_identical(plot_composition(d$fit, matrix, summary = "median")$data, baseline$data)
  constant <- plot_composition(d$fit, d$proportions, genes = "constant", summary = "mean")
  expect_true(all(is.na(constant$data$rho)))
  expect_true(all(constant$data$n_genes == 0L))
  single <- plot_composition(d$fit, d$proportions, genes = "a", summary = "median")
  original <- attr(single, "correlations")
  expect_equal(single$data$rho, original$rho)
  expect_equal(single$data$n_genes, as.integer(is.finite(original$rho)))
})

test_that("grouped composition summaries summarize correlations of group means", {
  skip_if_not_installed("ggplot2")
  d <- composition_summary_fit()
  groups <- data.frame(spot_id = paste0("spot", c(6, 2, 4, 1, 5, 3)),
    group_id = c("three", "one", "two", "one", "three", "two"))
  proportions <- data.frame(group_id = c("two", "three", "one"), tumor = c(.8, .4, .1))
  p <- plot_composition(d$fit, proportions, genes = c("a", "b", "c"),
    components = "total_he", groups = groups, summary = "mean")
  expected <- vapply(c("a", "b", "c"), function(gene) {
    spatial <- d$fit$spatial[[gene]]
    x <- spatial$total_he[match(groups$spot_id, spatial$spot_id)]
    means <- tapply(x, groups$group_id, mean)
    stats::cor(means[proportions$group_id], proportions$tumor, method = "spearman")
  }, numeric(1))
  expect_equal(p$data$rho, mean(expected), tolerance = 1e-12)
  expect_equal(p$data$n_genes, 3L)
  expect_equal(attr(p, "correlations")$rho, unname(expected), tolerance = 1e-12)
  expect_equal(attr(p, "correlations")$n, rep(3L, 3))
})

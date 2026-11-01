component_summary_fit <- function() {
  values <- list(
    a = list(expression = c(1, 10, 100), total_he = c(1, 2, 3), spatial_beyond_he = c(-3, -2, -1)),
    b = list(expression = c(10, 2, 20), total_he = c(10, -5, 4), spatial_beyond_he = c(6, 4, 0)),
    c = list(expression = c(3, 30, 4), total_he = c(6, 0, 20), spatial_beyond_he = c(-10, 1, 2))
  )
  spatial <- lapply(values, function(value) {
    data.frame(spot_id = c("s1", "s2", "s3"), x = c(1, 2, 3), y = c(4, 6, 5),
      value, slide_wide_he = c(.1, .2, .3), location_specific_he = value$total_he - c(.1, .2, .3),
      total_he_q = c(.05, .2, NA_real_), spatial_beyond_he_q = c(.1, .03, .02))
  })
  spatial$b <- spatial$b[c(3, 1, 2), ]
  spatial$c <- spatial$c[c(2, 3, 1), ]
  structure(list(statistics = data.frame(gene = names(spatial), r2_slide_wide_he = .1,
    r2_location_specific_he = .2, r2_spatial_beyond_he = .3), spatial = spatial,
    status = data.frame(gene = names(spatial), status = "complete", message = ""),
    settings = list()), class = "morphoanova_fit")
}

test_that("multi-gene component maps take spot-matched medians", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  fit <- component_summary_fit()
  before <- fit
  p <- plot_components(fit, gene = c("a", "b", "c"), family = "sans",
    components = c("expression", "total_he", "spatial_beyond_he"))
  expected <- data.frame(spot_id = c("s1", "s2", "s3"), x = c(1, 2, 3), y = c(4, 6, 5),
    expression = c(3, 10, 20), total_he = c(6, 0, 4), spatial_beyond_he = c(-3, 1, 0))
  expect_identical(attr(p, "values"), expected)
  expect_identical(p[[1]]$data$value, expected$expression)
  expect_identical(p[[2]]$data$value, expected$total_he)
  expect_identical(p[[3]]$data$value, expected$spatial_beyond_he)
  expect_equal(p[[2]]$scales$get_scales("colour")$limits, c(-6, 6))
  expect_equal(p[[3]]$scales$get_scales("colour")$limits, c(-3, 3))
  expect_equal(p[[1]]$labels$subtitle, "Median across 3 genes")
  expect_identical(fit, before)

  shared <- plot_components(fit, gene = c("a", "b", "c"), shared_scale = TRUE,
    components = c("spatial_beyond_he", "expression", "total_he"), family = "sans")
  expect_identical(shared[[2]]$data$value, expected$expression)
  expect_identical(shared[[2]]$scales$get_scales("colour")$name, "Expression")
  expect_equal(shared[[1]]$scales$get_scales("colour")$limits, c(-6, 6))
  expect_equal(shared[[3]]$scales$get_scales("colour")$limits, c(-6, 6))
})

test_that("component summaries read disk genes and retain unavailable estimates", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  fit <- component_summary_fit()
  fit$spatial$c$total_he[fit$spatial$c$spot_id == "s2"] <- NA_real_
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  files <- stats::setNames(file.path(directory, paste0(names(fit$spatial), ".rds")), names(fit$spatial))
  for (gene in names(files)) {
    saveRDS(list(gene = gene, spatial = fit$spatial[[gene]],
      status = list(status = "complete")), files[[gene]])
  }
  disk <- fit
  disk$spatial <- files
  memory_plot <- plot_components(fit, gene = c("a", "b", "c"), components = "total_he", family = "sans")
  disk_plot <- plot_components(disk, gene = c("a", "b", "c"), components = "total_he", family = "sans")
  expect_identical(attr(memory_plot, "values"), attr(disk_plot, "values"))
  expect_equal(attr(disk_plot, "values")$total_he, c(6, NA_real_, 4))
})

test_that("single-gene component selections preserve data and variance labels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  fit <- component_summary_fit()
  p <- plot_components(fit, gene = "a", components = c("total_he", "expression"), family = "sans")
  expect_identical(attr(p, "values"), fit$spatial$a[c("spot_id", "x", "y", "total_he", "expression")])
  expect_equal(p[[1]]$labels$subtitle, "30.0% variance explained")
  expect_equal(p[[2]]$labels$subtitle, "Observed")
  expression_only <- plot_components(fit, gene = "a", components = "expression", shared_scale = TRUE, family = "sans")
  expect_identical(expression_only[[1]]$data$value, fit$spatial$a$expression)
  significance <- plot_components(fit, gene = "a", type = "significance",
    components = "total_he", family = "sans")
  expect_equal(attr(significance, "values")$total_he_significance, c(-log10(.05), NA_real_, NA_real_))
})

test_that("pointwise mechanism maps use both q-values and retain missing results", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  q_he <- c(.05, .050001, .05, .2, NA_real_, .01, Inf)
  q_beyond <- c(.2, .05, .05, .3, .01, NA_real_, .01)
  spatial <- data.frame(spot_id = letters[1:7], x = 1:7, y = 7:1,
    total_he_q = q_he, spatial_beyond_he_q = q_beyond)
  fit <- structure(list(statistics = data.frame(gene = "a", F21 = NA_real_, F32 = 0, F43 = 0),
    spatial = list(a = spatial), status = data.frame(gene = "a", status = "partial", message = "")),
    class = "morphoanova_fit")
  before <- fit
  p <- plot_components(fit, type = "mechanism", family = "sans")
  plotted <- attr(p, "values")
  expect_equal(as.character(plotted$mechanism), c("H&E only", "Beyond H&E only", "Both", "Neither", NA, NA, NA))
  expect_identical(plotted[names(spatial)], spatial)
  expect_identical(p[[1]]$data, plotted)
  expect_identical(fit, before)
  tighter <- plot_components(fit, type = "mechanism", alpha = .01, flip_y = TRUE, family = "sans")
  expect_equal(as.character(attr(tighter, "values")$mechanism)[1:4], rep("Neither", 4))
  expect_equal(ggplot2::ggplot_build(tighter[[1]])$data[[1]]$y, -(7:1))
})

test_that("component selections and unmatched spatial tables are rejected", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  fit <- component_summary_fit()
  for (type in c("significance", "mechanism")) {
    expect_error(plot_components(fit, gene = c("a", "b"), type = type), "single gene")
  }
  expect_error(plot_components(fit, gene = "a", type = "mechanism", components = "total_he"), "both")
  expect_error(plot_components(fit, gene = c("a", "a")), "distinct")
  expect_error(plot_components(fit, gene = "unknown"), "fitted genes")
  expect_error(plot_components(fit, gene = "a", components = "residual"), "components from")
  expect_error(plot_components(fit, gene = "a", components = character()), "components from")
  changed <- fit
  changed$spatial$b$x[1] <- 50
  expect_error(plot_components(changed, gene = c("a", "b")), "matching spot IDs and coordinates")
  changed <- fit
  changed$spatial$b$spot_id[1] <- "missing"
  expect_error(plot_components(changed, gene = c("a", "b")), "matching spot IDs and coordinates")
})

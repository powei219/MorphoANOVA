scatter_fit <- function(statistics) {
  for (test in c("21", "32", "43")) {
    statistics[[paste0("p", test)]] <- ifelse(is.finite(statistics[[paste0("F", test)]]), .1, NA_real_)
    statistics[[paste0("q", test)]] <- ifelse(is.finite(statistics[[paste0("F", test)]]), .2, NA_real_)
  }
  status <- ifelse(rowSums(is.finite(as.matrix(statistics[c("F21", "F32", "F43")]))) == 3L,
    "complete", "partial")
  structure(list(statistics = statistics,
    status = data.frame(gene = statistics$gene, status = status)), class = "morphoanova_fit")
}

test_that("scatter uses all available genes and component-specific label gates", {
  skip_if_not_installed("ggplot2")
  statistics <- data.frame(gene = c("A", "B", "C", "D"),
    F21 = c(1, 2, 3, 4), F32 = c(100, 90, 80, 70), F43 = c(1, 2, 10, 9),
    r2_location_specific_he = c(.01, .03, .2, .2), r2_spatial_beyond_he = c(.1, .1, .01, .04))
  fit <- scatter_fit(statistics)
  p <- plot_scatter(fit, top_n = 1, family = "sans")
  expect_setequal(p$layers[[2L]]$data$gene, c("B", "D"))
  expect_identical(p$layers[[2L]]$data$label_colour, c("#FF7F00", "#E41A1C"))
  expect_identical(p$layers[[3L]]$data$gene, p$layers[[2L]]$data$gene)
  points <- ggplot2::ggplot_build(p)$data
  expect_equal(points[[2L]]$x, log10(c(90, 70) + 1))
  expect_equal(points[[2L]]$y, log10(c(2, 9) + 1))
  expect_length(plot_scatter(fit, top_n = 0)$layers, 1L)
  expect_error(plot_scatter(fit, top_n = -1), "non-negative")
})

test_that("scatter distinguishes component membership and permits explicit labels", {
  skip_if_not_installed("ggplot2")
  statistics <- data.frame(gene = c("local", "spatial", "both", "manual"),
    F21 = 1:4, F32 = c(100, 20, 90, 200), F43 = c(20, 100, 90, 200),
    r2_location_specific_he = c(.1, .1, .1, .001),
    r2_spatial_beyond_he = c(.1, .1, .1, .001))
  fit <- scatter_fit(statistics)
  before <- fit
  automatic <- plot_scatter(fit, top_n = 2, family = "sans")
  expect_identical(automatic$layers[[2L]]$data$gene, c("local", "spatial", "both"))
  expect_identical(automatic$layers[[2L]]$data$label_colour, c("#FF7F00", "#E41A1C", "#E41A1C"))
  explicit <- plot_scatter(fit, top_n = 2, labels = statistics$gene, family = "sans")
  expect_identical(explicit$layers[[2L]]$data$gene, statistics$gene)
  expect_identical(explicit$layers[[2L]]$data$label_colour, c("#FF7F00", "#E41A1C", "#E41A1C", "#444444"))
  expect_identical(fit, before)
})

test_that("component maps have separate symmetric scales unless sharing is requested", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  spatial <- data.frame(spot_id = letters[1:4], x = 1:4, y = c(1, 3, 2, 4),
    expression = 0:3, slide_wide_he = c(-1, -.5, .5, 1),
    location_specific_he = c(-3, -1, 1, 3), spatial_beyond_he = c(-6, -2, 2, 6))
  fit <- structure(list(statistics = data.frame(gene = "gene", r2_slide_wide_he = .1,
    r2_location_specific_he = .2, r2_spatial_beyond_he = .3),
    spatial = list(gene = spatial),
    status = data.frame(gene = "gene", status = "complete", message = ""),
    settings = list()), class = "morphoanova_fit")
  separate <- plot_components(fit, family = "sans")
  shared <- plot_components(fit, shared_scale = TRUE, family = "sans")
  for (i in 2:4) {
    expected_bound <- c(1, 3, 6)[i - 1L]
    scale <- separate[[i]]$scales$get_scales("colour")
    expect_equal(unname(scale$limits), c(-expected_bound, expected_bound))
    expect_identical(unname(scale$map(c(-expected_bound, 0, expected_bound))), c("#3AB370", "#EAE7CC", "#FD1593"))
    expect_equal(unname(shared[[i]]$scales$get_scales("colour")$limits), c(-6, 6))
    expect_identical(separate[[i]]$data$value, spatial[[c("slide_wide_he", "location_specific_he", "spatial_beyond_he")[i - 1L]]])
  }
  expect_identical(separate$patches$layout$guides, "keep")
  expect_identical(separate[[1L]]$data$value, spatial$expression)
  expect_equal(separate[[1L]]$layers[[1L]]$aes_params$stroke, 0)
  expect_error(plot_components(fit, shared_scale = NA), "shared_scale")
})

test_that("unavailable scatter coordinates do not promote lower-ranked labels", {
  skip_if_not_installed("ggplot2")
  statistics <- data.frame(gene = c("A", "B", "C", "D"),
    F21 = 1:4, F32 = c(100, 90, 80, NA), F43 = c(NA, 8, 9, 100),
    r2_location_specific_he = rep(.1, 4), r2_spatial_beyond_he = rep(.1, 4))
  fit <- scatter_fit(statistics)
  expect_warning(p <- plot_scatter(fit, top_n = 1), "2 genes have unavailable")
  expect_setequal(p$data$gene, c("B", "C"))
  expect_length(p$layers, 1L)
})

test_that("scatter uses the same named thresholds and rankings as rank_genes", {
  skip_if_not_installed("ggplot2")
  statistics <- data.frame(gene = c("local", "spatial", "both", "below", "failed"),
    F21 = 1:5, F32 = c(100, 20, 90, 200, 500), F43 = c(20, 100, 90, 200, 500),
    r2_location_specific_he = c(.1, .1, .1, .001, .5),
    r2_spatial_beyond_he = c(.02, .1, .05, .001, .5))
  fit <- scatter_fit(statistics)
  fit$status$status[5L] <- "failed"
  top_n <- c(spatial_beyond_he = 2, location_specific_he = 1)
  min_r2 <- c(location_specific_he = .03, spatial_beyond_he = .05)
  ranked <- rank_genes(fit, top_n = top_n, min_r2 = min_r2)
  p <- plot_scatter(fit, top_n = top_n, min_r2 = min_r2, family = "sans")
  expect_setequal(p$layers[[2L]]$data$gene,
    c(ranked$location_specific_he$gene, ranked$spatial_beyond_he$gene))
  expect_identical(p$layers[[2L]]$data$label_colour, c("#FF7F00", "#E41A1C", "#E41A1C"))
  expect_equal(p$data$log_F32, log10(statistics$F32 + 1))
  expect_equal(p$data$log_F43, log10(statistics$F43 + 1))
  expect_identical(p$data$F21, statistics$F21)
  expect_identical(p$data$gene, statistics$gene)
})

test_that("significance maps preserve the BH threshold boundary and display cap", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  spatial <- data.frame(spot_id = letters[1:5], x = 1:5, y = c(1, 3, 2, 4, 5),
    total_he_q = c(.05, .050001, NA_real_, 1e-25, 0),
    spatial_beyond_he_q = c(1, .01, .2, .03, NA_real_))
  fit <- structure(list(statistics = data.frame(gene = "gene"),
    spatial = list(gene = spatial),
    status = data.frame(gene = "gene", status = "complete", message = ""),
    settings = list()), class = "morphoanova_fit")
  before <- fit
  p <- plot_components(fit, type = "significance", family = "sans")
  expect_s3_class(p, "patchwork")
  expect_equal(p$patches$plots[[1L]]$data$value[1L], -log10(.05))
  expect_true(all(is.na(p$patches$plots[[1L]]$data$value[2:3])))
  expect_equal(p[[1L]]$data$value[4:5], c(20, 20))
  expect_equal(p[[1L]]$scales$get_scales("colour")$limits, c(0, 20))
  expect_identical(p[[1L]]$scales$get_scales("colour")$labels, c("0", "10", "20+"))
  expect_identical(p[[1L]]$scales$get_scales("colour")$na.value, "grey90")
  expect_identical(p$patches$layout$guides, "collect")
  expect_identical(p[[1L]]$guides$guides$colour$params$direction, "horizontal")
  capped <- plot_components(fit, type = "significance", significance_cap = 8, family = "sans")
  expect_equal(capped[[1L]]$data$value[4:5], c(8, 8))
  expect_identical(capped[[1L]]$scales$get_scales("colour")$labels, c("0", "4", "8+"))
  expect_identical(fit, before)
  expect_error(plot_components(fit, alpha = 0), "alpha")
  expect_error(plot_components(fit, significance_cap = 0), "significance_cap")
  expect_error(plot_components(fit, significance_cap = Inf), "significance_cap")
})

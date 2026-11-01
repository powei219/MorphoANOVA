composition_example <- function() {
  ids <- paste0("spot", 1:8)
  x <- c(2, 4, 1, 8, 5, 3, 7, 6)
  spatial <- list(
    geneB = data.frame(spot_id = ids, expression = x * 2, total_he = x,
      spatial_beyond_he = -x),
    geneA = data.frame(spot_id = ids, expression = -x, total_he = -x,
      spatial_beyond_he = x))
  fit <- structure(list(statistics = data.frame(gene = names(spatial)),
    spatial = spatial,
    status = data.frame(gene = names(spatial), status = "complete", message = ""),
    settings = list()), class = "morphoanova_fit")
  proportions <- data.frame(spot_id = ids, tumor = c(.2, .6, .1, .8, NA, .4, .7, .3))
  proportions$immune <- 1 - proportions$tumor
  list(fit = fit, proportions = proportions)
}

test_that("composition heatmaps match independent signed Spearman correlations", {
  skip_if_not_installed("ggplot2")
  d <- composition_example()
  before <- d$fit
  p <- plot_composition(d$fit, d$proportions)
  expect_s3_class(p, "ggplot")
  expect_identical(d$fit, before)
  for (i in seq_len(nrow(p$data))) {
    row <- p$data[i, ]
    x <- d$fit$spatial[[row$gene]][[as.character(row$component)]]
    y <- d$proportions[[row$cell_type]]
    complete <- is.finite(x) & is.finite(y)
    expect_equal(row$rho, stats::cor(x[complete], y[complete], method = "spearman"))
    expect_equal(row$n, sum(complete))
  }
  shuffled <- plot_composition(d$fit, d$proportions[c(4, 2, 8, 7, 1, 6, 3, 5), ])
  expect_identical(shuffled$data, p$data)
  matrix_proportions <- as.matrix(d$proportions[c("tumor", "immune")])
  rownames(matrix_proportions) <- d$proportions$spot_id
  expect_identical(plot_composition(d$fit, matrix_proportions)$data, p$data)
  single <- plot_composition(d$fit, d$proportions, genes = "geneA", components = "expression")
  expect_identical(unique(single$data$gene), "geneA")
  expect_identical(as.character(unique(single$data$component)), "expression")
  built <- ggplot2::ggplot_build(p)
  expect_identical(built$layout$panel_scales_x[[1]]$get_limits(), c("immune", "tumor"))
  expect_identical(built$layout$panel_scales_y[[1]]$get_limits(), c("geneA", "geneB"))
  expect_identical(levels(p$data$component), c("total_he", "spatial_beyond_he"))
  expect_identical(p$scales$get_scales("fill")$map(c(-1, 0, 1)), c("#2166AC", "#FFFFFF", "#B2182B"))
})

test_that("composition correlations report insufficient and constant inputs", {
  skip_if_not_installed("ggplot2")
  d <- composition_example()
  d$proportions$constant <- .5
  d$proportions$incomplete <- c(.2, .4, NA, Inf, NA, NA, NA, NA)
  p <- plot_composition(d$fit, d$proportions, genes = "geneA")
  expect_true(all(is.na(p$data$rho[p$data$cell_type %in% c("constant", "incomplete")])))
  expect_equal(p$data$n[p$data$cell_type == "incomplete"], c(2, 2))
  d$fit$spatial$geneA$total_he <- 1
  p <- plot_composition(d$fit, d$proportions, genes = "geneA", components = "total_he")
  expect_true(all(is.na(p$data$rho)))
})

test_that("composition groups use matched group means and disk results", {
  skip_if_not_installed("ggplot2")
  d <- composition_example()
  groups <- data.frame(spot_id = paste0("spot", c(8, 6, 4, 2, 7, 5, 3, 1)),
    group_id = rep(c("four", "three", "two", "one"), 2))
  proportions <- data.frame(group_id = c("two", "four", "one", "three"), tumor = c(.1, .6, .3, .8))
  p <- plot_composition(d$fit, proportions, genes = "geneB", components = "total_he", groups = groups)
  x <- d$fit$spatial$geneB$total_he
  means <- c(one = mean(x[1:2]), two = mean(x[3:4]), three = mean(x[5:6]), four = mean(x[7:8]))
  expected <- stats::cor(means[proportions$group_id], proportions$tumor, method = "spearman")
  expect_equal(p$data$rho, expected)
  expect_identical(p$data$n, 4L)
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(list(gene = "geneB", spatial = d$fit$spatial$geneB,
    status = d$fit$status[1, ]), file)
  disk <- d$fit
  disk$spatial <- NULL
  attr(disk, "internal") <- list(files = c(geneB = file, geneA = "not-read.rds"))
  expect_identical(plot_composition(disk, proportions, genes = "geneB", components = "total_he", groups = groups)$data, p$data)
  disk$spatial <- c(geneB = file, geneA = "not-read.rds")
  expect_identical(plot_composition(disk, proportions, genes = "geneB", components = "total_he", groups = groups)$data, p$data)
})

test_that("composition inputs require unambiguous genes and location IDs", {
  skip_if_not_installed("ggplot2")
  d <- composition_example()
  expect_error(plot_composition(d$fit, d$proportions, genes = "missing"), "gene names")
  expect_error(plot_composition(d$fit, d$proportions, genes = c("geneA", "geneA")), "unique gene")
  expect_error(plot_composition(d$fit, d$proportions, components = "residual"), "Choose")
  expect_error(plot_composition(d$fit, d$proportions[-1, ]), "Every fitted location")
  duplicate <- d$proportions
  duplicate$spot_id[2] <- duplicate$spot_id[1]
  expect_error(plot_composition(d$fit, duplicate), "IDs must be unique")
  duplicate$spot_id[2] <- NA_character_
  expect_error(plot_composition(d$fit, duplicate), "complete")
  expect_error(plot_composition(d$fit, matrix(.5, 8, 2)), "row names")
  matrix_ids <- matrix(.5, 8, 2, dimnames = list(d$proportions$spot_id, c("spot_id", "tumor")))
  expect_error(plot_composition(d$fit, matrix_ids), "not an ID column")
  groups <- data.frame(spot_id = c("spot1", "missing", "spot3"), group_id = 1:3)
  proportions <- data.frame(group_id = 1:3, tumor = c(.1, .2, .3))
  expect_error(plot_composition(d$fit, proportions, groups = groups), "grouped location")
  groups$spot_id[2] <- groups$spot_id[1]
  expect_error(plot_composition(d$fit, proportions, groups = groups), "unique spot_id")
})

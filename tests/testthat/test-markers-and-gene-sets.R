marker_measurements <- function() {
  counts <- rbind(higher = c(8, 9, 8, 9, 0, 1, 0, 1, 10),
    lower = c(0, 1, 0, 1, 8, 9, 8, 9, 0),
    ties = c(2, 2, 3, 3, 1, 1, 2, 2, 4),
    discordant = c(10, 10, 10, 10, 0, 0, 0, 1000, 0),
    zero = rep(0, 9), constant = rep(5, 9), unannotated_only = c(rep(0, 8), 99))
  colnames(counts) <- paste0("spot", 1:9)
  groups <- data.frame(spot_id = colnames(counts),
    group = c(rep("invasive", 4), "immune", "stroma", "immune", "stroma", "unannotated"))
  list(counts = counts, expression = log1p(counts), groups = groups)
}

test_that("markers match tie-corrected Wilcoxon inference without continuity correction", {
  d <- marker_measurements()
  result <- find_markers(d$expression, d$groups, "invasive", c("immune", "stroma"))
  expected <- lapply(seq_len(nrow(d$expression)), function(i) {
    x <- d$expression[i, 1:4]
    y <- d$expression[i, 5:8]
    u <- sum(outer(x, y, ">")) + .5 * sum(outer(x, y, "=="))
    ranks <- rank(c(x, y), ties.method = "average")
    variance <- length(x) * length(y) / (length(ranks) * (length(ranks) - 1)) *
      sum((ranks - mean(ranks))^2)
    z <- if (variance > 0) (u - length(x) * length(y) / 2) / sqrt(variance) else NA_real_
    data.frame(gene = rownames(d$expression)[i], U = u, z = z,
      p = if (variance > 0) 2 * stats::pnorm(-abs(z)) else NaN,
      log2_fold_change = log2((mean(d$counts[i, 1:4]) + 1) / (mean(d$counts[i, 5:8]) + 1)),
      detection_fraction = mean(d$counts[i, ] > 0))
  })
  expected <- do.call(rbind, expected)
  expected$q <- stats::p.adjust(expected$p, "BH")
  expected <- expected[match(result$gene, expected$gene), ]
  rownames(expected) <- NULL
  expect_equal(result[names(expected)], expected, tolerance = 1e-12)
  expect_true(all(is.na(result$z[result$gene %in% c("zero", "constant", "unannotated_only")])))
  expect_true(all(is.na(result$p[result$gene %in% c("zero", "constant", "unannotated_only")])))
  expect_true(result$z[result$gene == "discordant"] > 0)
  expect_true(result$log2_fold_change[result$gene == "discordant"] < 0)
  expect_true(is.na(result$positive_rank[result$gene == "discordant"]))
  positive <- result$gene[is.finite(result$z) & result$z > 0 & result$log2_fold_change > 0]
  expect_equal(result$positive_rank, match(result$gene, positive))
})

test_that("marker results are invariant to ID order and sparse storage", {
  d <- marker_measurements()
  before <- d$expression
  expected <- find_markers(d$expression, d$groups, "invasive", c("immune", "stroma"))
  reordered <- find_markers(d$expression[, c(9, 2, 4, 6, 8, 1, 3, 5, 7)],
    d$groups[9:1, ], "invasive", c("stroma", "immune"))
  expect_equal(reordered, expected, tolerance = 1e-12)
  expect_identical(d$expression, before)
  raw <- find_markers(d$counts, d$groups, "invasive", c("immune", "stroma"), log_transformed = FALSE)
  expect_equal(raw, expected, tolerance = 1e-12)
  skip_if_not_installed("Matrix")
  sparse <- find_markers(Matrix::Matrix(d$expression, sparse = TRUE), d$groups,
    "invasive", c("immune", "stroma"))
  expect_equal(sparse, expected, tolerance = 1e-12)
})

test_that("marker comparisons require unambiguous IDs and measurements", {
  d <- marker_measurements()
  expect_error(find_markers(d$expression, d$groups[-1, ], "invasive", "immune"), "every expression column")
  duplicated <- d$groups
  duplicated$spot_id[1] <- duplicated$spot_id[2]
  expect_error(find_markers(d$expression, duplicated, "invasive", "immune"), "spot_id")
  expect_error(find_markers(d$expression, d$groups, "invasive", "invasive"), "distinct reference")
  negative <- d$expression
  negative[1, 1] <- -1
  expect_error(find_markers(negative, d$groups, "invasive", "immune"), "finite non-negative")
  negative[1, 1] <- NA_real_
  expect_error(find_markers(negative, d$groups, "invasive", "immune"), "finite non-negative")
})

gene_set_fits <- function() {
  values <- list(
    A1 = list(F32 = c(12, 14, NA, -1, 0, 1, 4, 8, Inf, 2), F43 = c(1, 2, 4, 8, 1, 2, 4, 8, NA, Inf)),
    G2 = list(F32 = c(1, 1, 2, 3, 2, 4, 4, 8, NA, Inf), F43 = c(20, 10, 12, 18, 1, 1, 2, 3, -1, NA)),
    H1 = list(F32 = c(0, 1, 2, 2, 0, 1, 2, 2, NA, Inf), F43 = c(0, 0, NA, 1, 2, 3, 4, 4, NA, -1))
  )
  lapply(values, function(x) structure(list(statistics = data.frame(gene = paste0("g", 1:10),
    F21 = 1:10, x, q32 = 1, q43 = 1)), class = "morphoanova_fit"))
}

test_that("gene-set tests use continuity correction and BH across six comparisons", {
  skip_if_not_installed("ggplot2")
  fits <- gene_set_fits()
  markers <- paste0("g", 1:4)
  p <- plot_gene_sets(fits, list(markers = markers))
  tests <- attr(p, "tests")
  expect_equal(nrow(tests), 6L)
  expected_p <- numeric(nrow(tests))
  for (i in seq_len(nrow(tests))) {
    row <- tests[i, ]
    column <- if (row$component == "location_specific_he") "F32" else "F43"
    statistics <- fits[[row$section]]$statistics
    f <- statistics[[column]]
    selected <- statistics$gene %in% markers
    valid <- is.finite(f) & f >= 0
    x <- f[valid & selected]
    y <- f[valid & !selected]
    u <- sum(outer(x, y, ">")) + .5 * sum(outer(x, y, "=="))
    ranks <- rank(c(x, y), ties.method = "average")
    variance <- length(x) * length(y) / (length(ranks) * (length(ranks) - 1)) *
      sum((ranks - mean(ranks))^2)
    difference <- u - length(x) * length(y) / 2
    z <- (difference - .5 * sign(difference)) / sqrt(variance)
    expected_p[i] <- 2 * stats::pnorm(-abs(z))
    expect_equal(row$W, u)
    expect_equal(row$probability_superiority, u / (length(x) * length(y)))
    expect_equal(row$rank_biserial, 2 * u / (length(x) * length(y)) - 1)
    expect_equal(c(row$n_selected, row$n_remaining), c(length(x), length(y)))
    expect_equal(c(row$excluded_selected, row$excluded_remaining),
      c(sum(selected & !valid), sum(!selected & !valid)))
    drawn <- p$data[p$data$section == row$section & p$data$component == row$component, ]
    expect_identical(drawn$gene, statistics$gene[valid])
    expect_equal(drawn$log10_F_plus_1, log10(f[valid] + 1))
  }
  expect_equal(tests$p, expected_p, tolerance = 1e-12)
  expect_equal(tests$q, stats::p.adjust(expected_p, "BH"), tolerance = 1e-12)
  expect_true(any(p$data$F == 0))
})

test_that("gene-set BH includes all supplied sets and missing tests stay missing", {
  skip_if_not_installed("ggplot2")
  fits <- gene_set_fits()
  sets <- list(first = paste0("g", 1:4), second = c("g3", "g5", "g7"))
  p <- plot_gene_sets(fits, sets)
  tests <- attr(p, "tests")
  expect_equal(nrow(tests), 12L)
  expected <- vapply(seq_len(nrow(tests)), function(i) {
    row <- tests[i, ]
    f <- fits[[row$section]]$statistics
    column <- if (row$component == "location_specific_he") "F32" else "F43"
    selected <- f$gene %in% sets[[row$gene_set]]
    values <- f[[column]]
    keep <- is.finite(values) & values >= 0
    stats::wilcox.test(values[selected & keep], values[!selected & keep], exact = FALSE, correct = TRUE)$p.value
  }, numeric(1))
  expect_equal(tests$q, stats::p.adjust(expected, "BH"), tolerance = 1e-12)
  missing <- plot_gene_sets(fits$A1, list(absent = "not_analyzed"))
  expect_true(all(is.na(attr(missing, "tests")$p)))
  expect_true(all(is.na(attr(missing, "tests")$q)))
  expect_true(all(attr(missing, "tests")$n_selected == 0L))
  fits$A1$statistics$F32 <- 1
  constant <- plot_gene_sets(fits$A1, list(markers = paste0("g", 1:4)), components = "location_specific_he")
  expect_true(is.na(attr(constant, "tests")$p))
  expect_true(is.na(attr(constant, "tests")$q))
})

test_that("gene-set effect-size filters use each raw component increment for both groups", {
  skip_if_not_installed("ggplot2")
  fits <- gene_set_fits()
  for (section in names(fits)) {
    fits[[section]]$statistics$r2_location_specific_he <- c(.03, .029, .5, .02, .04, .01, NA, .5, .04, -.01)
    fits[[section]]$statistics$r2_spatial_beyond_he <- c(.1, 0, .03, .9, .02, .03, .2, NA, .02, .4)
  }
  markers <- paste0("g", 1:4)
  for (cutoffs in list(.03, c(spatial_beyond_he = .1, location_specific_he = -.01))) {
    p <- plot_gene_sets(fits, list(markers = markers), min_r2 = cutoffs)
    tests <- attr(p, "tests")
    expected_p <- numeric(6)
    for (i in seq_len(nrow(tests))) {
      row <- tests[i, ]
      statistics <- fits[[row$section]]$statistics
      f <- statistics[[if (row$component == "location_specific_he") "F32" else "F43"]]
      r2 <- statistics[[paste0("r2_", row$component)]]
      threshold <- if (length(cutoffs) == 1L) cutoffs else unname(cutoffs[row$component])
      selected <- statistics$gene %in% markers
      valid_f <- is.finite(f) & f >= 0
      pass_r2 <- is.finite(r2) & r2 >= threshold
      keep <- valid_f & pass_r2
      x <- f[keep & selected]
      y <- f[keep & !selected]
      expected_p[i] <- stats::wilcox.test(x, y, exact = FALSE, correct = TRUE)$p.value
      drawn <- p$data[p$data$section == row$section & p$data$component == row$component, ]
      expect_identical(drawn$gene, statistics$gene[keep])
      expect_equal(row$min_r2, threshold)
      expect_equal(row$excluded_selected_invalid_f, sum(selected & !valid_f))
      expect_equal(row$excluded_remaining_invalid_f, sum(!selected & !valid_f))
      expect_equal(row$excluded_selected_r2, sum(selected & valid_f & !pass_r2))
      expect_equal(row$excluded_remaining_r2, sum(!selected & valid_f & !pass_r2))
      expect_equal(row$excluded_selected, row$excluded_selected_invalid_f + row$excluded_selected_r2)
      expect_equal(row$excluded_remaining, row$excluded_remaining_invalid_f + row$excluded_remaining_r2)
      expect_equal(row$n_selected + row$excluded_selected, row$n_selected_total)
      expect_equal(row$n_remaining + row$excluded_remaining, row$n_remaining_total)
    }
    expect_equal(tests$p, expected_p, tolerance = 1e-12)
    expect_equal(tests$q, stats::p.adjust(expected_p, "BH"), tolerance = 1e-12)
  }
  expect_error(plot_gene_sets(fits, list(markers = markers), min_r2 = NA_real_), "finite numbers")
  expect_error(plot_gene_sets(fits, list(markers = markers), min_r2 = c(.03, .1)), "named vector")
  expect_error(plot_gene_sets(fits, list(markers = markers), min_r2 = c(F32 = .03, F43 = .1)), "names must match")
})

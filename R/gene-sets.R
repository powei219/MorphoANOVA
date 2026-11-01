#' Compare component statistics for selected genes
#'
#' Compares each gene set with the remaining analyzed genes using a two-sided
#' Wilcoxon rank-sum test (normal approximation with continuity correction).
#' BH adjustment is applied jointly across the supplied fits, gene sets, and
#' components. Only finite, non-negative F statistics enter the comparison;
#' excluded counts are reported in the attached tests table. An optional
#' component-specific adjusted-R-squared threshold is applied equally to the
#' selected and remaining genes. No p-value or q-value filter is applied.
#'
#' @param fits A fitted MorphoANOVA result or a named list of results.
#' @param gene_sets A named list of gene-name vectors.
#' @param components Components to compare.
#' @param family Font family.
#' @param min_r2 Minimum raw adjusted-R-squared increment, or `NULL` for no
#'   effect-size filter. Supply one finite number for all components, or a
#'   named vector with one threshold for each selected component.
#' @return A ggplot of log10(F+1) for selected and remaining genes. Its `data`
#'   contains the plotted values; `attr(plot, "tests")` contains comparison
#'   statistics, p-values, BH-adjusted p-values, and included/excluded counts.
#'   Exclusions are split into invalid F statistics and otherwise valid F
#'   statistics failing the R-squared filter. The table also reports probability
#'   of superiority, `W / (n_selected * n_remaining)`, and the rank-biserial
#'   correlation, twice that probability minus one.
#' @export
plot_gene_sets <- function(fits, gene_sets,
                           components = c("location_specific_he", "spatial_beyond_he"),
                           family = "sans", min_r2 = NULL) {
  .plot_dependencies("ggplot2")
  if (inherits(fits, "morphoanova_fit")) fits <- list(Section = fits)
  if (!is.list(fits) || !length(fits) || is.null(names(fits)) ||
      anyNA(names(fits)) || any(!nzchar(names(fits))) || anyDuplicated(names(fits)) ||
      !all(vapply(fits, inherits, logical(1), "morphoanova_fit"))) {
    stop("fits must be a fit or a named list of fits.", call. = FALSE)
  }
  if (!is.list(gene_sets) || !length(gene_sets) || is.null(names(gene_sets)) ||
      anyNA(names(gene_sets)) || any(!nzchar(names(gene_sets))) || anyDuplicated(names(gene_sets)) ||
      !all(vapply(gene_sets, function(x) is.character(x) && length(x) > 0L &&
        !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x), logical(1)))) {
    stop("gene_sets must be a named list of unique gene-name vectors.", call. = FALSE)
  }
  test_columns <- c(slide_wide_he = "F21", location_specific_he = "F32", spatial_beyond_he = "F43")
  if (!is.character(components) || !length(components) || anyNA(components) ||
      anyDuplicated(components) || !all(components %in% names(test_columns))) {
    stop("Choose valid, unique component names.", call. = FALSE)
  }
  if (!is.null(min_r2)) {
    if (!is.numeric(min_r2) || !length(min_r2) || any(!is.finite(min_r2))) {
      stop("min_r2 must contain finite numbers.", call. = FALSE)
    }
    if (is.null(names(min_r2))) {
      if (length(min_r2) != 1L) stop("min_r2 must be a scalar or a named vector matching components.", call. = FALSE)
      min_r2 <- stats::setNames(rep(min_r2, length(components)), components)
    } else {
      if (anyNA(names(min_r2)) || anyDuplicated(names(min_r2)) ||
          length(min_r2) != length(components) || !setequal(names(min_r2), components)) {
        stop("min_r2 names must match components exactly.", call. = FALSE)
      }
      min_r2 <- min_r2[components]
    }
  }
  values <- tests <- list()
  for (section in names(fits)) for (set in names(gene_sets)) for (component in components) {
    statistics <- fits[[section]]$statistics
    required <- c("gene", test_columns[[component]], if (!is.null(min_r2)) paste0("r2_", component))
    if (!is.data.frame(statistics) || !all(required %in% names(statistics)) ||
        anyNA(statistics$gene) || anyDuplicated(statistics$gene)) {
      stop("fit$statistics must contain unique genes and the selected component columns.", call. = FALSE)
    }
    f <- statistics[[test_columns[[component]]]]
    selected <- statistics$gene %in% gene_sets[[set]]
    valid_f <- is.finite(f) & f >= 0
    valid_r2 <- rep(TRUE, nrow(statistics))
    if (!is.null(min_r2)) {
      r2 <- statistics[[paste0("r2_", component)]]
      valid_r2 <- is.finite(r2) & r2 >= min_r2[[component]]
    }
    valid <- valid_f & valid_r2
    x <- f[valid & selected]; y <- f[valid & !selected]
    test <- if (length(x) && length(y)) suppressWarnings(stats::wilcox.test(
      x, y, exact = FALSE, correct = TRUE)) else NULL
    i <- length(tests) + 1L
    tests[[i]] <- data.frame(section = section, gene_set = set, component = component,
      min_r2 = if (is.null(min_r2)) NA_real_ else min_r2[[component]],
      n_selected_total = sum(selected), n_remaining_total = sum(!selected),
      n_selected = length(x), n_remaining = length(y),
      excluded_selected = sum(selected & !valid), excluded_remaining = sum(!selected & !valid),
      excluded_selected_invalid_f = sum(selected & !valid_f),
      excluded_remaining_invalid_f = sum(!selected & !valid_f),
      excluded_selected_r2 = sum(selected & valid_f & !valid_r2),
      excluded_remaining_r2 = sum(!selected & valid_f & !valid_r2),
      W = if (is.null(test)) NA_real_ else unname(test$statistic),
      p = if (is.null(test)) NA_real_ else test$p.value,
      probability_superiority = if (is.null(test)) NA_real_ else unname(test$statistic) / (length(x) * length(y)))
    values[[i]] <- data.frame(section = rep(section, sum(valid)), gene_set = rep(set, sum(valid)),
      component = rep(component, sum(valid)), gene = statistics$gene[valid],
      group = ifelse(selected[valid], "Selected genes", "Remaining genes"), F = f[valid],
      log10_F_plus_1 = log10(f[valid] + 1))
  }
  tests <- do.call(rbind, tests)
  tests$q <- stats::p.adjust(tests$p, "BH")
  tests$rank_biserial <- 2 * tests$probability_superiority - 1
  values <- do.call(rbind, values)
  values$group <- factor(values$group, levels = c("Selected genes", "Remaining genes"))
  labels <- c(slide_wide_he = "Slide-wide H&E effect",
    location_specific_he = "Location-specific H&E effect", spatial_beyond_he = "Spatial variation beyond H&E")
  plot <- ggplot2::ggplot(values, ggplot2::aes(x = group, y = log10_F_plus_1, fill = group)) +
    ggplot2::geom_boxplot(width = .6, outlier.size = .5) +
    ggplot2::facet_grid(gene_set + section ~ component, scales = "free_y",
      labeller = ggplot2::labeller(component = ggplot2::as_labeller(labels))) +
    ggplot2::scale_fill_manual(values = c("#E6A0C4", "#9DCCB3")) +
    ggplot2::labs(x = NULL, y = "log10(F + 1)") +
    ggplot2::theme_classic(base_family = family) +
    ggplot2::theme(legend.position = "none")
  attr(plot, "tests") <- tests
  plot
}

utils::globalVariables(c("group", "log10_F_plus_1", "gene_set", "section"))

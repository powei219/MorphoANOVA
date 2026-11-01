#' Plot fitted spatial components
#'
#' The default four panels show expression, slide-wide H&E effect,
#' location-specific H&E effect, and spatial variation beyond H&E. The three
#' effects have separate, zero-centered color scales by default. Significance
#' maps test whether the total H&E contribution or spatial variation beyond
#' H&E equals zero at each location, using pointwise BH-adjusted p-values.
#' Total H&E is `slide_wide_he + location_specific_he`. BH correction is across
#' locations within each gene and component. Values above the display cap share
#' its darkest blue; stored p-values and q-values are unchanged.
#' When several genes are selected, component maps show their median at each
#' location. Genes are matched by spot ID before taking medians. Pointwise
#' significance and mechanism maps are available for a single gene. Mechanism
#' maps classify locations as H&E only, beyond H&E only, both, or neither using
#' the two pointwise BH-adjusted p-values and `alpha`, without a gene-level F
#' gate. Locations with either adjusted p-value unavailable remain unclassified.
#'
#' @param fit A result from [fit_morphoanova()].
#' @param gene Gene or vector of genes to display. May be omitted for a
#'   single-gene result. Multiple genes produce median component maps.
#' @param type `"components"`, `"significance"`, or `"mechanism"`.
#' @param point_size Spatial point size.
#' @param alpha Pointwise BH significance threshold; includes equality.
#' @param flip_y Reverse the vertical axis for image-row coordinates.
#' @param family Plot font family.
#' @param shared_scale Use the same symmetric color range for all selected
#'   component maps. Expression always has its own scale.
#' @param significance_cap Maximum displayed \code{-log10(BH q)} in significance
#'   maps. The default is 20; the colorbar labels the upper endpoint \code{20+}.
#' @param components Components to display, in panel order. For component maps,
#'   choose from `"expression"`, `"slide_wide_he"`, `"location_specific_he"`,
#'   `"total_he"`, and `"spatial_beyond_he"`. The default is expression and the
#'   three separate components. For significance maps, choose `"total_he"`
#'   and/or `"spatial_beyond_he"`; both are shown by default. Mechanism maps
#'   require both components.
#' @return A patchwork plot, suitable for `ggplot2::ggsave()`. Its `values`
#'   attribute contains the spot IDs, coordinates, and values shown in each
#'   panel. For significance maps, the `<component>_significance` columns hold
#'   the displayed negative log10 BH-adjusted p-values after thresholding and
#'   capping. Mechanism-map tables contain both adjusted p-value columns and
#'   the `mechanism` classification. Retrieve this table with `attr(plot, "values")`.
#' @export
plot_components <- function(fit, gene = NULL, type = c("components", "significance", "mechanism"),
                            point_size = 0.8, alpha = 0.05, flip_y = FALSE,
                            family = "Arial", shared_scale = FALSE,
                            significance_cap = 20, components = NULL) {
  .plot_dependencies(c("ggplot2", "patchwork"))
  type <- match.arg(type)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha > 1) {
    stop("alpha must be in (0, 1].", call. = FALSE)
  }
  if (!is.logical(shared_scale) || length(shared_scale) != 1L || is.na(shared_scale)) {
    stop("shared_scale must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(significance_cap) || length(significance_cap) != 1L ||
      !is.finite(significance_cap) || significance_cap <= 0) {
    stop("significance_cap must be a positive finite number.", call. = FALSE)
  }
  if (is.null(gene)) {
    if (nrow(fit$statistics) != 1L) stop("Choose one or more `gene` values from this multi-gene fit.", call. = FALSE)
    gene <- fit$statistics$gene[[1L]]
  }
  if (!is.character(gene) || !length(gene) || anyNA(gene) || anyDuplicated(gene) ||
      !all(gene %in% fit$statistics$gene)) {
    stop("`gene` must name one or more distinct fitted genes.", call. = FALSE)
  }
  if (type != "components" && length(gene) != 1L) {
    stop("Significance and mechanism maps require a single gene; pointwise p-values are not pooled across genes.", call. = FALSE)
  }
  labels <- c(expression = "Gene expression", slide_wide_he = "Slide-wide H&E effect",
    location_specific_he = "Location-specific\nH&E effect", total_he = "Total H&E effect",
    spatial_beyond_he = "Spatial variation\nbeyond H&E")
  allowed <- if (type == "components") names(labels) else c("total_he", "spatial_beyond_he")
  if (is.null(components)) {
    components <- if (type == "components") {
      c("expression", "slide_wide_he", "location_specific_he", "spatial_beyond_he")
    } else allowed
  }
  if (!is.character(components) || !length(components) || anyNA(components) ||
      anyDuplicated(components) || !all(components %in% allowed)) {
    stop("Choose distinct components from: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
  if (type == "mechanism" && !setequal(components, allowed)) {
    stop("Mechanism maps require both total_he and spatial_beyond_he.", call. = FALSE)
  }
  columns <- if (type == "components") components else paste0(components, "_q")
  spatial <- .gene_spatial(fit, gene = gene[1L])
  if (!all(columns %in% names(spatial))) stop("Requested component values are missing from the fit.", call. = FALSE)
  if (length(gene) > 1L) {
    reference <- spatial[c("spot_id", "x", "y")]
    estimates <- stats::setNames(lapply(columns, function(column) {
      matrix(NA_real_, nrow(reference), length(gene))
    }), columns)
    for (i in seq_along(gene)) {
      values <- if (i == 1L) spatial else .gene_spatial(fit, gene = gene[i])
      order <- match(reference$spot_id, values$spot_id)
      if (nrow(values) != nrow(reference) || anyNA(order) || anyDuplicated(values$spot_id) ||
          !isTRUE(all.equal(unname(as.matrix(values[order, c("x", "y")])),
                            unname(as.matrix(reference[c("x", "y")])), tolerance = 0))) {
        stop("Gene spatial tables must have matching spot IDs and coordinates.", call. = FALSE)
      }
      if (!all(columns %in% names(values))) stop("Requested component values are missing from the fit.", call. = FALSE)
      for (column in columns) estimates[[column]][, i] <- values[[column]][order]
    }
    spatial <- reference
    for (column in columns) spatial[[column]] <- apply(estimates[[column]], 1L, stats::median)
    rm(estimates, values)
  }
  plotted <- spatial[c("spot_id", "x", "y")]
  blues <- c("#C7DFF0", "#94BFDC", "#609CC7", "#3D73AA", "#2B477A")
  panels <- list()
  if (type == "components") {
    if (length(gene) == 1L) {
      statistics <- fit$statistics[match(gene, fit$statistics$gene), , drop = FALSE]
      increments <- c(slide_wide_he = statistics$r2_slide_wide_he,
        location_specific_he = statistics$r2_location_specific_he,
        spatial_beyond_he = statistics$r2_spatial_beyond_he,
        total_he = statistics$r2_slide_wide_he + statistics$r2_location_specific_he)
      subtitles <- c(expression = "Observed", stats::setNames(
        sprintf("%.1f%% variance explained", 100 * increments), names(increments)))
    } else {
      subtitles <- stats::setNames(rep(paste("Median across", length(gene), "genes"), length(columns)), columns)
    }
    effect_columns <- setdiff(columns, "expression")
    effect_limits <- vapply(effect_columns, function(column) max(abs(spatial[[column]]), na.rm = TRUE), numeric(1))
    if (shared_scale && length(effect_limits)) effect_limits[] <- max(effect_limits)
    effect_limits[!is.finite(effect_limits) | effect_limits == 0] <- 1
    for (i in seq_along(columns)) {
      spatial$value <- spatial[[columns[i]]]
      plotted[[columns[i]]] <- spatial$value
      p <- ggplot2::ggplot(spatial, ggplot2::aes(x, y, colour = value)) +
        ggplot2::geom_point(size = point_size, shape = 16, stroke = 0) +
        ggplot2::coord_equal() + ggplot2::labs(title = labels[[columns[i]]], subtitle = subtitles[[columns[i]]]) +
        ggplot2::theme_void(base_family = family) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
          plot.subtitle = ggplot2::element_text(size = 8),
          legend.title = ggplot2::element_text(size = 8),
          legend.text = ggplot2::element_text(size = 7), legend.position = "right")
      if (columns[i] == "expression") {
        p <- p + ggplot2::scale_colour_gradientn(colours = blues, name = "Expression")
      } else {
        bound <- effect_limits[[columns[i]]]
        p <- p + ggplot2::scale_colour_gradient2(low = "#3AB370", mid = "#EAE7CC", high = "#FD1593",
          limits = c(-bound, bound), midpoint = 0, name = "Value")
      }
      if (flip_y) p <- p + ggplot2::scale_y_reverse()
      p <- p + ggplot2::guides(colour = ggplot2::guide_colourbar(
        title.position = "top", barwidth = grid::unit(2, "mm"), barheight = grid::unit(18, "mm")))
      panels[[i]] <- p
    }
  } else if (type == "significance") {
    labels[["spatial_beyond_he"]] <- "Spatial variation beyond H&E"
    breaks <- c(0, significance_cap / 2, significance_cap)
    break_labels <- c("0", format(significance_cap / 2, trim = TRUE), paste0(format(significance_cap, trim = TRUE), "+"))
    for (i in seq_along(columns)) {
      q <- spatial[[columns[i]]]
      spatial$value <- ifelse(is.finite(q) & q <= alpha,
        pmin(-log10(pmax(q, .Machine$double.xmin)), significance_cap), NA_real_)
      plotted[[paste0(components[i], "_significance")]] <- spatial$value
      p <- ggplot2::ggplot(spatial, ggplot2::aes(x, y, colour = value)) +
        ggplot2::geom_point(size = point_size, shape = 16, stroke = 0) +
        ggplot2::scale_colour_gradientn(colours = blues, limits = c(0, significance_cap),
          breaks = breaks, labels = break_labels, na.value = "grey90", name = expression(-log[10]("BH q"))) +
        ggplot2::coord_equal() + ggplot2::labs(title = labels[[components[i]]]) +
        ggplot2::theme_void(base_family = family) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
          legend.title = ggplot2::element_text(size = 8),
          legend.text = ggplot2::element_text(size = 7), legend.position = "right")
      if (flip_y) p <- p + ggplot2::scale_y_reverse()
      p <- p + ggplot2::guides(colour = ggplot2::guide_colourbar(
        direction = "horizontal", title.position = "left",
        barwidth = grid::unit(35, "mm"), barheight = grid::unit(2, "mm")))
      panels[[i]] <- p
    }
  } else {
    plotted$total_he_q <- spatial$total_he_q
    plotted$spatial_beyond_he_q <- spatial$spatial_beyond_he_q
    available <- is.finite(plotted$total_he_q) & is.finite(plotted$spatial_beyond_he_q)
    he_significant <- plotted$total_he_q <= alpha
    beyond_significant <- plotted$spatial_beyond_he_q <= alpha
    mechanism <- rep(NA_character_, nrow(plotted))
    mechanism[available] <- "Neither"
    mechanism[available & he_significant] <- "H&E only"
    mechanism[available & beyond_significant] <- "Beyond H&E only"
    mechanism[available & he_significant & beyond_significant] <- "Both"
    colours <- c("H&E only" = "#3AB370", "Beyond H&E only" = "#FD1593",
      "Both" = "#8B5AA9", "Neither" = "grey85")
    plotted$mechanism <- factor(mechanism, levels = names(colours))
    p <- ggplot2::ggplot(plotted, ggplot2::aes(x, y, colour = mechanism)) +
      ggplot2::geom_point(size = point_size, shape = 16, stroke = 0, show.legend = TRUE) +
      ggplot2::scale_colour_manual(values = colours, drop = FALSE,
        na.value = "grey50", name = NULL,
        labels = function(x) ifelse(is.na(x), "Unavailable", x)) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) +
      ggplot2::coord_equal() + ggplot2::labs(title = "Pointwise component classification") +
      ggplot2::theme_void(base_family = family) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 7), legend.position = "bottom")
    if (flip_y) p <- p + ggplot2::scale_y_reverse()
    panels[[1L]] <- p
  }
  p <- patchwork::wrap_plots(panels, nrow = 1,
    guides = if (type == "significance") "collect" else "keep") +
    patchwork::plot_annotation(title = if (length(gene) == 1L) gene else paste("Median of", length(gene), "genes"),
      theme = ggplot2::theme(text = ggplot2::element_text(family = family)))
  if (type == "significance") p <- p & ggplot2::theme(legend.position = "bottom")
  attr(p, "values") <- plotted
  p
}

#' Plot component evidence across genes
#'
#' Labels are selected
#' separately from the F32 and F43 rankings after each corresponding raw
#' adjusted-R-squared gain passes `min_r2`, using [rank_genes()]. Selected F32 genes are
#' marked orange and F43 genes red, with red taking precedence for overlap.
#' @param fit A result from [fit_morphoanova()].
#' @param top_n Number of labels from each component ranking. Use zero for none.
#'   A scalar applies to both components; a named vector supplies
#'   `location_specific_he` and `spatial_beyond_he` separately.
#' @param min_r2 Minimum corresponding adjusted-R-squared increment for labels.
#'   A scalar applies to both components; a named vector supplies
#'   `location_specific_he` and `spatial_beyond_he` separately.
#' @param labels Optional gene names to label instead of selecting top genes.
#'   Names outside the selected component rankings are marked in neutral gray.
#' @param point_size Scatter point size.
#' @param family Plot font family.
#' @return A ggplot, suitable for `ggplot2::ggsave()`.
#' @export
plot_scatter <- function(fit, top_n = 3, min_r2 = 0.03, labels = NULL,
                         point_size = 0.5, family = "Arial") {
  .plot_dependencies("ggplot2")
  if (!inherits(fit, "morphoanova_fit")) stop("fit must come from fit_morphoanova().", call. = FALSE)
  ranked <- rank_genes(fit, top_n = top_n, min_r2 = min_r2)
  d <- fit$statistics
  selected32 <- ranked$location_specific_he$gene
  selected43 <- ranked$spatial_beyond_he$gene
  if (is.null(labels)) {
    labels <- unique(c(selected32, selected43))
  }
  keep <- is.finite(d$F32) & is.finite(d$F43) & d$F32 > -1 & d$F43 > -1
  if (any(!keep)) warning(sum(!keep), " genes have unavailable log10(F + 1) coordinates.", call. = FALSE)
  d <- d[keep, , drop = FALSE]
  d$log_F32 <- log10(d$F32 + 1)
  d$log_F43 <- log10(d$F43 + 1)
  marked <- d[d$gene %in% labels, , drop = FALSE]
  marked$label_colour <- ifelse(marked$gene %in% selected43, "#E41A1C",
    ifelse(marked$gene %in% selected32, "#FF7F00", "#444444"))
  p <- ggplot2::ggplot(d, ggplot2::aes(log_F32, log_F43, colour = F21)) +
    ggplot2::geom_point(size = point_size, alpha = 0.6, stroke = 0) +
    ggplot2::scale_colour_viridis_c(name = "Slide-wide H&E effect (F21)", na.value = "grey70") +
    ggplot2::labs(x = bquote(atop("Location-specific H&E effect", log[10](F[32] + 1))),
      y = bquote(atop("Spatial variation beyond H&E", log[10](F[43] + 1)))) +
    ggplot2::theme_classic(base_family = family, base_size = 9) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 7),
      legend.box.spacing = grid::unit(1.5, "mm")) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(title.position = "top",
      barwidth = grid::unit(2, "mm"), barheight = grid::unit(18, "mm")))
  if (nrow(marked)) {
    p <- p + ggplot2::geom_point(data = marked, ggplot2::aes(log_F32, log_F43),
      inherit.aes = FALSE, shape = 21, fill = NA, colour = marked$label_colour,
      size = max(2, point_size + 1), stroke = 0.4, show.legend = FALSE)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(data = marked, ggplot2::aes(log_F32, log_F43, label = gene),
        inherit.aes = FALSE, colour = marked$label_colour, family = family, size = 2.5,
        seed = 1, box.padding = 0.3, point.padding = 0.25,
        min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf, show.legend = FALSE)
    } else {
      p <- p + ggplot2::geom_text(data = marked, ggplot2::aes(log_F32, log_F43, label = gene),
        inherit.aes = FALSE, colour = marked$label_colour, family = family, size = 2.5,
        vjust = -0.5, check_overlap = TRUE, show.legend = FALSE)
    }
  }
  p
}

.plot_dependencies <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install package '", package, "' to draw this plot.", call. = FALSE)
    }
  }
}

utils::globalVariables(c("x", "y", "value", "mechanism", "log_F32", "log_F43", "F21", "gene"))

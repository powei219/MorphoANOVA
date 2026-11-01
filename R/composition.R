#' Plot associations with tissue or cell-type proportions
#'
#' Examines how the total H&E-associated effect and additional spatial
#' variation beyond the selected H&E feature relate to tissue composition.
#' Spearman correlations with cell-type proportions are
#' calculated for each gene and displayed as a heatmap, with spots matched
#' by ID.
#'
#' @param fit A result from [fit_morphoanova()].
#' @param proportions A data frame containing `spot_id` and numeric proportion
#'   columns, one per tissue category or cell type. A numeric matrix with
#'   location IDs as row names is also accepted. Rows may be in any order;
#'   every analyzed location must be present.
#' @param genes Gene names to display. `NULL` uses all genes in `fit`.
#' @param components Components to correlate: `"total_he"`,
#'   `"spatial_beyond_he"`, and optionally `"expression"`.
#' @param groups Optional data frame with unique `spot_id` and a `group_id`
#'   for each included location. Components are averaged within each group
#'   before correlation. In this case `proportions` needs `group_id` instead
#'   of `spot_id`, or group IDs as matrix row names. Locations not listed in
#'   `groups` are excluded.
#' @param family Plot font family.
#' @return A ggplot. Its `data` table contains `gene`, `component`, `cell_type`,
#'   `n` (number of finite pairs), and signed Spearman `rho`.
#' @details Proportions are used as supplied.
#' @examples
#' \dontrun{
#' plot_composition(fit, ovarian$proportions, genes = c("IGHM", "PRG4"))
#' }
#' @export
plot_composition <- function(fit, proportions, genes = NULL,
                             components = c("total_he", "spatial_beyond_he"),
                             groups = NULL, family = "sans") {
  .plot_dependencies("ggplot2")
  correlations <- .component_correlations(fit, proportions, genes, components, groups)
  correlations$component <- factor(correlations$component, levels = components)
  labels <- c(expression = "Gene expression", total_he = "Total H&E effect",
              spatial_beyond_he = "Spatial variation beyond H&E")
  ggplot2::ggplot(correlations, ggplot2::aes(x = cell_type, y = gene, fill = rho)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~ component, labeller = ggplot2::as_labeller(labels)) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                 limits = c(-1, 1), name = "Spearman rho") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_family = family) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

.component_correlations <- function(fit, proportions, genes, components, groups) {
  if (!inherits(fit, "morphoanova_fit")) stop("fit must come from fit_morphoanova().", call. = FALSE)
  if (is.null(genes)) genes <- fit$statistics$gene
  if (!is.character(genes) || !length(genes) || anyNA(genes) ||
      any(!nzchar(genes)) || anyDuplicated(genes) ||
      !all(genes %in% fit$statistics$gene)) {
    stop("genes must contain unique gene names present in fit.", call. = FALSE)
  }
  if (!is.character(components) || !length(components) || anyNA(components) ||
      anyDuplicated(components) ||
      !all(components %in% c("expression", "total_he", "spatial_beyond_he"))) {
    stop("Choose expression, total_he, and/or spatial_beyond_he as components.", call. = FALSE)
  }
  id_column <- if (is.null(groups)) "spot_id" else "group_id"
  if (is.matrix(proportions)) {
    ids <- rownames(proportions)
    if (is.null(ids)) stop("A proportions matrix needs location IDs as row names.", call. = FALSE)
    if (id_column %in% colnames(proportions)) {
      stop("A proportions matrix must use row names for IDs, not an ID column.", call. = FALSE)
    }
    proportions <- data.frame(proportions, check.names = FALSE)
    proportions[[id_column]] <- ids
  }
  if (!is.data.frame(proportions) || !id_column %in% names(proportions)) {
    stop("proportions needs a ", id_column, " column.", call. = FALSE)
  }
  if (anyNA(names(proportions)) || any(!nzchar(names(proportions))) ||
      anyDuplicated(names(proportions))) {
    stop("proportions column names must be unique and complete.", call. = FALSE)
  }
  ids <- as.character(proportions[[id_column]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Proportion IDs must be unique and complete.", call. = FALSE)
  }
  cell_types <- setdiff(names(proportions), id_column)
  if (!length(cell_types) || !all(vapply(proportions[cell_types], is.numeric, logical(1)))) {
    stop("Supply numeric cell-type columns in proportions.", call. = FALSE)
  }
  if (!is.null(groups)) {
    if (!is.data.frame(groups) || !all(c("spot_id", "group_id") %in% names(groups))) {
      stop("groups needs spot_id and group_id columns.", call. = FALSE)
    }
    groups <- groups[c("spot_id", "group_id")]
    groups[] <- lapply(groups, as.character)
    if (!nrow(groups) || anyNA(groups) || any(!nzchar(as.matrix(groups))) ||
        anyDuplicated(groups$spot_id)) {
      stop("groups needs unique spot_id and complete group_id columns.", call. = FALSE)
    }
  }

  result <- list()
  for (gene in genes) {
    values <- .gene_spatial(fit, gene = gene)
    component_values <- values[components]
    location_ids <- as.character(values$spot_id)
    if (!is.null(groups)) {
      selected <- match(groups$spot_id, location_ids)
      if (anyNA(selected)) stop("Every grouped location must occur in the fit.", call. = FALSE)
      component_values <- component_values[selected, , drop = FALSE]
      component_values <- stats::aggregate(component_values, list(group_id = groups$group_id), mean)
      location_ids <- as.character(component_values$group_id)
    }
    matched <- match(location_ids, ids)
    if (anyNA(matched)) stop("Every fitted location or group needs cell-type proportions.", call. = FALSE)
    for (component in components) {
      for (cell_type in cell_types) {
        x <- component_values[[component]]
        y <- proportions[[cell_type]][matched]
        keep <- is.finite(x) & is.finite(y)
        n <- sum(keep)
        rho <- NA_real_
        if (n >= 3L && length(unique(x[keep])) > 1L && length(unique(y[keep])) > 1L) {
          rho <- stats::cor(x[keep], y[keep], method = "spearman")
        }
        result[[length(result) + 1L]] <- data.frame(gene = gene, component = component,
          cell_type = cell_type, n = n, rho = rho, stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, result)
}

utils::globalVariables(c("cell_type", "rho", "component"))

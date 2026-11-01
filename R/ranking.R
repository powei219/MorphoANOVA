#' Rank genes by component evidence
#'
#' For each component, retain genes whose adjusted-R-squared increment is at
#' least `min_r2`, then rank them by decreasing F statistic. Ties are ordered
#' by gene name. The default components use F32 and F43; `slide_wide_he` uses
#' F21. The first `top_n` genes are returned after filtering and ordering.
#'
#' @param fit A result from [fit_morphoanova()].
#' @param top_n Maximum number of genes to return per component. Use one
#'   non-negative whole number for all components, or a named vector with
#'   one value for every selected component.
#' @param min_r2 Minimum adjusted-R-squared increment, expressed as a fraction:
#'   `0.03` means 3%. Use one finite number for all components, or a named
#'   vector with one value for every selected component.
#' @param components Components to rank: `location_specific_he`,
#'   `spatial_beyond_he`, and optionally `slide_wide_he`.
#' @return An ordinary list with one data frame per component, containing
#'   `gene`, `rank`, `F`, `p`, `q`, and `r2`. `universe` contains all genes in
#'   `fit$statistics`, before ranking filters, and provides the default
#'   background for [plot_enrichment()]. It retains genes with unavailable
#'   component tests or failed fits. `settings` records the components and
#'   thresholds. Ranked tables include complete or partial fits when their
#'   component F statistic, p-value, and R-squared increment are available.
#' @export
rank_genes <- function(fit, top_n = 100, min_r2 = 0.03,
                       components = c("location_specific_he", "spatial_beyond_he")) {
  if (!inherits(fit, "morphoanova_fit")) stop("fit must come from fit_morphoanova().", call. = FALSE)
  tests <- c(slide_wide_he = "21", location_specific_he = "32", spatial_beyond_he = "43")
  if (!is.character(components) || !length(components) || anyNA(components) ||
      anyDuplicated(components) || !all(components %in% names(tests))) {
    stop("Choose unique components from slide_wide_he, location_specific_he, and spatial_beyond_he.", call. = FALSE)
  }

  parameters <- list(top_n = top_n, min_r2 = min_r2)
  for (parameter in names(parameters)) {
    value <- parameters[[parameter]]
    if (!is.numeric(value) || !length(value) || any(!is.finite(value))) {
      stop(parameter, " must contain finite numbers.", call. = FALSE)
    }
    if (is.null(names(value))) {
      if (length(value) != 1L) stop(parameter, " must be a scalar or a named vector matching components.", call. = FALSE)
      value <- stats::setNames(rep(value, length(components)), components)
    } else {
      if (anyNA(names(value)) || anyDuplicated(names(value)) ||
          length(value) != length(components) || !setequal(names(value), components)) {
        stop(parameter, " names must match components exactly.", call. = FALSE)
      }
      value <- value[components]
    }
    if (parameter == "top_n" && any(value < 0 | value != floor(value))) {
      stop("top_n must contain non-negative whole numbers.", call. = FALSE)
    }
    parameters[[parameter]] <- value
  }

  statistics <- fit$statistics
  required <- c("gene", paste0("F", tests[components]), paste0("p", tests[components]),
    paste0("q", tests[components]), paste0("r2_", components))
  if (!is.data.frame(statistics) || !all(required %in% names(statistics))) {
    stop("fit$statistics is missing the selected component columns.", call. = FALSE)
  }
  if (anyNA(statistics$gene) || anyDuplicated(statistics$gene)) {
    stop("fit$statistics gene names must be unique and complete.", call. = FALSE)
  }
  status <- fit$status
  if (!is.data.frame(status) || !all(c("gene", "status") %in% names(status)) ||
      anyNA(status$gene) || anyDuplicated(status$gene) ||
      anyNA(match(statistics$gene, status$gene))) {
    stop("fit$status must contain one status for every gene in fit$statistics.", call. = FALSE)
  }
  fitted_genes <- status$status[match(statistics$gene, status$gene)] %in% c("complete", "partial")
  result <- list()
  for (component in components) {
    test <- tests[[component]]
    F <- statistics[[paste0("F", test)]]
    p <- statistics[[paste0("p", test)]]
    r2 <- statistics[[paste0("r2_", component)]]
    valid <- fitted_genes & is.finite(F) & is.finite(p) & p >= 0 & p <= 1 & is.finite(r2)
    selected <- which(valid & r2 >= parameters$min_r2[[component]])
    selected <- selected[order(-F[selected], statistics$gene[selected])]
    selected <- utils::head(selected, parameters$top_n[[component]])
    result[[component]] <- data.frame(gene = as.character(statistics$gene[selected]),
      rank = seq_along(selected), F = F[selected], p = p[selected],
      q = statistics[[paste0("q", test)]][selected], r2 = r2[selected],
      row.names = NULL, stringsAsFactors = FALSE)
  }
  result$universe <- as.character(statistics$gene)
  result$settings <- c(list(components = components), parameters)
  result
}

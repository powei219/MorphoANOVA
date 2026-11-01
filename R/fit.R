.empty_statistics <- function(gene) {
  out <- data.frame(gene = gene, stringsAsFactors = FALSE)
  columns <- c("F21", "p21", "q21", "F32", "p32", "q32", "F43", "p43", "q43",
               "r2_slide_wide_he", "r2_location_specific_he", "r2_spatial_beyond_he",
               "r2_full", "r2_unexplained_variance")
  for (column in columns) out[[column]] <- NA_real_
  out
}

.sequential_test <- function(reduced, full) {
  result <- list(F = NA_real_, p = NA_real_, df1 = NA_real_, df2 = NA_real_,
                 raw_F = NA_real_, dispersion = NA_real_,
                 deviance_difference = NA_real_,
                 residual_df = c(reduced = NA_real_, full = NA_real_), reason = "")
  comparison <- tryCatch(mgcv::anova.gam(reduced, full, test = "F"), error = identity)
  if (inherits(comparison, "error")) {
    result$reason <- paste("GAM ANOVA failed:", conditionMessage(comparison))
    return(result)
  }
  if (!is.data.frame(comparison) || nrow(comparison) != 2L ||
      !all(c("Resid. Df", "Df", "Deviance") %in% names(comparison))) {
    result$reason <- "GAM ANOVA did not return a two-model comparison"
    return(result)
  }
  result$residual_df[] <- comparison[["Resid. Df"]]
  result$df1 <- comparison[["Df"]][2L]
  result$deviance_difference <- comparison[["Deviance"]][2L]

  # Save the ANOVA scale and degrees of freedom alongside its returned F and p.
  big <- order(result$residual_df)[1L]
  bigmodel <- list(reduced, full)[[big]]
  result$dispersion <- summary(bigmodel)$dispersion
  family_dispersion <- bigmodel$family$dispersion
  result$df2 <- if (!is.null(family_dispersion)) {
    if (is.na(family_dispersion)) result$residual_df[big] else Inf
  } else if (isTRUE(result$dispersion == 1)) Inf else min(result$residual_df)
  result$df2 <- unname(result$df2)
  if (is.finite(result$df1) && result$df1 != 0 &&
      is.finite(result$dispersion) && result$dispersion > 0) {
    result$raw_F <- result$deviance_difference / result$df1 / result$dispersion
  }
  if (all(c("F", "Pr(>F)") %in% names(comparison))) {
    result$F <- comparison[["F"]][2L]
    result$p <- comparison[["Pr(>F)"]][2L]
  }
  if (is.na(result$F) || is.na(result$p)) {
    result$reason <- if (isTRUE(result$df1 == 0)) {
      "GAM ANOVA F test unavailable: zero difference in test degrees of freedom"
    } else if (is.finite(result$raw_F) && result$raw_F < 0) {
      "GAM ANOVA F test unavailable: negative algebraic F statistic"
    } else "GAM ANOVA F test unavailable"
  }
  result
}

.component_se <- function(design, covariance) {
  variance <- rowSums((design %*% covariance) * design)
  tolerance <- 1e-10 * max(1, max(abs(variance)))
  if (any(!is.finite(variance)) || any(variance < -tolerance)) {
    stop("Invalid component covariance.", call. = FALSE)
  }
  sqrt(pmax(variance, 0))
}

.pointwise_p <- function(estimate, se) {
  p <- rep(NA_real_, length(estimate))
  valid <- is.finite(estimate) & is.finite(se) & se > 0
  p[valid] <- 2 * stats::pnorm(-abs(estimate[valid] / se[valid]))
  p
}

.fit_one_gene <- function(expression, spots, h, gene, k) {
  y <- as.numeric(expression[gene, , drop = TRUE])
  if (!is.finite(stats::var(y)) || stats::var(y) <= 0) {
    stop("Gene expression has no variation.", call. = FALSE)
  }
  d <- data.frame(y = y, x = spots$x, y_coord = spots$y, h = h)

  m1 <- mgcv::gam(y ~ 1, data = d, family = stats::gaussian(), method = "REML")
  m2 <- mgcv::gam(y ~ h, data = d, family = stats::gaussian(), method = "REML")
  m3 <- mgcv::gam(y ~ s(x, y_coord, by = h, k = k, bs = "tp"),
                  data = d, family = stats::gaussian(), method = "REML")
  m4 <- mgcv::gam(y ~ s(x, y_coord, by = h, k = k, bs = "tp") +
                    s(x, y_coord, k = k, bs = "tp"),
                  data = d, family = stats::gaussian(), method = "REML")
  models <- list(M1 = m1, M2 = m2, M3 = m3, M4 = m4)
  if (any(!vapply(models, function(model) isTRUE(model$converged), logical(1)))) {
    stop("A nested GAM did not converge.", call. = FALSE)
  }
  n <- length(y)
  diagnostics <- data.frame(
    model = names(models),
    rss = vapply(models, function(model) sum(stats::residuals(model, type = "response")^2), numeric(1)),
    edf = vapply(models, function(model) sum(model$edf), numeric(1)),
    stringsAsFactors = FALSE
  )
  diagnostics$residual_df <- n - diagnostics$edf
  diagnostics$adj_r2 <- 1 - (diagnostics$rss / diagnostics$residual_df) / stats::var(y)
  diagnostics$adj_r2[1] <- 0
  comparisons <- list(
    `21` = .sequential_test(m1, m2),
    `32` = .sequential_test(m2, m3),
    `43` = .sequential_test(m3, m4)
  )
  statistics <- .empty_statistics(gene)
  for (test in names(comparisons)) {
    statistics[[paste0("F", test)]] <- comparisons[[test]]$F
    statistics[[paste0("p", test)]] <- comparisons[[test]]$p
  }
  statistics$r2_slide_wide_he <- diagnostics$adj_r2[2]
  statistics$r2_location_specific_he <- diagnostics$adj_r2[3] - diagnostics$adj_r2[2]
  statistics$r2_spatial_beyond_he <- diagnostics$adj_r2[4] - diagnostics$adj_r2[3]
  statistics$r2_full <- diagnostics$adj_r2[4]
  statistics$r2_unexplained_variance <- 1 - diagnostics$adj_r2[4]

  # The numeric-by smooth includes the constant coefficient direction.
  # Split that coefficient surface after fitting, preserving the model space.
  he_smooth <- m4$smooth[[which(vapply(m4$smooth, function(smooth) identical(smooth$by, "h"), logical(1)))]]
  spatial_smooth <- m4$smooth[[which(vapply(m4$smooth, function(smooth) identical(smooth$by, "NA"), logical(1)))]]
  he_index <- seq.int(he_smooth$first.para, he_smooth$last.para)
  spatial_index <- seq.int(spatial_smooth$first.para, spatial_smooth$last.para)
  he_coef <- stats::coef(m4)[he_index]
  spatial_coef <- stats::coef(m4)[spatial_index]
  he_cov <- m4$Vp[he_index, he_index, drop = FALSE]
  spatial_cov <- m4$Vp[spatial_index, spatial_index, drop = FALSE]

  coefficient <- spatial_raw <- numeric(n)
  mean_spatial_basis <- numeric(length(spatial_index))
  blocks <- split(seq_len(n), ceiling(seq_len(n) / 2048L))
  for (index in blocks) {
    newdata <- d[index, , drop = FALSE]
    newdata$h <- 1
    coefficient[index] <- as.numeric(mgcv::PredictMat(he_smooth, newdata) %*% he_coef)
    basis <- mgcv::PredictMat(spatial_smooth, newdata)
    spatial_raw[index] <- as.numeric(basis %*% spatial_coef)
    mean_spatial_basis <- mean_spatial_basis + colSums(basis) / n
  }
  slide_wide_coefficient <- mean(coefficient)
  location_specific_coefficient <- coefficient - slide_wide_coefficient
  slide_wide_he <- slide_wide_coefficient * h
  location_specific_he <- location_specific_coefficient * h
  total_he <- slide_wide_he + location_specific_he
  spatial_beyond_he <- spatial_raw - mean(spatial_raw)
  total_he_se <- spatial_beyond_he_se <- numeric(n)
  for (index in blocks) {
    newdata <- d[index, , drop = FALSE]
    basis <- mgcv::PredictMat(he_smooth, newdata)
    total_he_se[index] <- .component_se(basis, he_cov)
    basis <- sweep(mgcv::PredictMat(spatial_smooth, newdata), 2, mean_spatial_basis, "-")
    spatial_beyond_he_se[index] <- .component_se(basis, spatial_cov)
  }
  total_he_p <- .pointwise_p(total_he, total_he_se)
  spatial_beyond_he_p <- .pointwise_p(spatial_beyond_he, spatial_beyond_he_se)
  fitted_values <- as.numeric(stats::fitted(m4))
  baseline <- unname(stats::coef(m4)["(Intercept)"]) + mean(spatial_raw)
  reconstruction_error <- max(abs(fitted_values - baseline - total_he - spatial_beyond_he))
  if (reconstruction_error > 1e-7 * max(1, max(abs(fitted_values)))) {
    stop("Component maps do not reconstruct M4 fitted expression.", call. = FALSE)
  }
  spatial <- data.frame(
    spot_id = spots$spot_id, x = spots$x, y = spots$y, expression = y,
    fitted = fitted_values, he = h,
    slide_wide_he = slide_wide_he, location_specific_he = location_specific_he,
    total_he = total_he, spatial_beyond_he = spatial_beyond_he,
    total_he_se = total_he_se, total_he_p = total_he_p,
    total_he_q = stats::p.adjust(total_he_p, method = "BH"),
    spatial_beyond_he_se = spatial_beyond_he_se, spatial_beyond_he_p = spatial_beyond_he_p,
    spatial_beyond_he_q = stats::p.adjust(spatial_beyond_he_p, method = "BH"),
    stringsAsFactors = FALSE
  )
  reasons <- vapply(comparisons, `[[`, character(1), "reason")
  reasons <- if (any(nzchar(reasons))) {
    paste(paste0(names(reasons)[nzchar(reasons)], ": ", reasons[nzchar(reasons)]), collapse = "; ")
  } else ""
  status <- data.frame(gene = gene, status = if (nzchar(reasons)) "partial" else "complete",
                       message = reasons, warnings = "", stringsAsFactors = FALSE)
  list(gene = gene, statistics = statistics, spatial = spatial,
       model_diagnostics = diagnostics, status = status,
       diagnostics = list(n = n, baseline = baseline,
                          slide_wide_coefficient = slide_wide_coefficient,
                          mean_location_specific_coefficient = mean(location_specific_coefficient),
                          component_reconstruction_error = reconstruction_error,
                          comparisons = comparisons, pointwise_covariance = "mgcv conditional Vp"))
}

.atomic_save_rds <- function(object, file) {
  temporary <- tempfile(pattern = ".writing-", tmpdir = dirname(file))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = TRUE)
  if (!file.rename(temporary, file)) stop("Could not write ", file, ".", call. = FALSE)
  invisible(file)
}

.fit_gene_job <- function(index, state) {
  gene <- state$genes[index]
  file <- if (is.null(state$output_dir)) NULL else file.path(state$output_dir, "genes", sprintf("gene_%06d.rds", index))
  if (state$resume && !is.null(file) && file.exists(file)) {
    previous <- tryCatch(readRDS(file), error = function(e) NULL)
    if (!is.null(previous) && identical(previous$run_fingerprint, state$fingerprint) &&
        identical(previous$gene, gene) && isTRUE(previous$status$status %in% c("complete", "partial"))) {
      previous$spatial <- NULL
      previous$file <- file
      return(previous)
    }
  }
  warnings <- character()
  result <- tryCatch(withCallingHandlers(
    .fit_one_gene(state$expression, state$spots, state$h, gene, state$k),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      list(gene = gene, statistics = .empty_statistics(gene), spatial = NULL,
           model_diagnostics = NULL, diagnostics = NULL,
           status = data.frame(gene = gene, status = "failed", message = conditionMessage(e), warnings = "", stringsAsFactors = FALSE))
    })
  result$status$warnings <- paste(unique(warnings), collapse = "; ")
  result$run_fingerprint <- state$fingerprint
  if (!is.null(file)) {
    .atomic_save_rds(result, file)
    result$spatial <- NULL
    result$file <- file
  }
  result
}

#' Fit MorphoANOVA to one or more genes
#'
#' Fits four nested Gaussian spatial models with REML and summarizes sequential
#' F tests, adjusted R-squared increments, and spatial component maps. Expression
#' is used as supplied; normalize and log-transform it before calling this
#' function. The selected H&E feature is mean-centered across locations.
#'
#' @param expression Numeric genes-by-locations matrix, including a sparse
#'   \code{Matrix}, with unique row and column names.
#' @param spots A data frame with \code{spot_id}, \code{x}, \code{y}, and the
#'   selected H&E feature. IDs must match expression column names; row order may
#'   differ. Each call analyzes one tissue section.
#' @param genes Gene names to fit. \code{NULL} uses every row of expression.
#' @param feature Name of the numeric morphology column in \code{spots}.
#' @param k Thin-plate spline basis dimension, used for both the H&E coefficient
#'   surface and the additional spatial smooth. The default is 30. This sets
#'   their maximum flexibility; REML selects the smoothing penalties and hence
#'   their effective degrees of freedom. Larger values allow finer spatial
#'   patterns and increase computation.
#' @param cores Number of parallel workers. The default is one.
#' @param output_dir Optional directory for per-gene results. With \code{NULL},
#'   spatial results stay in memory. With a directory, \code{spatial} contains
#'   named per-gene file paths; plotting functions read only the requested genes.
#' @param resume Resume a matching run in \code{output_dir}. Changed inputs,
#'   source code, model settings, or software versions require a new directory.
#' @return A \code{morphoanova_fit} list containing four entries:
#'   \code{statistics}, \code{spatial}, \code{status}, and \code{settings}.
#'   \code{spatial[[gene]]} is an ordinary data frame with the already-computed
#'   location-level results for that gene. Use \code{head()} to preview it.
#'   In disk mode, \code{spatial[[gene]]} is the saved file path. Read its
#'   location-level table with \code{readRDS(fit$spatial[[gene]])$spatial}.
#'   \code{settings} records the main model settings, centering,
#'   package version, and numbers of genes and locations. File paths, model
#'   diagnostics, and full provenance are retained in an internal attribute.
#'   \code{F21}, \code{F32}, and \code{F43} compare M2 versus M1, M3 versus M2,
#'   and M4 versus M3, respectively. \code{p21}, \code{p32}, and \code{p43}
#'   are unadjusted p-values. \code{q21}, \code{q32}, and \code{q43} are
#'   Benjamini-Hochberg (BH)-adjusted p-values across the genes analyzed in the
#'   current fit, separately for each comparison. Component contributions are
#'   raw adjusted R-squared increments:
#'   \code{r2_slide_wide_he} (M2 minus M1), \code{r2_location_specific_he}
#'   (M3 minus M2), and \code{r2_spatial_beyond_he} (M4 minus M3).
#'   \code{r2_full} is M4 adjusted R-squared; \code{r2_unexplained_variance} is
#'   \code{1 - r2_full}. Negative increments are retained.
#'   Use adjusted p-values to assess statistical evidence for each component,
#'   and incremental adjusted R-squared values to summarize its contribution
#'   to model fit.
#'   \code{statistics} is the summary table containing the completed across-gene
#'   BH adjustments. Failed genes retain missing statistics, a \code{NULL}
#'   spatial entry in memory mode, and an explanatory status.
#' @details Sequential F statistics and p-values are taken from pairwise
#' \code{mgcv::anova.gam(reduced, full, test = "F")} comparisons. The test
#' uses corrected effective degrees of freedom and the estimated dispersion
#' of the model with the smaller corrected residual degrees of freedom.
#' The fitted degrees of freedom, \code{sum(model$edf)}, remain the definition
#' used for adjusted R-squared and the model diagnostics.
#'
#' Unavailable comparisons retain \code{NA}, with \code{partial} gene status
#' and an explanatory message. Their fitted spatial components remain
#' available. Internal comparison diagnostics retain the signed test degrees
#' of freedom, denominator degrees of freedom, corrected residual degrees of
#' freedom, dispersion, deviance difference, and algebraic F statistic.
#' The algebraic statistic is diagnostic only when ANOVA returns \code{NA}.
#' 
#' In spatial tables, \code{fitted} contains the full model's fitted expression.
#' \code{total_he = slide_wide_he + location_specific_he} combines the two
#' fitted H&E contributions. \code{spatial_beyond_he} is M4's additional
#' spatial smooth, centered across locations.
#'
#' Pointwise tests assess whether a fitted component contribution is zero at
#' each location: total H&E effect or spatial variation beyond H&E. They are
#' two-sided normal Wald tests, using the estimate divided by its standard
#' error from M4's conditional coefficient covariance (\code{Vp}). BH
#' correction is across locations, separately within each gene and component.
#' Structurally zero effects with zero standard error have missing p-values.
#' @examples
#' \dontrun{
#' # Download ovarian_full.rds and tissue_hires_image.png from GitHub tutorial/data.
#' ovarian <- readRDS("ovarian_full.rds")
#' he <- extract_he_feature("tissue_hires_image.png", ovarian$coordinates,
#'                          channel = "red", patch_size = c(15, 11), flip_y = TRUE)
#' spots <- he$spots
#' fit <- fit_morphoanova(ovarian$expression, spots, genes = "IGHM", k = 30)
#' fit$statistics
#' head(fit$spatial[["IGHM"]])
#' }
#' @export
fit_morphoanova <- function(expression, spots, genes = NULL, feature = "he",
                            k = 30, cores = 1, output_dir = NULL, resume = FALSE) {
  input <- .prepare_fit_input(expression, spots, genes, feature, k, cores)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) {
    stop("`resume` must be TRUE or FALSE.", call. = FALSE)
  }
  if (resume && is.null(output_dir)) stop("`resume` requires `output_dir`.", call. = FALSE)
  helper_names <- c(".validate_spots", ".prepare_fit_input", ".empty_statistics",
                    ".sequential_test", ".component_se", ".pointwise_p",
                    ".fit_one_gene", ".atomic_save_rds", ".fit_gene_job")
  source_definition <- lapply(c(helper_names, "fit_morphoanova"), function(name) {
    fun <- get(name, envir = environment(fit_morphoanova))
    list(formals = deparse(formals(fun)), body = deparse(body(fun)))
  })
  settings <- list(feature = feature, feature_center = input$center, k = as.integer(k),
                   family = "gaussian", spline = "tp", method = "REML", edf = "sum(model$edf)",
                   adjusted_r2_edf = "sum(model$edf)",
                   sequential_test = "mgcv::anova.gam(reduced, full, test = 'F')",
                   test_definition_version = "mgcv_anova_gam_F_v1",
                   centering = "H&E mean across all locations in this tissue section",
                   models = c(M1 = "y ~ 1", M2 = "y ~ h", M3 = "y ~ s(x, y_coord, by = h)",
                              M4 = "y ~ s(x, y_coord, by = h) + s(x, y_coord)"),
                   version = tryCatch(as.character(utils::packageVersion("MorphoANOVA")), error = function(e) "development"),
                   R_version = as.character(getRversion()), mgcv_version = as.character(utils::packageVersion("mgcv")),
                   source_hash = digest::digest(source_definition, algo = "sha256"))
  fingerprint <- digest::digest(list(expression = input$expression,
                                    locations = input$spots[c("spot_id", "x", "y", feature)],
                                    settings = settings), algo = "sha256")
  settings$run_fingerprint <- fingerprint
  settings$genes <- input$genes
  settings$n_locations <- nrow(input$spots)
  if (!is.null(output_dir)) {
    if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
      stop("`output_dir` must be one directory path.", call. = FALSE)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    output_dir <- normalizePath(output_dir, mustWork = TRUE)
    manifest_file <- file.path(output_dir, "run.rds")
    if (file.exists(manifest_file)) {
      previous <- readRDS(manifest_file)
      if (!resume) stop("Output directory contains a run. Use `resume = TRUE` or a new directory.", call. = FALSE)
      if (!identical(previous$run_fingerprint, fingerprint)) {
        stop("Cannot resume: input data, model settings, source, or software versions changed.", call. = FALSE)
      }
    } else {
      if (length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
        stop("Use an empty output directory for a new run.", call. = FALSE)
      }
      .atomic_save_rds(settings, manifest_file)
    }
    dir.create(file.path(output_dir, "genes"), showWarnings = FALSE)
  }
  state <- list(expression = input$expression, spots = input$spots[c("spot_id", "x", "y")],
                h = input$h, genes = input$genes, k = as.integer(k),
                output_dir = output_dir, resume = resume, fingerprint = fingerprint)
  cores <- min(as.integer(cores), length(input$genes))
  if (cores == 1L) {
    results <- lapply(seq_along(input$genes), .fit_gene_job, state = state)
  } else {
    cluster <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
    if (inherits(input$expression, "Matrix")) {
      parallel::clusterCall(cluster, function() { loadNamespace("Matrix"); NULL })
    }
    parallel::clusterExport(cluster, helper_names, envir = environment(fit_morphoanova))
    parallel::clusterCall(cluster, function(value) assign(".morphoanova_state", value, envir = .GlobalEnv), state)
    worker <- function(index) {
      .fit_gene_job(index, .morphoanova_state)
    }
    environment(worker) <- .GlobalEnv
    results <- parallel::parLapply(cluster, seq_along(input$genes), worker)
  }
  names(results) <- input$genes
  statistics <- do.call(rbind, lapply(results, `[[`, "statistics"))
  rownames(statistics) <- NULL
  for (test in c("21", "32", "43")) {
    statistics[[paste0("q", test)]] <- stats::p.adjust(statistics[[paste0("p", test)]], method = "BH")
  }
  status <- do.call(rbind, lapply(results, `[[`, "status"))
  rownames(status) <- NULL
  files <- if (is.null(output_dir)) NULL else vapply(results, `[[`, character(1), "file")
  spatial <- if (is.null(output_dir)) lapply(results, `[[`, "spatial") else files
  diagnostics <- lapply(results, function(gene_result) {
    gene_result[c("model_diagnostics", "diagnostics")]
  })
  main_settings <- settings[c("feature", "feature_center", "k", "family", "spline",
                              "method", "version", "n_locations")]
  main_settings$n_genes <- length(input$genes)
  result <- structure(list(statistics = statistics, spatial = spatial,
                           status = status, settings = main_settings),
                      class = "morphoanova_fit")
  attr(result, "internal") <- list(files = files, diagnostics = diagnostics,
                                    provenance = settings)
  if (!is.null(output_dir)) .atomic_save_rds(result, file.path(output_dir, "fit.rds"))
  if (any(status$status != "complete")) {
    warning(sum(status$status == "failed"), " gene(s) failed; ",
            sum(status$status == "partial"), " gene(s) have unavailable tests. See `fit$status`.", call. = FALSE)
  }
  result
}

# Read one gene for component maps and composition plots.
.gene_spatial <- function(object, gene = NULL) {
  genes <- object$statistics$gene
  if (is.null(gene)) {
    if (length(genes) != 1L) stop("Choose one `gene` from this multi-gene fit.", call. = FALSE)
    gene <- genes[1]
  }
  if (!is.character(gene) || length(gene) != 1L || is.na(gene) || !gene %in% genes) {
    stop("`gene` must name one fitted gene.", call. = FALSE)
  }
  status <- object$status[match(gene, object$status$gene), , drop = FALSE]
  if (identical(status$status, "failed")) stop("Gene fit failed: ", status$message, call. = FALSE)
  if (is.list(object$spatial)) return(object$spatial[[gene]])

  # Retain access to results saved before spatial tables were moved to the top level.
  if (is.character(object$spatial)) {
    result <- readRDS(object$spatial[[gene]])
  } else if (!is.null(object$per_gene)) {
    result <- object$per_gene[[gene]]
  } else {
    files <- attr(object, "internal", exact = TRUE)$files
    if (is.null(files)) files <- object$files
    result <- readRDS(files[[gene]])
  }
  if (!identical(result$gene, gene)) stop("Stored result does not match the requested gene.", call. = FALSE)
  if (identical(result$status$status, "failed")) stop("Gene fit failed: ", result$status$message, call. = FALSE)
  result$spatial
}

#' @export
print.morphoanova_fit <- function(x, ...) {
  cat("MorphoANOVA fit:", nrow(x$statistics), "gene(s),", x$settings$n_locations, "locations\n")
  cat(sum(x$status$status == "complete"), "complete,",
      sum(x$status$status == "partial"), "with unavailable tests,",
      sum(x$status$status == "failed"), "failed\n")
  invisible(x)
}

utils::globalVariables(".morphoanova_state")

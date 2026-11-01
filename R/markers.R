#' Compare expression between annotated tissue regions
#'
#' Tests each gene using a two-sided Wilcoxon rank-sum test with a normal
#' approximation, tie correction, and no continuity correction. BH adjustment
#' includes all genes supplied to this function. Positive z scores indicate
#' higher expression in `group` than in `reference`.
#'
#' @param expression Numeric genes-by-locations matrix with row and column names.
#' @param groups Data frame with `spot_id` and `group`, matched by ID.
#' @param group Tissue label to compare against the reference labels.
#' @param reference One or more reference tissue labels.
#' @param log_transformed Whether expression contains natural-log(1+x) values.
#'   When TRUE, undo this transformation to calculate mean expression ratios.
#' @return A data frame sorted by decreasing z score, with gene, U statistic,
#'   tie-adjusted z score, unadjusted p-value, BH-adjusted p-value, log2 ratio
#'   of group means with a pseudocount of one, detection fraction, and positive
#'   rank. Positive rank includes genes with both z and log2 ratio above zero;
#'   it does not apply a significance cutoff. Subset `q < 0.05` for evidence
#'   of differential expression.
#' @export
find_markers <- function(expression, groups, group, reference, log_transformed = TRUE) {
  if (!(is.matrix(expression) || inherits(expression, "Matrix")) ||
      is.null(rownames(expression)) || is.null(colnames(expression)) ||
      anyNA(rownames(expression)) || anyNA(colnames(expression)) ||
      anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression)) ||
      !nrow(expression) || !ncol(expression)) {
    stop("expression must have unique gene and spot names.", call. = FALSE)
  }
  if (!is.data.frame(groups) || !all(c("spot_id", "group") %in% names(groups)) ||
      anyNA(groups$spot_id) || anyDuplicated(groups$spot_id) ||
      anyNA(match(colnames(expression), groups$spot_id))) {
    stop("groups must contain one spot_id and group for every expression column.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1L || is.na(group) ||
      !is.character(reference) || !length(reference) || anyNA(reference) ||
      group %in% reference || !is.logical(log_transformed) ||
      length(log_transformed) != 1L || is.na(log_transformed)) {
    stop("Choose one group and distinct reference labels; log_transformed must be TRUE or FALSE.", call. = FALSE)
  }
  labels <- as.character(groups$group[match(colnames(expression), groups$spot_id)])
  in_group <- which(labels == group)
  in_reference <- which(labels %in% reference)
  n1 <- length(in_group); n0 <- length(in_reference)
  if (n1 < 2L || n0 < 2L) stop("Both comparison groups need at least two locations.", call. = FALSE)
  rows <- vector("list", nrow(expression))
  for (i in seq_len(nrow(expression))) {
    values <- as.numeric(expression[i, ])
    if (any(!is.finite(values)) || any(values < 0)) {
      stop("expression must contain finite non-negative measurements.", call. = FALSE)
    }
    x <- values[in_group]; y <- values[in_reference]
    test <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE, correct = FALSE))
    u <- unname(test$statistic)
    ties <- as.numeric(table(c(x, y)))
    variance_u <- n1 * n0 / 12 *
      (n1 + n0 + 1 - sum(ties^3 - ties) / ((n1 + n0) * (n1 + n0 - 1)))
    z <- if (variance_u > 0) (u - n1 * n0 / 2) / sqrt(variance_u) else NA_real_
    original_x <- if (log_transformed) expm1(x) else x
    original_y <- if (log_transformed) expm1(y) else y
    rows[[i]] <- data.frame(gene = rownames(expression)[i], U = u, z = z,
      p = test$p.value, log2_fold_change = log2((mean(original_x) + 1) /
      (mean(original_y) + 1)), detection_fraction = mean(values > 0))
  }
  result <- do.call(rbind, rows)
  result$q <- stats::p.adjust(result$p, "BH")
  result <- result[order(-result$z, result$gene, na.last = TRUE), ]
  positive <- is.finite(result$z) & result$z > 0 & result$log2_fold_change > 0
  result$positive_rank <- match(result$gene, result$gene[positive])
  rownames(result) <- NULL
  result
}

.validate_spots <- function(spots, required = c("spot_id", "x", "y")) {
  if (!is.data.frame(spots)) stop("`spots` must be a data.frame.", call. = FALSE)
  if (anyDuplicated(names(spots))) stop("`spots` has duplicated column names.", call. = FALSE)
  missing <- setdiff(required, names(spots))
  if (length(missing)) {
    stop("Missing spot columns: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  ids <- as.character(spots$spot_id)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("`spots$spot_id` must contain unique, non-empty IDs.", call. = FALSE)
  }
  spots$spot_id <- ids
  for (column in setdiff(required, "spot_id")) {
    if (!is.numeric(spots[[column]]) || any(!is.finite(spots[[column]]))) {
      stop("`spots$", column, "` must contain finite numeric values.", call. = FALSE)
    }
  }
  spots
}

.prepare_fit_input <- function(expression, spots, genes, feature, k, cores) {
  if (!(is.matrix(expression) || inherits(expression, "Matrix"))) {
    stop("`expression` must be a genes-by-locations matrix.", call. = FALSE)
  }
  if (inherits(expression, "Matrix") && !"x" %in% methods::slotNames(expression)) {
    stop("`expression` must contain finite numeric values.", call. = FALSE)
  }
  values <- if (inherits(expression, "Matrix")) expression@x else expression
  if (!is.numeric(values) || is.complex(values) || any(!is.finite(values))) {
    stop("`expression` must contain finite numeric values.", call. = FALSE)
  }
  for (margin in c("genes", "locations")) {
    ids <- if (margin == "genes") rownames(expression) else colnames(expression)
    if (is.null(ids) || !length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("Expression ", margin, " must have unique, non-empty names.", call. = FALSE)
    }
  }
  if (!is.character(feature) || length(feature) != 1L || is.na(feature) || !nzchar(feature)) {
    stop("`feature` must name one numeric spot column.", call. = FALSE)
  }
  if (feature %in% c("spot_id", "x", "y")) {
    stop("`feature` must name a morphology column, separate from IDs and coordinates.", call. = FALSE)
  }
  spots <- .validate_spots(spots, c("spot_id", "x", "y", feature))
  missing <- setdiff(colnames(expression), spots$spot_id)
  extra <- setdiff(spots$spot_id, colnames(expression))
  if (length(missing) || length(extra)) {
    stop("Expression and spots must contain the same location IDs (",
         length(missing), " missing from spots; ", length(extra),
         " extra in spots). Subset both inputs explicitly.", call. = FALSE)
  }
  spots <- spots[match(colnames(expression), spots$spot_id), , drop = FALSE]
  rownames(spots) <- spots$spot_id
  if (is.null(genes)) genes <- rownames(expression)
  if (!is.character(genes) || !length(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("`genes` must contain unique gene names.", call. = FALSE)
  }
  missing <- setdiff(genes, rownames(expression))
  if (length(missing)) stop("Genes absent from expression: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  if (length(k) != 1L || !is.numeric(k) || !is.finite(k) || k < 4 || k != as.integer(k)) {
    stop("`k` must be an integer of at least 4.", call. = FALSE)
  }
  if (length(cores) != 1L || !is.numeric(cores) || !is.finite(cores) || cores < 1 || cores != as.integer(cores)) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  if (nrow(unique(spots[c("x", "y")])) < k || nrow(spots) < 2 * k) {
    stop("Too few locations for `k`; use a smaller basis dimension.", call. = FALSE)
  }
  h <- as.numeric(spots[[feature]])
  center <- mean(h)
  h <- h - center
  if (!any(h != 0)) stop("The selected H&E feature has no variation.", call. = FALSE)
  list(expression = expression[genes, , drop = FALSE], spots = spots,
       genes = genes, h = h, center = center)
}

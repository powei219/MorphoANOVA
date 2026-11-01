#' Extract one H&E feature at each spatial location
#'
#' Calculates mean image intensity in a rectangular patch around each location,
#' then subtracts the mean across the supplied locations.
#' The returned spot table can be passed directly to [fit_morphoanova()].
#' Coordinates used for extraction are one-based image pixels: x is the column
#' and y is the row. Model coordinates in `spots$x` and `spots$y` are preserved.
#' Use pixel coordinates aligned to the image's scale and orientation after
#' any rotation selected with `rotate`.
#' If model coordinates are in micrometers, retain separate image-coordinate
#' columns and select them with `pixel_columns`.
#'
#' @param image PNG/JPEG filename or numeric height-by-width-by-channel RGB array.
#'   Intensities may be in `[0, 1]` or `[0, 255]`; patch means use the `[0, 1]`
#'   scale before centering.
#' @param spots Data frame containing unique `spot_id`, model coordinates `x`,
#'   `y`, and the image coordinate columns selected by `pixel_columns`.
#' @param channel Image channel: `"rgb"` averages all three channels.
#' @param patch_size Patch width and height in pixels, or one square-patch size.
#' @param pixel_columns Names of the image x and y columns in `spots`.
#' @param flip_y Reflect rounded image y coordinates using `height - y + 1`,
#'   with the height measured after rotation.
#' @param rotate Rotate the image counterclockwise by 0, 90, 180 or 270 degrees
#'   before matching pixel coordinates. The supplied pixel coordinates refer
#'   to the rotated image; model coordinates are unchanged.
#' @return A list with `spots` (the input table plus the mean-centered `he`), `preview`, and
#'   `settings`. With ggplot2 and patchwork installed, the preview shows the
#'   aligned H&E image on the left and the mean-centered H&E feature on the
#'   right, using the paper's green-cream-pink scale. The preview uses the
#'   returned `spots$he` values directly. `settings$feature_center` records the
#'   subtracted mean, and `settings$intensity_scale` describes the `[0, 1]`
#'   intensity scale before centering. `settings$rotate`, `settings$width` and
#'   `settings$height` record the rotation and resulting image dimensions.
#'   Extraction uses the original image resolution. The histology preview is
#'   limited to 2,000 pixels on its longest edge, keeping the original field
#'   of view and coordinate scale.
#'   Patches at image edges use their available pixels; their IDs are recorded
#'   in settings.
#' @export
#' @examples
#' image <- array(seq(0, 1, length.out = 20 * 20 * 3), c(20, 20, 3))
#' spots <- data.frame(spot_id = c("a", "b"), x = c(5, 15), y = c(5, 15))
#' he <- extract_he_feature(image, spots, channel = "green", patch_size = 3)
#' he$spots
extract_he_feature <- function(image, spots, channel = c("rgb", "red", "green", "blue"),
                               patch_size, pixel_columns = c("x", "y"),
                               flip_y = FALSE, rotate = 0) {
  channel <- match.arg(channel)
  if (!is.data.frame(spots) || !nrow(spots) ||
      !all(c("spot_id", "x", "y") %in% names(spots))) {
    stop("spots must contain spot_id, x and y.", call. = FALSE)
  }
  ids <- as.character(spots$spot_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("spot_id must contain unique, non-missing IDs.", call. = FALSE)
  }
  if (length(pixel_columns) != 2L || anyNA(pixel_columns) ||
      anyDuplicated(pixel_columns) || !all(pixel_columns %in% names(spots))) {
    stop("pixel_columns must name the image x and y columns.", call. = FALSE)
  }
  for (column in unique(c("x", "y", pixel_columns))) {
    if (!is.numeric(spots[[column]]) || any(!is.finite(spots[[column]]))) {
      stop("Coordinate column must be finite numeric: ", column, call. = FALSE)
    }
  }
  if (!is.numeric(patch_size) || !length(patch_size) || length(patch_size) > 2L ||
      any(!is.finite(patch_size)) || any(patch_size < 1 | patch_size != round(patch_size))) {
    stop("patch_size must be one or two positive whole numbers.", call. = FALSE)
  }
  patch_size <- rep(as.integer(patch_size), length.out = 2L)
  if (!is.logical(flip_y) || length(flip_y) != 1L || is.na(flip_y)) {
    stop("flip_y must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(rotate) || length(rotate) != 1L || is.na(rotate) ||
      !rotate %in% c(0, 90, 180, 270)) {
    stop("rotate must be 0, 90, 180 or 270 degrees counterclockwise.", call. = FALSE)
  }

  source <- if (is.character(image) && length(image) == 1L) image else "array"
  native <- FALSE
  if (is.character(image) && length(image) == 1L) {
    extension <- tolower(tools::file_ext(image))
    reader <- if (extension == "png") "png" else if (extension %in% c("jpg", "jpeg")) "jpeg" else ""
    if (!nzchar(reader)) stop("Use a PNG/JPEG filename or an RGB array.", call. = FALSE)
    if (!requireNamespace(reader, quietly = TRUE)) stop("Install package '", reader, "' to read this image.", call. = FALSE)
    native <- reader == "jpeg"
    image <- if (native) jpeg::readJPEG(image, native = TRUE) else png::readPNG(image)
  }
  intensity_scale <- 1
  if (!native) {
    if (!is.array(image) || !is.numeric(image) || length(dim(image)) != 3L ||
        any(dim(image) == 0L) || dim(image)[3L] < 3L) {
      stop("image must be a finite RGB array with intensities in [0, 1] or [0, 255].", call. = FALSE)
    }
    # Validate numeric images in blocks without copying the whole RGB array.
    pixels_per_channel <- as.double(dim(image)[1L]) * dim(image)[2L]
    rgb_max <- 0
    for (j in seq_len(dim(image)[3L])) {
      offset <- (j - 1) * pixels_per_channel
      for (start in seq(1, pixels_per_channel, by = 1e6)) {
        values <- image[offset + seq.int(start, min(start + 1e6 - 1, pixels_per_channel))]
        if (any(!is.finite(values)) || any(values < 0 | values > 255)) {
          stop("image must be a finite RGB array with intensities in [0, 1] or [0, 255].", call. = FALSE)
        }
        if (j <= 3L) rgb_max <- max(rgb_max, max(values))
      }
    }
    if (rgb_max > 1) intensity_scale <- 255
    rm(values)
  }
  if (rotate != 0) {
    source_height <- dim(image)[1L]
    source_width <- dim(image)[2L]
    rotated_height <- if (rotate == 180) source_height else source_width
    rotated_width <- if (rotate == 180) source_width else source_height
    if (native) {
      # Rotate packed row-major pixels in blocks without expanding them to RGB.
      rotated <- integer(length(image))
      for (start in seq(1, length(image), by = 1e6)) {
        index <- seq.int(start, min(start + 1e6 - 1, length(image)))
        row <- (index - 1) %/% rotated_width + 1
        column <- (index - 1) %% rotated_width + 1
        if (rotate == 90) {
          source_row <- column
          source_column <- source_width - row + 1
        } else if (rotate == 180) {
          source_row <- source_height - row + 1
          source_column <- source_width - column + 1
        } else {
          source_row <- source_height - column + 1
          source_column <- row
        }
        rotated[index] <- image[(source_row - 1) * as.double(source_width) + source_column]
      }
      rotated <- structure(rotated, dim = c(rotated_height, rotated_width),
                           class = "nativeRaster")
      rm(index, row, column, source_row, source_column)
    } else {
      rotated <- array(0, c(rotated_height, rotated_width, dim(image)[3L]))
      for (j in seq_len(dim(image)[3L])) {
        plane <- matrix(image[, , j], nrow = source_height, ncol = source_width)
        if (rotate == 90) {
          rotated[, , j] <- t(plane[, source_width:1, drop = FALSE])
        } else if (rotate == 180) {
          rotated[, , j] <- plane[source_height:1, source_width:1, drop = FALSE]
        } else {
          rotated[, , j] <- t(plane[source_height:1, , drop = FALSE])
        }
      }
      rm(plane)
    }
    image <- rotated
    rm(rotated)
  }
  height <- dim(image)[1L]
  width <- dim(image)[2L]
  x <- round(spots[[pixel_columns[1L]]])
  y <- round(spots[[pixel_columns[2L]]])
  if (flip_y) y <- height - y + 1
  outside <- x < 1 | x > width | y < 1 | y > height
  if (any(outside)) {
    stop("Spot centers outside the image: ", paste(utils::head(ids[outside], 5L), collapse = ", "), call. = FALSE)
  }
  left <- floor((patch_size - 1L) / 2L)
  right <- patch_size - 1L - left
  x1 <- pmax(1, x - left[1L]); x2 <- pmin(width, x + right[1L])
  y1 <- pmax(1, y - left[2L]); y2 <- pmin(height, y + right[2L])
  area <- (x2 - x1 + 1) * (y2 - y1 + 1)
  clipped <- area < prod(patch_size)
  if (any(clipped)) warning(sum(clipped), " edge patches use the available image pixels.", call. = FALSE)

  channels <- if (channel == "rgb") 1:3 else match(channel, c("red", "green", "blue"))
  feature <- numeric(nrow(spots))
  for (i in seq_len(nrow(spots))) {
    rows <- seq.int(y1[i], y2[i])
    columns <- seq.int(x1[i], x2[i])
    if (native) {
      # nativeRaster stores packed RGBA pixels in row-major order.
      index <- as.vector(outer(columns, (rows - 1) * as.double(width), "+"))
      pixels <- image[index]
      for (j in channels) {
        feature[i] <- feature[i] +
          mean(bitwAnd(bitwShiftR(pixels, 8L * (j - 1L)), 255L)) / 255 / length(channels)
      }
    } else {
      feature[i] <- mean(image[rows, columns, channels, drop = FALSE]) / intensity_scale
    }
  }
  feature_center <- mean(feature)
  feature <- feature - feature_center
  spots$spot_id <- ids
  spots$he <- feature
  preview <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("patchwork", quietly = TRUE)) {
    preview_scale <- min(1, 2000 / max(width, height))
    preview_width <- max(1L, round(width * preview_scale))
    preview_height <- max(1L, round(height * preview_scale))
    columns <- round(seq(1, width, length.out = preview_width))
    rows <- round(seq(1, height, length.out = preview_height))
    if (native) {
      index <- as.vector(outer(columns, (rows - 1) * as.double(width), "+"))
      thumbnail <- structure(image[index], dim = c(preview_height, preview_width),
                             class = "nativeRaster")
      rm(index, pixels)
    } else {
      thumbnail <- grDevices::as.raster(image[rows, columns, 1:3, drop = FALSE] / intensity_scale)
    }
    # Plot mappings retain this environment; release the full image first.
    rm(image)
    mapped <- data.frame(pixel_x = x, pixel_y = height - y + 1,
                         value = feature)
    bound <- max(abs(mapped$value))
    frame <- ggplot2::coord_fixed(xlim = c(0.5, width + 0.5),
      ylim = c(0.5, height + 0.5), expand = FALSE)
    style <- ggplot2::theme_void(base_family = "Arial") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0.5),
        plot.margin = ggplot2::margin(5, 5, 5, 5))
    histology <- ggplot2::ggplot() +
      ggplot2::annotation_raster(thumbnail, xmin = 0.5,
        xmax = width + 0.5, ymin = 0.5, ymax = height + 0.5, interpolate = FALSE) +
      frame + style + ggplot2::labs(title = "Histology")
    feature_map <- ggplot2::ggplot(mapped, ggplot2::aes(pixel_x, pixel_y)) +
      ggplot2::geom_point(ggplot2::aes(colour = value), size = 1.2, stroke = 0) +
      ggplot2::scale_colour_gradient2(low = "#3AB370", mid = "#EAE7CC", high = "#FD1593",
        midpoint = 0, limits = c(-bound, bound), name = "Value",
        breaks = if (bound == 0) 0 else ggplot2::waiver()) +
      frame + style + ggplot2::labs(title = "Extracted feature")
    preview <- patchwork::wrap_plots(histology, feature_map, nrow = 1,
      widths = c(1, 1), guides = "collect")
  }
  list(spots = spots, preview = preview, settings = list(
    image = source, width = width, height = height, channel = channel,
    patch_size = stats::setNames(patch_size, c("width", "height")),
    pixel_columns = pixel_columns, flip_y = flip_y, rotate = rotate,
    coordinate_origin = "one_based",
    clipped_spot_ids = ids[clipped], intensity_scale = c(0, 1),
    feature_center = feature_center))
}

utils::globalVariables(c("pixel_x", "pixel_y", "value"))

# Pixel IDs for a three-row, five-column image, listed in display order.
rotation_pixels <- list(
  "0" = rbind(1:5, 6:10, 11:15),
  "90" = rbind(c(5, 10, 15), c(4, 9, 14), c(3, 8, 13), c(2, 7, 12), c(1, 6, 11)),
  "180" = rbind(15:11, 10:6, 5:1),
  "270" = rbind(c(11, 6, 1), c(12, 7, 2), c(13, 8, 3), c(14, 9, 4), c(15, 10, 5))
)

rotation_fixture <- function() {
  ids <- rotation_pixels[["0"]]
  image <- array(0, c(3, 5, 3))
  image[, , 1] <- ids / 16
  image[, , 2] <- (ids * 7) %% 17 / 16
  image[, , 3] <- (ids * 11) %% 19 / 18
  image
}

explicit_rotation <- function(image, angle) {
  ids <- rotation_pixels[[as.character(angle)]]
  rotated <- array(0, c(dim(ids), 3))
  for (j in 1:3) {
    pixels <- as.vector(t(image[, , j]))
    rotated[, , j] <- matrix(pixels[ids], nrow = nrow(ids), ncol = ncol(ids))
  }
  rotated
}

test_that("array rotations use coordinates in the rotated image", {
  image <- rotation_fixture()
  for (angle in c(0, 90, 180, 270)) {
    expected <- explicit_rotation(image, angle)
    height <- dim(expected)[1]
    width <- dim(expected)[2]
    coordinates <- expand.grid(px = seq_len(width), py = seq_len(height))
    spots <- data.frame(spot_id = paste0("s", seq_len(nrow(coordinates))),
      x = coordinates$px * 100, y = coordinates$py * 200, coordinates)
    for (flip in c(FALSE, TRUE)) {
      rows <- if (flip) height - spots$py + 1 else spots$py
      for (channel in c("red", "green", "blue", "rgb")) {
        selected <- if (channel == "rgb") 1:3 else match(channel, c("red", "green", "blue"))
        raw <- vapply(seq_len(nrow(spots)), function(i) {
          mean(expected[rows[i], spots$px[i], selected])
        }, numeric(1))
        he <- extract_he_feature(image, spots, channel = channel, patch_size = 1,
          pixel_columns = c("px", "py"), flip_y = flip, rotate = angle)
        expect_equal(he$spots$he, raw - mean(raw), tolerance = 1e-12)
        expect_equal(he$settings$feature_center, mean(raw), tolerance = 1e-12)
        expect_identical(he$spots[names(spots)], spots)
        expect_equal(c(he$settings$height, he$settings$width), c(height, width))
        expect_equal(he$settings$rotate, angle)
      }
    }
    interior <- spots[spots$px > 1 & spots$px < width, ]
    actual <- extract_he_feature(image * 255, interior, channel = "red",
      patch_size = c(3, 1), pixel_columns = c("px", "py"), flip_y = TRUE, rotate = angle)
    reference <- extract_he_feature(expected, interior, channel = "red",
      patch_size = c(3, 1), pixel_columns = c("px", "py"), flip_y = TRUE)
    expect_equal(actual$spots, reference$spots, tolerance = 1e-12)
    expect_equal(actual$settings$feature_center, reference$settings$feature_center,
      tolerance = 1e-12)
  }
})

test_that("JPEG rotations preserve native row-major pixels and rectangular patches", {
  skip_if_not_installed("jpeg")
  file <- tempfile(fileext = ".jpg")
  on.exit(unlink(file), add = TRUE)
  jpeg::writeJPEG(rotation_fixture(), file, quality = 1)
  decoded <- jpeg::readJPEG(file)
  for (angle in c(0, 90, 180, 270)) {
    expected <- explicit_rotation(decoded, angle)
    height <- dim(expected)[1]
    width <- dim(expected)[2]
    coordinates <- expand.grid(px = 2:(width - 1), py = seq_len(height))
    spots <- data.frame(spot_id = paste0("s", seq_len(nrow(coordinates))),
      x = coordinates$px * 100, y = coordinates$py * 200, coordinates)
    for (flip in c(FALSE, TRUE)) {
      for (channel in c("red", "green", "blue", "rgb")) {
        he <- extract_he_feature(file, spots, channel = channel, patch_size = c(3, 1),
          pixel_columns = c("px", "py"), flip_y = flip, rotate = angle)
        reference <- extract_he_feature(expected, spots, channel = channel,
          patch_size = c(3, 1), pixel_columns = c("px", "py"), flip_y = flip)
        expect_equal(he$spots, reference$spots, tolerance = 1e-12)
        expect_equal(he$settings$feature_center, reference$settings$feature_center,
          tolerance = 1e-12)
      }
    }
  }
})

test_that("rotated previews align image pixels and spot locations", {
  skip_if_not_installed("jpeg")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  file <- tempfile(fileext = ".jpg")
  on.exit(unlink(file), add = TRUE)
  image <- rotation_fixture()
  jpeg::writeJPEG(image, file, quality = 1)
  decoded <- jpeg::readJPEG(file)
  for (native in c(FALSE, TRUE)) {
    for (angle in c(90, 180, 270)) {
      expected <- explicit_rotation(if (native) decoded else image, angle)
      height <- dim(expected)[1]
      width <- dim(expected)[2]
      spots <- data.frame(spot_id = c("a", "b", "c"), x = c(100, 200, 300),
        y = c(1000, 2000, 3000), px = c(1, width, 2), py = c(1, height, 2))
      he <- extract_he_feature(if (native) file else image, spots, channel = "red",
        patch_size = 1, pixel_columns = c("px", "py"), flip_y = TRUE, rotate = angle)
      histology <- he$preview[[1]]
      feature_map <- he$preview[[2]]
      raster <- histology$layers[[1]]$geom_params$raster
      expect_equal(dim(raster), c(height, width))
      if (native) {
        expect_s3_class(raster, "nativeRaster")
        rgb <- vapply(0:2, function(j) {
          bitwAnd(bitwShiftR(as.integer(raster), 8L * j), 255L) / 255
        }, numeric(height * width))
      } else {
        rgb <- t(grDevices::col2rgb(as.vector(raster))) / 255
      }
      expected_rgb <- vapply(1:3, function(j) as.vector(t(expected[, , j])),
        numeric(height * width))
      expect_equal(unname(rgb), expected_rgb, tolerance = 1 / 255)
      expect_equal(feature_map$data$pixel_x, spots$px)
      expect_equal(feature_map$data$pixel_y, spots$py)
      expect_identical(feature_map$data$value, he$spots$he)
      image_frame <- ggplot2::ggplot_build(histology)$layout$panel_params[[1]]
      feature_frame <- ggplot2::ggplot_build(feature_map)$layout$panel_params[[1]]
      expect_equal(image_frame$x.range, c(0.5, width + 0.5))
      expect_equal(image_frame$y.range, c(0.5, height + 0.5))
      expect_equal(feature_frame$x.range, image_frame$x.range)
      expect_equal(feature_frame$y.range, image_frame$y.range)
    }
  }
})

test_that("rotation rejects unsupported angles", {
  image <- rotation_fixture()
  spots <- data.frame(spot_id = "a", x = 2, y = 2)
  for (angle in list(-90, 45, 360, NA_real_, Inf, c(90, 180), "90", TRUE, numeric())) {
    expect_error(extract_he_feature(image, spots, patch_size = 1, rotate = angle),
      "rotate must be")
  }
})

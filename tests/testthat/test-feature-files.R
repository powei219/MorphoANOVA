test_that("JPEG patch means preserve channels, orientation, and edge pixels", {
  skip_if_not_installed("jpeg")
  image <- array(0, c(37, 61, 3))
  image[, , 1] <- outer(seq_len(37), seq_len(61), function(y, x) (3 * y + x) %% 47 / 46)
  image[, , 2] <- outer(seq_len(37), seq_len(61), function(y, x) (y + 5 * x) %% 53 / 52)
  image[, , 3] <- outer(seq_len(37), seq_len(61), function(y, x) (7 * y + 2 * x) %% 59 / 58)
  file <- tempfile(fileext = ".jpg")
  on.exit(unlink(file), add = TRUE)
  jpeg::writeJPEG(image, file, quality = 1)
  decoded <- jpeg::readJPEG(file)
  spots <- data.frame(spot_id = c("a", "b", "c", "d"),
    x = c(100, 200, 300, 400), y = c(10, 20, 30, 40),
    pixel_x = c(1, 17, 43, 61), pixel_y = c(1, 23, 9, 37))

  for (flip in c(FALSE, TRUE)) {
    row <- if (flip) 38 - spots$pixel_y else spots$pixel_y
    for (channel in c("red", "green", "blue", "rgb")) {
      selected <- if (channel == "rgb") 1:3 else match(channel, c("red", "green", "blue"))
      expected <- vapply(seq_len(nrow(spots)), function(i) {
        # A 4-by-2 patch includes x-1 through x+2, and y through y+1.
        columns <- max(1, spots$pixel_x[i] - 1):min(61, spots$pixel_x[i] + 2)
        rows <- row[i]:min(37, row[i] + 1)
        mean(decoded[rows, columns, selected, drop = FALSE])
      }, numeric(1))
      expect_warning(he <- extract_he_feature(file, spots, channel = channel,
        patch_size = c(4, 2), pixel_columns = c("pixel_x", "pixel_y"),
        flip_y = flip), "edge patches")
      expect_equal(he$spots$he, expected - mean(expected), tolerance = 1e-12)
      expect_equal(he$settings$feature_center, mean(expected), tolerance = 1e-12)
    }
  }
  expect_equal(he$spots[names(spots)], spots)
})

test_that("large-image previews are bounded without changing feature resolution", {
  skip_if_not_installed("jpeg")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  image <- array(0, c(31, 2107, 3))
  image[, , 1] <- outer(seq_len(31), seq_len(2107), function(y, x) (x + 11 * y) %% 251 / 250)
  image[, , 2] <- 0.2
  image[, , 3] <- outer(seq_len(31), seq_len(2107), function(y, x) (3 * x + y) %% 241 / 240)
  file <- tempfile(fileext = ".jpeg")
  on.exit(unlink(file), add = TRUE)
  jpeg::writeJPEG(image, file, quality = 1)
  decoded <- jpeg::readJPEG(file)
  spots <- data.frame(spot_id = c("a", "b", "c", "d"),
    x = c(1, 2, 1031, 2107), y = c(1, 2, 17, 31))
  he <- extract_he_feature(file, spots, channel = "red", patch_size = 1)
  expected <- decoded[cbind(spots$y, spots$x, 1)]
  expect_equal(he$spots$he, expected - mean(expected), tolerance = 1e-12)
  expect_equal(he$settings$feature_center, mean(expected), tolerance = 1e-12)
  histology <- he$preview[[1]]
  feature_map <- he$preview[[2]]
  raster <- histology$layers[[1]]$geom_params$raster
  expect_true(all(dim(raster) <= 2000L))
  expect_equal(nrow(feature_map$data), nrow(spots))
  expect_equal(feature_map$data$pixel_x, spots$x)
  expect_equal(feature_map$data$pixel_y, 32 - spots$y)
  expect_identical(feature_map$data$value, he$spots$he)
  image_frame <- ggplot2::ggplot_build(histology)$layout$panel_params[[1]]
  feature_frame <- ggplot2::ggplot_build(feature_map)$layout$panel_params[[1]]
  expect_equal(image_frame$x.range, c(0.5, 2107.5))
  expect_equal(image_frame$y.range, c(0.5, 31.5))
  expect_equal(feature_frame$x.range, image_frame$x.range)
  expect_equal(feature_frame$y.range, image_frame$y.range)

  # The corners detect transposed or upside-down native-raster previews.
  if (inherits(raster, "nativeRaster")) {
    index <- c(1L, ncol(raster), length(raster) - ncol(raster) + 1L, length(raster))
    preview_red <- bitwAnd(as.integer(raster[index]), 255L) / 255
  } else {
    preview_red <- grDevices::col2rgb(c(raster[1, 1], raster[1, ncol(raster)],
      raster[nrow(raster), 1], raster[nrow(raster), ncol(raster)]))[1, ] / 255
  }
  expect_equal(unname(preview_red), decoded[cbind(c(1, 1, 31, 31),
    c(1, 2107, 1, 2107), 1)], tolerance = 1 / 255)
})

test_that("PNG extraction retains 16-bit channel precision", {
  skip_if_not_installed("png")
  # The fixture contains these six RGB pixels, row by row, at 16-bit precision.
  pixels <- rbind(c(1, 257, 1025), c(8193, 32769, 65534), c(65535, 1, 4097),
                  c(17, 511, 2049), c(16385, 49153, 60001), c(12345, 23456, 34567)) / 65535
  spots <- data.frame(spot_id = letters[1:6], x = rep(1:3, 2), y = rep(1:2, each = 3))
  file <- test_path("fixtures", "rgb16.png")
  for (channel in c("red", "green", "blue", "rgb")) {
    selected <- if (channel == "rgb") 1:3 else match(channel, c("red", "green", "blue"))
    he <- extract_he_feature(file, spots, channel = channel, patch_size = 1)
    expected <- rowMeans(pixels[, selected, drop = FALSE])
    expect_equal(he$spots$he, expected - mean(expected), tolerance = 1e-12)
    expect_equal(he$settings$feature_center, mean(expected), tolerance = 1e-12)
  }
})

test_that("array extraction retains non-8-bit values and ignores alpha", {
  image <- array(seq(0.001, 0.999, length.out = 8 * 11 * 4), c(8, 11, 4))
  image[, , 4] <- 0.4
  spots <- data.frame(spot_id = c("a", "b"), x = c(3, 8), y = c(2, 6))
  expected <- c(mean(image[1:3, 2:4, 1:3]), mean(image[5:7, 7:9, 1:3]))
  for (scale in c(1, 255)) {
    he <- extract_he_feature(image * scale, spots, channel = "rgb", patch_size = 3)
    expect_equal(he$spots$he, expected - mean(expected), tolerance = 1e-12)
    expect_equal(he$settings$feature_center, mean(expected), tolerance = 1e-12)
  }
  image[, , 4] <- 1
  he <- extract_he_feature(image, spots, channel = "rgb", patch_size = 3)
  expect_equal(he$spots$he, expected - mean(expected), tolerance = 1e-12)
  expect_equal(he$settings$feature_center, mean(expected), tolerance = 1e-12)
})

test_that("patch features equal independently centered pixel means", {
  image <- array(seq(0, 1, length.out = 12 * 15 * 3), c(12, 15, 3))
  spots <- data.frame(spot_id = c("s2", "s1"), x = c(800, 400), y = c(300, 100), px = c(8, 4), py = c(7, 4))
  he <- extract_he_feature(image, spots, channel = "green", patch_size = c(5, 3), pixel_columns = c("px", "py"))
  expected <- c(mean(image[6:8, 6:10, 2]), mean(image[3:5, 2:6, 2]))
  expect_equal(he$spots$he, expected - mean(expected), tolerance = 1e-12)
  expect_equal(he$settings$feature_center, mean(expected), tolerance = 1e-12)
  expect_equal(he$settings$intensity_scale, c(0, 1))
  expect_equal(he$spots[c("spot_id", "x", "y")], spots[c("spot_id", "x", "y")])
  rgb <- extract_he_feature(image, spots, patch_size = c(4, 2), pixel_columns = c("px", "py"))
  expected <- c(mean(image[7:8, 7:10, ]), mean(image[4:5, 3:6, ]))
  expect_equal(rgb$spots$he, expected - mean(expected), tolerance = 1e-12)
  expect_equal(rgb$settings$feature_center, mean(expected), tolerance = 1e-12)
  flip <- extract_he_feature(image * 255, spots, channel = "red", patch_size = 1,
    pixel_columns = c("px", "py"), flip_y = TRUE)
  expected <- image[cbind(c(6, 9), c(8, 4), 1)]
  expect_equal(flip$spots$he, expected - mean(expected), tolerance = 1e-12)
  expect_equal(flip$settings$feature_center, mean(expected), tolerance = 1e-12)
})

test_that("image edges and bad coordinate inputs are explicit", {
  image <- array(0.5, c(6, 6, 3))
  spots <- data.frame(spot_id = "a", x = 1, y = 1)
  expect_warning(he <- extract_he_feature(image, spots, patch_size = 3), "edge patches")
  expect_equal(he$settings$clipped_spot_ids, "a")
  expect_equal(he$spots$he, 0)
  expect_equal(he$settings$feature_center, 0.5)
  expect_error(extract_he_feature(image, transform(spots, x = 20), patch_size = 3), "outside")
  expect_error(extract_he_feature(image, rbind(spots, spots), patch_size = 1), "unique")
  expect_error(extract_he_feature(image, spots, patch_size = 1.5), "whole")
  expect_equal(extract_he_feature(array(1, c(1, 1, 3)), spots, patch_size = 1)$spots$he, 0)
})

test_that("paired preview preserves the image frame and uses returned features directly", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  image <- array(seq(0, 1, length.out = 12 * 15 * 3), c(12, 15, 3))
  spots <- data.frame(spot_id = c("a", "b"), x = c(800, 400), y = c(300, 100),
                      px = c(8, 4), py = c(7, 4))
  he <- extract_he_feature(image, spots, channel = "red", patch_size = 1,
                           pixel_columns = c("px", "py"), flip_y = TRUE)
  raw <- image[cbind(c(6, 9), c(8, 4), 1)]
  expect_equal(he$spots$he, raw - mean(raw), tolerance = 1e-12)
  expect_equal(he$settings$feature_center, mean(raw), tolerance = 1e-12)
  expect_s3_class(he$preview, "patchwork")
  histology <- he$preview[[1]]
  feature_map <- he$preview[[2]]
  expect_equal(length(histology$layers), 1L)
  expect_s3_class(histology$layers[[1]]$geom, "GeomRasterAnn")
  expect_s3_class(feature_map$layers[[1]]$geom, "GeomPoint")
  expect_equal(feature_map$data$pixel_x, c(8, 4))
  expect_equal(feature_map$data$pixel_y, c(7, 4))
  expect_identical(feature_map$data$value, he$spots$he)
  image_frame <- ggplot2::ggplot_build(histology)$layout$panel_params[[1]]
  feature_frame <- ggplot2::ggplot_build(feature_map)$layout$panel_params[[1]]
  expect_equal(image_frame$x.range, c(0.5, 15.5))
  expect_equal(image_frame$y.range, c(0.5, 12.5))
  expect_equal(feature_frame$x.range, image_frame$x.range)
  expect_equal(feature_frame$y.range, image_frame$y.range)
  colour <- feature_map$scales$get_scales("colour")
  expect_equal(colour$limits, c(-max(abs(raw - mean(raw))), max(abs(raw - mean(raw)))), tolerance = 1e-12)
  expect_equal(colour$map(c(-1, 0, 1) * max(abs(feature_map$data$value))),
                c("#3AB370", "#EAE7CC", "#FD1593"))
})

test_that("constant and single-spot features and their displays are zero", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  image <- array(0.5, c(6, 6, 3))
  spots <- data.frame(spot_id = c("a", "b"), x = c(2, 5), y = c(2, 5))
  for (n in 1:2) {
    he <- extract_he_feature(image, spots[seq_len(n), ], patch_size = 1)
    expect_equal(he$spots$he, rep(0, n))
    expect_equal(he$settings$feature_center, 0.5)
    expect_identical(he$preview[[2]]$data$value, he$spots$he)
    display <- ggplot2::ggplot_build(he$preview[[2]])$data[[1]]
    expect_equal(display$colour, rep("#EAE7CC", n))
  }
})

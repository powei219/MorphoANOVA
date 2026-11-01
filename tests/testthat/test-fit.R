core_example <- function() {
  set.seed(41)
  n <- 180
  spots <- data.frame(spot_id = paste0("s", seq_len(n)), x = runif(n),
                      y = runif(n), he = rnorm(n))
  expression <- rbind(
    geneA = spots$he * (1 + sin(2 * pi * spots$x)) + 2 * cos(2 * pi * spots$y) + rnorm(n),
    geneB = spots$he + rnorm(n)
  )
  colnames(expression) <- spots$spot_id
  list(expression = expression, spots = spots)
}

test_that("spatial tables are directly accessible and IDs align explicitly", {
  d <- core_example()
  fit <- fit_morphoanova(d$expression, d$spots, k = 8)
  expect_null(getS3method("fitted", "morphoanova_fit", optional = TRUE))
  expect_named(fit, c("statistics", "spatial", "status", "settings"))
  expect_named(fit$spatial, rownames(d$expression))
  expect_identical(class(fit$spatial$geneA), "data.frame")
  expect_identical(fit$spatial$geneA, .gene_spatial(fit, "geneA"))
  expect_named(fit$spatial$geneA, c("spot_id", "x", "y", "expression", "fitted",
    "he", "slide_wide_he", "location_specific_he", "total_he",
    "spatial_beyond_he", "total_he_se", "total_he_p", "total_he_q",
    "spatial_beyond_he_se", "spatial_beyond_he_p", "spatial_beyond_he_q"))
  expect_equal(fit$settings$n_genes, nrow(d$expression))
  expect_equal(fit$settings$n_locations, nrow(d$spots))
  expect_true(fit$settings$center)
  expect_equal(fit$settings$feature_center, mean(d$spots$he))
  expect_false(any(c("source_hash", "run_fingerprint", "genes") %in% names(fit$settings)))
  expect_true(nzchar(attr(fit, "internal", exact = TRUE)$provenance$run_fingerprint))
  reordered <- fit_morphoanova(d$expression, d$spots[rev(seq_len(nrow(d$spots))), ], k = 8)
  expect_equal(fit$statistics, reordered$statistics, tolerance = 1e-10)
  expect_equal(fit$spatial$geneA, reordered$spatial$geneA, tolerance = 1e-10)
  expect_error(fit_morphoanova(d$expression, d$spots[-1, ], k = 8), "same location IDs")
  duplicate <- d$spots
  duplicate$spot_id[2] <- duplicate$spot_id[1]
  expect_error(fit_morphoanova(d$expression, duplicate, k = 8), "unique")
  expect_error(fit_morphoanova(d$expression, d$spots, genes = "absent", k = 8), "absent")
  expect_error(fit_morphoanova(unname(d$expression), d$spots, k = 8), "names")
  constant <- d$spots
  constant$he <- 5
  expect_error(fit_morphoanova(d$expression, constant, k = 8), "no variation")
  expect_error(fit_morphoanova(d$expression, constant, k = 8, center = FALSE), "no variation")
  for (invalid in list(NA, NULL, 1, "TRUE", c(TRUE, FALSE))) {
    expect_error(fit_morphoanova(d$expression, d$spots, k = 8, center = invalid), "center")
  }
})

test_that("automatic centering matches a manually centered feature", {
  d <- core_example()
  automatic <- fit_morphoanova(d$expression, d$spots, k = 8)
  centered <- d$spots
  centered$he <- centered$he - mean(centered$he)
  manual <- fit_morphoanova(d$expression, centered, k = 8, center = FALSE)
  expect_equal(manual$statistics, automatic$statistics, tolerance = 1e-10)
  expect_equal(manual$spatial, automatic$spatial, tolerance = 1e-10)
  expect_false(manual$settings$center)
  expect_equal(manual$settings$feature_center, 0)
})

test_that("an uncentered feature is fitted as supplied", {
  d <- core_example()
  d$spots$he <- d$spots$he + 3
  fit <- fit_morphoanova(d$expression, d$spots, genes = "geneA", k = 8, center = FALSE)
  input <- data.frame(y = as.numeric(d$expression["geneA", ]),
                       x = d$spots$x, z = d$spots$y, h = d$spots$he)
  m1 <- mgcv::gam(y ~ 1, data = input, family = gaussian(), method = "REML")
  m2 <- mgcv::gam(y ~ h, data = input, family = gaussian(), method = "REML")
  m3 <- mgcv::gam(y ~ s(x, z, by = h, k = 8), data = input, family = gaussian(), method = "REML")
  m4 <- mgcv::gam(y ~ s(x, z, by = h, k = 8) + s(x, z, k = 8),
                  data = input, family = gaussian(), method = "REML")
  models <- list(m1, m2, m3, m4)
  comparisons <- lapply(seq_len(3), function(i) mgcv::anova.gam(models[[i]], models[[i + 1]], test = "F"))
  expect_equal(unname(unlist(fit$statistics[c("F21", "F32", "F43")])),
    vapply(comparisons, function(a) a[["F"]][2], numeric(1)), tolerance = 1e-8)
  expect_equal(unname(unlist(fit$statistics[c("p21", "p32", "p43")])),
    vapply(comparisons, function(a) a[["Pr(>F)"]][2], numeric(1)), tolerance = 1e-8)
  expect_equal(fit$spatial$geneA$he, d$spots$he)
  expect_gt(abs(mean(fit$spatial$geneA$he)), 2)
  expect_equal(fit$spatial$geneA$fitted, as.numeric(fitted(m4)), tolerance = 1e-9)
  expect_false(fit$settings$center)
  expect_equal(fit$settings$feature_center, 0)
})

test_that("F tests agree with GAM ANOVA while R2 retains fitted EDF", {
  d <- core_example()
  fit <- fit_morphoanova(d$expression, d$spots, genes = "geneA", k = 8)
  input <- data.frame(y = as.numeric(d$expression["geneA", ]),
                       x = d$spots$x, z = d$spots$y,
                       h = d$spots$he - mean(d$spots$he))
  m1 <- mgcv::gam(y ~ 1, data = input, method = "REML")
  m2 <- mgcv::gam(y ~ h, data = input, method = "REML")
  m3 <- mgcv::gam(y ~ s(x, z, by = h, k = 8), data = input, method = "REML")
  m4 <- mgcv::gam(y ~ s(x, z, by = h, k = 8) + s(x, z, k = 8), data = input, method = "REML")
  models <- list(m1, m2, m3, m4)
  rss <- vapply(models, function(m) sum(residuals(m)^2), numeric(1))
  edf <- vapply(models, function(m) sum(m$edf), numeric(1))
  n <- nrow(input)
  comparisons <- lapply(seq_len(3), function(i) stats::anova(models[[i]], models[[i + 1]], test = "F"))
  statistic <- vapply(comparisons, function(a) a[["F"]][2], numeric(1))
  p <- vapply(comparisons, function(a) a[["Pr(>F)"]][2], numeric(1))
  manual_edf_f <- diff(-rss) / diff(edf) / (rss[-1] / (n - edf[-1]))
  expect_gt(max(abs(statistic - manual_edf_f)), 1)
  adjusted <- 1 - rss / (n - edf) / stats::var(input$y)
  adjusted[1] <- 0
  expect_equal(unname(unlist(fit$statistics[c("F21", "F32", "F43")])), statistic, tolerance = 1e-8)
  expect_equal(unname(unlist(fit$statistics[c("p21", "p32", "p43")])), p, tolerance = 1e-8)
  expect_equal(unname(unlist(fit$statistics[c("r2_slide_wide_he", "r2_location_specific_he", "r2_spatial_beyond_he")])), diff(adjusted), tolerance = 1e-10)
  expect_equal(fit$statistics$r2_full, adjusted[4], tolerance = 1e-10)
  expect_equal(fit$statistics$r2_unexplained_variance, 1 - adjusted[4], tolerance = 1e-10)
  expect_named(fit$statistics, c("gene", "F21", "p21", "q21", "F32", "p32", "q32",
    "F43", "p43", "q43", "r2_slide_wide_he", "r2_location_specific_he",
    "r2_spatial_beyond_he", "r2_full", "r2_unexplained_variance"))
  expect_equal(attr(fit, "internal", exact = TRUE)$diagnostics$geneA$model_diagnostics$edf, edf, tolerance = 1e-10)
  expect_equal(attr(fit, "internal", exact = TRUE)$diagnostics$geneA$model_diagnostics$rss, rss, tolerance = 1e-10)
  provenance <- attr(fit, "internal", exact = TRUE)$provenance
  expect_identical(provenance$sequential_test, "mgcv::anova.gam(reduced, full, test = 'F')")
  expect_identical(provenance$adjusted_r2_edf, "sum(model$edf)")
  expect_identical(provenance$test_definition_version, "mgcv_anova_gam_F_v1")
  for (i in seq_along(comparisons)) {
    saved <- attr(fit, "internal", exact = TRUE)$diagnostics$geneA$diagnostics$comparisons[[i]]
    residual_df <- vapply(models[c(i, i + 1)], function(m) {
      correction <- if (is.null(m$edf2)) 0 else sum(m$edf2) - sum(m$edf)
      n - sum(m$edf1) - correction
    }, numeric(1))
    dispersion <- summary(models[[c(i, i + 1)[which.min(residual_df)]]])$dispersion
    expect_equal(unname(saved$residual_df), residual_df, tolerance = 1e-9)
    expect_equal(saved$df1, -diff(residual_df), tolerance = 1e-9)
    expect_equal(saved$df2, min(residual_df), tolerance = 1e-9)
    expect_equal(saved$dispersion, dispersion, tolerance = 1e-9)
    expect_equal(saved$raw_F, statistic[i], tolerance = 1e-9)
  }

  map <- fit$spatial[["geneA"]]
  xp <- predict(m4, type = "lpmatrix")
  he_index <- seq.int(m4$smooth[[1]]$first.para, m4$smooth[[1]]$last.para)
  spatial_index <- seq.int(m4$smooth[[2]]$first.para, m4$smooth[[2]]$last.para)
  he_basis <- xp[, he_index, drop = FALSE]
  spatial_basis <- scale(xp[, spatial_index, drop = FALSE], center = TRUE, scale = FALSE)
  expected_he_se <- sqrt(rowSums((he_basis %*% m4$Vp[he_index, he_index]) * he_basis))
  expected_spatial_se <- sqrt(rowSums((spatial_basis %*% m4$Vp[spatial_index, spatial_index]) * spatial_basis))
  expect_equal(map$total_he_se, unname(expected_he_se), tolerance = 1e-9)
  expect_equal(map$spatial_beyond_he_se, unname(expected_spatial_se), tolerance = 1e-9)
  expect_equal(map$fitted, as.numeric(fitted(m4)), tolerance = 1e-9)
})

test_that("M4 components retain centering, reconstruction and raw variance increments", {
  d <- core_example()
  fit <- fit_morphoanova(d$expression, d$spots, k = 8)
  for (gene in fit$statistics$gene) {
    map <- fit$spatial[[gene]]
    diagnostics <- attr(fit, "internal", exact = TRUE)$diagnostics[[gene]]$diagnostics
    expect_equal(mean(map$he), 0, tolerance = 1e-14)
    expect_equal(mean(map$spatial_beyond_he), 0, tolerance = 1e-12)
    expect_equal(mean(map$location_specific_he / map$he), 0, tolerance = 1e-12)
    expect_equal(map$total_he, map$slide_wide_he + map$location_specific_he, tolerance = 1e-12)
    expect_equal(map$fitted, diagnostics$baseline + map$total_he + map$spatial_beyond_he, tolerance = 1e-9)
    expect_equal(map$total_he_q, p.adjust(map$total_he_p, "BH"))
    expect_equal(map$spatial_beyond_he_q, p.adjust(map$spatial_beyond_he_p, "BH"))
  }
  for (test in c("21", "32", "43")) {
    expect_equal(fit$statistics[[paste0("q", test)]], p.adjust(fit$statistics[[paste0("p", test)]], "BH"))
  }
  gene_b <- fit$statistics[fit$statistics$gene == "geneB", ]
  expect_lt(gene_b$r2_location_specific_he, 0)
  expect_equal(sum(unlist(gene_b[c("r2_slide_wide_he", "r2_location_specific_he", "r2_spatial_beyond_he", "r2_unexplained_variance")])), 1, tolerance = 1e-12)
})

test_that("structural zero pointwise tests are preserved", {
  d <- core_example()
  d$spots$he <- c(0, seq_len(89), -seq_len(89), 0)
  fit <- fit_morphoanova(d$expression, d$spots, genes = "geneA", k = 8)
  map <- fit$spatial[["geneA"]]
  zero <- which(map$he == 0)
  expect_length(zero, 2)
  expect_equal(map$total_he[zero], c(0, 0))
  expect_equal(map$total_he_se[zero], c(0, 0))
  expect_true(all(is.na(map$total_he_p[zero])))
  expect_true(all(is.na(map$total_he_q[zero])))
})

test_that("ANOVA preserves signed df and keeps unavailable tests distinct from nulls", {
  compare <- get(".sequential_test", envir = environment(fit_morphoanova))
  d <- core_example()
  input <- data.frame(y = as.numeric(d$expression["geneA", ]), h = d$spots$he)
  m1 <- mgcv::gam(y ~ 1, data = input, method = "REML")
  m2 <- mgcv::gam(y ~ h, data = input, method = "REML")
  forward <- compare(m1, m2)
  reversed <- compare(m2, m1)
  expect_lt(reversed$df1, 0)
  expect_equal(reversed$F, forward$F, tolerance = 1e-12)
  expect_equal(reversed$p, forward$p, tolerance = 1e-12)
  expect_equal(reversed$dispersion, m2$sig2, tolerance = 1e-12)
  expect_identical(reversed$reason, "")

  invalid <- compare(m2, m2)
  expect_true(is.na(invalid$F))
  expect_true(is.na(invalid$p))
  expect_match(invalid$reason, "zero difference")

  # Isolate ANOVA's negative-F handling without relying on a near-singular fit.
  worse <- m2
  worse$deviance <- m1$deviance + 10
  negative <- compare(m1, worse)
  expect_lt(negative$raw_F, 0)
  expect_equal(negative$deviance_difference, -10, tolerance = 1e-12)
  expect_true(is.na(negative$F))
  expect_true(is.na(negative$p))
  expect_match(negative$reason, "negative algebraic F")
  short <- mgcv::gam(y ~ h, data = input[-1, ], method = "REML")
  bad <- compare(m1, short)
  expect_true(is.na(bad$F))
  expect_true(is.na(bad$p))
  expect_match(bad$reason, "same size of dataset")
})

test_that("unavailable ANOVA comparisons retain maps and partial gene status", {
  compare <- get(".sequential_test", envir = environment(fit_morphoanova))
  local_mocked_bindings(.sequential_test = function(reduced, full) {
    if (length(full$smooth) == 1L) compare(full, full) else compare(reduced, full)
  })
  d <- core_example()
  expect_warning(fit <- fit_morphoanova(d$expression, d$spots, k = 8), "2 gene\\(s\\) have unavailable tests")
  expect_true(all(fit$status$status == "partial"))
  expect_true(all(grepl("32:.*zero difference", fit$status$message)))
  expect_true(all(is.na(fit$statistics$F32)))
  expect_true(all(is.na(fit$statistics$p32)))
  expect_true(all(is.na(fit$statistics$q32)))
  expect_true(all(is.finite(fit$statistics$F21)))
  expect_true(all(is.finite(fit$statistics$F43)))
  expect_equal(nrow(.gene_spatial(fit, "geneA")), nrow(d$spots))
  expect_true(all(is.finite(.gene_spatial(fit, "geneA")$total_he)))
})

test_that("failed genes keep their status and a named empty spatial entry", {
  d <- core_example()
  expression <- rbind(d$expression, constant = 1)
  expect_warning(fit <- fit_morphoanova(expression, d$spots, k = 8), "1 gene\\(s\\) failed")
  expect_identical(fit$status$status, c("complete", "complete", "failed"))
  expect_true(all(is.na(fit$statistics[3, -1])))
  expect_match(fit$status$message[3], "no variation")
  expect_named(fit$spatial, rownames(expression))
  expect_null(fit$spatial$constant)
  expect_error(.gene_spatial(fit, "constant"), "fit failed")
  expect_error(.gene_spatial(fit), "Choose one")
})

test_that("disk results, compatible resume, and stale-input rejection work", {
  d <- core_example()
  directory <- tempfile("morphoanova-test-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  memory <- fit_morphoanova(d$expression, d$spots, k = 8)
  disk <- fit_morphoanova(d$expression, d$spots, k = 8, output_dir = directory)
  expect_identical(names(disk), names(memory))
  expect_type(disk$spatial, "character")
  expect_named(disk$spatial, rownames(d$expression))
  expect_true(all(file.exists(disk$spatial)))
  expect_identical(readRDS(disk$spatial[["geneA"]])$spatial, memory$spatial$geneA)
  expect_named(attr(disk, "internal", exact = TRUE)$files, rownames(d$expression))
  expect_equal(memory$statistics, disk$statistics, tolerance = 1e-10)
  expect_equal(.gene_spatial(memory, "geneA"), .gene_spatial(disk, "geneA"), tolerance = 1e-10)
  before <- unname(tools::md5sum(attr(disk, "internal", exact = TRUE)$files))
  resumed <- fit_morphoanova(d$expression, d$spots, k = 8, output_dir = directory, resume = TRUE)
  expect_identical(unname(tools::md5sum(attr(resumed, "internal", exact = TRUE)$files)), before)
  expect_equal(resumed$statistics, disk$statistics)
  changed <- d$expression
  changed[1, 1] <- changed[1, 1] + 0.1
  expect_error(fit_morphoanova(changed, d$spots, k = 8, output_dir = directory, resume = TRUE), "Cannot resume")
  expect_error(fit_morphoanova(d$expression, d$spots, k = 9, output_dir = directory, resume = TRUE), "Cannot resume")
  expect_error(fit_morphoanova(d$expression, d$spots, k = 8, output_dir = directory), "contains a run")

  file <- file.path(directory, "memory-fit.rds")
  saveRDS(memory, file)
  restored <- readRDS(file)
  expect_identical(.gene_spatial(restored, "geneA"), memory$spatial$geneA)
  expect_identical(.gene_spatial(readRDS(file.path(directory, "fit.rds")), "geneB"),
    memory$spatial$geneB)

  # Plotting can still read fits saved with previous result layouts.
  legacy <- structure(list(statistics = memory$statistics, status = memory$status,
    settings = attr(memory, "internal", exact = TRUE)$provenance,
    per_gene = lapply(rownames(d$expression), function(gene) {
      list(gene = gene, spatial = memory$spatial[[gene]],
        status = memory$status[memory$status$gene == gene, , drop = FALSE])
    }), files = NULL), class = "morphoanova_fit")
  names(legacy$per_gene) <- rownames(d$expression)
  names(legacy$statistics)[names(legacy$statistics) == "r2_unexplained_variance"] <-
    "residual_fraction"
  saveRDS(legacy, file)
  expect_identical(.gene_spatial(readRDS(file), "geneA"), memory$spatial$geneA)
  legacy$per_gene <- NULL
  legacy$files <- attr(disk, "internal", exact = TRUE)$files
  saveRDS(legacy, file)
  expect_identical(.gene_spatial(readRDS(file), "geneB"), memory$spatial$geneB)

  # Access to one disk-backed gene must not load unrelated gene files.
  files <- attr(disk, "internal", exact = TRUE)$files
  unlink(files[["geneB"]])
  expect_identical(.gene_spatial(disk, "geneA"), memory$spatial$geneA)
  stored <- readRDS(files[["geneA"]])
  stored$gene <- "geneB"
  saveRDS(stored, files[["geneA"]])
  expect_error(.gene_spatial(disk, "geneA"), "does not match")
  stored$gene <- "geneA"
  stored$status$status <- "failed"
  stored$status$message <- "stored gene failure"
  saveRDS(stored, files[["geneA"]])
  expect_error(.gene_spatial(disk, "geneA"), "stored gene failure")
})

test_that("changing centering or the F test invalidates existing checkpoints", {
  d <- core_example()
  d$spots$he <- d$spots$he - mean(d$spots$he)
  directory <- tempfile("morphoanova-test-definition-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  saved <- fit_morphoanova(d$expression, d$spots, genes = "geneA", k = 8, output_dir = directory)
  files <- c(file.path(directory, "run.rds"), attr(saved, "internal", exact = TRUE)$files)
  before <- unname(tools::md5sum(files))
  expect_equal(saved$settings$feature_center, 0, tolerance = 1e-14)
  expect_error(fit_morphoanova(d$expression, d$spots, genes = "geneA", k = 8,
    output_dir = directory, resume = TRUE, center = FALSE), "Cannot resume")
  expect_identical(unname(tools::md5sum(files)), before)
  local_mocked_bindings(.sequential_test = function(reduced, full) {
    stop("Changed test implementation must not execute when resuming an old run")
  })
  expect_error(fit_morphoanova(d$expression, d$spots, genes = "geneA", k = 8,
    output_dir = directory, resume = TRUE), "Cannot resume.*source")
  expect_identical(unname(tools::md5sum(files)), before)
})

test_that("uncentered memory, disk, parallel and sparse results agree", {
  skip_on_cran()
  d <- core_example()
  d$spots$he <- d$spots$he + 3
  serial <- fit_morphoanova(d$expression, d$spots, k = 8, center = FALSE)
  parallel <- fit_morphoanova(d$expression, d$spots, k = 8, cores = 2, center = FALSE)
  expect_equal(serial$statistics, parallel$statistics, tolerance = 1e-10)
  expect_equal(.gene_spatial(serial, "geneB"), .gene_spatial(parallel, "geneB"), tolerance = 1e-10)
  skip_if_not_installed("Matrix")
  sparse <- fit_morphoanova(Matrix::Matrix(d$expression, sparse = TRUE), d$spots, k = 8, center = FALSE)
  expect_equal(serial$statistics, sparse$statistics, tolerance = 1e-10)
  dense <- fit_morphoanova(Matrix::Matrix(d$expression, sparse = FALSE), d$spots, k = 8, center = FALSE)
  expect_equal(serial$statistics, dense$statistics, tolerance = 1e-10)
  sparse_expression <- Matrix::Matrix(d$expression, sparse = TRUE)
  sparse_parallel <- fit_morphoanova(sparse_expression, d$spots, k = 8, cores = 2, center = FALSE)
  expect_true(all(sparse_parallel$status$status == "complete"))
  expect_equal(serial$statistics, sparse_parallel$statistics, tolerance = 1e-10)
  expect_equal(.gene_spatial(serial, "geneA"), .gene_spatial(sparse_parallel, "geneA"), tolerance = 1e-10)
  dense_parallel <- fit_morphoanova(Matrix::Matrix(d$expression, sparse = FALSE), d$spots, k = 8, cores = 2, center = FALSE)
  expect_true(all(dense_parallel$status$status == "complete"))
  expect_equal(serial$statistics, dense_parallel$statistics, tolerance = 1e-10)
  directory <- tempfile("sparse-parallel-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  disk <- fit_morphoanova(sparse_expression, d$spots, k = 8, cores = 2, output_dir = directory, center = FALSE)
  expect_true(all(disk$status$status == "complete"))
  expect_equal(serial$statistics, disk$statistics, tolerance = 1e-10)
  expect_equal(.gene_spatial(serial, "geneB"), .gene_spatial(disk, "geneB"), tolerance = 1e-10)
  before <- unname(tools::md5sum(attr(disk, "internal", exact = TRUE)$files))
  resumed <- fit_morphoanova(sparse_expression, d$spots, k = 8, cores = 2,
    output_dir = directory, resume = TRUE, center = FALSE)
  expect_identical(unname(tools::md5sum(attr(resumed, "internal", exact = TRUE)$files)), before)
  expect_equal(serial$statistics, resumed$statistics, tolerance = 1e-10)
})

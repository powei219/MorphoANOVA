enrichment_example <- function() {
  symbols <- c("TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "ACTB", "GAPDH", "STAT3")
  entrez <- c("7157", "672", "675", "1956", "3845", "60", "2597", "6774")
  statistics <- data.frame(gene = c(symbols, "unmapped_gene", "failed_gene"),
    F21 = 10, p21 = 0.01, q21 = 0.02, F32 = 10:1, p32 = 0.01, q32 = 0.02,
    F43 = 1:10, p43 = 0.01, q43 = 0.02,
    r2_slide_wide_he = 0.1, r2_location_specific_he = 0.1, r2_spatial_beyond_he = 0.1)
  status <- data.frame(gene = statistics$gene, status = c(rep("complete", 9), "failed"))
  fit <- structure(list(statistics = statistics, status = status), class = "morphoanova_fit")
  cache <- list(
    TERM2GENE = data.frame(term = rep(c("path_a", "path_b", "path_c"), each = 4),
      gene = c(entrez[1:4], entrez[5:8], entrez[c(1, 2, 7, 8)])),
    TERM2NAME = data.frame(term = c("path_a", "path_b", "path_c"),
      name = c("First gene set", "Second gene set", "Overlapping gene set")),
    provenance = list(source = "Deterministic test gene sets"))
  list(fit = fit, cache = cache, symbols = symbols, entrez = entrez)
}

test_that("enrichment uses tested background and agrees with an independent test", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  d <- enrichment_example()
  ranked <- rank_genes(d$fit, top_n = 2, components = "location_specific_he")
  original <- ranked
  result <- suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "KEGG", pvalue_cutoff = 1, qvalue_cutoff = 1,
    min_gene_set_size = 1, kegg_cache = d$cache))
  component <- result$location_specific_he
  expect_s4_class(component$KEGG, "enrichResult")
  expect_null(component$GO)
  expect_setequal(component$genes, d$entrez[1:2])
  expect_setequal(component$universe, d$entrez)
  expect_setequal(component$background$KEGG, d$entrez)
  expect_identical(result$unmapped$location_specific_he$universe, c("unmapped_gene", "failed_gene"))
  expect_false("failed_gene" %in% result$mapping$gene)
  expect_equal(result$summary$n_universe, 10)
  expect_equal(result$summary$n_annotation_background, 8)
  table <- as.data.frame(component$KEGG)
  expect_setequal(table$ID, c("path_a", "path_c"))
  expect_equal(table$pvalue, rep(stats::phyper(1, 4, 4, 2, lower.tail = FALSE), 2))
  expect_equal(table$p.adjust, stats::p.adjust(table$pvalue, "BH"))
  expect_identical(ranked, original)
  expect_equal(result$settings$background_source, "rank_genes")
  expect_equal(result$settings$kegg_provenance, d$cache$provenance)
  expect_match(result$settings$annotation_version, "^[0-9]+\\.")
})

test_that("unavailable component tests do not change the background", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  d <- enrichment_example()
  d$fit$statistics$p43[1:4] <- NA_real_
  ranked <- rank_genes(d$fit, top_n = 2)
  # The unmapped symbol ranks first for F43; keep it visible in the mapping report.
  result <- suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "KEGG", pvalue_cutoff = 1, qvalue_cutoff = 1,
    min_gene_set_size = 1, kegg_cache = d$cache))
  expect_setequal(result$location_specific_he$universe, d$entrez)
  expect_setequal(result$spatial_beyond_he$universe, d$entrez)
  expect_equal(result$unmapped$spatial_beyond_he$genes, "unmapped_gene")
  expect_equal(result$summary$n_mapped, c(2, 1))

  ranked_ids <- ranked
  ranked_ids$settings$components <- "location_specific_he"
  ranked_ids$spatial_beyond_he <- NULL
  ranked_ids$location_specific_he$gene <- d$entrez[1:2]
  ranked_ids$universe <- list(location_specific_he = d$entrez)
  converted <- suppressMessages(.enrich_ranked_genes(ranked_ids, org.Hs.eg.db::org.Hs.eg.db,
    key_type = "ENTREZID", databases = "KEGG", pvalue_cutoff = 1, qvalue_cutoff = 1,
    min_gene_set_size = 1, kegg_cache = d$cache))
  expect_equal(as.data.frame(converted$location_specific_he$KEGG),
    as.data.frame(result$location_specific_he$KEGG))
  expect_equal(converted$mapping$gene, converted$mapping$ENTREZID)
  expect_length(converted$unmapped$location_specific_he$universe, 0)
})

test_that("GO and KEGG use explicit database cutoffs and standard result objects", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  d <- enrichment_example()
  ranked <- rank_genes(d$fit, top_n = 2, components = "location_specific_he")
  cache_path <- tempfile(fileext = ".rds")
  on.exit(unlink(cache_path))
  saveRDS(d$cache, cache_path)
  result <- suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    pvalue_cutoff = c(KEGG = 1, GO = 1), qvalue_cutoff = c(GO = 1, KEGG = 0.5),
    min_gene_set_size = 1, kegg_cache = cache_path, simplify_go = FALSE))
  go <- result$location_specific_he$GO
  expect_s4_class(go, "enrichResult")
  expect_equal(go@ontology, "BP")
  expect_true(go@readable)
  expect_equal(go@pAdjustMethod, "BH")
  expect_equal(go@pvalueCutoff, 1)
  expect_equal(go@qvalueCutoff, 1)
  expect_equal(result$location_specific_he$KEGG@qvalueCutoff, 0.5)
  expect_equal(result$settings$qvalue_cutoff, c(GO = 1, KEGG = 0.5))
  expect_equal(result$settings$kegg_source, basename(cache_path))
  expect_equal(result$settings$kegg_cache_sha256, digest::digest(file = cache_path, algo = "sha256"))
  expect_setequal(result$location_specific_he$background$GO, go@universe)
})

test_that("empty selections are distinct from unavailable identifier mappings", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  d <- enrichment_example()
  ranked <- rank_genes(d$fit, top_n = 0, components = "location_specific_he")
  result <- suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "KEGG", min_gene_set_size = 1, kegg_cache = d$cache))
  expect_null(result$location_specific_he$KEGG)
  expect_equal(result$summary$status, "no_selected_genes")
  expect_equal(result$summary$n_mapped_universe, 8)
  expect_true(is.na(result$summary$n_annotation_background))
  ranked$location_specific_he <- data.frame(gene = "unmapped_gene")
  expect_warning(empty <- suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "KEGG", kegg_cache = d$cache)), "no_mapped_genes")
  expect_equal(empty$summary$status, "no_mapped_genes")
  expect_equal(empty$unmapped$location_specific_he$genes, "unmapped_gene")
  ranked$universe <- list(location_specific_he = "unmapped_gene")
  expect_error(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db, databases = "GO"),
    "No universe genes mapped")
})

test_that("backgrounds, cutoffs and cached annotations are checked before enrichment", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  d <- enrichment_example()
  ranked <- rank_genes(d$fit, top_n = 2, components = "location_specific_he")
  expect_error(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    universe = "ACTB", databases = "GO"), "Every selected gene")
  expect_error(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    key_type = "invalid"), "key_type")
  expect_error(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    pvalue_cutoff = c(0.1, 0.2)), "named")
  expect_error(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    qvalue_cutoff = NA_real_), "qvalue_cutoff")
  bad_cache <- d$cache
  bad_cache$TERM2GENE$gene[1] <- "hsa:7157"
  expect_error(suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "KEGG", kegg_cache = bad_cache)), "unprefixed Entrez")
  bad_cache <- d$cache
  bad_cache$TERM2NAME <- bad_cache$TERM2NAME[-1, ]
  expect_error(suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "KEGG", kegg_cache = bad_cache)), "one description")
})


test_that("annotation backgrounds and GO simplification match independent analyses", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  d <- enrichment_example()
  ranked <- rank_genes(d$fit, top_n = 2, components = "location_specific_he")
  result <- suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "KEGG", universe = "annotation", pvalue_cutoff = 1, qvalue_cutoff = 1,
    min_gene_set_size = 1, kegg_cache = d$cache))
  independent <- clusterProfiler::enricher(gene = d$entrez[1:2],
    TERM2GENE = d$cache$TERM2GENE, TERM2NAME = d$cache$TERM2NAME,
    pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
  expect_equal(as.data.frame(result$location_specific_he$KEGG), as.data.frame(independent))
  expect_identical(result$settings$background_source, "annotation")
  expect_true(is.na(result$summary$n_universe))
  expect_null(result$location_specific_he$universe)
  expect_setequal(result$location_specific_he$background$KEGG, d$entrez)

  go_terms <- data.frame(ID = c("GO:0006915", "GO:0097190", "GO:0042981"),
    Description = c("apoptotic process", "apoptotic signaling pathway", "regulation of apoptotic process"),
    GeneRatio = "2/2", BgRatio = "4/8", pvalue = c(.001, .01, .02),
    p.adjust = c(.003, .015, .02), qvalue = c(.003, .015, .02),
    geneID = "7157/672", Count = 2)
  raw <- methods::new("enrichResult", result = go_terms, ontology = "BP",
    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", universe = d$entrez)
  testthat::local_mocked_bindings(enrichGO = function(...) raw, .package = "clusterProfiler")
  result <- suppressMessages(.enrich_ranked_genes(ranked, org.Hs.eg.db::org.Hs.eg.db,
    databases = "GO", pvalue_cutoff = 1, qvalue_cutoff = 1, min_gene_set_size = 1))
  simplified <- clusterProfiler::simplify(raw, cutoff = .7,
    by = "p.adjust", select_fun = min, measure = "Wang")
  expect_equal(result$location_specific_he$GO_full, raw)
  expect_equal(as.data.frame(result$location_specific_he$GO), as.data.frame(simplified))
  expect_gte(nrow(as.data.frame(raw)), nrow(as.data.frame(simplified)))
  expect_identical(result$settings$simplify_go, "location_specific_he")
  expect_identical(result$settings$similarity_measure, "Wang")
  both <- rank_genes(d$fit, top_n = 2)
  selective <- suppressMessages(.enrich_ranked_genes(both, org.Hs.eg.db::org.Hs.eg.db,
    databases = "GO", simplify_go = "location_specific_he"))
  expect_equal(as.data.frame(selective$location_specific_he$GO), as.data.frame(simplified))
  expect_equal(selective$spatial_beyond_he$GO, raw)
})

test_that("ranked genes can be analyzed once and replotted from saved results", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  d <- enrichment_example()
  ranked <- rank_genes(d$fit, top_n = 2, components = "location_specific_he")
  p <- suppressMessages(plot_enrichment(ranked, databases = "KEGG", top_n = 1,
    pvalue_cutoff = 1, qvalue_cutoff = 1, min_gene_set_size = 1, kegg_cache = d$cache))
  result <- attr(p, "enrichment")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 1L)
  expect_equal(nrow(as.data.frame(result$location_specific_he$KEGG)), 2L)
  expect_equal(result$summary$n_universe, nrow(d$fit$statistics))
  expect_setequal(result$location_specific_he$background$KEGG, d$entrez)
  saved <- plot_enrichment(result, top_n = 2)
  reused <- plot_enrichment(p, top_n = 2)
  expect_identical(saved$data, reused$data)
  expect_equal(nrow(saved$data), 2L)
  expect_identical(attr(saved, "enrichment"), result)
  expect_identical(attr(reused, "enrichment"), result)
  expect_error(plot_enrichment(result, universe = "annotation"), "ranked genes")
  expect_error(plot_enrichment(ranked, organism = "mmu"), "Supply OrgDb")
})

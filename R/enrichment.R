# Internal calculation used by plot_enrichment().
.enrich_ranked_genes <- function(ranked, OrgDb, organism = "hsa", key_type = "SYMBOL",
                               ontology = "BP", databases = c("GO", "KEGG"),
                               pvalue_cutoff = 0.05, qvalue_cutoff = 0.2,
                               min_gene_set_size = 10, max_gene_set_size = 500,
                               universe = NULL, kegg_cache = NULL,
                               simplify_go = TRUE, similarity_cutoff = 0.7) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Install clusterProfiler and AnnotationDbi to run enrichment.", call. = FALSE)
  }
  if (!is.list(ranked) || !is.list(ranked$settings)) {
    stop("ranked must be a result from rank_genes().", call. = FALSE)
  }
  components <- ranked$settings$components
  if (!is.character(components) || !length(components) ||
      anyNA(components) || anyDuplicated(components) ||
      !all(components %in% c("slide_wide_he", "location_specific_he", "spatial_beyond_he")) ||
      !all(components %in% names(ranked))) {
    stop("ranked must be a result from rank_genes().", call. = FALSE)
  }
  if (!inherits(OrgDb, "OrgDb")) stop("OrgDb must be an AnnotationDbi OrgDb object.", call. = FALSE)
  if (!is.character(key_type) || length(key_type) != 1L || is.na(key_type) ||
      !key_type %in% AnnotationDbi::keytypes(OrgDb)) {
    stop("Choose key_type from AnnotationDbi::keytypes(OrgDb).", call. = FALSE)
  }
  ontology <- match.arg(ontology, c("BP", "MF", "CC", "ALL"))
  if (!is.character(databases) || !length(databases) || anyNA(databases) ||
      anyDuplicated(databases) || !all(databases %in% c("GO", "KEGG"))) {
    stop("databases must contain GO, KEGG, or both, without duplicates.", call. = FALSE)
  }
  if (!is.character(organism) || length(organism) != 1L || is.na(organism) || !nzchar(organism)) {
    stop("organism must be a KEGG organism code.", call. = FALSE)
  }
  p_cutoff <- .enrichment_cutoffs(pvalue_cutoff, databases, "pvalue_cutoff")
  q_cutoff <- .enrichment_cutoffs(qvalue_cutoff, databases, "qvalue_cutoff")
  sizes <- c(min_gene_set_size, max_gene_set_size)
  if (!is.numeric(min_gene_set_size) || length(min_gene_set_size) != 1L ||
      !is.numeric(max_gene_set_size) || length(max_gene_set_size) != 1L ||
      any(!is.finite(sizes)) || any(sizes < 1) || any(sizes != floor(sizes)) ||
      min_gene_set_size > max_gene_set_size) {
    stop("Gene-set sizes must be positive whole numbers with minimum <= maximum.", call. = FALSE)
  }

  if (!is.numeric(similarity_cutoff) || length(similarity_cutoff) != 1L ||
      !is.finite(similarity_cutoff) || similarity_cutoff < 0 || similarity_cutoff > 1) {
    stop("similarity_cutoff must be one number in [0, 1].", call. = FALSE)
  }
  if (identical(simplify_go, TRUE)) simplify_go <- components
  else if (identical(simplify_go, FALSE)) simplify_go <- character()
  if (!is.character(simplify_go) || anyNA(simplify_go) || anyDuplicated(simplify_go) ||
      !all(simplify_go %in% components)) {
    stop("simplify_go must be TRUE, FALSE, or names of ranked components.", call. = FALSE)
  }

  annotation_background <- identical(universe, "annotation")
  background_source <- if (annotation_background) "annotation" else if (is.null(universe)) "rank_genes" else "provided"
  if (is.null(universe)) universe <- ranked$universe
  if (annotation_background) universe <- stats::setNames(lapply(components,
    function(component) ranked[[component]]$gene), components)
  if (is.character(universe)) universe <- stats::setNames(rep(list(universe), length(components)), components)
  if (!is.list(universe) || is.null(names(universe)) || anyDuplicated(names(universe)) ||
      !setequal(names(universe), components)) {
    stop("universe must be a character vector or a list named for each ranked component.", call. = FALSE)
  }
  selected <- list()
  for (component in components) {
    if (!is.data.frame(ranked[[component]]) || !is.character(ranked[[component]]$gene) ||
        anyNA(ranked[[component]]$gene) || any(!nzchar(ranked[[component]]$gene))) {
      stop("Each ranked component needs a gene column of input identifiers.", call. = FALSE)
    }
    selected[[component]] <- unique(ranked[[component]]$gene)
    background <- universe[[component]]
    if (!is.character(background) || anyNA(background) || any(!nzchar(background))) {
      stop("Each universe must contain complete character gene identifiers.", call. = FALSE)
    }
    universe[[component]] <- unique(background)
    if (length(setdiff(selected[[component]], background))) {
      stop("Every selected gene must belong to its component's universe: ", component, ".", call. = FALSE)
    }
  }

  all_ids <- unique(unlist(universe[components], use.names = FALSE))
  valid_ids <- intersect(all_ids, AnnotationDbi::keys(OrgDb, keytype = key_type))
  mapping <- data.frame(gene = character(), ENTREZID = character())
  if (length(valid_ids)) {
    mapped <- AnnotationDbi::select(OrgDb, keys = valid_ids, keytype = key_type, columns = "ENTREZID")
    mapping <- unique(data.frame(gene = as.character(mapped[[key_type]]),
                                 ENTREZID = as.character(mapped$ENTREZID)))
    mapping <- mapping[!is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID), , drop = FALSE]
    rownames(mapping) <- NULL
  }
  if (length(all_ids) && !nrow(mapping)) {
    stop("No universe genes mapped to Entrez IDs; check key_type and OrgDb.", call. = FALSE)
  }

  cache_source <- "online KEGG via clusterProfiler"
  cache_sha256 <- NULL
  cache_provenance <- NULL
  if (!is.null(kegg_cache)) {
    if (!"KEGG" %in% databases) stop("kegg_cache requires KEGG in databases.", call. = FALSE)
    if (is.character(kegg_cache)) {
      if (length(kegg_cache) != 1L || is.na(kegg_cache) || !file.exists(kegg_cache)) {
        stop("kegg_cache must name an existing RDS file.", call. = FALSE)
      }
      cache_source <- basename(kegg_cache)
      cache_sha256 <- digest::digest(file = kegg_cache, algo = "sha256")
      kegg_cache <- readRDS(kegg_cache)
    } else {
      cache_source <- "provided annotation tables"
      cache_sha256 <- digest::digest(kegg_cache, algo = "sha256")
    }
    if (!is.list(kegg_cache) || !all(c("TERM2GENE", "TERM2NAME") %in% names(kegg_cache))) {
      stop("kegg_cache needs TERM2GENE and TERM2NAME data frames.", call. = FALSE)
    }
    for (table in c("TERM2GENE", "TERM2NAME")) {
      d <- kegg_cache[[table]]
      if (!is.data.frame(d) || ncol(d) != 2L || !nrow(d)) {
        stop(table, " must be a nonempty two-column data frame.", call. = FALSE)
      }
      d[] <- lapply(d, as.character)
      if (anyNA(d) || any(!nzchar(as.matrix(d)))) {
        stop(table, " must not contain missing or empty values.", call. = FALSE)
      }
      kegg_cache[[table]] <- unique(d)
    }
    if (any(!grepl("^[0-9]+$", kegg_cache$TERM2GENE[[2L]]))) {
      stop("TERM2GENE must use unprefixed Entrez IDs in its second column.", call. = FALSE)
    }
    if (anyDuplicated(kegg_cache$TERM2NAME[[1L]]) ||
        length(setdiff(kegg_cache$TERM2GENE[[1L]], kegg_cache$TERM2NAME[[1L]]))) {
      stop("TERM2NAME needs one description for every TERM2GENE term.", call. = FALSE)
    }
    cache_provenance <- kegg_cache$provenance
  }

  out <- list()
  unmapped <- list()
  summary_rows <- list()
  for (component in components) {
    background <- if (annotation_background) NULL else unique(mapping$ENTREZID[mapping$gene %in% universe[[component]]])
    genes <- unique(mapping$ENTREZID[mapping$gene %in% selected[[component]]])
    unmapped[[component]] <- list(genes = setdiff(selected[[component]], mapping$gene),
                                  universe = if (annotation_background) character() else setdiff(universe[[component]], mapping$gene))
    result <- list(GO = NULL, KEGG = NULL, GO_full = NULL)
    actual_background <- list(GO = character(), KEGG = character())
    status <- "complete"
    if (!length(selected[[component]])) status <- "no_selected_genes"
    else if (!annotation_background && !length(universe[[component]])) status <- "no_tested_genes"
    else if (!annotation_background && !length(background)) status <- "no_mapped_background"
    else if (!length(genes)) status <- "no_mapped_genes"
    if (status %in% c("no_mapped_background", "no_mapped_genes")) {
      warning(component, ": ", status, "; see the returned unmapped identifiers.", call. = FALSE)
    }
    for (database in databases) {
      analysis <- NULL
      if (status == "complete") {
        if (database == "GO") {
          analysis <- clusterProfiler::enrichGO(gene = genes, universe = background,
            OrgDb = OrgDb, keyType = "ENTREZID", ont = ontology, pAdjustMethod = "BH",
            pvalueCutoff = p_cutoff[[database]], qvalueCutoff = q_cutoff[[database]],
            minGSSize = min_gene_set_size, maxGSSize = max_gene_set_size, readable = TRUE)
        } else if (is.null(kegg_cache)) {
          analysis <- clusterProfiler::enrichKEGG(gene = genes, universe = background,
            organism = organism, keyType = "ncbi-geneid", pAdjustMethod = "BH",
            pvalueCutoff = p_cutoff[[database]], qvalueCutoff = q_cutoff[[database]],
            minGSSize = min_gene_set_size, maxGSSize = max_gene_set_size)
        } else {
          analysis <- clusterProfiler::enricher(gene = genes, universe = background,
            TERM2GENE = kegg_cache$TERM2GENE, TERM2NAME = kegg_cache$TERM2NAME,
            pAdjustMethod = "BH", pvalueCutoff = p_cutoff[[database]],
            qvalueCutoff = q_cutoff[[database]], minGSSize = min_gene_set_size,
            maxGSSize = max_gene_set_size)
        }
      }
      if (database == "GO") {
        result["GO_full"] <- list(analysis)
        if (!is.null(analysis) && nrow(as.data.frame(analysis)) && component %in% simplify_go) {
          analysis <- clusterProfiler::simplify(analysis, cutoff = similarity_cutoff,
            by = "p.adjust", select_fun = min, measure = "Wang")
        }
      }
      result[database] <- list(analysis)
      if (!is.null(analysis)) actual_background[[database]] <- as.character(analysis@universe)
      n_terms <- nrow(as.data.frame(analysis))
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(component = component,
        database = database, n_selected = length(selected[[component]]), n_mapped = length(genes),
        n_universe = if (annotation_background) NA_integer_ else length(universe[[component]]),
        n_mapped_universe = if (annotation_background) NA_integer_ else length(background),
        n_annotation_background = if (is.null(analysis)) NA_integer_ else length(actual_background[[database]]),
        n_terms = n_terms,
        status = if (status == "complete" && !n_terms) "no_terms" else status)
    }
    out[[component]] <- c(result, list(genes = genes, universe = background, background = actual_background))
  }
  annotation_package <- OrgDb$packageName
  annotation_version <- if (length(annotation_package) == 1L &&
    requireNamespace(annotation_package, quietly = TRUE)) as.character(utils::packageVersion(annotation_package)) else NA_character_
  out$summary <- do.call(rbind, summary_rows)
  out$mapping <- mapping
  out$unmapped <- unmapped
  out$settings <- list(components = components, databases = databases, organism = organism,
    key_type = key_type, ontology = ontology, adjustment = "BH", pvalue_cutoff = p_cutoff,
    qvalue_cutoff = q_cutoff, min_gene_set_size = min_gene_set_size,
    max_gene_set_size = max_gene_set_size, background_source = background_source,
    simplify_go = simplify_go, similarity_cutoff = similarity_cutoff,
    similarity_measure = "Wang",
    annotation_package = annotation_package, annotation_version = annotation_version,
    annotation_metadata = AnnotationDbi::metadata(OrgDb),
    clusterProfiler_version = as.character(utils::packageVersion("clusterProfiler")),
    kegg_source = if ("KEGG" %in% databases) cache_source else NULL,
    kegg_cache_sha256 = cache_sha256, kegg_provenance = cache_provenance,
    analyzed_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  out
}

.enrichment_cutoffs <- function(value, databases, name) {
  if (!is.numeric(value) || !length(value) || any(!is.finite(value)) || any(value < 0 | value > 1)) {
    stop(name, " must contain values in [0, 1].", call. = FALSE)
  }
  if (length(value) == 1L && is.null(names(value))) return(stats::setNames(rep(value, length(databases)), databases))
  if (is.null(names(value)) || anyNA(names(value)) || anyDuplicated(names(value)) ||
      !setequal(names(value), databases)) {
    stop(name, " must be a single value or a vector named for each requested database.", call. = FALSE)
  }
  value[databases]
}

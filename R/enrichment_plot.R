#' Analyze and plot GO and KEGG enrichment
#'
#' Tests genes selected by [rank_genes()] for over-representation and draws
#' the leading terms for each component and database. By default, the
#' background contains all genes analyzed in the fit, before ranking filters.
#' Supplying a saved enrichment result or a previous plot redraws the results
#' without repeating the analysis.
#'
#' @param enrichment A result from [rank_genes()], a plot returned by this
#'   function, or its saved enrichment results (see Value).
#' @param components Components to display. `NULL` uses all components in
#'   the input, in their original order.
#' @param databases Databases to analyze or display: `"GO"`, `"KEGG"`, or both.
#'   `NULL` analyzes both for ranked genes, or uses the saved databases.
#' @param top_n Maximum number of terms per component and database. Terms
#'   are ordered by increasing BH-adjusted p-value, then raw p-value when
#'   available, and finally term ID. This changes the display, not the genes
#'   or tests used for enrichment.
#' @param family Plot font family.
#' @param OrgDb An `AnnotationDbi` organism database. For human genes
#'   (`organism = "hsa"`), `NULL` uses `org.Hs.eg.db`. Supply the matching
#'   database for other organisms. Used only when analyzing ranked genes.
#' @param organism KEGG organism code, such as `"hsa"` for human.
#' @param key_type Input identifier type supported by `OrgDb`, such as
#'   `"SYMBOL"` or `"ENSEMBL"`.
#' @param ontology GO ontology: `"BP"`, `"MF"`, `"CC"`, or `"ALL"`.
#' @param universe Background gene identifiers. `NULL` uses all genes saved
#'   in `ranked$universe`. Supply a character vector, a component-named list,
#'   or `"annotation"` to use each database's full annotation background.
#'   An explicit gene background must include the selected genes.
#' @param simplify_go Reduce redundant GO terms using Wang semantic similarity,
#'   retaining the term with the smallest adjusted p-value. Use `TRUE` for all
#'   components, `FALSE` for none, or component names for selected components.
#'   KEGG terms are not simplified.
#' @param similarity_cutoff GO semantic-similarity cutoff, between 0 and 1.
#' @param pvalue_cutoff Cutoff for raw and BH-adjusted enrichment p-values.
#'   Use one value or a vector named for the requested databases.
#' @param qvalue_cutoff Additional cutoff on clusterProfiler's estimated
#'   q-values. Use one value or a database-named vector, such as
#'   `c(GO = 0.1, KEGG = 0.2)`.
#' @param min_gene_set_size,max_gene_set_size Minimum and maximum gene-set
#'   size within the annotation background.
#' @param kegg_cache Optional list, or path to an RDS list, with two-column
#'   `TERM2GENE` and `TERM2NAME` tables. They contain term/Entrez-ID and
#'   term/description pairs, respectively. Entrez IDs have no organism prefix.
#'   `NULL` retrieves KEGG annotations through clusterProfiler.
#' @return A ggplot. Its `data` contains the displayed terms, original
#'   `p.adjust` and `GeneRatio` values, and numeric `gene_ratio` and `log_p`.
#'   `attr(plot, "enrichment")` contains all enrichment results, before the
#'   `top_n` display limit: component-specific `GO` and `KEGG` enrichResult
#'   objects, unsimplified `GO_full`, identifier mapping, actual backgrounds,
#'   and analysis settings. This list can be supplied to `plot_enrichment()`.
#' @details
#' The horizontal axis shows gene ratio, point size shows gene count, and
#' color shows `-log10(BH-adjusted p-value)`. Empty panels are labeled
#' `No enriched terms`. Zero adjusted p-values use the smallest positive
#' normalized double when computing colors; the stored p-values are unchanged.
#'
#' Analysis requires `clusterProfiler` and `AnnotationDbi`; default human
#' annotations also require `org.Hs.eg.db`. Gene identifiers are mapped to
#' Entrez IDs, retaining all one-to-many mappings and removing duplicates.
#' GO and KEGG use the BH procedure separately for each component. The actual
#' background is restricted to identifiers represented in each database.
#' Analysis settings apply to ranked inputs; saved results are plotted as-is.
#' @examples
#' \dontrun{
#' ranked_full <- rank_genes(fit_full, top_n = 100, min_r2 = 0.03)
#' p <- plot_enrichment(ranked_full, top_n = 3)
#' enrichment <- attr(p, "enrichment")
#' plot_enrichment(enrichment, top_n = 5)
#' }
#' @export
plot_enrichment <- function(enrichment, components = NULL, databases = NULL,
                            top_n = 5, family = "sans", OrgDb = NULL,
                            organism = "hsa", key_type = "SYMBOL", ontology = "BP",
                            universe = NULL, simplify_go = TRUE, similarity_cutoff = 0.7,
                            pvalue_cutoff = 0.05, qvalue_cutoff = 0.2,
                            min_gene_set_size = 10, max_gene_set_size = 500,
                            kegg_cache = NULL) {
  .plot_dependencies("ggplot2")
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
      top_n < 1 || top_n != floor(top_n)) {
    stop("top_n must be a positive whole number.", call. = FALSE)
  }
  saved <- attr(enrichment, "enrichment", exact = TRUE)
  if (!is.null(saved)) enrichment <- saved
  if (!is.list(enrichment) || !is.list(enrichment$settings)) {
    stop("Supply ranked genes, saved enrichment results, or an enrichment plot.", call. = FALSE)
  }
  labels <- c(slide_wide_he = "Slide-wide H&E effect",
              location_specific_he = "Location-specific H&E effect",
              spatial_beyond_he = "Spatial variation beyond H&E")
  if (is.null(components)) components <- enrichment$settings$components
  if (!is.character(components) || !length(components) || anyNA(components) ||
      anyDuplicated(components) || !all(components %in% names(labels)) ||
      !all(components %in% enrichment$settings$components)) {
    stop("Choose components present in the input.", call. = FALSE)
  }
  ranked_input <- is.null(enrichment$settings$databases)
  if (is.null(databases)) {
    databases <- if (ranked_input) c("GO", "KEGG") else enrichment$settings$databases
  }
  if (!is.character(databases) || !length(databases) || anyNA(databases) ||
      anyDuplicated(databases) || !all(databases %in% c("GO", "KEGG")) ||
      (!ranked_input && !all(databases %in% enrichment$settings$databases))) {
    stop("Choose GO and/or KEGG from the input.", call. = FALSE)
  }
  if (ranked_input) {
    if (is.null(OrgDb)) {
      if (!identical(organism, "hsa")) {
        stop("Supply OrgDb for the selected organism.", call. = FALSE)
      }
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Install org.Hs.eg.db for human enrichment, or supply OrgDb.", call. = FALSE)
      }
      OrgDb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    }
    enrichment <- .enrich_ranked_genes(enrichment, OrgDb = OrgDb,
      organism = organism, key_type = key_type, ontology = ontology, databases = databases,
      pvalue_cutoff = pvalue_cutoff, qvalue_cutoff = qvalue_cutoff,
      min_gene_set_size = min_gene_set_size, max_gene_set_size = max_gene_set_size,
      universe = universe, kegg_cache = kegg_cache,
      simplify_go = simplify_go, similarity_cutoff = similarity_cutoff)
  } else {
    analysis_arguments <- c("OrgDb", "organism", "key_type", "ontology", "universe",
      "simplify_go", "similarity_cutoff", "pvalue_cutoff", "qvalue_cutoff",
      "min_gene_set_size", "max_gene_set_size", "kegg_cache")
    if (any(analysis_arguments %in% names(match.call()))) {
      stop("Supply ranked genes to change enrichment analysis settings.", call. = FALSE)
    }
  }

  terms <- list()
  empty <- character()
  panels <- as.vector(outer(databases, labels[components], paste, sep = "\n"))
  for (database in databases) {
    for (component in components) {
      panel <- paste(database, labels[[component]], sep = "\n")
      if (!is.list(enrichment[[component]]) || !database %in% names(enrichment[[component]])) {
        stop("The enrichment result is missing ", component, " / ", database, ".", call. = FALSE)
      }
      d <- as.data.frame(enrichment[[component]][[database]])
      if (!nrow(d)) {
        empty <- c(empty, panel)
        next
      }
      columns <- c("ID", "Description", "GeneRatio", "Count", "p.adjust")
      if (!all(columns %in% names(d))) stop("Enrichment terms are missing plot columns.", call. = FALSE)
      if (!"pvalue" %in% names(d)) d$pvalue <- NA_real_
      d <- d[, c(columns, "pvalue"), drop = FALSE]
      d$ID <- as.character(d$ID)
      d$Description <- as.character(d$Description)
      if (anyNA(d$ID) || any(!nzchar(d$ID)) || anyDuplicated(d$ID) ||
          anyNA(d$Description) || !is.numeric(d$Count) ||
          any(!is.finite(d$Count) | d$Count < 0) || !is.numeric(d$p.adjust) ||
          any(!is.finite(d$p.adjust) | d$p.adjust < 0 | d$p.adjust > 1)) {
        stop("Enrichment terms need unique IDs, descriptions, counts and adjusted p-values.", call. = FALSE)
      }
      ratios <- strsplit(as.character(d$GeneRatio), "/", fixed = TRUE)
      d$gene_ratio <- vapply(ratios, function(value) {
        numbers <- suppressWarnings(as.numeric(value))
        if (length(numbers) != 2L || any(!is.finite(numbers)) ||
            numbers[1L] < 0 || numbers[2L] <= 0 || numbers[1L] > numbers[2L]) return(NA_real_)
        numbers[1L] / numbers[2L]
      }, numeric(1))
      if (anyNA(d$gene_ratio)) stop("GeneRatio must contain ratios such as '5/100'.", call. = FALSE)
      raw_p <- d$pvalue
      if (!is.numeric(raw_p) || any(!is.na(raw_p) & (!is.finite(raw_p) | raw_p < 0 | raw_p > 1))) {
        stop("Available raw enrichment p-values must be finite numbers in [0, 1].", call. = FALSE)
      }
      d <- utils::head(d[order(d$p.adjust, raw_p, d$ID), , drop = FALSE], top_n)
      d$component <- component
      d$database <- database
      d$panel <- panel
      terms[[length(terms) + 1L]] <- d
    }
  }
  if (length(terms)) {
    d <- do.call(rbind, terms)
    rownames(d) <- NULL
  } else {
    d <- data.frame(ID = character(), Description = character(), GeneRatio = character(),
      Count = numeric(), p.adjust = numeric(), pvalue = numeric(), gene_ratio = numeric(),
      component = character(), database = character(), panel = character())
  }
  d$log_p <- -log10(pmax(d$p.adjust, .Machine$double.xmin))
  # Give each term its own axis entry, even when descriptions repeat across panels.
  keys <- as.character(seq_len(nrow(d)))
  term_labels <- stats::setNames(vapply(d$Description, function(value) {
    paste(strwrap(value, width = 38), collapse = "\n")
  }, character(1)), keys)
  d$term <- factor(keys, levels = rev(keys))
  d$panel <- factor(d$panel, levels = panels)
  database_labels <- c(GO = "GO", KEGG = "KEGG Pathways")
  ontology_labels <- c(BP = "GO Biological Process", MF = "GO Molecular Function",
                       CC = "GO Cellular Component", ALL = "GO All Ontologies")
  ontology <- enrichment$settings$ontology
  if (length(ontology) == 1L && ontology %in% names(ontology_labels)) {
    database_labels["GO"] <- ontology_labels[[ontology]]
  }
  panel_labels <- stats::setNames(as.vector(outer(database_labels[databases],
    labels[components], paste, sep = "\n")), panels)
  for (panel in names(panel_labels)) {
    title_parts <- strsplit(panel_labels[[panel]], "\n", fixed = TRUE)[[1L]]
    panel_labels[[panel]] <- paste(unlist(lapply(title_parts, strwrap, width = 22)), collapse = "\n")
  }
  empty_panels <- data.frame(panel = factor(empty, levels = panels),
    term = rep("empty", length(empty)),
    gene_ratio = rep(if (nrow(d)) mean(range(d$gene_ratio)) else 0.5, length(empty)))
  term_labels <- c(term_labels, empty = "")

  # Use the manuscript's log-spaced palette; keep legend ticks in log units.
  cutoff <- -log10(.05)
  maximum <- max(c(3, d$log_p))
  color_values <- c(0, .25 * cutoff, .6 * cutoff, cutoff,
                    cutoff + .45 * (maximum - cutoff), maximum) / maximum
  color_breaks <- sort(unique(c(0, cutoff, 3, maximum)))
  color_labels <- vapply(color_breaks, function(value) {
    if (abs(value - cutoff) < 1e-8) "1.30 (p = 0.05)"
    else if (abs(value - 3) < 1e-8) "3 (p = 0.001)"
    else format(signif(value, 3), trim = TRUE)
  }, character(1))
  count_max <- max(c(1, d$Count))
  count_breaks <- unique(round(pretty(c(0, count_max), n = 3)))
  count_breaks <- count_breaks[count_breaks > 0 & count_breaks <= count_max]

  p <- ggplot2::ggplot(d, ggplot2::aes(gene_ratio, term, size = Count, colour = log_p)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(data = empty_panels,
      ggplot2::aes(gene_ratio, term, label = "No enriched terms"),
      inherit.aes = FALSE, colour = "grey40", size = 3, family = family) +
    ggplot2::facet_wrap(~ panel, ncol = length(databases), scales = "free", drop = FALSE,
      labeller = ggplot2::as_labeller(panel_labels)) +
    ggplot2::scale_y_discrete(labels = term_labels, expand = ggplot2::expansion(add = .65)) +
    ggplot2::scale_x_continuous(n.breaks = 3, expand = ggplot2::expansion(mult = .2),
      limits = if (!nrow(d)) c(0, 1) else NULL) +
    ggplot2::scale_colour_gradientn(
      colours = c("grey90", "grey70", "skyblue3", "seagreen3", "gold2", "orange2", "purple3"),
      values = color_values, limits = c(0, maximum), breaks = color_breaks, labels = color_labels) +
    ggplot2::scale_size_continuous(range = c(0, 6), limits = c(0, count_max), breaks = count_breaks) +
    ggplot2::labs(x = "Gene Ratio", y = NULL, size = "Gene Count",
      colour = "-log10(BH-adjusted p-value)") +
    ggplot2::theme_minimal(base_family = family) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      strip.background = ggplot2::element_rect(fill = "grey90", colour = NA),
      strip.text = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 9),
      panel.spacing = grid::unit(4, "mm")) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(order = 1, title.position = "top",
      barheight = grid::unit(35, "mm"), barwidth = grid::unit(4, "mm")),
      size = ggplot2::guide_legend(order = 2, title.position = "top",
        override.aes = list(colour = "black")))
  attr(p, "enrichment") <- enrichment
  p
}

utils::globalVariables(c("gene_ratio", "term", "Count", "log_p"))

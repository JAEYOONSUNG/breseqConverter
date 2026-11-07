# Compact views are derived from canonical calls; the companion wide workbook
# remains the source of all original fields and unassigned evidence rows.
.breseq_pub_reference <- function(reference_file = NULL, gff_file = NULL) {
  labels <- character(); features <- list()
  if (!is.null(reference_file)) {
    if (!file.exists(reference_file)) stop("reference_file does not exist.")
    records <- .breseq_gb_read_records(reference_file)
    for (record in records) {
      label <- record$seq_id
      for (feature in record$features) {
        q <- feature$qualifiers
        get <- function(key) {
          at <- which(vapply(q, .breseq_gb_qualifier_key, "") == key)
          if (!length(at)) return("")
          text <- paste(trimws(q[[at[1]]]), collapse = " ")
          sub('"$', "", sub('^[^=]*="?', "", text))
        }
        if (feature$key == "source") {
          plasmid <- get("plasmid"); chromosome <- get("chromosome")
          if (nzchar(plasmid)) label <- plasmid else if (nzchar(chromosome)) label <- paste("Chromosome", chromosome)
        }
        locus <- get("locus_tag"); gene <- get("gene"); product <- get("product")
        if (nzchar(locus) || nzchar(gene)) features[[length(features) + 1L]] <-
          data.frame(seq_id = record$seq_id, locus = locus, gene = gene, product = product,
                     type = feature$key, stringsAsFactors = FALSE)
      }
      if (identical(label, record$seq_id) && any(grepl("chromosome", record$header, ignore.case = TRUE))) label <- "Chromosome"
      labels[record$seq_id] <- label
    }
  }
  if (!is.null(gff_file)) {
    if (!file.exists(gff_file)) stop("gff_file does not exist.")
    gff <- suppressMessages(parse_gff_file(gff_file))
    if (!is.null(gff)) for (i in seq_len(nrow(gff))) features[[length(features) + 1L]] <-
      data.frame(seq_id = .breseq_seq(gff$seqid[i]), locus = gff$locus_tag[i], gene = gff$gene[i],
        product = utils::URLdecode(gff$product[i]), type = gff$type[i], stringsAsFactors = FALSE)
  }
  data <- if (length(features)) unique(do.call(rbind, features)) else
    data.frame(seq_id = character(), locus = character(), gene = character(), product = character(), type = character())
  for (column in names(data)) data[[column]][is.na(data[[column]])] <- ""
  # CDS/RNA product qualifiers take precedence over product-less gene features.
  data <- data[order(!nzchar(data$product), data$type == "gene"), , drop = FALSE]
  if (anyDuplicated(unname(labels))) labels[duplicated(labels) | duplicated(labels, fromLast = TRUE)] <-
    names(labels)[duplicated(labels) | duplicated(labels, fromLast = TRUE)]
  list(labels = labels, features = data)
}

.breseq_pub_coalesce <- function(data, columns, fields, rows, samples, labels) {
  values <- character(); stages <- character()
  for (sample in samples) {
    text <- rep(NA_character_, length(rows))
    for (field in fields) {
      candidate <- .breseq_field(data, columns, field, sample)[rows]
      take <- !.breseq_present(text) & .breseq_present(candidate)
      text[take] <- candidate[take]
    }
    text <- unique(text[.breseq_present(text)])
    if (length(text)) { values <- c(values, text); stages <- c(stages, rep(unname(labels[sample]), length(text))) }
  }
  distinct <- unique(values)
  if (!length(distinct)) return("")
  if (length(distinct) == 1L) return(distinct)
  paste(vapply(distinct, function(value) paste0(paste(unique(stages[values == value]), collapse = "/"), ": ", value), ""), collapse = "\n")
}

.breseq_publication_tables <- function(merged_results, sample_order = NULL, sample_labels = NULL,
                                       reference_file = NULL, gff_file = NULL, plot_evidence = NULL,
                                       annotation_overrides = NULL) {
  if (!is.list(merged_results) || is.null(names(merged_results))) stop("merged_results must be a named list of canonical data frames.")
  kinds <- vapply(names(merged_results), .breseq_table_kind, "")
  if (anyDuplicated(kinds[kinds != "other"])) stop("Only one source table of each mutation/evidence kind is supported.")
  inputs <- merged_results[kinds != "other"]; kinds <- kinds[kinds != "other"]
  if (!length(inputs)) stop("No mutation or evidence tables were supplied.")
  if (!all(vapply(inputs, is.data.frame, logical(1)))) stop("Every source table must be a data frame.")
  detected <- unique(unlist(lapply(inputs, function(d) attr(.breseq_stage_columns(d), "sample_order")), use.names = FALSE))
  samples <- if (is.null(sample_order)) detected else sample_order
  if (!length(samples)) stop("No stage columns were found in the canonical tables.")
  columns <- lapply(inputs, .breseq_stage_columns, sample_order = samples)
  labels <- .breseq_validate_labels(sample_labels, samples)
  reference <- .breseq_pub_reference(reference_file, gff_file)
  replicon <- function(id) {
    mapped <- unname(reference$labels[id]); mapped[is.na(mapped)] <- id[is.na(mapped)]; mapped
  }
  seq_order <- unique(c(names(reference$labels), unlist(lapply(seq_along(inputs), function(i) .breseq_wide_seq(inputs[[i]], columns[[i]])))))
  baseline <- samples[toupper(samples) == "WT"]
  baseline <- baseline[vapply(baseline, function(sample) !any(vapply(seq_along(inputs), function(i) {
    cols <- columns[[i]]$index[columns[[i]]$sample == sample]
    any(vapply(inputs[[i]][cols], function(x) any(.breseq_present(x)), logical(1)))
  }, logical(1))), logical(1))]
  absent <- function(sample) if (sample %in% baseline) "Ref." else "\u2014"
  states_for <- function(d, name, mutation = NULL) .breseq_stage_status(d, name, samples, mutation)$states
  mutation_at <- which(kinds == "mutation")
  mutation <- if (length(mutation_at)) inputs[[mutation_at]] else NULL
  sheets <- list(); audit <- list(); mc_mapping <- list()
  colors <- stats::setNames(rep(c("#DEEBF7", "#DAE1EF", "#E4DFF0", "#EBE2F3"), length.out = length(samples)), samples)
  experimental <- samples[!samples %in% baseline]
  colors[experimental] <- rep(c("#DEEBF7", "#DAE1EF", "#E4DFF0", "#EBE2F3"), length.out = length(experimental))
  colors[baseline] <- "#F2F7FD"
  number <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  interval <- function(x) {
    b <- .breseq_bounds(x)
    if (anyNA(b)) return(as.character(x))
    paste(number(unique(b)), collapse = "\u2013")
  }
  make_sheet <- function(data, states, stage_columns, title, notes, widths, groups) {
    list(data = data, states = states, stage_columns = stage_columns, title = title, notes = notes,
         widths = widths, groups = groups)
  }
  for (table_index in seq_along(inputs)) {
    d <- inputs[[table_index]]; cols <- columns[[table_index]]; kind <- kinds[table_index]
    ids <- .breseq_wide_seq(d, cols)
    if (kind == "mutation") {
      state <- states_for(d, "Mutations")
      position <- suppressWarnings(as.numeric(.breseq_coord(d$position)))
      if (anyNA(position)) stop("Mutation positions must be numeric for the compact view.")
      order <- order(match(ids, seq_order), position)
      out <- data.frame(Replicon = replicon(ids), `Position (bp)` = position,
        `Gene(s)` = vapply(seq_len(nrow(d)), function(i) .breseq_pub_coalesce(d, cols, "gene", i, samples, labels), ""),
        `Product / region annotation` = vapply(seq_len(nrow(d)), function(i) .breseq_pub_coalesce(d, cols, c("description", "product"), i, samples, labels), ""),
        `Variant annotation` = vapply(seq_len(nrow(d)), function(i) .breseq_pub_coalesce(d, cols, "annotation", i, samples, labels), ""), check.names = FALSE)
      for (sample in samples) {
        value <- .breseq_field(d, cols, "mutation", sample)
        value[! .breseq_present(value)] <- absent(sample)
        out[[labels[[sample]]]] <- value
      }
      out <- out[order, , drop = FALSE]; state <- state[order, , drop = FALSE]
      rownames(out) <- NULL
      sheets[["Mutations"]] <- make_sheet(out, state, 5L + seq_along(samples),
        "Sequence variants by construction stage",
        c(paste0(sum(state != "not_reported"), " source calls at ", nrow(out), " loci. Ref. denotes the reference-assembly baseline."),
          "Pale stage colors mark every reported call; bold marks new, changed or reappearing calls.",
          "An em dash means no call was reported; it does not establish sequence restoration."),
        c(19, 15, 36, 78, 48, rep(22, length(samples))),
        list(list(label = "Genomic locus", cols = 1:2), list(label = "Annotation", cols = 3:5), list(label = "Reported calls by strain", cols = 5L + seq_along(samples))))
      audit$mutation_source_rows <- order
      audit$mutation_keys <- data.frame(output_row = seq_along(order), seq_id = ids[order],
                                        position = .breseq_coord(d$position[order]), stringsAsFactors = FALSE)
    } else if (kind == "missing") {
      state <- states_for(d, "Missing", mutation)
      records <- list()
      for (sample in samples) for (i in which(state[, sample] != "not_reported")) {
        end <- .breseq_field(d, cols, "end", sample)[i]
        a <- .breseq_bounds(d$start[i]); b <- .breseq_bounds(end)
        records[[length(records) + 1L]] <- data.frame(row = i, sample = sample, seq_id = ids[i],
          start = as.character(d$start[i]), end = end, core_start = a[2], core_end = b[1],
          gene = .breseq_field(d, cols, "gene", sample)[i],
          description = .breseq_field(d, cols, "description", sample)[i], state = state[i, sample],
          stringsAsFactors = FALSE)
      }
      r <- if (length(records)) do.call(rbind, records) else data.frame()
      groups <- list()
      if (nrow(r)) {
        parent <- seq_len(nrow(r))
        root <- function(i) { while (parent[i] != i) i <- parent[i]; i }
        common_locus <- function(a, b) {
          if (!.breseq_present(a) || !.breseq_present(b)) return(TRUE)
          loci <- function(x) regmatches(x, gregexpr("[A-Za-z0-9]+_RS[0-9]+", x))[[1]]
          left <- loci(a); right <- loci(b)
          if (length(left) && length(right)) length(intersect(left, right)) > 0L else identical(a, b)
        }
        if (nrow(r) > 1L) for (a in seq_len(nrow(r) - 1L)) for (b in seq.int(a + 1L, nrow(r))) {
          if (r$seq_id[a] == r$seq_id[b] && common_locus(r$gene[a], r$gene[b]) &&
              max(r$core_start[a], r$core_start[b]) <= min(r$core_end[a], r$core_end[b])) parent[root(b)] <- root(a)
        }
        groups <- split(seq_len(nrow(r)), vapply(seq_len(nrow(r)), root, 1L))
        groups <- groups[order(vapply(groups, function(g) match(r$seq_id[g[1]], seq_order), 1L),
                               vapply(groups, function(g) min(r$core_start[g]), 1))]
      }
      omitted <- vapply(groups, function(g) all(r$state[g] == "prior_del_region"), logical(1))
      audit$mc_records <- r; audit$mc_families <- groups; audit$omitted_mc_families <- groups[omitted]
      shown <- groups[!omitted]
      audit$mc_keys <- if (length(shown)) do.call(rbind, lapply(seq_along(shown), function(j) {
        g <- shown[[j]]
        unique(data.frame(output_row = j, seq_id = r$seq_id[g], position = .breseq_coord(r$start[g]), stringsAsFactors = FALSE))
      })) else data.frame(output_row = integer(), seq_id = character(), position = character())
      out <- data.frame(Replicon = character(), `Reported gene region` = character(),
                        `Gene / region annotation` = character(), check.names = FALSE)
      state_out <- matrix("not_reported", length(shown), length(samples), dimnames = list(NULL, samples))
      for (j in seq_along(shown)) {
        g <- shown[[j]]; source_rows <- unique(r$row[g]); seqid <- r$seq_id[g[1]]
        gene <- .breseq_pub_coalesce(d, cols, "gene", source_rows, samples, labels)
        raw_description <- unique(r$description[g][.breseq_present(r$description[g])])
        tokens <- unique(unlist(strsplit(raw_description, ",", fixed = TRUE), use.names = FALSE))
        annotated <- character()
        for (token in tokens) {
          clean <- trimws(gsub("\\[|\\]", "", token))
          candidates <- reference$features[reference$features$seq_id == seqid &
            (reference$features$locus == clean | reference$features$gene == clean), , drop = FALSE]
          # A single-gene MC record may already contain a product instead of a
          # gene-token list; resolve its gene field for the annotation audit.
          if (!nrow(candidates) && length(tokens) == 1L) {
            reported_genes <- unique(r$gene[g][.breseq_present(r$gene[g])])
            if (length(reported_genes) == 1L) {
              clean_gene <- gsub("\\[|\\]", "", reported_genes)
              candidates <- reference$features[reference$features$seq_id == seqid &
                (reference$features$locus == clean_gene | reference$features$gene == clean_gene), , drop = FALSE]
            }
          }
          if (length(unique(candidates$locus[nzchar(candidates$locus)])) > 1L) candidates <- candidates[0, , drop = FALSE]
          if (nrow(candidates) && any(nzchar(candidates$product))) {
            candidates <- candidates[nzchar(candidates$product), , drop = FALSE]
            product <- candidates$product[1]; locus <- candidates$locus[1]
            annotated <- c(annotated, product)
            mc_mapping[[length(mc_mapping) + 1L]] <- data.frame(family = j, seq_id = seqid,
              source_token = token, locus = locus, product = product, partial = grepl("\\[", token), stringsAsFactors = FALSE)
          } else annotated <- c(annotated, token)
        }
        products <- unique(annotated)
        annotation <- paste(vapply(products, function(p) {
          count <- sum(annotated == p)
          if (count > 1L) paste0(p, " (", count, " loci)") else p
        }, ""), collapse = "; ")
        row <- list(replicon(seqid), gene, annotation)
        for (sample in samples) {
          at <- g[r$sample[g] == sample]
          if (!length(at)) { row[[length(row) + 1L]] <- absent(sample); next }
          texts <- vapply(at, function(k) {
            size <- .breseq_field(d, cols, "size", sample)[r$row[k]]
            if (!.breseq_present(size)) size <- paste(r$core_end[k] - r$core_start[k] + 1,
              max(.breseq_bounds(r$end[k])) - min(.breseq_bounds(r$start[k])) + 1, sep = "\u2013")
            paste0("Start: ", interval(r$start[k]), "\nEnd: ", interval(r$end[k]), "\nLength: ", interval(size), " bp")
          }, "")
          row[[length(row) + 1L]] <- paste(unique(texts), collapse = "\n\n")
          priority <- c("new", "reappeared", "updated", "retained", "prior_del_region")
          state_out[j, sample] <- priority[match(TRUE, priority %in% r$state[at])]
        }
        names(row) <- c(names(out)[1:3], unname(labels))
        out[j, names(row)] <- row
      }
      if (!length(shown)) for (sample in samples) out[[labels[[sample]]]] <- character()
      retained_records <- if (length(shown)) length(unique(r$row[unlist(shown)])) else 0L
      sheets[["Missing coverage"]] <- make_sheet(out, state_out, 3L + seq_along(samples),
        "Unassigned missing-coverage evidence by construction stage",
        c(paste0(nrow(out), " region families summarize ", retained_records, " source MC records. ", sum(omitted), " prior-deletion-only families are omitted here; all ", nrow(d), " raw records remain in the companion wide workbook."),
          "Pale colors mark reported MC evidence. Bold marks first/reappearing regional evidence; boundary revisions and evidence within a prior deletion are regular weight.",
          "MC is unassigned evidence, not a finalized deletion. Products come from reference annotation when available; uncertain boundaries do not establish complete gene loss."),
        c(19, 39, 88, rep(37, length(samples))),
        list(list(label = "Reported region and annotation", cols = 1:3), list(label = "Reported missing-coverage evidence", cols = 3L + seq_along(samples))))
    } else if (kind == "junction") {
      keys <- if ("comparison_key" %in% names(d)) as.character(d$comparison_key) else attr(d, "breseq_keys")
      paired <- !is.null(keys)
      if (paired) {
        state <- states_for(d, "Junction evidence")
        groups <- split(seq_len(nrow(d)), sub("\\|\\|end[12]$", "", keys))
      } else {
        # Without oriented pair keys, retain every endpoint row. Never infer
        # pairs from row adjacency or a shared endpoint after an Excel reload.
        groups <- as.list(seq_len(nrow(d)))
        state <- matrix("not_reported", nrow(d), length(samples), dimnames = list(NULL, samples))
        for (i in seq_len(nrow(d))) {
          previous <- FALSE; seen <- FALSE
          for (sample in samples) {
            presence_field <- if (any(cols$field == "seq id" & cols$sample == sample)) "seq id" else "reads (cov)"
            present <- .breseq_present(.breseq_field(d, cols, presence_field, sample)[i])
            if (present) state[i, sample] <- if (previous) "retained" else if (seen) "reappeared" else "new"
            seen <- seen || present; previous <- present
          }
        }
      }
      coordinate <- suppressWarnings(as.numeric(gsub("[^0-9]", "", .breseq_coord(d$position))))
      if (anyNA(coordinate)) stop("Junction endpoints must have a reference coordinate.")
      groups <- lapply(groups, function(g) g[order(match(ids[g], seq_order), coordinate[g], .breseq_coord(d$position[g]))])
      groups <- groups[order(vapply(groups, function(g) match(ids[g[1]], seq_order), 1L),
                             vapply(groups, function(g) coordinate[g[1]], 1))]
      endpoint <- function(i) paste0(replicon(ids[i]), ": ", .breseq_text(d$position[i]))
      out <- data.frame(); state_out <- matrix("not_reported", length(groups), length(samples), dimnames = list(NULL, samples))
      for (j in seq_along(groups)) {
        g <- groups[[j]]
        context <- unique(vapply(g, function(i) .breseq_pub_coalesce(d, cols, "gene", i, samples, labels), ""))
        row <- if (paired) list(sprintf("J%02d", j), endpoint(g[1]), endpoint(g[2]), paste(context[nzchar(context)], collapse = "\n")) else
          list(sprintf("E%02d", j), endpoint(g[1]), paste(context[nzchar(context)], collapse = "\n"))
        for (sample in samples) row[[length(row) + 1L]] <- if (state[g[1], sample] == "not_reported") absent(sample) else "Detected"
        names(row) <- c(if (paired) c("Event", "Oriented endpoint 1", "Oriented endpoint 2", "Gene context") else
          c("Endpoint", "Oriented endpoint", "Gene context"), unname(labels))
        out[j, names(row)] <- row
        state_out[j, ] <- state[g[1], ]
      }
      fixed <- if (paired) 4L else 3L
      if (!length(groups)) {
        headers <- c(if (paired) c("Event", "Oriented endpoint 1", "Oriented endpoint 2", "Gene context") else
          c("Endpoint", "Oriented endpoint", "Gene context"), unname(labels))
        out <- as.data.frame(stats::setNames(rep(list(character()), length(headers)), headers), check.names = FALSE)
      }
      audit$junction_mode <- if (paired) "paired_events" else "unpaired_endpoints"
      audit$junction_source_rows <- groups
      sheets[["Junction evidence"]] <- make_sheet(out, state_out, fixed + seq_along(samples),
        "Unassigned junction evidence by construction stage",
        c(if (paired) paste0(nrow(out), " junction events retain both oriented endpoints from ", nrow(d), " source endpoint rows.") else
          paste0(nrow(out), " endpoint rows are shown separately because paired comparison keys are unavailable; no junction pairs were inferred."),
          "Pale colors mark reported evidence. Bold marks new or reappearing events (or endpoints in unpaired mode).",
          "Detected denotes unassigned junction evidence, not a finalized structural mutation."),
        c(12, rep(34, fixed - 2L), 56, rep(22, length(samples))),
        list(list(label = if (paired) "Oriented junction" else "Unpaired endpoint", cols = seq_len(fixed)), list(label = "Reported junction evidence", cols = fixed + seq_along(samples))))
    }
  }
  audit$mc_annotation_map <- if (length(mc_mapping)) do.call(rbind, mc_mapping) else data.frame()
  if (!is.null(plot_evidence) && !is.data.frame(plot_evidence))
    stop("plot_evidence must be an audited display-evidence data frame.")
  if (!is.null(plot_evidence) && "appended" %in% names(plot_evidence))
    plot_evidence <- plot_evidence[!is.na(plot_evidence$appended) & plot_evidence$appended, , drop = FALSE]
  if (!is.null(plot_evidence) && nrow(plot_evidence)) {
    required <- c("seq_id", "start", "end", "size", "evidence_origin", "source_start", "source_end", "criterion", "depth_ratio", "raw_final_DEL_call")
    if (!is.data.frame(plot_evidence) || !all(required %in% names(plot_evidence))) stop("plot_evidence must contain the audited display-evidence fields.")
    p <- plot_evidence
    sample <- if ("sample_id" %in% names(p)) as.character(p$sample_id) else as.character(p$sample)
    if (anyNA(sample) || any(!sample %in% samples)) stop("plot_evidence contains unknown sample IDs.")
    order <- order(match(sample, samples))
    p <- p[order, , drop = FALSE]; sample <- sample[order]
    gene <- label <- character(nrow(p))
    missing_at <- which(kinds == "missing")
    for (i in seq_len(nrow(p))) {
      if (p$evidence_origin[i] == "COVERAGE_RATIO") gene[i] <- replicon(.breseq_seq(p$seq_id[i])) else if (length(missing_at)) {
        d <- inputs[[missing_at]]; cols <- columns[[missing_at]]
        at <- which(.breseq_wide_seq(d, cols) == .breseq_seq(p$seq_id[i]) &
          .breseq_coord(d$start) == .breseq_coord(p$source_start[i]) &
          .breseq_coord(.breseq_field(d, cols, "end", sample[i])) == .breseq_coord(p$source_end[i]))
        if (length(at)) gene[i] <- .breseq_field(d, cols, "gene", sample[i])[at[1]]
      }
    }
    displayed <- paste0("\u0394", number(p$size), " bp")
    label <- trimws(paste(displayed, gene))
    if ("figure_label" %in% names(p)) label[.breseq_present(p$figure_label)] <- p$figure_label[.breseq_present(p$figure_label)]
    out <- data.frame(Strain = unname(labels[sample]), Replicon = as.character(p$seq_id),
      `Start (bp)` = as.numeric(p$start), `End (bp)` = as.numeric(p$end), `Displayed mutation` = displayed,
      `Gene / region` = gene, `Figure label` = label,
      `Evidence source` = ifelse(p$evidence_origin == "MC", "Unassigned missing coverage", "Replicon coverage ratio"),
      `Source interval (bp)` = paste(p$source_start, p$source_end, sep = " - "),
      `Evidence rule` = as.character(p$criterion), `Depth ratio` = as.numeric(p$depth_ratio),
      `Raw DEL call` = ifelse(toupper(as.character(p$raw_final_DEL_call)) %in% c("TRUE", "YES"), "Yes", "No"), check.names = FALSE)
    sheets[["Plot-derived evidence"]] <- make_sheet(out, matrix(character(), nrow(out), 0), integer(),
      "Additional DEL evidence displayed in the figure",
      c(paste0(nrow(out), " explicitly supplied display marks are listed separately from finalized mutation calls."),
        "These rows reproduce audited plotting evidence; no inferred DEL is added to the mutation table.",
        "MC boundaries and coverage-based replicon loss do not establish a reconstructed sequence or exact deletion breakpoints."),
      c(19, 20, 15, 15, 22, 40, 48, 30, 35, 45, 17, 17),
      list(list(label = "Figure mark coordinates", cols = 1:4), list(label = "Displayed change", cols = 5:7), list(label = "Evidence supporting the display", cols = 8:12)))
    audit$plot_evidence <- plot_evidence
    audit$plot_keys <- data.frame(output_row = seq_len(nrow(out)),
      seq_id = .breseq_seq(p$seq_id), position = .breseq_coord(p$start),
      stringsAsFactors = FALSE)
  }
  sheets <- sheets[intersect(c("Mutations", "Missing coverage", "Junction evidence", "Plot-derived evidence"), names(sheets))]
  if (!is.null(annotation_overrides)) {
    required <- c("sheet", "seq_id", "position", "field", "value")
    if (!is.data.frame(annotation_overrides) || !all(required %in% names(annotation_overrides))) {
      stop("annotation_overrides must contain sheet, seq_id, position, field and value columns.")
    }
    overrides <- annotation_overrides[required]
    for (field in required) overrides[[field]] <- as.character(overrides[[field]])
    if (anyNA(overrides) || any(!nzchar(overrides$sheet)) || any(!nzchar(overrides$field))) stop("Annotation override keys and values cannot be missing.")
    keys <- paste(overrides$sheet, .breseq_seq(overrides$seq_id), .breseq_coord(overrides$position), overrides$field, sep = "\r")
    if (anyDuplicated(keys)) stop("Duplicate annotation overrides are ambiguous.")
    overrides$original_value <- character(nrow(overrides))
    overrides$output_row <- integer(nrow(overrides))
    for (i in seq_len(nrow(overrides))) {
      name <- overrides$sheet[i]
      allowed <- switch(name, Mutations = c("Gene(s)", "Product / region annotation", "Variant annotation"),
        `Missing coverage` = c("Reported gene region", "Gene / region annotation"),
        `Plot-derived evidence` = c("Gene / region", "Figure label"), character())
      if (!name %in% names(sheets) || !overrides$field[i] %in% allowed) {
        stop("Annotation overrides may change annotation fields only, never calls, coordinates or evidence.")
      }
      lookup <- switch(name, Mutations = audit$mutation_keys,
        `Missing coverage` = audit$mc_keys, `Plot-derived evidence` = audit$plot_keys)
      at <- unique(lookup$output_row[lookup$seq_id == .breseq_seq(overrides$seq_id[i]) &
        lookup$position == .breseq_coord(overrides$position[i])])
      if (length(at) != 1L) stop("Annotation override must match exactly one compact row: ", name, " / ", overrides$seq_id[i], " / ", overrides$position[i])
      overrides$original_value[i] <- sheets[[name]]$data[[overrides$field[i]]][at]
      overrides$output_row[i] <- at
      sheets[[name]]$data[[overrides$field[i]]][at] <- overrides$value[i]
    }
    audit$annotation_overrides <- overrides
  }
  list(sheets = sheets, sample_order = samples, sample_labels = labels, stage_colors = colors,
       baseline_samples = baseline, audit = audit)
}

#' Save Compact Publication Mutation and Evidence Tables
#'
#' Creates an additional, organized workbook without changing the canonical wide
#' tables. Reference-only WT cells are labelled Ref.; absent calls are not carried
#' forward. MC families contained entirely in prior finalized deletion regions are
#' omitted from this compact view and retained in the companion wide source file.
#' @param merged_results Named list of canonical merged mutation/evidence tables.
#' @param output_file Output Excel path for the additional compact workbook.
#' @param sample_order Optional exact chronological sample IDs.
#' @param sample_labels Optional named sample-ID to display-label map.
#' @param reference_file Optional GenBank/FASTA reference for replicon labels and products.
#' @param gff_file Optional GFF annotation used to resolve gene products.
#' @param plot_evidence Optional audited plotting-evidence data frame, kept in a separate sheet.
#' @param annotation_overrides Optional data frame with sheet, seq_id, position,
#'   field and value columns. Applies explicit annotation aliases without changing calls.
#' @return Invisibly, the normalized compact tables and presentation metadata.
#' @export
save_publication_results_to_excel <- function(merged_results, output_file, sample_order = NULL,
                                               sample_labels = NULL, reference_file = NULL,
                                               gff_file = NULL, plot_evidence = NULL,
                                               annotation_overrides = NULL) {
  prepared <- .breseq_publication_tables(merged_results, sample_order, sample_labels,
                                        reference_file, gff_file, plot_evidence, annotation_overrides)
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10, fontColour = "#222222")
  body <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, fontColour = "#222222", fgFill = "#FFFFFF", wrapText = TRUE, valign = "center")
  header <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, fontColour = "#222222", textDecoration = "bold", fgFill = "#FFFFFF", wrapText = TRUE, valign = "center")
  band <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, fontColour = "#222222", textDecoration = "bold", fgFill = "#F5F5F5", wrapText = TRUE, halign = "center", valign = "center")
  for (name in names(prepared$sheets)) {
    sheet <- prepared$sheets[[name]]; data <- sheet$data; width <- ncol(data)
    openxlsx::addWorksheet(wb, name, gridLines = FALSE)
    for (row in 1:4) {
      openxlsx::writeData(wb, name, if (row == 1L) sheet$title else sheet$notes[row - 1L], startRow = row, colNames = FALSE)
      openxlsx::mergeCells(wb, name, cols = seq_len(width), rows = row)
    }
    for (group in sheet$groups) {
      openxlsx::writeData(wb, name, group$label, startRow = 6L, startCol = min(group$cols), colNames = FALSE)
      if (length(group$cols) > 1L) openxlsx::mergeCells(wb, name, cols = group$cols, rows = 6L)
    }
    if (nrow(data)) openxlsx::writeDataTable(wb, name, data, startRow = 7L, tableStyle = "TableStyleLight1", bandedRows = FALSE, withFilter = TRUE) else
      openxlsx::writeData(wb, name, data, startRow = 7L)
    openxlsx::addStyle(wb, name, body, rows = seq_len(nrow(data) + 7L), cols = seq_len(width), gridExpand = TRUE)
    openxlsx::addStyle(wb, name, header, rows = c(1L, 7L), cols = seq_len(width), gridExpand = TRUE, stack = TRUE)
    openxlsx::addStyle(wb, name, band, rows = 6L, cols = seq_len(width), gridExpand = TRUE, stack = TRUE)
    openxlsx::setColWidths(wb, name, cols = seq_len(width), widths = sheet$widths)
    openxlsx::setRowHeights(wb, name, rows = 1:7, heights = c(24, 34, 34, 34, 8, 24, 30))
    if (nrow(data)) {
      heights <- vapply(seq_len(nrow(data)), function(i) {
        lines <- vapply(seq_len(width), function(j) {
          text <- as.character(data[[j]][i]); if (is.na(text)) return(1)
          sum(pmax(1, ceiling(nchar(strsplit(text, "\n", fixed = TRUE)[[1]]) / max(1, sheet$widths[j] - 3))))
        }, 1)
        min(409.5, max(30, max(lines) * 14 + 6))
      }, 1)
      openxlsx::setRowHeights(wb, name, rows = seq_len(nrow(data)) + 7L, heights = heights)
      coordinate_columns <- which(names(data) %in% c("Position (bp)", "Start (bp)", "End (bp)"))
      if (length(coordinate_columns)) openxlsx::addStyle(wb, name,
        openxlsx::createStyle(numFmt = "#,##0"), rows = seq_len(nrow(data)) + 7L,
        cols = coordinate_columns, gridExpand = TRUE, stack = TRUE)
      ratio_column <- which(names(data) == "Depth ratio")
      if (length(ratio_column)) openxlsx::addStyle(wb, name,
        openxlsx::createStyle(numFmt = "0.000000"), rows = seq_len(nrow(data)) + 7L,
        cols = ratio_column, gridExpand = TRUE, stack = TRUE)
    }
    for (j in seq_along(sheet$stage_columns)) {
      sample <- prepared$sample_order[j]; column <- sheet$stage_columns[j]
      openxlsx::addStyle(wb, name, openxlsx::createStyle(fgFill = prepared$stage_colors[[sample]]),
        rows = 7L, cols = column, stack = TRUE)
      present <- which(sheet$states[, sample] != "not_reported")
      if (length(present)) openxlsx::addStyle(wb, name,
        openxlsx::createStyle(fgFill = prepared$stage_colors[[sample]]),
        rows = present + 7L, cols = column, gridExpand = TRUE, stack = TRUE)
      bold <- which(sheet$states[, sample] %in% c("new", "changed", "reappeared"))
      if (length(bold)) openxlsx::addStyle(wb, name,
        openxlsx::createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold"),
        rows = bold + 7L, cols = column, gridExpand = TRUE, stack = TRUE)
      absent_rows <- which(sheet$states[, sample] == "not_reported")
      if (length(absent_rows)) openxlsx::addStyle(wb, name,
        openxlsx::createStyle(fontName = "Arial", fontSize = 10, fontColour = "#999999",
          fgFill = if (sample %in% prepared$baseline_samples) prepared$stage_colors[[sample]] else "#FFFFFF"),
        rows = absent_rows + 7L, cols = column, gridExpand = TRUE, stack = TRUE)
    }
  }
  .breseq_save_table_workbook(wb, output_file)
  invisible(prepared)
}

#!/usr/bin/env Rscript

fail <- function(message) {
  stop(message, call. = FALSE)
}

load_required_packages <- function() {
  required_packages <- c("edgeR", "limma", "corral", "PCAtools", "scran", "WGCNA", "rtracklayer")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    fail(paste0(
      "Missing required R packages: ", paste(missing_packages, collapse = ", "),
      ". Install the RNA container dependencies before running this script."
    ))
  }
}

parse_cli <- function(arguments) {
  defaults <- list(
    n_geno_pcs = NULL,
    phenotype_pc_noise = NULL,
    connectivity_z = -3,
    logcpm_drop = 1,
    z_cutoffs_file = NULL,
    z_cutoff = NULL,
    threads = 1L
  )
  required <- c("counts", "genotype_covariates", "gencode", "out_dir")
  values <- defaults
  supplied <- character()

  if (length(arguments) %% 2L != 0L) {
    fail("Arguments must be supplied as --name value pairs")
  }
  for (index in seq.int(1L, length(arguments), by = 2L)) {
    flag <- arguments[[index]]
    raw_key <- sub("^--", "", flag)
    key <- gsub("-", "_", raw_key, fixed = TRUE)
    if (!startsWith(flag, "--") || !nzchar(raw_key)) {
      fail("Arguments must use --name value pairs")
    }
    if (!key %in% c(required, names(defaults))) {
      fail(paste0("Unknown argument --", raw_key))
    }
    if (key %in% supplied) {
      fail(paste0("Argument --", key, " was supplied more than once"))
    }
    values[[key]] <- arguments[[index + 1L]]
    supplied <- c(supplied, key)
  }
  missing <- required[vapply(required, function(key) is.null(values[[key]]), logical(1))]
  if (length(missing) > 0L) {
    fail(paste("Missing required arguments:", paste(paste0("--", missing), collapse = ", ")))
  }

  parse_integer <- function(value, name, minimum) {
    number <- suppressWarnings(as.integer(value))
    if (length(number) != 1L || is.na(number) || as.character(number) != value || number < minimum) {
      fail(paste0("--", name, " must be an integer >= ", minimum))
    }
    number
  }
  parse_number <- function(value, name) {
    number <- suppressWarnings(as.numeric(value))
    if (length(number) != 1L || !is.finite(number)) {
      fail(paste0("--", name, " must be a finite number"))
    }
    number
  }

  if (!is.null(values$n_geno_pcs)) {
    values$n_geno_pcs <- parse_integer(values$n_geno_pcs, "n-geno-pcs", 1L)
  }
  values$threads <- parse_integer(values$threads, "threads", 1L)
  values$connectivity_z <- parse_number(values$connectivity_z, "connectivity-z")
  values$logcpm_drop <- parse_number(values$logcpm_drop, "logcpm-drop")
  if (values$logcpm_drop < 0) {
    fail("--logcpm-drop must be non-negative")
  }
  if (!is.null(values$z_cutoff) && !is.null(values$z_cutoffs_file)) {
    fail("--z-cutoff and --z-cutoffs-file cannot be supplied together")
  }
  if (!is.null(values$z_cutoff)) {
    values$z_cutoff <- parse_number(values$z_cutoff, "z-cutoff")
    values$z_cutoffs <- values$z_cutoff
  } else if (is.null(values$z_cutoffs_file)) {
    values$z_cutoffs <- seq(-1, -10, by = -1)
  } else {
    if (!file.exists(values$z_cutoffs_file)) {
      fail(paste0("--z-cutoffs-file does not exist: ", values$z_cutoffs_file))
    }
    raw_cutoffs <- trimws(readLines(values$z_cutoffs_file, warn = FALSE))
    raw_cutoffs <- raw_cutoffs[nzchar(raw_cutoffs)]
    if (length(raw_cutoffs) == 0L) {
      fail("--z-cutoffs-file must contain at least one cutoff")
    }
    values$z_cutoffs <- vapply(
      raw_cutoffs,
      function(value) parse_number(value, "z-cutoff"),
      numeric(1)
    )
  }
  if (any(values$z_cutoffs >= 0) || anyDuplicated(values$z_cutoffs)) {
    fail("z-score cutoffs must be unique and strictly negative")
  }
  if (!is.null(values$phenotype_pc_noise)) {
    values$phenotype_pc_noise <- parse_number(values$phenotype_pc_noise, "phenotype-pc-noise")
    if (values$phenotype_pc_noise < 0) {
      fail("--phenotype-pc-noise must be non-negative")
    }
  }
  values
}

read_tsv_rows <- function(path, label) {
  if (!file.exists(path)) {
    fail(paste0(label, " file does not exist: ", path))
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(lines) < 2L) {
    fail(paste0(label, " TSV must contain a header and at least one data row"))
  }
  split_tsv_line <- function(line) {
    without_terminal_tabs <- sub("\t+$", "", line)
    terminal_tab_count <- nchar(line) - nchar(without_terminal_tabs)
    if (!nzchar(without_terminal_tabs)) {
      return(rep("", terminal_tab_count + 1L))
    }
    c(
      strsplit(without_terminal_tabs, "\t", fixed = TRUE)[[1L]],
      rep("", terminal_tab_count)
    )
  }
  fields <- lapply(lines, split_tsv_line)
  fields[[1L]][1L] <- sub("^\ufeff", "", fields[[1L]][1L])
  width <- length(fields[[1L]])
  if (width < 2L || any(vapply(fields, length, integer(1)) != width)) {
    fail(paste0(label, " TSV has inconsistent tab-delimited field counts"))
  }
  list(header = fields[[1L]], rows = fields[-1L])
}

numeric_matrix <- function(rows, row_ids, column_ids, label, require_nonnegative = FALSE) {
  values <- unlist(rows, use.names = FALSE)
  numbers <- suppressWarnings(as.numeric(values))
  if (length(numbers) != length(values) || any(!is.finite(numbers))) {
    fail(paste0(label, " contains nonnumeric or non-finite values"))
  }
  if (require_nonnegative && any(numbers < 0)) {
    fail(paste0(label, " contains a negative count"))
  }
  matrix(numbers, nrow = length(row_ids), ncol = length(column_ids),
         byrow = TRUE, dimnames = list(row_ids, column_ids))
}

read_counts <- function(path) {
  parsed <- read_tsv_rows(path, "Counts")
  if (parsed$header[[1L]] != "gene_id") {
    fail("Counts TSV header must begin with gene_id")
  }
  sample_ids <- parsed$header[-1L]
  if (any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    fail("Counts TSV sample IDs must be nonempty and unique")
  }
  gene_ids <- vapply(parsed$rows, `[[`, character(1), 1L)
  if (any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    fail("Counts TSV gene_id values must be nonempty and unique")
  }
  count_rows <- lapply(parsed$rows, function(row) row[-1L])
  numeric_matrix(count_rows, gene_ids, sample_ids, "Counts TSV", require_nonnegative = TRUE)
}

read_genotype_covariates <- function(path, sample_ids, n_geno_pcs) {
  parsed <- read_tsv_rows(path, "Genotype covariates")
  sample_column <- which(parsed$header == "sample_id")
  if (length(sample_column) != 1L || anyDuplicated(parsed$header)) {
    fail("Genotype-covariate TSV must contain one unique sample_id column")
  }
  pc_columns <- parsed$header[-sample_column]
  pc_numbers <- suppressWarnings(as.integer(sub("^(?:Genotype_PC|GENETICPC)", "", pc_columns)))
  if (any(is.na(pc_numbers)) || anyDuplicated(pc_numbers) ||
      !identical(sort(pc_numbers), seq_along(pc_numbers))) {
    fail("Genotype-covariate TSV columns after sample_id must be named Genotype_PC1... or GENETICPC1... without gaps")
  }
  pc_order <- order(pc_numbers)
  pc_columns <- pc_columns[pc_order]
  pc_numbers <- pc_numbers[pc_order]
  canonical_pc_columns <- paste0("Genotype_PC", pc_numbers)
  needed <- if (is.null(n_geno_pcs)) canonical_pc_columns else paste0("Genotype_PC", seq_len(n_geno_pcs))
  if (!is.null(n_geno_pcs) && n_geno_pcs > length(canonical_pc_columns)) {
    fail(paste0(
      "Requested ", n_geno_pcs, " genotype PCs but the covariate TSV contains only ",
      length(canonical_pc_columns)
    ))
  }
  covariate_samples <- vapply(parsed$rows, `[[`, character(1), sample_column)
  if (any(!nzchar(covariate_samples)) || anyDuplicated(covariate_samples)) {
    fail("Genotype-covariate sample_id values must be nonempty and unique")
  }
  if (!setequal(covariate_samples, sample_ids) || length(covariate_samples) != length(sample_ids)) {
    shared_samples <- sample_ids[sample_ids %in% covariate_samples]
    if (length(shared_samples) == 0L) {
      fail("Genotype-covariate sample_id values have no overlap with counts TSV samples")
    }
    warning(
      paste0(
        "Retaining ", length(shared_samples), " samples shared by counts and genotype covariates; ",
        length(setdiff(sample_ids, covariate_samples)), " counts samples and ",
        length(setdiff(covariate_samples, sample_ids)), " covariate samples were excluded"
      ),
      call. = FALSE
    )
  } else {
    shared_samples <- sample_ids
  }
  all_covariates <- numeric_matrix(
    lapply(parsed$rows, function(row) row[-sample_column]), covariate_samples, pc_columns,
    "Genotype covariates"
  )
  colnames(all_covariates) <- canonical_pc_columns
  all_covariates[shared_samples, needed, drop = FALSE]
}

parse_gff_attributes <- function(attributes) {
  pieces <- strsplit(attributes, ";", fixed = TRUE)[[1L]]
  pairs <- strsplit(pieces, "=", fixed = TRUE)
  keys <- vapply(pairs, function(pair) if (length(pair) > 0L) pair[[1L]] else "", character(1))
  values <- vapply(
    pairs,
    function(pair) if (length(pair) > 1L) paste(pair[-1L], collapse = "=") else "",
    character(1)
  )
  stats::setNames(utils::URLdecode(values), keys)
}

read_gene_metadata <- function(path) {
  if (!file.exists(path)) {
    fail(paste0("GENCODE file does not exist: ", path))
  }
  gff <- rtracklayer::readGFF(path)
  if (!"type" %in% names(gff) || !"gene_id" %in% names(gff)) {
    fail("GENCODE GFF/GFF3 must contain gene records with a gene_id attribute")
  }
  type_column <- if ("gene_type" %in% names(gff)) "gene_type" else if ("gene_biotype" %in% names(gff)) "gene_biotype" else NA_character_
  if (is.na(type_column)) {
    fail("GENCODE GFF/GFF3 must contain gene_type or gene_biotype attributes")
  }
  genes <- gff[gff$type == "gene" & gff[[type_column]] == "protein_coding", , drop = FALSE]
  if (nrow(genes) == 0L) {
    fail("GENCODE GFF/GFF3 contains no protein-coding gene records")
  }
  gene_id <- sub("\\..*$", "", as.character(genes$gene_id))
  symbol <- if ("gene_name" %in% names(genes)) as.character(genes$gene_name) else gene_id
  keep <- !duplicated(gene_id)
  data.frame(gene_nv = gene_id[keep], symbol = symbol[keep], stringsAsFactors = FALSE)
}

filter_expression <- function(counts) {
  poisson <- scran::modelGeneVarByPoisson(counts)
  top_poisson <- scran::getTopHVGs(poisson)
  cv2 <- scran::modelGeneCV2(counts)
  top_cv2 <- scran::getTopHVGs(cv2, var.field = "ratio", var.threshold = 1.0)
  selected <- intersect(top_poisson, top_cv2)
  if (length(selected) == 0L) {
    fail("Expression filtering retained no genes")
  }
  counts[selected, , drop = FALSE]
}

freeman_tukey_normalize <- function(counts) {
  normalized <- as.matrix(corral::corral_preproc(counts, rtype = "freemantukey"))
  colnames(normalized) <- colnames(counts)
  normalized
}

connectivity_kept_samples <- function(normalized_expression, threshold) {
  adjacency <- 0.5 + 0.5 * stats::cor(as.matrix(normalized_expression))
  connectivity <- WGCNA::fundamentalNetworkConcepts(adjacency)$Connectivity
  names(connectivity) <- colnames(normalized_expression)
  connectivity_sd <- stats::sd(connectivity)
  if (!is.finite(connectivity_sd) || connectivity_sd == 0) {
    return(colnames(normalized_expression))
  }
  connectivity_z <- (connectivity - mean(connectivity)) / connectivity_sd
  kept <- names(connectivity_z)[connectivity_z > threshold]
  if (length(kept) < 2L) {
    fail("Connectivity QC retained fewer than two samples")
  }
  kept
}

# PCAtools::chooseGavishDonoho expects the variance of the entry-level noise.
# The INT/unit-variance working convention uses one as the temporary null.
default_phenotype_pc_noise <- 1

strict_output_tag <- function(z_cutoff) {
  formatted_cutoff <- format(z_cutoff, trim = TRUE, scientific = FALSE, digits = 15)
  formatted_cutoff <- sub("[.]0+$", "", formatted_cutoff)
  paste0("z_", formatted_cutoff)
}

validate_covariate_design <- function(covariates, sample_ids, context) {
  if (!identical(rownames(covariates), sample_ids)) {
    fail(paste0(context, " covariates are not ordered to match the retained samples"))
  }
  if (ncol(covariates) == 0L || any(!nzchar(colnames(covariates))) || anyDuplicated(colnames(covariates))) {
    fail(paste0(context, " must contain uniquely named covariate columns"))
  }
  if (any(!is.finite(covariates))) {
    fail(paste0(context, " contains non-finite values"))
  }
  design <- cbind(Intercept = 1, covariates)
  design_rank <- qr(design)$rank
  design_columns <- colnames(covariates)
  if (design_rank != ncol(design)) {
    fail(sprintf(
      "%s is rank-deficient: %d samples, %d covariates (%s), design rank %d of %d",
      context, nrow(covariates), ncol(covariates), paste(design_columns, collapse = ", "),
      design_rank, ncol(design)
    ))
  }
  residual_degrees_freedom <- nrow(covariates) - design_rank
  if (residual_degrees_freedom < 1L) {
    fail(sprintf(
      "%s leaves no residual degrees of freedom: %d samples, design rank %d, covariates (%s)",
      context, nrow(covariates), design_rank, paste(design_columns, collapse = ", ")
    ))
  }
  list(
    columns = design_columns,
    design_rank = design_rank,
    residual_degrees_freedom = residual_degrees_freedom
  )
}

remove_covariates <- function(expression, covariates) {
  limma::removeBatchEffect(expression, covariates = covariates)
}

selected_pc_metadata <- function(
  selected_pcs,
  selected_raw,
  available_rank,
  noise,
  noise_source,
  n_genotype_pcs,
  phenotype_design,
  genotype_pc_columns,
  residual_design,
  n_samples_after_qc,
  n_genes_after_qc
) {
  data.frame(
    phenotype_pc_method = "PCAtools::chooseGavishDonoho",
    selected_phenotype_pcs = selected_pcs,
    gavish_donoho_raw = as.integer(selected_raw),
    available_rank = available_rank,
    noise_variance = noise,
    noise_source = noise_source,
    n_genotype_pcs = n_genotype_pcs,
    phenotype_pc_columns = paste(phenotype_design$columns, collapse = ","),
    genotype_pc_columns = paste(genotype_pc_columns, collapse = ","),
    residualization_covariate_columns = paste(residual_design$columns, collapse = ","),
    phenotype_design_rank = phenotype_design$design_rank,
    residualization_design_rank = residual_design$design_rank,
    residual_degrees_freedom = residual_design$residual_degrees_freedom,
    n_samples_after_qc = n_samples_after_qc,
    n_genes_after_qc = n_genes_after_qc,
    stringsAsFactors = FALSE
  )
}

long_expression <- function(expression, value_name) {
  data.frame(
    gene_id = rep(rownames(expression), times = ncol(expression)),
    subjectid = rep(colnames(expression), each = nrow(expression)),
    value = as.vector(expression),
    stringsAsFactors = FALSE
  ) |> stats::setNames(c("gene_id", "subjectid", value_name))
}

call_underliers <- function(expr_z_join, z_cutoff, logcpm_drop) {
  gene_means <- stats::ave(expr_z_join$expression_logcpm, expr_z_join$gene_id, FUN = mean)
  expr_z_join$gene_mean_logcpm <- gene_means
  underliers <- expr_z_join[
    expr_z_join$expression_logcpm < expr_z_join$gene_mean_logcpm - logcpm_drop,
    , drop = FALSE
  ]
  if (!is.null(z_cutoff)) {
    underliers <- underliers[underliers$expression_zscore < z_cutoff, , drop = FALSE]
  }
  underliers
}

prevalence_by_gene <- function(underliers, gene_metadata, n_subjects_tested) {
  if (nrow(underliers) == 0L) {
    return(data.frame(
      gene_id = character(), gene_nv = character(), symbol = character(),
      n_underlier_subjects = integer(), n_subjects_tested = integer(),
      rna_outlier_prevalence = numeric(), stringsAsFactors = FALSE
    ))
  }
  events <- unique(underliers[c("gene_id", "subjectid")])
  counts <- stats::aggregate(subjectid ~ gene_id, events, length)
  names(counts)[[2L]] <- "n_underlier_subjects"
  counts$gene_nv <- sub("\\..*$", "", counts$gene_id)
  counts$n_subjects_tested <- n_subjects_tested
  counts$rna_outlier_prevalence <- counts$n_underlier_subjects / n_subjects_tested
  result <- merge(counts, gene_metadata, by = "gene_nv", all.x = TRUE, sort = FALSE)
  result <- result[c(
    "gene_id", "gene_nv", "symbol", "n_underlier_subjects", "n_subjects_tested",
    "rna_outlier_prevalence"
  )]
  result[order(-result$rna_outlier_prevalence, result$gene_id), , drop = FALSE]
}

write_tsv <- function(data, path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.table(data, connection, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

run_pipeline <- function(options) {
  dir.create(options$out_dir, recursive = TRUE, showWarnings = FALSE)
  counts <- read_counts(options$counts)
  genotype_covariates <- read_genotype_covariates(
    options$genotype_covariates, colnames(counts), options$n_geno_pcs
  )
  counts <- counts[, rownames(genotype_covariates), drop = FALSE]
  gene_metadata <- read_gene_metadata(options$gencode)
  gene_nv <- sub("\\..*$", "", rownames(counts))
  coding <- match(gene_nv, gene_metadata$gene_nv)
  coding_counts <- counts[!is.na(coding), , drop = FALSE]
  if (nrow(coding_counts) == 0L) {
    fail("No count-matrix genes matched protein-coding GENCODE genes")
  }

  normalized_pass1 <- freeman_tukey_normalize(filter_expression(coding_counts))
  kept_samples <- connectivity_kept_samples(normalized_pass1, options$connectivity_z)
  qc_counts <- coding_counts[, kept_samples, drop = FALSE]
  qc_normalized <- freeman_tukey_normalize(filter_expression(qc_counts))
  genotype_covariates <- genotype_covariates[kept_samples, , drop = FALSE]

  pca_result <- PCAtools::pca(as.matrix(qc_normalized), scale = TRUE)
  variance_spectrum <- pca_result$sdev^2
  available_rank <- min(length(variance_spectrum), ncol(pca_result$rotated))
  if (available_rank < 1L) {
    fail("PCA returned no available phenotype-PC components")
  }
  noise_source <- if (is.null(options$phenotype_pc_noise)) "unit_variance_default" else "override"
  noise <- if (is.null(options$phenotype_pc_noise)) default_phenotype_pc_noise else options$phenotype_pc_noise
  selected_raw <- PCAtools::chooseGavishDonoho(
    as.matrix(qc_normalized), var.explained = variance_spectrum, noise = noise
  )
  if (length(selected_raw) != 1L || !is.finite(selected_raw) || selected_raw < 1L) {
    fail("PCAtools::chooseGavishDonoho selected zero phenotype PCs")
  }
  selected_pcs <- min(as.integer(selected_raw), available_rank)
  expression_pc_names <- paste0("PC", seq_len(selected_pcs))
  if (!all(expression_pc_names %in% colnames(pca_result$rotated))) {
    fail("PCA output does not contain the selected phenotype-PC score columns")
  }
  rotated_scores <- as.matrix(pca_result$rotated)
  storage.mode(rotated_scores) <- "double"
  expression_pcs <- rotated_scores[kept_samples, expression_pc_names, drop = FALSE]
  phenotype_design <- validate_covariate_design(
    expression_pcs, kept_samples, "Selected phenotype-PC adjustment"
  )

  counts_aligned <- qc_counts[rownames(qc_normalized), kept_samples, drop = FALSE]
  logcpm <- edgeR::cpm(edgeR::DGEList(counts_aligned), log = TRUE)
  logcpm_adjusted <- remove_covariates(logcpm, expression_pcs)
  covariates <- cbind(expression_pcs, genotype_covariates)
  residual_design <- validate_covariate_design(
    covariates, kept_samples, "Phenotype/genotype residualization"
  )
  residualized <- remove_covariates(qc_normalized, covariates)
  expression_zscore <- t(scale(t(residualized)))
  if (any(!is.finite(expression_zscore))) {
    fail("Residual expression z-scores contain non-finite values")
  }

  expr_z_join <- merge(
    long_expression(expression_zscore, "expression_zscore"),
    long_expression(logcpm_adjusted, "expression_logcpm"),
    by = c("gene_id", "subjectid"), sort = FALSE
  )
  expr_z_join <- expr_z_join[order(expr_z_join$gene_id, expr_z_join$subjectid), , drop = FALSE]
  output_path <- function(name) file.path(options$out_dir, name)
  write_tsv(expr_z_join, output_path("expr_z_join.tsv.gz"))

  metadata <- selected_pc_metadata(
    selected_pcs,
    selected_raw,
    available_rank,
    noise,
    noise_source,
    ncol(genotype_covariates),
    phenotype_design,
    colnames(genotype_covariates),
    residual_design,
    length(kept_samples),
    nrow(qc_normalized)
  )
  write_tsv(metadata, output_path("selected_phenotype_pcs.tsv"))

  definitions <- c(
    list(list(tag = "haplo", z_cutoff = NULL)),
    lapply(options$z_cutoffs, function(z_cutoff) {
      list(tag = strict_output_tag(z_cutoff), z_cutoff = z_cutoff)
    })
  )
  for (definition in definitions) {
    underliers <- call_underliers(expr_z_join, definition$z_cutoff, options$logcpm_drop)
    write_tsv(underliers, output_path(paste0("underliers_", definition$tag, ".tsv.gz")))
    prevalence <- prevalence_by_gene(underliers, gene_metadata, length(kept_samples))
    write_tsv(
      prevalence,
      output_path(paste0("rna_outlier_prevalence_per_gene_", definition$tag, ".tsv"))
    )
  }
}

run_main <- function() {
  options <- parse_cli(commandArgs(trailingOnly = TRUE))
  load_required_packages()
  run_pipeline(options)
}

if (!identical(Sys.getenv("TOPMED_RNA_UNDERLIER_NO_MAIN"), "1")) {
  tryCatch(
    run_main(),
    error = function(error) {
      message("run_rna_underlier.R: error: ", conditionMessage(error))
      quit(status = 2L)
    }
  )
}

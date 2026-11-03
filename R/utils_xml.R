## R/utils_xml.R
## EvoTempo: BEAST XML Auto-detector and Date Randomizer
## Supports both BEAST 1.x and BEAST 2.x syntax.
## CSV clusters: clustered permutation adapted from Murray et al. (2016),
## perm.rand(); https://doi.org/10.1111/2041-210X.12466.
## Mclust clusters: independent outside-cluster sampling adapted from
## TipDatingBeast (Rieux & Khatchikian, 2017; doi:10.1111/1755-0998.12603).

generate_randomized_xmls <- function(xml_path, out_dir, base_name, nreps = 20, 
                                     use_clusters = FALSE, clusters = NULL,
                                     cluster_method = "murray") {
  xml_lines <- readLines(xml_path)
  is_beast2 <- any(grepl("<trait.*traitname=\"date(-forward)?\"", xml_lines))
  is_beast1 <- any(grepl("<date value=", xml_lines))
  
  if (!is_beast1 && !is_beast2) {
    stop("Unrecognized XML format: Neither BEAST 1 '<date value=...>' nor BEAST 2 '<trait ... traitname=\"date\">' tags were found.")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is_beast2) {
    message("BEAST 2 format detected. Initializing BEAST 2 parsing engine...")
    process_beast2_xml(xml_lines, out_dir, base_name, nreps, use_clusters, clusters, cluster_method)
  } else {
    message("BEAST 1 format detected. Initializing BEAST 1 parsing engine...")
    process_beast1_xml(xml_lines, out_dir, base_name, nreps, use_clusters, clusters, cluster_method)
  }
}

validate_drt_clusters <- function(original_dates, clusters) {
  if (length(clusters) != length(original_dates) || anyNA(clusters) ||
      any(!nzchar(trimws(as.character(clusters)))) || length(unique(clusters)) < 2) {
    stop("Cluster-based DRT requires a non-missing cluster for every date and at least two clusters.")
  }
}

# TipDatingBeast RandomCluster core: each taxon draws independently from the
# original observations outside its own cluster. Pools are not depleted.
get_shuffled_dates <- function(original_dates, use_clusters = FALSE, clusters = NULL) {
  if (!use_clusters) {
    return(original_dates[sample.int(length(original_dates))])
  } else {
    validate_drt_clusters(original_dates, clusters)
    new_dates <- numeric(length(original_dates))   
    for (i in seq_along(original_dates)) {
      c_id <- clusters[i]
      valid_pool <- original_dates[clusters != c_id]      
      if (length(valid_pool) == 0) {
        new_dates[i] <- original_dates[i]
      } else {
        new_dates[i] <- valid_pool[sample.int(length(valid_pool), size = 1)]
      }
    }
    return(new_dates)
  }
}

# Adapted from Murray's perm.rand(): one permutation maps all destination
# clusters to source clusters; original dates are sampled within each source.
# For single-date clusters this retains a common date within each cluster.
# Differences from the reference: normalize labels, unambiguous cache keys,
# deduplicate actual date vectors as well as mappings, and fail explicitly
# rather than returning NA when new permutations cannot be found.
get_murray_date_permutations <- function(original_dates, clusters, nreps) {
  validate_drt_clusters(original_dates, clusters)
  clusters <- match(clusters, unique(clusters))
  n_clusters <- length(unique(clusters))
  pools <- split(original_dates, clusters)
  single_date <- all(vapply(pools, function(x) length(unique(x)) == 1L, logical(1)))
  log_total <- lfactorial(n_clusters)
  if (single_date) {
    cluster_dates <- vapply(pools, function(x) x[1], numeric(1))
    log_total <- log_total - sum(lfactorial(table(match(cluster_dates, unique(cluster_dates)))))
  }
  # For mixed-date clusters, the number of source-cluster mappings is an
  # upper bound; retain Murray's rule of accepting each mapping at most once.
  n_generate <- if (log_total < log(nreps + 1)) {
    min(nreps, max(0, round(exp(log_total)) - 1))
  } else nreps
  if (n_generate < 1) stop("No non-original clustered date permutations are available. Check the cluster dates.")

  date_key <- function(x) paste(sprintf("%a", x), collapse = ",")
  seen_maps <- new.env(hash = TRUE, parent = emptyenv())
  seen_dates <- new.env(hash = TRUE, parent = emptyenv())
  assign(paste(seq_len(n_clusters), collapse = ","), TRUE, envir = seen_maps)
  assign(date_key(original_dates), TRUE, envir = seen_dates)
  result <- matrix(NA_real_, nrow = length(original_dates), ncol = n_generate)
  accepted <- 0L
  attempts <- 0L
  max_attempts <- max(10000, n_generate * 1000)
  while (accepted < n_generate) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop("Unable to find the requested number of distinct clustered date permutations. Request fewer replicates or check the cluster dates.")
    }
    mapping <- sample.int(n_clusters)
    map_key <- paste(mapping, collapse = ",")
    source_clusters <- mapping[clusters]
    new_dates <- vapply(source_clusters, function(id) {
      pool <- pools[[id]]
      pool[sample.int(length(pool), size = 1)]
    }, numeric(1))
    if (exists(map_key, envir = seen_maps, inherits = FALSE)) next
    key <- date_key(new_dates)
    if (exists(key, envir = seen_dates, inherits = FALSE)) next
    assign(map_key, TRUE, envir = seen_maps)
    assign(key, TRUE, envir = seen_dates)
    accepted <- accepted + 1L
    result[, accepted] <- new_dates
  }
  result
}

make_drt_date_replicates <- function(original_dates, nreps, use_clusters = FALSE,
                                   clusters = NULL, cluster_method = "murray") {
  if (!is.numeric(original_dates) || length(original_dates) < 2 || any(!is.finite(original_dates))) {
    stop("The XML must contain at least two finite numeric sampling dates.")
  }
  if (length(nreps) != 1 || !is.finite(nreps) || nreps < 1 || nreps != floor(nreps)) {
    stop("The number of randomizations must be a positive integer.")
  }
  if (use_clusters) {
    cluster_method <- match.arg(cluster_method, c("murray", "tipdatingbeast"))
    validate_drt_clusters(original_dates, clusters)
    if (cluster_method == "murray") return(get_murray_date_permutations(original_dates, clusters, nreps))
  }
  replicate(nreps, get_shuffled_dates(original_dates, use_clusters, clusters))
}

process_beast2_xml <- function(xml_lines, out_dir, base_name, nreps, use_clusters, clusters,
                               cluster_method = "murray") {
  trait_line_idx <- grep("<trait.*traitname=\"date(-forward)?\"", xml_lines)
  trait_line <- xml_lines[trait_line_idx]
  value_str <- sub(".*value=\"([^\"]+)\".*", "\\1", trait_line)
  pairs <- strsplit(value_str, ",")[[1]]
  pairs <- trimws(pairs)
  taxa_names <- sapply(strsplit(pairs, "="), `[`, 1)
  original_dates <- as.numeric(sapply(strsplit(pairs, "="), `[`, 2))
  date_sets <- make_drt_date_replicates(original_dates, nreps, use_clusters, clusters, cluster_method)
  for (i in seq_len(ncol(date_sets))) {
    shuffled_dates <- date_sets[, i]
    new_pairs <- paste(taxa_names, shuffled_dates, sep = "=")
    new_value_str <- paste(new_pairs, collapse = ",\n\t\t\t\t")
    new_trait_line <- sub(value_str, new_value_str, trait_line, fixed = TRUE)
    new_xml <- xml_lines
    new_xml[trait_line_idx] <- new_trait_line
    new_xml <- gsub("fileName=\"[^\"]+\\.log\"", sprintf("fileName=\"%s_%02d.log\"", base_name, i), new_xml)
    new_xml <- gsub("fileName=\"[^\"]+\\.trees\"", sprintf("fileName=\"%s_%02d.trees\"", base_name, i), new_xml)
    
    out_file <- file.path(out_dir, sprintf("%s_%02d.xml", base_name, i))
    writeLines(new_xml, out_file)
  }
  invisible(ncol(date_sets))
}

process_beast1_xml <- function(xml_lines, out_dir, base_name, nreps, use_clusters, clusters,
                               cluster_method = "murray") {
  datelines <- grep("<date value=", xml_lines)
  original_dates <- as.numeric(
    gsub('\" direction.*', '', 
         gsub('.*<date value=\"', '', xml_lines[datelines]))
  )
  date_sets <- make_drt_date_replicates(original_dates, nreps, use_clusters, clusters, cluster_method)
  for (i in seq_len(ncol(date_sets))) {
    shuffled_dates <- date_sets[, i]
    new_xml <- xml_lines
    for (j in seq_along(datelines)) {
      line_idx <- datelines[j]
      old_line <- new_xml[line_idx]
      new_date_str <- as.character(shuffled_dates[j])
      pattern <- 'value="[^"]*"'
      replacement <- sprintf("value=\"%s\"", new_date_str)    
      new_xml[line_idx] <- sub(pattern, replacement, old_line)
    }
    new_xml <- gsub("fileName=\"[^\"]+\\.log\"", sprintf("fileName=\"%s_%02d.log\"", base_name, i), new_xml)
    new_xml <- gsub("fileName=\"[^\"]+\\.trees\"", sprintf("fileName=\"%s_%02d.trees\"", base_name, i), new_xml)  
    out_file <- file.path(out_dir, sprintf("%s_%02d.xml", base_name, i))
    writeLines(new_xml, out_file)
  }
  invisible(ncol(date_sets))
}

generate_drt_xmls <- function(xml_path, out_dir, n_reps, prefix) {
  generate_randomized_xmls(
    xml_path = xml_path, 
    out_dir = out_dir, 
    base_name = prefix, 
    nreps = n_reps,
    use_clusters = FALSE, 
    clusters = NULL
  )
}

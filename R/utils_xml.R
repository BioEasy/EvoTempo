## R/utils_xml.R
## EvoTempo: BEAST XML Auto-detector and Date Randomizer
## Supports both BEAST 1.x and BEAST 2.x syntax.
## Cluster-based exclusive random sampling logic adapted from TipDatingBeast 
## (Rieux & Khatchikian, 2017; https://doi.org/10.1111/1755-0998.12603).

generate_randomized_xmls <- function(xml_path, out_dir, base_name, nreps = 20, 
                                     use_clusters = FALSE, clusters = NULL) {
  xml_lines <- readLines(xml_path)
  is_beast2 <- any(grepl("<trait.*traitname=\"date(-forward)?\"", xml_lines))
  is_beast1 <- any(grepl("<date value=", xml_lines))
  
  if (!is_beast1 && !is_beast2) {
    stop("Unrecognized XML format: Neither BEAST 1 '<date value=...>' nor BEAST 2 '<trait ... traitname=\"date\">' tags were found.")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is_beast2) {
    message("BEAST 2 format detected. Initializing BEAST 2 parsing engine...")
    process_beast2_xml(xml_lines, out_dir, base_name, nreps, use_clusters, clusters)
  } else {
    message("BEAST 1 format detected. Initializing BEAST 1 parsing engine...")
    process_beast1_xml(xml_lines, out_dir, base_name, nreps, use_clusters, clusters)
  }
}

get_shuffled_dates <- function(original_dates, use_clusters = FALSE, clusters = NULL) {
  if (!use_clusters || is.null(clusters)) {
    return(sample(original_dates))
  } else {
    if (length(clusters) != length(original_dates) || anyNA(clusters) ||
        length(unique(clusters)) < 2) {
      stop("Cluster-based DRT requires a non-missing cluster for every date and at least two clusters.")
    }
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

process_beast2_xml <- function(xml_lines, out_dir, base_name, nreps, use_clusters, clusters) {
  trait_line_idx <- grep("<trait.*traitname=\"date(-forward)?\"", xml_lines)
  trait_line <- xml_lines[trait_line_idx]
  value_str <- sub(".*value=\"([^\"]+)\".*", "\\1", trait_line)
  pairs <- strsplit(value_str, ",")[[1]]
  pairs <- trimws(pairs)
  taxa_names <- sapply(strsplit(pairs, "="), `[`, 1)
  original_dates <- as.numeric(sapply(strsplit(pairs, "="), `[`, 2))
  for (i in 1:nreps) {
    shuffled_dates <- get_shuffled_dates(original_dates, use_clusters, clusters)
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
}

process_beast1_xml <- function(xml_lines, out_dir, base_name, nreps, use_clusters, clusters) {
  datelines <- grep("<date value=", xml_lines)
  original_dates <- as.numeric(
    gsub('\" direction.*', '', 
         gsub('.*<date value=\"', '', xml_lines[datelines]))
  )
  for (i in 1:nreps) {
    shuffled_dates <- get_shuffled_dates(original_dates, use_clusters, clusters)    
    new_xml <- xml_lines
    for (j in seq_along(datelines)) {
      line_idx <- datelines[j]
      old_line <- new_xml[line_idx]
      old_date_str <- as.character(original_dates[j])
      new_date_str <- as.character(shuffled_dates[j])
      pattern <- sprintf("value=\"%s\"", old_date_str)
      replacement <- sprintf("value=\"%s\"", new_date_str)    
      new_xml[line_idx] <- sub(pattern, replacement, old_line, fixed = TRUE)
    }
    new_xml <- gsub("fileName=\"[^\"]+\\.log\"", sprintf("fileName=\"%s_%02d.log\"", base_name, i), new_xml)
    new_xml <- gsub("fileName=\"[^\"]+\\.trees\"", sprintf("fileName=\"%s_%02d.trees\"", base_name, i), new_xml)  
    out_file <- file.path(out_dir, sprintf("%s_%02d.xml", base_name, i))
    writeLines(new_xml, out_file)
  }
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

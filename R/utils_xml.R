## R/utils_xml.R
## EvoTempo: BEAST XML Auto-detector and Date Randomizer
## Supports both BEAST 1.x and BEAST 2.x syntax.
## Cluster-based exclusive random sampling logic adapted from TipDatingBeast 
## (Rieux & Khatchikian, 2017; https://doi.org/10.1111/1755-0998.12603).

# ---------------------------------------------------------
# 1. Core randomizer generator (cross-version support)
# ---------------------------------------------------------
generate_randomized_xmls <- function(xml_path, out_dir, base_name, nreps = 20, 
                                     use_clusters = FALSE, clusters = NULL) {
  
  # Read all lines of the XML file
  xml_lines <- readLines(xml_path)
  
  # Automatic version detection
  is_beast2 <- any(grepl("<trait.*traitname=\"date(-forward)?\"", xml_lines))
  is_beast1 <- any(grepl("<date value=", xml_lines))
  
  if (!is_beast1 && !is_beast2) {
    stop("Unrecognized XML format: Neither BEAST 1 '<date value=...>' nor BEAST 2 '<trait ... traitname=\"date\">' tags were found.")
  }
  
  # Create output directory
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Route to the corresponding parsing engine
  if (is_beast2) {
    message("BEAST 2 format detected. Initializing BEAST 2 parsing engine...")
    process_beast2_xml(xml_lines, out_dir, base_name, nreps, use_clusters, clusters)
  } else {
    message("BEAST 1 format detected. Initializing BEAST 1 parsing engine...")
    process_beast1_xml(xml_lines, out_dir, base_name, nreps, use_clusters, clusters)
  }
}

# ---------------------------------------------------------
# 2. Internal permutation algorithm: Generate shuffled date vectors
# ---------------------------------------------------------
get_shuffled_dates <- function(original_dates, use_clusters = FALSE, clusters = NULL) {
  if (!use_clusters || is.null(clusters)) {
    # Standard mode: completely random shuffle of all sequence dates
    return(sample(original_dates))
  } else {
    # Cluster mode (Exclusive Random Sampling, adapted from Rieux et al.)
    new_dates <- numeric(length(original_dates))
    
    for (i in seq_along(original_dates)) {
      c_id <- clusters[i]
      
      # Pool of all dates excluding the current sequence's own cluster
      valid_pool <- original_dates[clusters != c_id]
      
      if (length(valid_pool) == 0) {
        # Fallback if only one cluster exists (prevents crashing)
        new_dates[i] <- original_dates[i]
      } else {
        # Randomly sample ONE date from the valid external pool
        new_dates[i] <- sample(valid_pool, size = 1)
      }
    }
    return(new_dates)
  }
}

# ---------------------------------------------------------
# 3. BEAST 2 Parsing and Injection Engine
# ---------------------------------------------------------
process_beast2_xml <- function(xml_lines, out_dir, base_name, nreps, use_clusters, clusters) {
  # Locate the trait line
  trait_line_idx <- grep("<trait.*traitname=\"date(-forward)?\"", xml_lines)
  trait_line <- xml_lines[trait_line_idx]
  
  # Extract the "taxon=date,taxon=date" long string from the value attribute
  value_str <- sub(".*value=\"([^\"]+)\".*", "\\1", trait_line)
  
  # Split into key-value pairs
  pairs <- strsplit(value_str, ",")[[1]]
  pairs <- trimws(pairs)
  
  taxa_names <- sapply(strsplit(pairs, "="), `[`, 1)
  original_dates <- as.numeric(sapply(strsplit(pairs, "="), `[`, 2))
  
  # Loop to generate the specified number of randomized files
  for (i in 1:nreps) {
    # Core: Get the shuffled dates
    shuffled_dates <- get_shuffled_dates(original_dates, use_clusters, clusters)
    
    # Reassemble the string
    new_pairs <- paste(taxa_names, shuffled_dates, sep = "=")
    new_value_str <- paste(new_pairs, collapse = ",\n\t\t\t\t")
    
    # Replace original line content (use fixed = TRUE to avoid regex escape errors)
    new_trait_line <- sub(value_str, new_value_str, trait_line, fixed = TRUE)
    
    # Inject new line and save
    new_xml <- xml_lines
    new_xml[trait_line_idx] <- new_trait_line
    
    # Automatically modify the output names of log and trees files in BEAST 2
    new_xml <- gsub("fileName=\"[^\"]+\\.log\"", sprintf("fileName=\"%s_%02d.log\"", base_name, i), new_xml)
    new_xml <- gsub("fileName=\"[^\"]+\\.trees\"", sprintf("fileName=\"%s_%02d.trees\"", base_name, i), new_xml)
    
    out_file <- file.path(out_dir, sprintf("%s_%02d.xml", base_name, i))
    writeLines(new_xml, out_file)
  }
}

# ---------------------------------------------------------
# 4. BEAST 1 Parsing and Injection Engine
# ---------------------------------------------------------
process_beast1_xml <- function(xml_lines, out_dir, base_name, nreps, use_clusters, clusters) {
  # Locate all lines containing date settings
  datelines <- grep("<date value=", xml_lines)
  
  # Extract original dates
  original_dates <- as.numeric(
    gsub('\" direction.*', '', 
         gsub('.*<date value=\"', '', xml_lines[datelines]))
  )
  
  # Loop to generate randomized files
  for (i in 1:nreps) {
    shuffled_dates <- get_shuffled_dates(original_dates, use_clusters, clusters)
    
    new_xml <- xml_lines
    
    # Replace date values line by line
    for (j in seq_along(datelines)) {
      line_idx <- datelines[j]
      old_line <- new_xml[line_idx]
      old_date_str <- as.character(original_dates[j])
      new_date_str <- as.character(shuffled_dates[j])
      
      # Exactly replace the numerical value inside value=""
      pattern <- sprintf("value=\"%s\"", old_date_str)
      replacement <- sprintf("value=\"%s\"", new_date_str)
      
      new_xml[line_idx] <- sub(pattern, replacement, old_line, fixed = TRUE)
    }
    
    # Automatically modify the output names of log and trees files in BEAST 1
    new_xml <- gsub("fileName=\"[^\"]+\\.log\"", sprintf("fileName=\"%s_%02d.log\"", base_name, i), new_xml)
    new_xml <- gsub("fileName=\"[^\"]+\\.trees\"", sprintf("fileName=\"%s_%02d.trees\"", base_name, i), new_xml)
    
    out_file <- file.path(out_dir, sprintf("%s_%02d.xml", base_name, i))
    writeLines(new_xml, out_file)
  }
}

# ---------------------------------------------------------
# 5. UI interface bridge
# ---------------------------------------------------------
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
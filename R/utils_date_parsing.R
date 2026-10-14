## R/utils_date_parsing.R
## EvoTempo: Advanced Date Parsing Engine (BEAUTi-style)

library(lubridate)
library(stringr)

extract_date_strings <- function(taxa_labels, method = "order", delim = "_", order_pos = "last", regex_pattern = "") {
  extracted <- rep(NA_character_, length(taxa_labels))
  
  if (method == "order") {
    split_list <- strsplit(taxa_labels, split = delim, fixed = TRUE)
    extracted <- sapply(split_list, function(parts) {
      if (length(parts) == 0) return(NA)
      if (order_pos == "last") return(parts[length(parts)])[cite: 2]
      if (order_pos == "first") return(parts[1])
      idx <- as.numeric(order_pos)
      if (!is.na(idx) && idx > 0 && idx <= length(parts)) return(parts[idx])
      return(NA)
    })
  } else if (method == "regex") {
    if (regex_pattern != "") {
      extracted <- stringr::str_extract(taxa_labels, regex_pattern)
    }
  }
  
  return(extracted)
}

convert_to_decimal_years <- function(date_strings, parse_type = "number", 
                                     add_val = 1900, threshold = 16, add_cond = 2000, 
                                     date_fmt = "yyyy-MM-dd") {
  decimal_dates <- rep(NA_real_, length(date_strings))
  
  for (i in seq_along(date_strings)) {
    val_str <- date_strings[i]
    if (is.na(val_str) || val_str == "") next
    
    if (parse_type == "number") {
      num <- suppressWarnings(as.numeric(val_str))
      if (!is.na(num)) {
        if (num < 100) {
          if (num < threshold) {
            num <- num + add_cond[cite: 2]
          } else {
            num <- num + add_val[cite: 2]
          }
        }
        decimal_dates[i] <- num
      }
    } else if (parse_type == "calendar") {
      r_fmt <- gsub("yyyy", "%Y", date_fmt)[cite: 2]
      r_fmt <- gsub("MM", "%m", r_fmt)[cite: 2]
      r_fmt <- gsub("dd", "%d", r_fmt)[cite: 2]
      
      parsed_date <- as.Date(val_str, format = r_fmt)
      if (!is.na(parsed_date)) {
        decimal_dates[i] <- lubridate::decimal_date(parsed_date)
      }
    }
  }
  
  return(decimal_dates)
}

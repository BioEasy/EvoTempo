## R/utils_date_parsing.R
## EvoTempo: Advanced Date Parsing Engine (BEAUTi-style)

library(lubridate)
library(stringr)

extract_date_strings <- function(taxa_labels, method = "order", delim = "_", order_pos = "last", regex_pattern = "") {
  taxa_clean <- gsub("^['\"]|['\"]$", "", trimws(as.character(taxa_labels)))
  extracted <- rep(NA_character_, length(taxa_clean))
  
  if (method == "order") {
    if (!nzchar(delim)) return(extracted)
    
    split_list <- strsplit(taxa_clean, split = delim, fixed = TRUE)
    extracted <- sapply(split_list, function(parts) {
      parts <- trimws(parts)
      n_parts <- length(parts)
      if (n_parts == 0) return(NA_character_)
      
      pos_str <- tolower(trimws(as.character(order_pos)))
      cand <- NA_character_
      
      if (pos_str == "last") {
        cand <- parts[n_parts]
      } else if (pos_str == "first") {
        cand <- parts[1]
      } else {
        idx <- suppressWarnings(as.numeric(pos_str))
        if (!is.na(idx)) {
          if (idx > 0 && idx <= n_parts) cand <- parts[idx]
          else if (idx < 0 && (n_parts + idx + 1) >= 1) cand <- parts[n_parts + idx + 1]
        }
      }

      is_date_like <- function(x) {
        !is.na(x) && grepl("^[0-9]{4}([-/.][0-9]{1,2}([-/.][0-9]{1,2})?|\\.[0-9]+)?$", x)
      }
      if (!is_date_like(cand)) {
        hits <- parts[sapply(parts, is_date_like)]
        if (length(hits) > 0) {
          cand <- hits[length(hits)]
        }
      }
      
      return(cand)
    })
  } else if (method == "regex") {
    if (nzchar(regex_pattern)) {
      extracted <- stringr::str_extract(taxa_clean, regex_pattern)
    }
  }
  
  return(trimws(extracted))
}

convert_to_decimal_years <- function(date_strings, parse_type = "number", 
                                     add_val = 1900, threshold = 16, add_cond = 2000, 
                                     date_fmt = "yyyy-MM-dd", enable_2digit = FALSE) {
  decimal_dates <- rep(NA_real_, length(date_strings))
  r_fmt <- if (!is.null(date_fmt) && nzchar(date_fmt)) date_fmt else "yyyy-MM-dd"
  r_fmt <- gsub("yyyy", "%Y", r_fmt, ignore.case = TRUE)
  r_fmt <- gsub("yy", "%y", r_fmt, ignore.case = TRUE)
  r_fmt <- gsub("MM", "%m", r_fmt)
  r_fmt <- gsub("dd", "%d", r_fmt, ignore.case = TRUE)
  
  parse_as_calendar <- function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_real_)
    d_obj <- tryCatch(
      suppressWarnings(as.Date(s, format = r_fmt)),
      error = function(e) as.Date(NA)
    )
    if (is.na(d_obj)) {
      d_obj <- tryCatch(
        suppressWarnings(as.Date(s, tryFormats = c(
          "%Y-%m-%d", "%Y/%m/%d", "%d/%m/%Y", "%m/%d/%Y", "%Y.%m.%d", "%Y-%m", "%Y/%m"
        ))),
        error = function(e) as.Date(NA)
      )
    }
    if (is.na(d_obj) && requireNamespace("lubridate", quietly = TRUE)) {
      d_parsed <- suppressWarnings(lubridate::parse_date_time(s, orders = c("ymd", "Ymd", "dmy", "mdy", "ym", "Ym"), quiet = TRUE))
      if (!is.na(d_parsed)) d_obj <- as.Date(d_parsed)
    }
    
    if (!is.na(d_obj)) {
      if (requireNamespace("lubridate", quietly = TRUE)) {
        return(lubridate::decimal_date(d_obj))
      } else {
        yr <- as.numeric(format(d_obj, "%Y"))
        doy <- as.numeric(format(d_obj, "%j"))
        days_in_yr <- ifelse((yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0), 366, 365)
        return(yr + (doy - 1) / days_in_yr)
      }
    }

    if (grepl("^[0-9]{4}[-/][0-9]{1,2}$", s)) {
      parts <- as.numeric(strsplit(s, "[-/]")[[1]])
      if (parts[2] >= 1 && parts[2] <= 12) {
        return(parts[1] + (parts[2] - 0.5) / 12)
      }
    }
    return(NA_real_)
  }
  
  for (i in seq_along(date_strings)) {
    val_str <- gsub("^['\"]|['\"]$", "", trimws(as.character(date_strings[i])))
    if (is.na(val_str) || val_str == "") next
    
    if (parse_type == "number") {
      if (grepl("^[0-9]+(\\.[0-9]+)?$", val_str)) {
        num <- suppressWarnings(as.numeric(val_str))
        if (!is.na(num)) {
          if (isTRUE(enable_2digit) && num < 100) {
            if (num < threshold) {
              num <- num + add_cond
            } else {
              num <- num + add_val
            }
          }
          decimal_dates[i] <- num
        }
      } else {
        decimal_dates[i] <- parse_as_calendar(val_str)
      }
    } else if (parse_type == "calendar") {
      cal_val <- parse_as_calendar(val_str)
      if (!is.na(cal_val)) {
        decimal_dates[i] <- cal_val
      } else if (grepl("^[0-9]+(\\.[0-9]+)?$", val_str)) {
        decimal_dates[i] <- as.numeric(val_str)
      }
    }
  }
  
  return(decimal_dates)
}

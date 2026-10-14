## R/utils_date_parsing.R
## EvoTempo: Advanced Date Parsing Engine (BEAUTi-style)

library(lubridate)
library(stringr)

# ---------------------------------------------------------
# 1. 提取时间字符串 (Extraction)
# ---------------------------------------------------------
extract_date_strings <- function(taxa_labels, method = "order", delim = "_", order_pos = "last", regex_pattern = "") {
  extracted <- rep(NA_character_, length(taxa_labels))
  
  if (method == "order") {
    # 按分隔符提取
    split_list <- strsplit(taxa_labels, split = delim, fixed = TRUE)
    extracted <- sapply(split_list, function(parts) {
      if (length(parts) == 0) return(NA)
      if (order_pos == "last") return(parts[length(parts)])[cite: 2]
      if (order_pos == "first") return(parts[1])
      # 处理指定数字索引
      idx <- as.numeric(order_pos)
      if (!is.na(idx) && idx > 0 && idx <= length(parts)) return(parts[idx])
      return(NA)
    })
  } else if (method == "regex") {
    # 按正则表达式提取[cite: 2]
    if (regex_pattern != "") {
      extracted <- stringr::str_extract(taxa_labels, regex_pattern)
    }
  }
  
  return(extracted)
}

# ---------------------------------------------------------
# 2. 格式转换与清洗 (Transformation)
# ---------------------------------------------------------
convert_to_decimal_years <- function(date_strings, parse_type = "number", 
                                     add_val = 1900, threshold = 16, add_cond = 2000, 
                                     date_fmt = "yyyy-MM-dd") {
  decimal_dates <- rep(NA_real_, length(date_strings))
  
  for (i in seq_along(date_strings)) {
    val_str <- date_strings[i]
    if (is.na(val_str) || val_str == "") next
    
    if (parse_type == "number") {
      # 纯数字解析[cite: 2]
      num <- suppressWarnings(as.numeric(val_str))
      if (!is.na(num)) {
        # 两位数年份智能补全逻辑[cite: 2]
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
      # 日历日期解析 (使用 lubridate 转换为十进制年份)[cite: 2]
      # 将 Java 风格的 yyyy-MM-dd 转为 R 的 %Y-%m-%d
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
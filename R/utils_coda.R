## R/utils_coda.R
## EvoTempo: MCMC Log parsing and Date Randomization Test evaluation

library(coda)
get_param_summary <- function(log_file, param_name = "clock.rate", burnin_frac = 0.1) {
  dat <- read.table(log_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
  if (!(param_name %in% colnames(dat))) {
    stop(sprintf("Parameter '%s' not found in the log file.", param_name))
  }
  start_idx <- floor(nrow(dat) * burnin_frac) + 1
  if (start_idx >= nrow(dat)) {
    stop("Burn-in fraction is too large, resulting in empty data.")
  }  
  param_trace <- dat[start_idx:nrow(dat), param_name]
  param_mcmc <- as.mcmc(param_trace)
  hpd_95 <- HPDinterval(param_mcmc, prob = 0.95)
  
  return(list(
    mean = mean(param_mcmc),
    lower = hpd_95[1, "lower"],
    upper = hpd_95[1, "upper"]
  ))
}

evaluate_drt <- function(true_log_path, rand_log_paths, param_name = "clock.rate", burnin_frac = 0.1) {
  true_stats <- get_param_summary(true_log_path, param_name, burnin_frac)
  rand_stats_list <- lapply(rand_log_paths, function(path) {
    get_param_summary(path, param_name, burnin_frac)
  })

  cr1_passed <- logical(length(rand_stats_list))
  cr2_passed <- logical(length(rand_stats_list))
  
  for (i in seq_along(rand_stats_list)) {
    r_stat <- rand_stats_list[[i]]
    cr1_passed[i] <- (true_stats$mean < r_stat$lower) || (true_stats$mean > r_stat$upper)
    cr2_passed[i] <- (true_stats$upper < r_stat$lower) || (true_stats$lower > r_stat$upper)
  }
  final_cr1 <- all(cr1_passed)
  final_cr2 <- all(cr2_passed)
  
  return(list(
    true_stats = true_stats,
    rand_stats = rand_stats_list,
    pass_cr1 = final_cr1,
    pass_cr2 = final_cr2
  ))
}

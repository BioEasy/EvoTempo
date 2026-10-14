## R/utils_coda.R
## EvoTempo: MCMC Log parsing and Date Randomization Test evaluation

library(coda)

# ---------------------------------------------------------
# 1. 核心解析函数：读取单个 BEAST Log 并提取统计量
# ---------------------------------------------------------
get_param_summary <- function(log_file, param_name = "clock.rate", burnin_frac = 0.1) {
  # 读取 BEAST 生成的 tab 分隔日志文件 (跳过以 # 开头的注释行)
  dat <- read.table(log_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
  
  # 检查目标参数是否存在 (BEAST 1 常用 clock.rate, BEAST 2 常用 ucldMean 等)
  if (!(param_name %in% colnames(dat))) {
    stop(sprintf("Parameter '%s' not found in the log file.", param_name))
  }
  
  # 剔除 Burn-in 数据 (默认前 10%)
  start_idx <- floor(nrow(dat) * burnin_frac) + 1
  if (start_idx >= nrow(dat)) {
    stop("Burn-in fraction is too large, resulting in empty data.")
  }
  
  param_trace <- dat[start_idx:nrow(dat), param_name]
  
  # 转换为 MCMC 对象并利用 coda 计算 95% HPD 区间
  param_mcmc <- as.mcmc(param_trace)
  hpd_95 <- HPDinterval(param_mcmc, prob = 0.95)
  
  return(list(
    mean = mean(param_mcmc),
    lower = hpd_95[1, "lower"],
    upper = hpd_95[1, "upper"]
  ))
}

# ---------------------------------------------------------
# 2. DRT 判定引擎：比对真实数据与随机化数据池
# ---------------------------------------------------------
evaluate_drt <- function(true_log_path, rand_log_paths, param_name = "clock.rate", burnin_frac = 0.1) {
  
  # 提取真实数据的均值与 95% HPD
  true_stats <- get_param_summary(true_log_path, param_name, burnin_frac)
  
  # 遍历提取所有随机化数据集的统计量
  rand_stats_list <- lapply(rand_log_paths, function(path) {
    get_param_summary(path, param_name, burnin_frac)
  })
  
  # 初始化逻辑判定向量
  cr1_passed <- logical(length(rand_stats_list))
  cr2_passed <- logical(length(rand_stats_list))
  
  for (i in seq_along(rand_stats_list)) {
    r_stat <- rand_stats_list[[i]]
    
    # CR1 判定：真实数据的 Mean 必须严格落在当前随机化数据的 95% HPD 之外
    cr1_passed[i] <- (true_stats$mean < r_stat$lower) || (true_stats$mean > r_stat$upper)
    
    # CR2 判定：真实数据的 95% HPD 与当前随机化数据的 95% HPD 完全无交集
    cr2_passed[i] <- (true_stats$upper < r_stat$lower) || (true_stats$lower > r_stat$upper)
  }
  
  # 汇总结果：必须战胜 *所有* 随机化数据集才算 Pass
  final_cr1 <- all(cr1_passed)
  final_cr2 <- all(cr2_passed)
  
  return(list(
    true_stats = true_stats,
    rand_stats = rand_stats_list,
    pass_cr1 = final_cr1,
    pass_cr2 = final_cr2
  ))
}
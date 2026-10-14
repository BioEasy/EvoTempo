## R/utils_temporal.R
## Core algorithms adapted from Murray et al. (https://doi.org/10.1111/2041-210X.12466)
## EvoTempo: Underlying statistical and phylogenetic algorithms

library(ape)
library(ade4)
library(seqinr)

# ==============================================================================
# 模块 1：数据读取与 Mantel 检验
# ==============================================================================
# 1. 读取序列/树并提取时间标签
read_sequence_data <- function(filepath, format = "fasta", delim = "_") {
if (format == "fasta") {
    # 1. 读取 FASTA 序列
    seqs <- ape::read.dna(filepath, format = "fasta")
    names_list <- rownames(seqs)
    
    # 2. 计算距离矩阵 (开启 pairwise.deletion 忽略空位)
    dist_mat <- ape::dist.dna(seqs, model = "TN93", pairwise.deletion = TRUE)
    
    # 3. 【核心修复】智能建树判定
    # 如果距离矩阵中出现无法计算的 NA，自动启用 njs() 容错算法
    if (any(is.na(dist_mat))) {
      phy <- ape::njs(dist_mat)
    } else {
      phy <- ape::nj(dist_mat)
    }
    
    # 返回统一的数据结构
    return(list(
      names = names_list,
      phy = phy,
      g = as.matrix(dist_mat)
    )) 
  } else if (format == "newick" || format == "nexus") {
    d <- tryCatch(read.tree(filepath), error = function(e) {
      stop("Failed to parse Newick tree. Please ensure it is a valid Newick format and matches the selected dropdown format.")
    })
    g <- cophenetic(d)
    phy <- d
    
  } else {
    stop("Unsupported file format selected.")
  }
  
  if(phy$Nnode < (length(phy$tip.label) - 1)) {
    phy <- multi2di(phy)
  }
  
  names <- rownames(g)
  x <- paste0(".*\\", date_sep)
  dates <- as.numeric(gsub(x, "", names))
  
  return(list(g = g, dates = dates, names = names, phy = phy))
}
  
# 2. 运行 Mantel Confounding Test[cite: 1]
run_mantel_test <- function(filepath, format = "fasta", date_sep = "_", nreps = 1000) {
  parsed_data <- read_sequence_data(filepath, format, date_sep)
  g_matrix <- parsed_data$g
  dates <- parsed_data$dates
  
  m_test <- suppressWarnings(mantel.rtest(as.dist(g_matrix), dist(dates), nrepet = nreps))
  p_val <- m_test$pval
  r_val <- m_test$obs
  
  dist_df <- data.frame(
    Temporal_Distance = as.vector(dist(dates)),
    Genetic_Distance = as.vector(as.dist(g_matrix))
  )
  
  is_confounded <- p_val < 0.05
  
  return(list(
    p_value = p_val,
    r_value = r_val,
    is_confounded = is_confounded,
    plot_data = dist_df,
    phylo_tree = parsed_data$phy,
    dates = dates
  ))
}

# ==============================================================================
# 模块 2：树聚类算法 (Cluster-based methods)
# ==============================================================================

# 提取单系群分支[cite: 1]
get.clades <- function(phy) {
  phy <- reorder(phy)
  ue <- unique(phy$edge[, 1])	
  fm <- vector("list", phy$Nnode)
  root.lab <- length(phy$tip.label) + 1 
  for(i in 1:length(phy$tip.label)) {
    z <- i
    while(z != root.lab) {
      new.z <- phy$edge[which(phy$edge[, 2] == z), 1]
      m <- match(new.z, ue)
      fm[[m]] <- c(fm[[m]], phy$tip.label[i])
      z <- new.z
    }
  }
  return(fm)
}

# 检查特定分支中包含多少个不同的采样时间[cite: 1]
checkndates <- function(clade, dates, noms) {
  length(unique(dates[match(clade, noms)]))
}

# 核心：自动将系统发育树划分为聚类 (Clusters)[cite: 1]
get_clusters <- function(phy, dates, auto.round.dates = FALSE, rounded.dates = NULL, 
                         do.plot = FALSE, stat = 'nve-rms', reroot = TRUE) {
  
  if(!is.rooted(phy) || reroot) {	
    phy.allroots <- get.all.roots(phy)	
    phy <- get.best.root(phy.allroots, dates, stat)
  }
  
  if(auto.round.dates) { dates <- round(dates + 0.000000000000000001) }
  cl <- 1:length(dates)
  
  gc <- get.clades(phy)
  da <- if(!is.null(rounded.dates)) rounded.dates else dates
  
  single.date.clades <- which(unlist(lapply(gc, FUN = checkndates, dates = da, noms = phy$tip.label)) == 1)
  
  if(length(single.date.clades) > 0) {
    for(i in 1:length(single.date.clades)) {
      m <- match(gc[[single.date.clades[i]]], phy$tip.label)		
      cl[m] <- cl[m[1]]
    }
  }
  
  cl <- match(cl, unique(cl))
  return(cl)
}

# ==============================================================================
# 模块 3：最佳定根与 RTT 回归分析 (Root-to-Tip Optimization)
# ==============================================================================

prepare.phy <- function(phy, node, nt) {
  
  # 1. 严格区分：定根在外部枝(Tip) 还是 内部枝(Internal Node)
  if (node <= nt) {
    # 外部节点必须传给 outgroup 参数
    rooted_tree <- root(phy, outgroup = node, resolve.root = TRUE)
  } else {
    # 内部节点传给 node 参数
    rooted_tree <- root(phy, node = node, resolve.root = TRUE)
  }
  
  # 2. 精准定位父节点为根节点 (nt + 1) 的那两条分支（完美规避 Segfault 闪退）
  root_edges <- which(rooted_tree$edge[, 1] == (nt + 1))
  
  # 3. 提取这两条分支的总长度，并平分赋值 (加入 na.rm 防止意外 NA 传递)
  b <- sum(rooted_tree$edge.length[root_edges], na.rm = TRUE)
  rooted_tree$edge.length[root_edges] <- b / 2
  
  return(rooted_tree)
}

get.all.roots <- function(phy) {
  nt <- length(phy$tip.label) 
  n <- nt + phy$Nnode
  phy.all.roots <- vector('list', n)
  for(i in c(1:nt, (nt+2):n)) {
    phy.all.roots[[i]] <- prepare.phy(phy, i, nt)
  }
  return(phy.all.roots[-(nt + 1)])
}

get.altered.tree <- function(phy, p, root.node) {
  w <- which(phy$edge[, 1] == root.node)
  s <- phy$edge.length[w[1]] + phy$edge.length[w[2]]
  phy$edge.length[w[1]] <- p * s
  phy$edge.length[w[2]] <- (1 - p) * s
  return(phy)
}

myfun <- function(p, phy, dates, stat = 'nve-rms', nt) {
  phy <- get.altered.tree(phy, p, nt + 1)
  div <- dist.nodes(phy)[nt + 1, 1:nt]
  mylm <- lm(div ~ dates)
  
  if (stat == 'r') {
    s <- summary(mylm)
    return(sqrt(s$r.squared) * sign(s$coefficients[2, 1]))
  } else if (stat == 'r2') {
    return(summary(mylm)$r.squared)
  } else if (stat == 'nve-rms') {
    return(-1e5 * (sum((mean(residuals(mylm)) - residuals(mylm))^2) / mylm$df.residual))
  }
}

optfun <- function(phy, dates, nt, stat = 'nve-rms') {
  ans <- optim(runif(1), fn = myfun, method = 'L-BFGS', lower = 0, upper = 1, 
               control = list(fnscale = -1), phy = phy, dates = dates, stat = stat, nt = nt)
  return(c(ans$par, ans$value))
}

get.best.root <- function(phy.all.roots, dates, stat = 'nve-rms') {
  nt <- length(phy.all.roots[[1]]$tip.label)
  p_list <- lapply(phy.all.roots, FUN = optfun, dates = dates, stat = stat, nt = nt)	
  p_mat <- matrix(unlist(p_list), length(phy.all.roots), 2, byrow = TRUE)
  w <- which.max(p_mat[, 2])
  return(get.altered.tree(phy.all.roots[[w]], p_mat[w, 1], nt + 1))
}

# 核心暴露函数：供 Shiny RTT 模块调用
run_rtt_regression <- function(phy, dates, stat = 'nve-rms', reroot = TRUE) {
  
  if (!is.rooted(phy) || reroot) {
    phy.allroots <- get.all.roots(phy)	
    phy <- get.best.root(phy.allroots, dates, stat)
  }
  
  nt <- length(phy$tip.label)
  root_node <- nt + 1
  distances <- dist.nodes(phy)[root_node, 1:nt]
  
  mylm <- lm(distances ~ dates)
  slm <- summary(mylm)
  
  r <- sqrt(slm$r.squared) * sign(slm$coefficients[2, 1])
  slope <- slm$coefficients[2, 1]
  intercepty <- slm$coefficients[1, 1]
  tmrca <- -intercepty / slope
  rms <- sum((mean(residuals(mylm)) - residuals(mylm))^2) / mylm$df.residual
  
  df_plot <- data.frame(
    Taxon = phy$tip.label,
    Date = dates,
    Distance = distances,
    Residuals = residuals(mylm)
  )
  
  return(list(
    phylo_rooted = phy,
    df_plot = df_plot,
    stats = list(
      r_value = r,
      r_squared = slm$r.squared,
      clock_rate = slope,
      tmrca = tmrca,
      rms = rms
    ),
    lm_model = mylm
  ))
}

# ---------------------------------------------------------
# 新增：接收预处理矩阵和日期的 Mantel 检验函数
# ---------------------------------------------------------
run_mantel_from_data <- function(phy, g_mat, dates, nreps = 1000) {
  # 1. 确保对角线自身距离为 0
  diag(g_mat) <- 0
  
  # 2. 建立时间距离矩阵
  t_dist <- as.matrix(dist(dates))
  rownames(t_dist) <- colnames(t_dist) <- phy$tip.label
  
  # 3. 提取矩阵下三角，展开为成对向量
  g_vec <- g_mat[lower.tri(g_mat)]
  t_vec <- t_dist[lower.tri(t_dist)]
  
  # 【核心修复】：找出所有非 NA、非 Inf 的有效序列配对
  valid_idx <- is.finite(g_vec) & is.finite(t_vec)
  g_vec_valid <- g_vec[valid_idx]
  t_vec_valid <- t_vec[valid_idx]
  
  if (length(g_vec_valid) < 3) {
    stop("Not enough valid pairwise distances to run Mantel Test. Please check your alignment.")
  }
  
  # 4. 只针对有效配对计算相关系数
  obs_r <- cor(g_vec_valid, t_vec_valid, method = "pearson")
  
  # 5. 置换检验 (Permutation) 也要同步过滤
  null_r <- numeric(nreps)
  for (i in 1:nreps) {
    shuff_dates <- sample(dates)
    shuff_t_dist <- as.matrix(dist(shuff_dates))
    shuff_t_vec <- shuff_t_dist[lower.tri(shuff_t_dist)]
    
    # 用相同的 valid_idx 提取洗牌后的时间距离
    shuff_t_vec_valid <- shuff_t_vec[valid_idx]
    
    null_r[i] <- cor(g_vec_valid, shuff_t_vec_valid, method = "pearson")
  }
  
  # 6. 计算 p-value
  p_val <- sum(null_r >= obs_r, na.rm = TRUE) / nreps
  
  # 7. 严格保证 is_confounded 输出必须是 TRUE 或 FALSE
  is_conf <- (!is.na(p_val) && p_val < 0.05)
  
  # 8. 组装只包含有效配对的绘图数据
  df_plot <- data.frame(
    Genetic_Distance = g_vec_valid,
    Temporal_Distance = t_vec_valid
  )
  
  list(
    phylo_tree = phy,
    dates = dates,
    r_value = obs_r,
    p_value = p_val,
    is_confounded = is_conf,
    plot_data = df_plot
  )
}
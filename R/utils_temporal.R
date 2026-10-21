## R/utils_temporal.R
## Core algorithms adapted from Murray et al. (https://doi.org/10.1111/2041-210X.12466)

library(ape)
library(ade4)
library(seqinr)


read_sequence_data <- function(filepath, format = "fasta", delim = "_") {
  if (format == "fasta") {
    seqs <- ape::read.dna(filepath, format = "fasta")
    names_list <- rownames(seqs)
    dist_mat <- ape::dist.dna(seqs, model = "TN93", pairwise.deletion = TRUE)
    if (any(is.na(dist_mat))) {
      phy <- ape::njs(dist_mat)
    } else {
      phy <- ape::nj(dist_mat)
    }
    
    g <- as.matrix(dist_mat)
    names <- names_list
    
  } else if (format == "newick" || format == "nexus") {
    lines <- readLines(filepath, warn = FALSE)
    lines_trim <- trimws(lines)
    lines_non_empty <- lines_trim[nzchar(lines_trim)]
    if (length(lines_non_empty) == 0) {
      stop("The uploaded file is empty.")
    }
    first_line <- lines_non_empty[1]
    
    phy <- NULL
    if (format == "nexus" || grepl("^#NEXUS", first_line, ignore.case = TRUE)) {
      phy <- tryCatch(
        ape::read.nexus(filepath),
        error = function(e) {
          nwk_text <- paste(lines_non_empty, collapse = "")
          ape::read.tree(text = nwk_text)
        }
      )
    } else {
      nwk_text <- paste(lines_non_empty, collapse = "")
      phy <- ape::read.tree(text = nwk_text)
    }
    
    if (inherits(phy, "multiPhylo")) {
      phy <- phy[[1]]
    }
    if (is.null(phy) || !inherits(phy, "phylo")) {
      stop("Failed to parse valid phylogenetic tree from file. Please ensure it is a valid Newick/NEXUS format.")
    }
    
    phy$tip.label <- gsub("^['\"]|['\"]$", "", trimws(phy$tip.label))
    
    if (is.null(phy$edge.length)) {
      stop("The uploaded tree does not contain branch lengths.")
    }
    
    g <- ape::cophenetic.phylo(phy)
    g <- g[phy$tip.label, phy$tip.label, drop = FALSE]
    names <- phy$tip.label
    
  } else {
    stop("Unsupported file format selected.")
  }
  
  if(phy$Nnode < (length(phy$tip.label) - 1)) {
    phy <- multi2di(phy)
  }
  
  dates <- suppressWarnings(as.numeric(sapply(names, function(nm) {
    m <- regmatches(nm, regexpr("[0-9]{4}(\\.[0-9]+)?(-[0-9]{2}-[0-9]{2})?", nm))
    if (length(m) > 0) {
      if (grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", m[1])) {
        d_obj <- as.Date(m[1])
        if (!is.na(d_obj) && requireNamespace("lubridate", quietly = TRUE)) {
          return(lubridate::decimal_date(d_obj))
        }
      }
      return(as.numeric(m[1]))
    }
    parts <- trimws(unlist(strsplit(nm, delim, fixed = TRUE)))
    cand <- parts[sapply(parts, function(x) grepl("^[0-9]{4}", x))]
    if (length(cand) > 0) {
      s <- cand[length(cand)]
      d_obj <- suppressWarnings(as.Date(s))
      if (!is.na(d_obj) && requireNamespace("lubridate", quietly = TRUE)) {
        return(lubridate::decimal_date(d_obj))
      }
      return(suppressWarnings(as.numeric(s)))
    }
    return(NA_real_)
  })))
  
  return(list(g = g, dates = dates, names = names, phy = phy))
}
  
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

checkndates <- function(clade, dates, noms) {
  length(unique(dates[match(clade, noms)]))
}

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


prepare.phy <- function(phy, node, nt) {
  
  if (node <= nt) {
    rooted_tree <- root(phy, outgroup = node, resolve.root = TRUE)
  } else {
    rooted_tree <- root(phy, node = node, resolve.root = TRUE)
  }
  root_edges <- which(rooted_tree$edge[, 1] == (nt + 1))
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
  
  div <- ape::node.depth.edgelength(phy)[1:nt]
  
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

run_rtt_regression <- function(phy, dates, stat = 'nve-rms', reroot = TRUE) {
  
  if (!is.rooted(phy) || reroot) {
    phy.allroots <- get.all.roots(phy)  
    phy <- get.best.root(phy.allroots, dates, stat)
  }
  
  nt <- length(phy$tip.label)
  
  distances <- ape::node.depth.edgelength(phy)[1:nt]
  
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

run_mantel_from_data <- function(phy, g_mat, dates, nreps = 1000) {
  diag(g_mat) <- 0
  t_dist <- as.matrix(dist(dates))
  rownames(t_dist) <- colnames(t_dist) <- phy$tip.label
  g_vec <- g_mat[lower.tri(g_mat)]
  t_vec <- t_dist[lower.tri(t_dist)]
  valid_idx <- is.finite(g_vec) & is.finite(t_vec)
  g_vec_valid <- g_vec[valid_idx]
  t_vec_valid <- t_vec[valid_idx]
  
  if (length(g_vec_valid) < 3) {
    stop("Not enough valid pairwise distances to run Mantel Test. Please check your alignment.")
  }
  
  obs_r <- cor(g_vec_valid, t_vec_valid, method = "pearson")

  null_r <- numeric(nreps)
  for (i in 1:nreps) {
    shuff_dates <- sample(dates)
    shuff_t_dist <- as.matrix(dist(shuff_dates))
    shuff_t_vec <- shuff_t_dist[lower.tri(shuff_t_dist)]
    shuff_t_vec_valid <- shuff_t_vec[valid_idx]
    
    null_r[i] <- cor(g_vec_valid, shuff_t_vec_valid, method = "pearson")
  }
  p_val <- sum(null_r >= obs_r, na.rm = TRUE) / nreps
  is_conf <- (!is.na(p_val) && p_val < 0.05)
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

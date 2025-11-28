library(MASS)

detect_outliers_mad <- function(x, thresh = 3) {
  med <- median(x)
  mad_val <- mad(x, constant = 1.4826)
  z <- abs(x - med) / mad_val
  which(z > thresh)
}

one_rep <- function(n, rho, contam, thresh) {
  
  # 1. Generate clean data
  Sigma <- matrix(c(1, rho, rho, 1), 2)
  z <- mvrnorm(n, mu = c(0,0), Sigma = Sigma)
  x <- z[,1]
  y <- z[,2]
  
  # 2. Inject outliers
  outl_num <- floor(n * contam)
  if(outl_num > 0) {
    idx <- sample(1:n, outl_num)
    x[idx] <- rnorm(outl_num, mean = 8, sd = 1)
    y[idx] <- rnorm(outl_num, mean = 8, sd = 1)
  }
  
  # raw correlations & tests
  raw_cor_pearson  <- cor(x, y)
  raw_cor_spearman <- cor(x, y, method="spearman")
  raw_cor_kendall  <- cor(x, y, method="kendall")
  
  raw_p_pearson  <- cor.test(x, y)$p.value
  raw_p_spearman <- cor.test(x, y, method="spearman")$p.value
  raw_p_kendall  <- cor.test(x, y, method="kendall")$p.value
  
  # clean data
  out <- union(detect_outliers_mad(x, thresh), detect_outliers_mad(y, thresh))
  if(length(out) > 0) {
    x_clean <- x[-out]
    y_clean <- y[-out]
  } else {
    x_clean <- x
    y_clean <- y
  }
  
  clean_cor_pearson  <- cor(x_clean, y_clean)
  clean_cor_spearman <- cor(x_clean, y_clean, method="spearman")
  clean_cor_kendall  <- cor(x_clean, y_clean, method="kendall")
  
  clean_p_pearson  <- cor.test(x_clean, y_clean)$p.value
  clean_p_spearman <- cor.test(x_clean, y_clean, method="spearman")$p.value
  clean_p_kendall  <- cor.test(x_clean, y_clean, method="kendall")$p.value
  
  return(c(
    raw_cor_pearson = raw_cor_pearson,
    clean_cor_pearson = clean_cor_pearson,
    raw_cor_spearman = raw_cor_spearman,
    clean_cor_spearman = clean_cor_spearman,
    raw_cor_kendall = raw_cor_kendall,
    clean_cor_kendall = clean_cor_kendall,
    raw_p_pearson = raw_p_pearson,
    clean_p_pearson = clean_p_pearson,
    raw_p_spearman = raw_p_spearman,
    clean_p_spearman = clean_p_spearman,
    raw_p_kendall = raw_p_kendall,
    clean_p_kendall = clean_p_kendall
  ))
}

run_sim <- function(rho = 0, contam_levels = seq(0, 0.5, 0.05),
                    R = 1000, n = 200, thresh = 3) {
  
  results <- list()
  
  for (cont in contam_levels) {
    print(cont)
    sims <- replicate(R, one_rep(n, rho, cont, thresh))
    sims <- t(sims)
    
    results[[as.character(cont)]] <- list(
      mean_raw_cor_pearson = mean(sims[, "raw_cor_pearson"]),
      mean_clean_cor_pearson = mean(sims[, "clean_cor_pearson"]),
      mean_raw_cor_spearman = mean(sims[, "raw_cor_spearman"]),
      mean_clean_cor_spearman = mean(sims[, "clean_cor_spearman"]),
      mean_raw_cor_kendall = mean(sims[, "raw_cor_kendall"]),
      mean_clean_cor_kendall = mean(sims[, "clean_cor_kendall"]),
      type1_raw_pearson = mean(sims[, "raw_p_pearson"] < 0.05),
      type1_clean_pearson = mean(sims[, "clean_p_pearson"] < 0.05),
      type1_raw_spearman = mean(sims[, "raw_p_spearman"] < 0.05),
      type1_clean_spearman = mean(sims[, "clean_p_spearman"] < 0.05),
      type1_raw_kendall = mean(sims[, "raw_p_kendall"] < 0.05),
      type1_clean_kendall = mean(sims[, "clean_p_kendall"] < 0.05)
    )
  }
  
  return(results)
}

results <- run_sim()

rates <- names(results)
make_summary_table <- function(results) {
  
  rates <- names(results)
  
  out <- do.call(rbind, lapply(rates, function(r) {
    c(
      contamination = as.numeric(r),
      
      # Pearson
      mean_raw_cor_pearson   = results[[r]]$mean_raw_cor_pearson,
      mean_clean_cor_pearson = results[[r]]$mean_clean_cor_pearson,
      type1_raw_pearson      = results[[r]]$type1_raw_pearson,
      type1_clean_pearson    = results[[r]]$type1_clean_pearson,
      
      # Spearman
      mean_raw_cor_spearman   = results[[r]]$mean_raw_cor_spearman,
      mean_clean_cor_spearman = results[[r]]$mean_clean_cor_spearman,
      type1_raw_spearman      = results[[r]]$type1_raw_spearman,
      type1_clean_spearman    = results[[r]]$type1_clean_spearman,
      
      # Kendall
      mean_raw_cor_kendall   = results[[r]]$mean_raw_cor_kendall,
      mean_clean_cor_kendall = results[[r]]$mean_clean_cor_kendall,
      type1_raw_kendall      = results[[r]]$type1_raw_kendall,
      type1_clean_kendall    = results[[r]]$type1_clean_kendall
    )
  }))
  
  return(as.data.frame(out))
}

summary_table <- make_summary_table(results)


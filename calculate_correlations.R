library(dplyr)

# ---------------------------------------------------------------
# 1. Outlier detection function using MAD
# ---------------------------------------------------------------
detect_outliers_mad <- function(x, thresh = 3) {
  x <- x[!is.na(x)]   # remove NAs if any
  
  med <- median(x)
  mad_val <- mad(x, constant = 1.4826)
  
  # If MAD is NA or zero, no outliers can be detected
  if (is.na(mad_val) || mad_val == 0) {
    return(integer(0))
  }
  
  z <- abs(x - med) / mad_val
  which(z > thresh)
}

# ---------------------------------------------------------------
# 2. Load generated datasets
# ---------------------------------------------------------------
df <- read.csv("generated_corr_contam_datasets.csv")

# split by dataset (iteration × rho × contam)
datasets <- split(df, list(df$iteration, df$rho, df$contam), drop = TRUE)

total_datasets <- length(datasets)
cat("Total datasets to process:", total_datasets, "\n")

# ---------------------------------------------------------------
# 3. Loop through datasets and compute correlations
# ---------------------------------------------------------------
results_list <- vector("list", total_datasets)
counter <- 1

for (name in names(datasets)) {
  
  dat <- datasets[[name]]
  x <- dat$x
  y <- dat$y
  
  iteration <- dat$iteration[1]
  rho <- dat$rho[1]
  contam <- dat$contam[1]
  
  # ---- Print progress ----
  cat("Processing dataset", counter, "/", total_datasets,
      "| iteration:", iteration,
      "| rho:", rho,
      "| contam:", contam, "\n")
  
  # RAW correlations
  raw_corr_pearson  <- cor(x, y)
  raw_corr_spearman <- 2 * sin(cor(x, y, method = "spearman") * pi / 6)
  raw_corr_kendall  <- sin(cor(x, y, method = "kendall") * pi /2)
  
  raw_p_pearson  <- cor.test(x, y)$p.value
  raw_p_spearman <- cor.test(x, y, method = "spearman")$p.value
  raw_p_kendall  <- cor.test(x, y, method = "kendall")$p.value
  
  # Outliers
  out <- union(detect_outliers_mad(x), detect_outliers_mad(y))
  
  if (length(out) > 0) {
    x_clean <- x[-out]
    y_clean <- y[-out]
  } else {
    x_clean <- x
    y_clean <- y
  }
  
  # CLEAN correlations and correct spearman and kendall
  clean_corr_pearson  <- cor(x_clean, y_clean)
  clean_corr_spearman <- 2 * sin(cor(x_clean, y_clean, method = "spearman") * pi / 6)
  clean_corr_kendall  <- sin(cor(x_clean, y_clean, method = "kendall") * pi / 2)
  
  clean_p_pearson  <- cor.test(x_clean, y_clean)$p.value
  clean_p_spearman <- cor.test(x_clean, y_clean, method = "spearman")$p.value
  clean_p_kendall  <- cor.test(x_clean, y_clean, method = "kendall")$p.value

  
  # Store results
  results_list[[counter]] <- data.frame(
    iteration = iteration,
    rho = rho,
    contam = contam,
    raw_corr_pearson,
    raw_corr_spearman,
    raw_corr_kendall,
    clean_corr_pearson,
    clean_corr_spearman,
    clean_corr_kendall,
    raw_p_pearson,
    raw_p_spearman,
    raw_p_kendall,
    clean_p_pearson,
    clean_p_spearman,
    clean_p_kendall
  )
  
  counter <- counter + 1
}




# ---------------------------------------------------------------
# 5. Combine and save results
# ---------------------------------------------------------------
results_df <- bind_rows(results_list)
write.csv(results_df, "correlation_results.csv", row.names = FALSE)

cat("\nFinished! Saved correlation_results.csv\n")

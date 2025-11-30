library(MASS)
library(dplyr)

# ---------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------
n <- 500   # number of observations per dataset

corr_levels <- c(0, 0.5, 0.6, 0.7, 0.8, 0.9)
cont_levels <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25,
                 0.30, 0.35, 0.40, 0.45, 0.50, 0.55)

total_iter <- 500

total_datasets <- total_iter * length(corr_levels) * length(cont_levels)
cat("Expected total datasets:", total_datasets, "\n\n")

# ---------------------------------------------------------------
# STORAGE FOR ALL GENERATED DATASETS
# ---------------------------------------------------------------
results_list <- vector("list", total_datasets)
counter <- 1

# ---------------------------------------------------------------
# MAIN GENERATION LOOP
# ---------------------------------------------------------------
for (i in 1:total_iter) {
  for (rho in corr_levels) {
    for (contam in cont_levels) {
      
      # ---- Print progress ----
      cat("Generating dataset",
          counter, "/", total_datasets,
          "| iteration:", i,
          "| rho:", rho,
          "| contam:", contam, "\n")
      
      # ---------------------------------------------------------
      # 1. Generate clean bivariate normal data
      # ---------------------------------------------------------
      Sigma <- matrix(c(1, rho, rho, 1), 2)
      z <- mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
      
      x <- z[, 1]
      y <- z[, 2]
      
      # ---------------------------------------------------------
      # 2. Inject contamination
      # ---------------------------------------------------------
      outl_num <- floor(n * contam)
      if (outl_num > 0) {
        idx <- sample(1:n, outl_num)
        x[idx] <- rnorm(outl_num, 8, 1)
        y[idx] <- rnorm(outl_num, 8, 1)
      }
      
      # ---------------------------------------------------------
      # 3. Save dataset in tidy format
      # ---------------------------------------------------------
      results_list[[counter]] <- data.frame(
        dataset_id = counter,
        iteration = i,
        rho = rho,
        contam = contam,
        n = n,
        obs = 1:n,
        x = x,
        y = y
      )
      
      counter <- counter + 1
    }
  }
}

# ---------------------------------------------------------------
# COMBINE ALL DATASETS INTO ONE DATAFRAME
# ---------------------------------------------------------------
results_df <- bind_rows(results_list)

# Check dataset count before writing
cat("\nGenerated datasets:", 
    results_df %>% summarise(datasets = n_distinct(dataset_id)) %>% pull(datasets),
    "\nExpected:", total_datasets, "\n")

# ---------------------------------------------------------------
# SAVE AS CSV
# ---------------------------------------------------------------
write.csv(results_df, "generated_corr_contam_datasets.csv", row.names = FALSE)

cat("\n✔ Generation complete!\n")
cat("✔ File saved as: generated_corr_contam_datasets.csv\n")

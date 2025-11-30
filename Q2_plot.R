library(dplyr)
library(ggplot2)

# ---- STEP 1: Read CSV ----
df_raw <- read.csv("results_flat.csv")

# ---- STEP 2: Parse n and contamination level from key ----
df <- df_raw %>%
  mutate(
    n    = as.numeric(sub("^n([0-9]+)_.*", "\\1", key)),
    cont = as.numeric(sub(".*cont", "", key)),
    raw  = type1_raw_pearson,
    clean = type1_clean_pearson
  ) %>%
  select(n, cont, raw, clean)

# ---- STEP 3: Keep only rows with n = 1000 ----
df_1000 <- df %>% filter(n == 1000)

# ---- STEP 4: Convert to long format for ggplot ----
df_long <- df_1000 %>%
  pivot_longer(cols = c(raw, clean),
               names_to = "metric",
               values_to = "value")

# ---- STEP 5: Bar plot ----
ggplot(df_long, aes(x = factor(cont), y = value, fill = metric)) +
  geom_col(position = position_dodge()) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  labs(
    title = "Type I Error at n = 1000 for Each Contamination Level",
    x = "Contamination Level",
    y = "Type I Error",
    fill = "Metric"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw()

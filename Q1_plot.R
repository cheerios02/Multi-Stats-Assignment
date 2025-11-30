library(dplyr)
library(tidyr)
library(ggplot2)

# ---- STEP 0: Read CSV with Option 1 format ----
df_raw <- read.csv("results_flat.csv")

# ---- STEP 1: Parse n and contamination level from `key` ----
df <- df_raw %>%
  mutate(
    n    = as.numeric(sub("^n([0-9]+)_.*", "\\1", key)),
    cont = as.numeric(sub(".*cont", "", key)),
    raw  = mean_cor_raw_pearson,
    clean = mean_cor_clean_pearson
  ) %>%
  select(n, cont, raw, clean)

# ---- STEP 2: Build wide-format tables (indexed by contamination) ----
table_raw <- df %>%
  select(cont, n, raw) %>%
  arrange(cont, n) %>%
  pivot_wider(
    names_from = n,
    values_from = raw
  )

table_clean <- df %>%
  select(cont, n, clean) %>%
  arrange(cont, n) %>%
  pivot_wider(
    names_from = n,
    values_from = clean
  )

# ---- STEP 3: Plot one graph per contamination level ----
ggplot(df, aes(x = n)) +
  geom_line(aes(y = raw, color = "Raw")) +
  geom_point(aes(y = raw, color = "Raw")) +
  geom_line(aes(y = clean, color = "Clean")) +
  geom_point(aes(y = clean, color = "Clean")) +
  facet_wrap(~ cont, scales = "free_y") +
  labs(
    title = "Type I Error vs n for Each Contamination Level",
    x = "n",
    y = "Type I Error",
    color = "Metric"
  ) +
  theme_bw()

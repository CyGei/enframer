library(tidyr)
library(arrow)
library(dplyr)

df <- read_parquet("data/RSV_ALL.parquet") |>
  transmute(
    value = as.numeric(value),
    row = as.integer(fluseason_week),
    col = factor(location_code, levels = sort(unique(location_code))),
    dim1 = as.integer(fluseason),
    dim2 = factor(datasetH1, levels = sort(unique(datasetH1))),
    dim3 = factor(datasetH2, levels = sort(unique(datasetH2))),
    dim4 = as.character(sample)
  )
col_levels <- levels(df$col)

df_complete <- df |>
  group_by(dim1, dim2, dim3, dim4) |>
  tidyr::complete(
    row = 1:length(unique(df$row)),
    col = col_levels,
    fill = list(value = NA_real_)
  ) |>
  ungroup()
df_complete


df_padded <- df_complete |>
  group_by(across(-row)) |>
  arrange(as.integer(row)) |>
  mutate(
    # Identify first and last non-NA week
    first_obs = min(row[!is.na(value)], na.rm = TRUE),
    last_obs = max(row[!is.na(value)], na.rm = TRUE),
    # External gaps → 0
    value = case_when(
      as.integer(row) < first_obs ~ 0,
      as.integer(row) > last_obs ~ 0,
      TRUE ~ value
    ),
    # Internal gaps → LOCF
    value = zoo::na.locf(value, na.rm = FALSE)
  ) |>
  ungroup()
df_padded

df_padded

install.packages("dlookr")
library(dlookr)
library(tidyr)
library(arrow)
library(dplyr)

df <- read_parquet("data/RSV_ALL.parquet") |>
  transmute(
    value = as.numeric(value),
    row = factor(
      fluseason_week,
      levels = as.character(sort(unique(fluseason_week)))
    ),
    col = factor(location_code, levels = sort(unique(location_code))),
    dim1 = as.integer(fluseason),
    dim2 = factor(datasetH1, levels = sort(unique(datasetH1))),
    dim3 = factor(datasetH2, levels = sort(unique(datasetH2))),
    dim4 = as.character(sample)
  )
row_levels <- levels(df$row)
col_levels <- levels(df$col)

install.packages("missRanger")
library(missRanger)

df_imputed <- df |>
  group_by(dim1, dim2, dim3, dim4) |>
  tidyr::complete(
    row = row_levels,
    col = col_levels,
    fill = list(value = NA_real_)
  ) |>
  ungroup() |>
  missRanger(
    formula = value ~ .,
    num.trees = 100,
    seed = 123
  )

library(ggplot2)
df_imputed |>
  filter(
    dim1 == 2023,
    dim2 == "RSV-Net",
    dim3 == "RSV-Net",
    dim4 == "1"
  ) |>
  ggplot(aes(x = row, y = col)) +
  geom_tile(aes(fill = value))

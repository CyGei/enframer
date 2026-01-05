# ------------------------------------
#           Data Load
# ------------------------------------
library(purrr)
library(tidytable)
library(arrow)
library(ggplot2)

df <- read_parquet("data/RSV_ALL.parquet") |>
  as_tidytable() |>
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
# Rows are weeks, Columns are locations
# dup_check <- df |>
#  janitor::get_dupes(dim1, dim2, dim3, dim4, row, col)

# ------------------------------------
#           Scaling Distribution
# ------------------------------------
# The scaling distribution is derived from the peak weekly sums of the
# "RSV_SMH" dataset across all seasons, datasets, and samples.
# This distribution is then used to scale other datasets to have
# similar peak weekly sums, ensuring comparability across datasets.

scaling_dist <- df |>
  filter(dim2 == "RSV_SMH") |>
  summarise(
    weekly_sum = sum(value, na.rm = TRUE),
    .by = c("dim1", "dim3", "dim4", "row")
  ) |>
  slice_max(
    weekly_sum,
    n = 1,
    with_ties = FALSE,
    .by = c("dim1", "dim3", "dim4")
  ) |>
  pull(weekly_sum)
hist(scaling_dist)

framelist <- df |>
  group_by(dim1, dim2, dim3, dim4) |>
  summarise(
    mat = list({
      m <- matrix(
        NA_real_,
        nrow = length(row_levels),
        ncol = length(col_levels),
        dimnames = list(row_levels, col_levels)
      )
      row <- factor(row, levels = row_levels)
      col <- factor(col, levels = col_levels)
      m[cbind(as.integer(row), as.integer(col))] <- value
      m
    }),
    .groups = "drop"
  ) |>
  mutate(
    mat = map(
      mat,
      ~ {
        peak <- max(rowSums(.x, na.rm = TRUE), na.rm = TRUE)
        if (is.na(peak) || peak == 0) {
          return(.x)
        }
        target_max <- sample(scaling_dist, size = 1)
        scaling_factor <- target_max / peak
        .x * scaling_factor
      }
    )
  )

# ------------------------------------
#           Frame visualisation
# ------------------------------------
plot_frame <- function(frame) {
  plot_df <- frame |>
    as.data.frame() |>
    tibble::rownames_to_column("row") |>
    pivot_longer(-row, names_to = "col", values_to = "value") |>
    mutate(
      row = factor(row, levels = row_levels),
      col = factor(col, levels = col_levels)
    )

  ggplot(plot_df, aes(x = row, y = col, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal()
}

pdf("all_frames.pdf", width = 8, height = 6)

framelist |>
  group_by(dim2) |>
  filter(
    dim2 != "RSV_SMH" |
      row_number() %in% sample.int(n(), size = min(10, n()))
  ) |>
  ungroup() |>
  pwalk(function(dim1, dim2, dim3, dim4, mat, mat_scaled) {
    p <- plot_frame(mat) +
      ggtitle(paste(
        "Year:",
        dim1,
        "\n",
        dim2,
        dim3,
        dim4
      ))
    print(p)
  })

dev.off()


# ------------------------------------
#           Filter Methods
# ------------------------------------
# Rows are weeks, Columns are locations

# Define your filters
filters <- list(
  # Each column must have at least X % of data observed
  col_prop_all = function(mat, prop = 0.5) {
    all(apply(mat, 2, function(v) mean(!is.na(v))) >= prop)
  },
  # At least one column must have at least X % of data observed
  col_prop_any = function(mat, prop = 0.5) {
    any(apply(mat, 2, function(v) mean(!is.na(v))) >= prop)
  },
  # Must observe at least X consecutive weeks per column
  col_rle = function(mat, len = 10) {
    all(apply(mat, 2, function(v) {
      rle_v <- rle(!is.na(v))
      any(rle_v$values & rle_v$lengths >= len)
    }))
  }
)

framelist_filtered <- map_dfc(
  filters,
  ~ map_lgl(framelist$mat, .x)
) |>
  (\(x) bind_cols(framelist, x))()


# ------------------------------------
#           Imputation Methods
# ------------------------------------
library(softImpute)
library(imputeTS)
library(missForest)
library(mice)

imputers <- list(
  zeros = function(mat) {
    mat[is.na(mat)] <- 0
    mat
  },
  linear = function(mat) {
    apply(mat, 2, imputeTS::na_interpolation, option = "linear")
  },
  kalman = function(mat) {
    apply(mat, 2, imputeTS::na_kalman) |> as.matrix()
  },
  softimpute = function(mat) {
    fit <- softImpute::softImpute(mat, rank.max = 10)
    softImpute::complete(mat, fit)
  },
  missforest = function(mat) {
    missForest::missForest(mat)$ximp
  }
)

x = framelist_filtered |> filter(col_rle) |> slice_sample(n = 10)

framelist_imputed <- map_dfc(names(imputers), function(name) {
  mat_list <- map(framelist_filtered$mat, imputers[[name]])
  tibble::tibble(!!name := mat_list)
}) |>
  (\(x) bind_cols(framelist_filtered, x))()

plot_frame <- function(frame) {
  plot_df <- frame |>
    as.data.frame() |>
    tibble::rownames_to_column("row") |>
    pivot_longer(-row, names_to = "col", values_to = "value") |>
    mutate(
      row = factor(row, levels = row_levels),
      col = factor(col, levels = col_levels)
    )

  ggplot(plot_df, aes(x = row, y = col, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "midnightblue",
      mid = "white",
      high = "firebrick",
      midpoint = 0,
    ) +
    theme_classic()
}

framelist_imputed |>
  slice(1) |>
  pluck("mat", 1) |>
  plot_frame() +
  ggtitle("Original")

map(
  names(imputers),
  function(name) {
    framelist_imputed |>
      slice(1) |>
      pluck(name, 1) |>
      plot_frame() +
      ggtitle(paste("Imputation method:", name))
  }
)

# ------------------------------------
#           Data Load
# ------------------------------------
library(purrr)
library(tibble)
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
#       Imputation Pipeline
# ------------------------------------
library(imputeTS)
library(softImpute)
library(missForest)
library(mice)

# --------- dianostics ---------
frame_diagnostics <- function(mat) {
  tibble(
    n_obs_total = sum(!is.na(mat)),
    min_obs_per_col = min(colSums(!is.na(mat))),
    max_run_per_col = min(apply(mat, 2, function(v) {
      r <- rle(!is.na(v))
      if (any(r$values)) max(r$lengths[r$values]) else 0
    })),
    prop_obs_min_col = min(colMeans(!is.na(mat)))
  )
}

# --------- Imputer Definitions ---------
nonneg <- function(mat) {
  mat[mat < 0] <- 0
  mat
}

imputers <- list(
  zeros = list(
    requires = function(d) TRUE,
    impute = function(mat) {
      mat[is.na(mat)] <- 0
      mat
    }
  ),
  linear = list(
    requires = function(d) d$min_obs_per_col >= 2,
    impute = function(mat) {
      apply(mat, 2, na_interpolation, option = "linear") |> nonneg()
    }
  ),
  kalman = list(
    requires = function(d) d$min_obs_per_col >= 2,
    impute = function(mat) apply(mat, 2, na_kalman) |> nonneg()
  ),
  softimpute = list(
    requires = function(d) d$n_obs_total > 0,
    impute = function(mat) {
      fit <- softImpute::softImpute(mat, rank.max = 10)
      softImpute::complete(mat, fit) |> nonneg()
    }
  ),
  missforest = list(
    requires = function(d) d$min_obs_per_col >= 1,
    impute = function(mat) missForest::missForest(mat)$ximp |> nonneg()
  )
)

activate_imputers <- function(d) {
  map_lgl(imputers, \(imp) imp$requires(d))
}


apply_imputers <- function(mat, active) {
  imap(imputers, \(imp, name) {
    if (!active[[name]]) {
      return(NULL)
    }
    imp$impute(mat)
  })
}


impute_framelist <- function(framelist) {
  framelist |>
    mutate(
      diagnostics = map(mat, frame_diagnostics),
      active_imputers = map(diagnostics, activate_imputers),
      imputed = map2(mat, active_imputers, apply_imputers)
    )
}

# --------- usage ---------
framelist_imputed <- impute_framelist(framelist)

framelist_imputed |>
  slice(n = 1) |>
  pluck("imputed", 1) |>
  pluck("softimpute") |>
  plot_frame()


framelist_imputed |>
  slice(3113) |>
  pluck("imputed", 1) |>
  imap(\(mat, method) {
    # plot only non-NULL imputations
    if (is.null(mat)) {
      return(NULL)
    } else {
      plot_frame(mat) +
        ggtitle(paste("Imputation:", method))
    }
  }) |>
  compact() |>
  patchwork::wrap_plots(guides = 'collect')

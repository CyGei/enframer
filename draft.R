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
    x = factor(
      fluseason_week,
      levels = as.character(sort(unique(fluseason_week)))
    ),
    y = factor(location_code, levels = sort(unique(location_code))),
    dim1 = as.integer(fluseason),
    dim2 = factor(datasetH1, levels = sort(unique(datasetH1))),
    dim3 = factor(datasetH2, levels = sort(unique(datasetH2))),
    dim4 = as.character(sample)
  )
x_levels <- levels(df$x)
y_levels <- levels(df$y)
# Rows (x) are weeks, Columns (y) are locations
# dup_check <- df |>
#  janitor::get_dupes(dim1, dim2, dim3, dim4, x, y)

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
    .by = c("dim1", "dim3", "dim4", "x")
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
        nrow = length(x_levels),
        ncol = length(y_levels),
        dimnames = list(x_levels, y_levels)
      )
      x <- factor(x, levels = x_levels)
      y <- factor(y, levels = y_levels)
      m[cbind(as.integer(x), as.integer(y))] <- value
      m
    }),
    .groups = "drop"
  ) |>
  mutate(
    mat_scaled = map(
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
    tibble::rownames_to_column("x") |>
    pivot_longer(-x, names_to = "y", values_to = "value") |>
    mutate(
      x = factor(x, levels = x_levels),
      y = factor(y, levels = y_levels)
    )

  ggplot(plot_df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal()
}

pdf("all_frames.pdf", width = 8, height = 6)

framelist |>
  filter(dim2 != "RSV_SMH") |>
  pwalk(function(dim1, dim2, dim3, dim4, mat, mat_scaled) {
    p <- plot_frame(mat_scaled) +
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
#        Imputation Methods
# ------------------------------------
# See imputator.R for full evaluation framework

library(softImpute)
library(imputeTS)
library(missForest)
library(mice)

imputelist <- framelist |>
  slice_sample(n = 10) |>
  mutate(
    softImpute = map(
      mat,
      ~ {
        # Uses correlation across locations (columns) AND time (rows)
        fits <- softImpute(.x, rank.max = 5, type = "als", trace = FALSE)
        imputed_mat <- complete(.x, fits)
        imputed_mat[imputed_mat < 0] <- 0
        imputed_mat
      }
    ),
    kalman = map(
      mat,
      ~ {
        # Treats every location independently
        apply(.x, 2, function(x) {
          if (all(is.na(x))) {
            return(x)
          }
          na_kalman(x)
        })
      }
    ),
    linear = map(
      mat,
      ~ t(apply(.x, 2, na_interpolation, option = "linear"))
    ),
    missForest = map(
      mat,
      ~ {
        if (sum(!is.na(.x)) / length(.x) > 0.1) {
          as.matrix(missForest(.x, verbose = FALSE)$ximp)
        } else {
          .x
        }
      }
    )
  )

library(purrr)
library(tidytable)
library(softImpute)
library(imputeTS)
library(missForest)

imputelist <- framelist |>
  slice_sample(n = 10) |>
  mutate(
    softImpute = map(
      mat,
      ~ {
        if (sum(!is.na(.x)) / length(.x) < 0.1) {
          return(.x)
        } # skip if too sparse
        fits <- softImpute(.x, rank.max = 5, type = "als", trace = FALSE)
        imputed <- complete(.x, fits)
        imputed[is.na(imputed)] <- 0
        imputed[imputed < 0] <- 0
        imputed
      }
    ),
    kalman = map(
      mat,
      ~ {
        res <- apply(.x, 2, function(x) {
          if (all(is.na(x))) {
            return(rep(0, length(x)))
          }
          na_kalman(x)
        })
        res[res < 0] <- 0
        res
      }
    ),
    linear = map(
      mat,
      ~ {
        res <- t(apply(.x, 2, function(x) {
          if (all(is.na(x))) {
            return(rep(0, length(x)))
          }
          na_interpolation(x, option = "linear")
        }))
        res[res < 0] <- 0
        res
      }
    ),
    missForest = map(
      mat,
      ~ {
        if (sum(!is.na(.x)) / length(.x) < 0.1) {
          return(.x)
        }
        res <- missForest(.x, verbose = FALSE)$ximp
        res[res < 0] <- 0
        res
      }
    )
  )


X <- framelist |>
  slice(1) |>
  pluck("mat", 1)

A = imputelist |>
  slice(1) |>
  pluck("softImpute", 1)
B = imputelist |>
  slice(1) |>
  pluck("missForest", 1)

B |>
  as.data.frame() |>
  tibble::rownames_to_column("x") |>
  pivot_longer(-x, names_to = "y", values_to = "value") |>
  mutate(
    x = factor(x, levels = rownames(X)),
    y = factor(y, levels = colnames(X))
  ) |>
  ggplot(aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

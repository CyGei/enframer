build_global_lookup <- function(all_datasets_df) {
  all_datasets_df |>
    group_by(col, dim1, dim2, dim4) %>%
    nest() %>%
    group_by(col) %>%
    nest(.key = "col_data")
}

x = build_global_lookup(df)
x[1, ]$col_data

library(zoo)
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
      apply(mat, 2, function(v) na.approx(v, maxgap = 2, rule = 2)) |> nonneg()
    }
  ),
  kalman = list(
    requires = function(d) d$min_obs_per_col >= 2,
    impute = function(mat) {
      apply(mat, 2, na_kalman, model = "StructTS") |> nonneg()
    }
  ),
  softimpute = list(
    requires = function(d) d$n_obs_total > 0,
    impute = function(mat) {
      fit <- softImpute::softImpute(mat, rank.max = 10, lambda = 0)
      softImpute::complete(mat, fit) |> nonneg()
    }
  ),
  missforest = list(
    requires = function(d) d$min_obs_per_col >= 1,
    impute = function(mat) {
      missForest::missForest(mat, verbose = FALSE)$ximp |> nonneg()
    }
  ),
  intelligent_fill = list(
    requires = function(d) d$n_obs_total > 0,
    impute = function(mat) {
      filled_mat <- mat

      for (col in 1:ncol(mat)) {
        series <- mat[, col]

        if (all(is.na(series))) {
          series <- rep(0, length(series))
          filled_mat[, col] <- series
          next
        }

        obs_idx <- which(!is.na(series))
        first_obs <- obs_idx[1]
        last_obs <- obs_idx[length(obs_idx)]
        if (first_obs > 1) {
          series[1:(first_obs - 1)] <- 0
        }
        if (last_obs < length(series)) {
          series[(last_obs + 1):length(series)] <- 0
        }

        # Use imputeTS::na_locf without na.rm argument
        series <- imputeTS::na_locf(series)
        series <- imputeTS::na_locf(series, option = "nocb") # nocb = next obs carried backward
        series <- pmax(series, 0)

        filled_mat[, col] <- series
      }

      filled_mat
    }
  )
)

# --------- usage ---------
framelist_imputed <- impute_framelist(framelist)

framelist_imputed |>
  slice(n = 1) |>
  pluck("imputed", 1) |>
  pluck("intelligent_fill") |>
  plot_frame()

framelist_imputed |>
  slice(3114) |>
  pluck("imputed", 1) |>
  imap(\(mat, method) {
    if (is.null(mat)) {
      return(NULL)
    }
    plot_frame(mat) + ggtitle(paste("Imputation:", method))
  }) |>
  compact() |>
  patchwork::wrap_plots(guides = 'collect')


build_location_lookup <- function(framelist) {
  framelist |>
    mutate(
      peak = map_dbl(mat, ~ max(rowSums(.x, na.rm = TRUE), na.rm = TRUE))
    ) |>
    filter(peak > 0) |>
    group_by(col) |>
    summarise(
      candidates = list(cur_data()),
      .groups = "drop"
    )
}
location_lookup <- build_location_lookup(framelist)


library(VIM)
mat <- framelist$mat[[1]]
mat
kNN(mat) |> as.matrix()
install.packages("simputation")
library(simputation)
mat |>
  impute_lm()
impute_cart(Species ~ .) |> # use all variables except 'Species' as predictor
  head(10)


library(tidyr)

df_complete <- df |>
  tidyr::complete(
    row = row_levels,
    col = col_levels,
    dim1,
    dim2,
    dim3,
    dim4,
    fill = list(value = NA_real_)
  )

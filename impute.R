library(purrr)
library(tidytable)
library(softImpute)
library(imputeTS)
library(missForest)

# Define core imputation methods (reduced set for conciseness)
imputation_methods <- list(
  softImpute = function(mat) {
    if (sum(!is.na(mat)) < 10) {
      return(mat)
    }
    fits <- softImpute(
      mat,
      rank.max = min(5, ncol(mat) - 1),
      type = "als",
      lambda = 0
    )
    imputed <- complete(mat, fits)
    imputed[imputed < 0] <- 0
    imputed
  },
  kalman = function(mat) {
    apply(mat, 2, function(col) {
      n_non_na <- sum(!is.na(col))
      if (n_non_na < 3) {
        if (n_non_na == 0) {
          return(rep(0, length(col)))
        } else {
          fill_val <- mean(col, na.rm = TRUE)
          col[is.na(col)] <- fill_val
          col
        }
      } else {
        imputed <- na_kalman(col, model = "auto.arima", smooth = TRUE)
        imputed[imputed < 0] <- 0
        imputed
      }
    })
  },
  spline = function(mat) {
    apply(mat, 2, function(col) {
      n_non_na <- sum(!is.na(col))
      if (n_non_na < 3) {
        if (n_non_na == 0) {
          return(rep(0, length(col)))
        } else {
          fill_val <- mean(col, na.rm = TRUE)
          col[is.na(col)] <- fill_val
          col
        }
      } else {
        imputed <- na_interpolation(col, option = "spline")
        imputed[imputed < 0] <- 0
        imputed
      }
    })
  },
  missForest = function(mat) {
    if (mean(!is.na(mat)) < 0.1) {
      return(mat)
    }
    imputed <- missForest(mat, verbose = FALSE, maxiter = 5)$ximp
    imputed[imputed < 0] <- 0
    imputed
  },
  mean_location = function(mat) {
    apply(mat, 2, function(col) {
      mean_val <- mean(col, na.rm = TRUE)
      col[is.na(col)] <- if (is.na(mean_val)) 0 else mean_val
      col
    })
  }
)

# Create validation set
create_validation_sets <- function(mat, missing_prop = 0.2, seed = 123) {
  set.seed(seed)
  observed_idx <- which(!is.na(mat))
  n_to_mask <- floor(length(observed_idx) * missing_prop)
  mask_idx <- sample(observed_idx, n_to_mask)
  true_values <- mat[mask_idx]
  mat_masked <- mat
  mat_masked[mask_idx] <- NA
  list(masked_mat = mat_masked, true_values = true_values, mask_idx = mask_idx)
}

# Metrics calculation
calculate_metrics <- function(true_values, imputed_values) {
  valid_idx <- !is.na(imputed_values)
  true_values <- true_values[valid_idx]
  imputed_values <- imputed_values[valid_idx]
  if (length(true_values) == 0) {
    return(list(rmse = NA, mae = NA, mape = NA, r2 = NA))
  }
  residuals <- true_values - imputed_values
  list(
    rmse = sqrt(mean(residuals^2)),
    mae = mean(abs(residuals)),
    mape = mean(abs(residuals) / (true_values + 1e-6)) * 100,
    r2 = cor(true_values, imputed_values)^2
  )
}

# Evaluate methods on samples
evaluate_imputation <- function(
  framelist,
  missing_props = c(0.1, 0.2, 0.3)
) {
  sampled <- framelist
  map_dfr(seq_len(nrow(sampled)), function(i) {
    map_dfr(missing_props, function(miss_prop) {
      val <- create_validation_sets(sampled$mat[[i]], miss_prop)
      imputed_mats <- map(imputation_methods, ~ .x(val$masked_mat))
      map_dfr(names(imputed_mats), function(method_name) {
        imputed_values <- imputed_mats[[method_name]][val$mask_idx]
        metrics <- calculate_metrics(val$true_values, imputed_values)
        tidytable(
          dim1 = sampled$dim1[i],
          dim2 = sampled$dim2[i],
          dim3 = sampled$dim3[i],
          dim4 = sampled$dim4[i],
          missing_prop = miss_prop,
          method = method_name,
          rmse = metrics$rmse,
          mae = metrics$mae,
          mape = metrics$mape,
          r2 = metrics$r2
        )
      })
    })
  })
}

# Run evaluation and summarize
set.seed(123)
evaluation_results <- evaluate_imputation(
  framelist |>
    filter(dim2 != "RSV_SMH")
)

evaluation_results |>
  pivot_longer(
    cols = c(rmse, mae, mape, r2),
    names_to = "metric",
    values_to = "value"
  ) |>
  ggplot(aes(x = as.factor(missing_prop), y = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(~metric, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "top") +
  labs(
    x = "Prop missing data",
    y = "Metric value",
    caption = "missing data represents the proportion of data removed for validation"
  )

summary_results <- evaluation_results |>
  summarise(
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_mape = mean(mape, na.rm = TRUE),
    mean_r2 = mean(r2, na.rm = TRUE),
    .by = c(method, missing_prop)
  ) |>
  arrange(missing_prop, mean_rmse)

print(summary_results)

# Select and apply best method (for 20% missing)
best_method <- summary_results |>
  filter(missing_prop == 0.2) |>
  slice_min(mean_rmse, n = 1) |>
  pull(method)

cat("Best performing method:", best_method, "\n")
cat("Applying", best_method, "to all matrices...\n")

framelist_imputed <- framelist |>
  mutate(
    mat_imputed = map(
      mat,
      ~ {
        m <- imputation_methods[[best_method]](.x)
        m[is.na(m)] <- 0
        m[m < 0] <- 0
        m
      }
    )
  )

cat(
  "Imputation complete! No NAs remaining:",
  !any(map_lgl(framelist_imputed$mat_imputed, ~ any(is.na(.x)))),
  "\n"
)
x = framelist$mat_scaled[[8]][, 40]
ggplot_na_distribution()
y = na_interpolation(x)
ggplot_na_imputations(x, y)

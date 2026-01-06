# ----------------------------
#        Diagnostics
# ----------------------------
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

# ----------------------------
#        Imputer Definitions & Conditions
# ----------------------------
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

# ----------------------------
#        Activate Imputers
# ----------------------------
activate_imputers <- function(d) {
  map_lgl(imputers, \(imp) imp$requires(d))
}

# ----------------------------
#        Apply Imputers
# ----------------------------
apply_imputers <- function(mat, active) {
  imap(imputers, \(imp, name) {
    if (!active[[name]]) {
      return(NULL)
    }
    imp$impute(mat)
  })
}

# ----------------------------
#        Main Preprocessing Pipeline
# ----------------------------
preprocess_framelist <- function(framelist) {
  framelist |>
    mutate(
      diagnostics = map(mat, frame_diagnostics),
      active_imputers = map(diagnostics, activate_imputers),
      imputed = map2(mat, active_imputers, apply_imputers)
    )
}

# ----------------------------
# Example usage
# ----------------------------
framelist_imputed <- preprocess_framelist(framelist)
# Check which imputers were applied for :
framelist_imputed |>
  filter(dim1 == 2024 & dim2 == "NHSN") |>
  pull(active_imputers)

framelist_imputed |>
  filter(dim1 == 2025 & dim2 == "RSV-Net") |>
  pull(active_imputers)

## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(cmdstanr, posterior, loo, dplyr, tibble)

source("scripts/model_load_data.R")

group_metadata <- tibble::tribble(
  ~group_code, ~group_id, ~group_label,
  "Rev", "nominal_extinction", "Nominal extinction",
  "RevFB", "reversal_with_feedback", "Reversal with feedback"
)

model_metadata <- tibble::tribble(
  ~model_code, ~model_id, ~model_label,
  "model_based", "model_based_revaluation", "Model-based revaluation",
  "model_free", "model_free_learning", "Model-free learning"
)

phase_lookup <- tibble(
  Epoch = seq_len(12),
  Phase = c(rep("Acquisition", 8), rep("Reversal", 4))
)

model_files <- list(
  model_based = "models/model_based_rev.stan",
  model_free = "models/model_free_rev.stan"
)

stan_models <- lapply(model_files, function(path) {
  cmdstan_model(path, force_recompile = TRUE)
})

sampling_args <- list(
  seed = 2025,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2500,
  adapt_delta = 0.95,
  refresh = 100
)

build_init <- function(group_data, model_name) {
  J <- group_data$J

  if (identical(model_name, "model_based")) {
    return(function(chain_id = 1) {
      list(
        mu = rnorm(6, 0, 0.5),
        sigma = abs(rnorm(6, 0, 0.2)) + 1e-3,
        L_phi = diag(2),
        Z_beta = matrix(rnorm(6 * J, 0, 0.3), nrow = 6, ncol = J)
      )
    })
  }

  function(chain_id = 1) {
    list(
      mu = rnorm(5, 0, 0.5),
      sigma = abs(rnorm(5, 0, 0.2)) + 1e-3,
      L_phi = diag(2),
      Z_beta = matrix(rnorm(5 * J, 0, 0.3), nrow = 5, ncol = J)
    )
  }
}

compute_loo_safe <- function(fit) {
  base_loo <- fit$loo(save_psis = TRUE)
  pareto_k <- loo::pareto_k_values(base_loo)

  if (!any(is.finite(pareto_k) & pareto_k > 0.7)) {
    return(list(loo = base_loo, method = "psis"))
  }

  mm_loo <- tryCatch(
    fit$loo(save_psis = TRUE, moment_match = TRUE),
    error = function(e) NULL
  )

  if (!is.null(mm_loo)) {
    return(list(loo = mm_loo, method = "psis+moment_match"))
  }

  list(loo = base_loo, method = "psis", warning = "moment_match_failed")
}

extract_loo_summary <- function(group_name, model_name, loo_result, method_name) {
  ks <- loo::pareto_k_values(loo_result)
  group_row <- group_metadata[group_metadata$group_code == group_name, ]
  model_row <- model_metadata[model_metadata$model_code == model_name, ]

  tibble(
    group = group_name,
    group_id = group_row$group_id,
    group_label = group_row$group_label,
    model = model_name,
    model_id = model_row$model_id,
    model_label = model_row$model_label,
    loo_method = method_name,
    elpd_loo = unname(loo_result$estimates["elpd_loo", "Estimate"]),
    se_elpd_loo = unname(loo_result$estimates["elpd_loo", "SE"]),
    p_loo = unname(loo_result$estimates["p_loo", "Estimate"]),
    looic = unname(loo_result$estimates["looic", "Estimate"]),
    n_bad_k = sum(is.finite(ks) & ks > 0.7),
    max_k = suppressWarnings(max(ks, na.rm = TRUE))
  )
}

target_groups <- names(dStan)
target_models <- names(stan_models)

fits <- list()
loos <- list()
summaries <- list()
comparisons <- list()
posterior_predictions <- list()
posterior_vmac <- list()

for (group_name in target_groups) {
  fits[[group_name]] <- list()
  loos[[group_name]] <- list()
  group_row <- group_metadata[group_metadata$group_code == group_name, ]

  for (model_name in target_models) {
    model_row <- model_metadata[model_metadata$model_code == model_name, ]
    message("Fitting ", group_row$group_label, " (", group_name, ") / ", model_row$model_label)

    fit <- do.call(
      stan_models[[model_name]]$sample,
      c(
        list(
          data = dStan[[group_name]],
          init = build_init(dStan[[group_name]], model_name)
        ),
        sampling_args
      )
    )

    loo_info <- compute_loo_safe(fit)

    fits[[group_name]][[model_name]] <- fit
    loos[[group_name]][[model_name]] <- loo_info$loo
    summaries[[paste(group_name, model_name, sep = "__")]] <-
      extract_loo_summary(group_name, model_name, loo_info$loo, loo_info$method)
    posterior_predictions[[paste(group_name, model_name, sep = "__")]] <-
      get_PPC_epochsY(
        fit = fit,
        group_size = dStan[[group_name]]$I,
        epoch_len = 48,
        credMass = 0.95
      ) %>%
      mutate(
        group = group_name,
        group_id = group_row$group_id,
        group_label = group_row$group_label,
        model = model_name,
        model_id = model_row$model_id,
        model_label = model_row$model_label,
        Singleton = ifelse(Distractor == "Distractor1", "High", "Low")
      ) %>%
      left_join(phase_lookup, by = "Epoch")
    posterior_vmac[[paste(group_name, model_name, sep = "__")]] <-
      get_PPC_diff_epochs(
        fit = fit,
        group_size = dStan[[group_name]]$I,
        epoch_len = 48,
        credMass = 0.95
      ) %>%
      mutate(
        group = group_name,
        group_id = group_row$group_id,
        group_label = group_row$group_label,
        model = model_name,
        model_id = model_row$model_id,
        model_label = model_row$model_label
      ) %>%
      left_join(phase_lookup, by = "Epoch")
  }

  if (all(c("model_based", "model_free") %in% names(loos[[group_name]]))) {
    comparisons[[group_name]] <- loo_compare(
      loos[[group_name]][["model_based"]],
      loos[[group_name]][["model_free"]]
    )
  }
}

loo_summary <- bind_rows(summaries)
posterior_predictions_df <- bind_rows(posterior_predictions)
posterior_vmac_df <- bind_rows(posterior_vmac)

comparison_summary <- bind_rows(lapply(names(comparisons), function(group_name) {
  comparison_df <- as.data.frame(comparisons[[group_name]])
  comparison_df$comparison_row <- rownames(comparison_df)
  rownames(comparison_df) <- NULL

  if (nrow(comparison_df) < 2) {
    stop(sprintf("Model comparison for group '%s' does not contain the expected two rows.", group_name))
  }

  free_row <- comparison_df[2, , drop = FALSE]
  group_row <- group_metadata[group_metadata$group_code == group_name, , drop = FALSE]

  delta_elpd <- -free_row$elpd_diff[[1]]
  se_diff <- free_row$se_diff[[1]]

  tibble(
    group = group_name,
    group_id = group_row$group_id,
    group_label = group_row$group_label,
    reference_model = "model_based",
    comparison_model = "model_free",
    comparison_row = free_row$comparison_row[[1]],
    elpd_diff = delta_elpd,
    se_diff = se_diff,
    p_superiority = stats::pnorm(delta_elpd / se_diff)
  )
}))

saveRDS(
  list(
    fits = fits,
    loos = loos,
    loo_summary = loo_summary,
    comparisons = comparisons,
    comparison_summary = comparison_summary,
    posterior_predictions = posterior_predictions_df,
    posterior_vmac = posterior_vmac_df
  ),
  file = "output/models/model_comparison_results.rds"
)

write.csv(
  loo_summary,
  file = "output/models/model_loo_summary.csv",
  row.names = FALSE
)

write.csv(
  comparison_summary,
  file = "output/models/model_comparison_summary.csv",
  row.names = FALSE
)

write.csv(
  posterior_predictions_df,
  file = "output/models/posterior_predictions_by_epoch.csv",
  row.names = FALSE
)

write.csv(
  posterior_vmac_df,
  file = "output/models/posterior_vmac_by_epoch.csv",
  row.names = FALSE
)

print(loo_summary)
print(comparisons)
print(comparison_summary)

## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

# Script for the model recovery of the model-free and model-based models (Supplementary Material)
# Run scripts/setup_model_recovery.R first. If the run is interrupted, running it again resumes it

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(cmdstanr, posterior, loo, dplyr, tidyr, purrr, tibble, ggplot2, ggdist, future, furrr)

source("scripts/model_load_data.R") # plyr (Rmisc) masks some dplyr verbs, hence dplyr::

# Settings ----
args <- commandArgs(trailingOnly = TRUE)
test_run <- identical(args, "test")
# Simulated datasets per group and generating model
N_iter <- if (length(args) >= 1 && !test_run) as.integer(args[1]) else 50
# Parallel jobs (4 cores each)
workers <- if (length(args) >= 2) as.integer(args[2]) else 8

output_dir <- "output/models/recovery"
if (test_run) output_dir <- file.path(output_dir, "test")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Same models as the main fits, without y_rep
model_files <- list(
  model_free = "models/recovery/model_free_rev_recovery.stan",
  model_based = "models/recovery/model_based_rev_recovery.stan"
)

# Same settings as the main fits, with fewer iterations
sampling_args <- list(
  seed = 2025,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  refresh = 0
)

# Same initial values as the main fits
build_init <- function(group_data, model_name) {
  K <- if (identical(model_name, "model_based")) 6 else 5
  J <- group_data$J

  function(chain_id = 1) {
    list(
      mu = rnorm(K, 0, 0.5),
      sigma = abs(rnorm(K, 0, 0.2)) + 1e-3,
      L_phi = diag(2),
      Z_beta = matrix(rnorm(K * J, 0, 0.3), nrow = K, ncol = J)
    )
  }
}

# Fitting ----
# Each parallel session compiles its own copy of the models
get_worker_models <- function() {
  if (!exists("worker_models", envir = globalenv())) {
    worker_dir <- file.path(tempdir(), "stan_models")
    dir.create(worker_dir, showWarnings = FALSE)
    assign("worker_models", lapply(model_files, function(path) {
      cmdstan_model(path, dir = worker_dir, force_recompile = TRUE)
    }), envir = globalenv())
  }
  get("worker_models", envir = globalenv())
}

# Population parameters of each posterior draw (mu, sigma, and lr_cor):
get_pop_draws <- function(draws) {
  draws <- unclass(as_draws_matrix(subset_draws(draws, variable = c("mu", "sigma", "L_phi"))))
  K <- sum(startsWith(colnames(draws), "mu["))
  list(
    mu = draws[, sprintf("mu[%d]", seq_len(K)), drop = FALSE],
    sigma = draws[, sprintf("sigma[%d]", seq_len(K)), drop = FALSE],
    lr_cor = draws[, "L_phi[2,1]"]
  )
}

# Function to fit one model and keep its PSIS-LOO, diagnostics, and population parameters:
fit_and_summarise <- function(stan_data, model_name, init_seed) {
  start_time <- Sys.time()
  fit_dir <- tempfile("model_recovery_fit_")
  dir.create(fit_dir)

  # Compiled before setting the seed (compiling uses R's RNG)
  model <- get_worker_models()[[model_name]]
  set.seed(init_seed, kind = "Mersenne-Twister")
  fit <- do.call(
    model$sample,
    c(
      list(
        data = stan_data,
        init = build_init(stan_data, model_name),
        output_dir = fit_dir
      ),
      sampling_args
    )
  )

  # PSIS-LOO without moment matching (only used in the main fits)
  loo_result <- fit$loo()
  pareto_k <- loo::pareto_k_values(loo_result)

  pop_draws <- get_pop_draws(fit$draws(variables = c("mu", "sigma", "L_phi")))

  summary <- list(
    loo = loo_result,
    n_div = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent),
    n_bad_k = sum(is.finite(pareto_k) & pareto_k > 0.7),
    minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
    pop_quantiles = lapply(pop_draws, function(x) {
      t(apply(as.matrix(x), 2, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE))
    })
  )

  rm(fit)
  gc()
  unlink(fit_dir, recursive = TRUE)
  summary
}

# Generating parameters ----
# Posterior draws of the fits to the real data
main_results <- "output/models/model_comparison_results.rds"
if (!file.exists(main_results)) {
  stop("Could not find '", main_results, "'. Run scripts/run_model_comparison.R first.")
}

main_fits <- readRDS(main_results)
generating_draws <- lapply(main_fits$param_draws, function(group_draws) {
  lapply(group_draws, get_pop_draws)
})

# Observed ELPD difference (model-based minus model-free) and z:
observed_diff <- main_fits$comparison_summary %>%
  transmute(
    group,
    elpd_diff_obs = ifelse(reference_model == "model_based", elpd_diff, -elpd_diff),
    z_obs = elpd_diff_obs / se_diff
  )
rm(main_fits)

# One simulated dataset ----
# Each job saves checkpoints, so an interrupted job resumes with the same data
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

run_recovery_job <- function(group, gen_model, iteration) {
  job_name <- sprintf("%s_%s_%d", group, gen_model, iteration)
  checkpoint <- function(step) file.path(checkpoint_dir, sprintf("%s_%s.rds", job_name, step))
  log_file <- file.path(output_dir, "model_recovery_log.txt")
  resumed <- file.exists(checkpoint("data"))
  cat(format(Sys.time()), if (resumed) "resume" else "start ", group, gen_model, iteration, "\n",
      file = log_file, append = TRUE)

  job_seed <- 100000 * match(group, names(dStan)) + 10000 * match(gen_model, names(model_files)) +
    iteration

  if (resumed) {
    simulated <- readRDS(checkpoint("data"))
  } else {
    # New participants from a random posterior draw (same RNG in all sessions)
    set.seed(job_seed, kind = "Mersenne-Twister")
    gen_draws <- generating_draws[[group]][[gen_model]]
    draw <- sample(nrow(gen_draws$mu), 1)
    pop_true <- list(mu = gen_draws$mu[draw, ], sigma = gen_draws$sigma[draw, ],
                     lr_cor = gen_draws$lr_cor[[draw]])
    beta <- simulate_subject_params(dStan[[group]]$J, pop_true$mu, pop_true$sigma,
                                    pop_true$lr_cor)
    # Reward of each distractor during reversal (none under nominal extinction)
    reversal_reward <- if (group == "RevFB") c(10, 500) / 500 else c(0, 0)
    simulated <- c(list(draw = draw, pop_true = pop_true),
                   simulate_omissions(dStan[[group]], beta, gen_model, reversal_reward))
    saveRDS(simulated, checkpoint("data"))
  }

  stan_data <- dStan[[group]]
  stan_data$Omission <- simulated$Omission
  stan_data$reward <- simulated$reward

  fit_summaries <- list()
  for (fit_model in names(model_files)) {
    if (file.exists(checkpoint(fit_model))) {
      fit_summaries[[fit_model]] <- readRDS(checkpoint(fit_model))
    } else {
      fit_summaries[[fit_model]] <- fit_and_summarise(
        stan_data, fit_model, init_seed = 10 * job_seed + match(fit_model, names(model_files))
      )
      saveRDS(fit_summaries[[fit_model]], checkpoint(fit_model))
    }
  }

  # Row 1 of loo_compare() is the winner
  loos <- lapply(fit_summaries, function(x) x$loo)
  comparison <- loo_compare(loos)

  result <- tibble(
    group = group,
    gen_model = gen_model,
    iteration = iteration,
    draw = simulated$draw,
    winner = rownames(comparison)[1],
    elpd_diff = comparison[2, "elpd_diff"],
    se_diff = comparison[2, "se_diff"],
    clear = abs(elpd_diff) > 2 * se_diff,
    elpd_model_free = loos$model_free$estimates["elpd_loo", "Estimate"],
    elpd_model_based = loos$model_based$estimates["elpd_loo", "Estimate"],
    div_model_free = fit_summaries$model_free$n_div,
    div_model_based = fit_summaries$model_based$n_div,
    bad_k_model_free = fit_summaries$model_free$n_bad_k,
    bad_k_model_based = fit_summaries$model_based$n_bad_k,
    minutes = sum(sapply(fit_summaries, function(x) x$minutes)), # Without interruptions
    # Kept to inspect parameter recovery
    pop_true = list(simulated$pop_true),
    pop_fit = list(lapply(fit_summaries, function(x) x$pop_quantiles))
  )

  saveRDS(result, file.path(output_dir, sprintf("model_recovery_%s.rds", job_name)))
  unlink(sapply(names(model_files), checkpoint))
  cat(format(Sys.time()), "end   ", group, gen_model, iteration, "\n",
      file = log_file, append = TRUE)

  result
}

# Test run ----
# One dataset, to time a job
if (test_run) {
  test_result <- run_recovery_job("RevFB", "model_based", 1)
  print(as.data.frame(test_result))
  quit(save = "no")
}

# Run ----
jobs <- expand_grid(
  group = names(dStan),
  gen_model = names(model_files),
  iteration = seq_len(N_iter)
)

# Finished datasets are skipped, and the rest run by iteration (unfinished ones first)
done <- list.files(output_dir, pattern = "^model_recovery_.*\\.rds$")
jobs <- jobs %>%
  filter(!sprintf("model_recovery_%s_%s_%d.rds", group, gen_model, iteration) %in% done) %>%
  dplyr::mutate(started = file.exists(file.path(
    checkpoint_dir, sprintf("%s_%s_%d_data.rds", group, gen_model, iteration)
  ))) %>%
  dplyr::arrange(desc(started), iteration, gen_model, group) %>%
  select(-started)

print(jobs)

plan(multisession, workers = workers)

# A failed job is logged and retried in the next run
future_pmap(
  jobs,
  function(group, gen_model, iteration) {
    tryCatch(run_recovery_job(group, gen_model, iteration), error = function(e) {
      cat(format(Sys.time()), "error ", group, gen_model, iteration, conditionMessage(e), "\n",
          file = file.path(output_dir, "model_recovery_log.txt"), append = TRUE)
      NULL
    })
  },
  .options = furrr_options(
    seed = TRUE,
    scheduling = Inf,
    packages = c("cmdstanr", "posterior", "loo", "dplyr", "tibble")
  )
)

plan(sequential)

# Results ----
model_recovery <- list.files(output_dir, pattern = "^model_recovery_.*\\.rds$",
                             full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  filter(iteration <= N_iter) %>%
  dplyr::mutate(
    # Within 2 SE: no clear difference
    outcome = case_when(
      !clear ~ "No clear difference",
      winner == "model_based" ~ "Model-based",
      TRUE ~ "Model-free"
    ),
    winner_div = ifelse(winner == "model_based", div_model_based, div_model_free),
    winner_bad_k = ifelse(winner == "model_based", bad_k_model_based, bad_k_model_free),
    elpd_mb_mf = elpd_model_based - elpd_model_free,
    z_mb_mf = elpd_mb_mf / se_diff
  )

saveRDS(model_recovery, file.path(output_dir, "recovery_results.rds"))

# Datasets in which the winning model had divergences or Pareto k > 0.7
model_recovery %>% filter(winner_div > 0 | winner_bad_k > 0)

# Expected: no clear difference with model-free data, and model-based wins with model-based data
recovery_summary <- model_recovery %>%
  dplyr::summarise(
    n = n(),
    model_based_clear = mean(outcome == "Model-based"),
    no_difference = mean(outcome == "No clear difference"),
    model_free_clear = mean(outcome == "Model-free"),
    .by = c(group, gen_model)
  )

print(recovery_summary)

# Confusion matrix ----
outcome_levels <- c("Model-free", "No clear difference", "Model-based")

confusion <- model_recovery %>%
  dplyr::mutate(
    group = factor(group, levels = group_metadata$group_code, labels = group_metadata$group_label),
    gen_model = factor(gen_model, levels = c("model_free", "model_based"),
                       labels = c("Model-free", "Model-based")),
    outcome = factor(outcome, levels = outcome_levels)
  ) %>%
  dplyr::count(group, gen_model, outcome) %>%
  complete(group, gen_model, outcome, fill = list(n = 0)) %>%
  dplyr::mutate(prop = n / sum(n), .by = c(group, gen_model))

print(confusion)

# Common style of the two recovery figures (sizes for a figure as wide as the text)
dir.create("output/figures", showWarnings = FALSE)
recovery_theme <- theme_Publication(text_size = 10, base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 11, colour = "black"),
    panel.spacing.x = unit(1.2, "lines"),
    plot.margin = margin(4, 8, 4, 4)
  )

# The counts are printed in each cell, so the fill (proportion) needs no legend
fig_recovery <- ggplot(confusion, aes(x = outcome, y = gen_model, fill = prop)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%d/%d", n, N_iter), colour = prop > 0.5), size = 4.6) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  scale_fill_gradient(low = "grey96", high = "#228BCB", limits = c(0, 1), guide = "none") +
  scale_x_discrete(labels = c("Model-free", "No clear\ndifference", "Model-based"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(~group) +
  labs(x = "Model selected by PSIS-LOO", y = "Generating model") +
  recovery_theme +
  theme(panel.grid.major = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/model_recovery.png", fig_recovery, width = 18, height = 7,
       units = "cm", dpi = 600, bg = "white")

# Model mimicry ----
# Proportion of simulated datasets with a z at least as large as the observed one:
mimicry <- model_recovery %>%
  left_join(observed_diff, by = "group") %>%
  dplyr::summarise(
    z_obs = first(z_obs),
    median_z = median(z_mb_mf),
    low = quantile(z_mb_mf, 0.025),
    high = quantile(z_mb_mf, 0.975),
    prop_at_least_obs = mean(z_mb_mf >= z_obs),
    .by = c(group, gen_model)
  )

print(mimicry)

label_group <- function(x) {
  factor(x, levels = group_metadata$group_code, labels = group_metadata$group_label)
}

fig_mimicry <- model_recovery %>%
  dplyr::mutate(
    group = label_group(group),
    gen_model = factor(gen_model, levels = c("model_free", "model_based"),
                       labels = c("Model-free", "Model-based"))
  ) %>%
  ggplot(aes(x = z_mb_mf, y = gen_model, fill = gen_model, colour = gen_model)) +
  # No clear difference (within 2 SE)
  annotate("rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_vline(data = dplyr::mutate(observed_diff, group = label_group(group)),
             aes(xintercept = z_obs), linetype = "dashed") +
  # Density and intervals above, one dot per dataset below
  stat_halfeye(height = 0.6, justification = -0.15, density = "unbounded", trim = FALSE,
               slab_alpha = 0.4, .width = c(0.66, 0.95), point_interval = median_qi) +
  stat_dots(side = "bottom", height = 0.6, justification = 1.15,
            binwidth = unit(1.2, "mm"), overflow = "compress", alpha = 0.6) +
  scale_fill_manual(values = c("Model-free" = "grey40", "Model-based" = "#228BCB"), guide = "none") +
  scale_colour_manual(values = c("Model-free" = "grey40", "Model-based" = "#228BCB"), guide = "none") +
  facet_wrap(~group) +
  labs(x = expression(bold(italic(z)[Delta * ELPD])), y = "Generating model") +
  recovery_theme

ggsave("output/figures/model_mimicry.png", fig_mimicry, width = 18, height = 8,
       units = "cm", dpi = 600, bg = "white")

# Population parameter recovery ----
# Only for inspection (mu on the parameter scale, sigma on the latent scale)
param_names <- list(
  model_free = c("intercept", "slope", "lr_acq", "lr_rev", "alpha0"),
  model_based = c("intercept", "slope", "lr_acq", "lr_rev", "rho", "alpha0")
)

to_parameter_scale <- function(x, parameter) {
  case_when(
    parameter == "intercept" ~ x,
    parameter == "slope" ~ 10 * phi_approx(x),
    TRUE ~ phi_approx(x)
  )
}

parameter_recovery <- model_recovery %>%
  dplyr::mutate(recovery = pmap(list(gen_model, pop_true, pop_fit), function(gen_model, pop_true, pop_fit) {
    map_dfr(c("mu", "sigma", "lr_cor"), function(pop_param) {
      fitted <- pop_fit[[gen_model]][[pop_param]]
      parameter <- if (pop_param == "lr_cor") "lr_acq_lr_rev" else param_names[[gen_model]]
      tibble(level = pop_param, parameter = parameter, true = pop_true[[pop_param]],
             low = fitted[, 1], est = fitted[, 2], high = fitted[, 3])
    })
  })) %>%
  select(group, gen_model, iteration, recovery) %>%
  unnest(recovery) %>%
  dplyr::mutate(
    covered = true >= low & true <= high,
    across(c(true, est), ~ if_else(level == "mu", to_parameter_scale(.x, parameter), .x))
  ) %>%
  dplyr::summarise(true_median = median(true), est_median = median(est),
            bias = median(est - true), coverage = mean(covered),
            .by = c(group, gen_model, level, parameter))

print(parameter_recovery, n = Inf)

# Instructed reversal parameter: population median, Phi(mu_rho), in the model-based datasets
rho_index <- which(param_names$model_based == "rho")
rho_recovery <- model_recovery %>%
  filter(gen_model == "model_based") %>%
  dplyr::mutate(
    true = map_dbl(pop_true, ~ phi_approx(.x$mu[[rho_index]])),
    low = map_dbl(pop_fit, ~ phi_approx(.x$model_based$mu[rho_index, 1])),
    est = map_dbl(pop_fit, ~ phi_approx(.x$model_based$mu[rho_index, 2])),
    high = map_dbl(pop_fit, ~ phi_approx(.x$model_based$mu[rho_index, 3]))
  ) %>%
  dplyr::summarise(n = n(), covered = sum(true >= low & true <= high),
                   r = cor(true, est), bias = median(est - true), .by = group)

print(rho_recovery)

# Phi(mu_rho) estimated by the model-based model in the model-free datasets (no instructed reversal)
rho_model_free_data <- model_recovery %>%
  filter(gen_model == "model_free") %>%
  dplyr::mutate(est = map_dbl(pop_fit, ~ phi_approx(.x$model_based$mu[rho_index, 2]))) %>%
  dplyr::summarise(n = n(), median_est = median(est), max_est = max(est), .by = group)

print(rho_model_free_data)

# Values reported in the Supplementary Material
saveRDS(
  list(
    n_iter = N_iter,
    confusion = confusion,
    mimicry = mimicry,
    max_z_model_free = model_recovery %>%
      filter(gen_model == "model_free") %>%
      dplyr::summarise(max_z = max(z_mb_mf), .by = group),
    parameter_recovery = parameter_recovery,
    rho_recovery = rho_recovery,
    rho_model_free_data = rho_model_free_data
  ),
  file.path(output_dir, "recovery_report.rds")
)

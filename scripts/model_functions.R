## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

create_epochs <- function(blocks, epoch = 2) {
  vapply(blocks, function(x, e = epoch) {
    ceiling(x / e)
  }, FUN.VALUE = numeric(1))
}

theme_Publication <- function(base_size = 12,
                              base_family = "sans",
                              text_size = 11) {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5),
      text = ggplot2::element_text(size = text_size),
      panel.background = ggplot2::element_rect(colour = NA),
      plot.background = ggplot2::element_rect(colour = NA),
      panel.border = ggplot2::element_rect(colour = NA),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      axis.title.y = ggplot2::element_text(angle = 90, vjust = 2),
      axis.title.x = ggplot2::element_text(vjust = -0.2),
      axis.text = ggplot2::element_text(),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.ticks = ggplot2::element_line(),
      panel.grid.major = ggplot2::element_line(colour = "#f0f0f0"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.2, "cm"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = ggplot2::element_rect(colour = NA, fill = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

compute_hdi_matrix <- function(x, cred_mass = 0.95) {
  alpha <- (1 - cred_mass) / 2
  t(apply(x, 2, function(vec) {
    stats::quantile(vec, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
  }))
}

reshape_draws_by_distractor <- function(pred_matrix, group_size) {
  n_draws <- nrow(pred_matrix)
  n_total <- ncol(pred_matrix) / 2
  n_subjects <- n_total / group_size

  if (n_subjects != floor(n_subjects)) {
    stop("group_size must divide ncol(pred_matrix)/2 exactly")
  }

  mat_1 <- pred_matrix[, 1:n_total, drop = FALSE]
  mat_2 <- pred_matrix[, (n_total + 1):(2 * n_total), drop = FALSE]

  list(
    n_draws = n_draws,
    group_size = group_size,
    n_subjects = n_subjects,
    distractor_1 = array(mat_1, dim = c(n_draws, group_size, n_subjects)),
    distractor_2 = array(mat_2, dim = c(n_draws, group_size, n_subjects))
  )
}

summarise_epoch_draws <- function(draw_array, epoch_len, cred_mass) {
  mean_by_trial <- apply(draw_array, c(1, 2), mean)
  n_draws <- nrow(mean_by_trial)
  n_epochs <- ceiling(ncol(mean_by_trial) / epoch_len)

  epoch_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_epochs)

  for (epoch_idx in seq_len(n_epochs)) {
    trial_idx <- ((epoch_idx - 1) * epoch_len + 1):min(epoch_idx * epoch_len, ncol(mean_by_trial))
    epoch_matrix[, epoch_idx] <- rowMeans(mean_by_trial[, trial_idx, drop = FALSE])
  }

  interval <- compute_hdi_matrix(epoch_matrix, cred_mass = cred_mass)

  data.frame(
    Epoch = seq_len(n_epochs),
    Mean = colMeans(epoch_matrix),
    Low = interval[, 1],
    High = interval[, 2]
  )
}

get_PPC_epochsV <- function(fit, group_size = 576, epoch_len = 24, credMass = 0.89) {
  preds <- fit$draws(variables = "V_rep", format = "matrix")
  reshaped <- reshape_draws_by_distractor(preds, group_size = group_size)

  out_1 <- summarise_epoch_draws(reshaped$distractor_1, epoch_len = epoch_len, cred_mass = credMass)
  out_2 <- summarise_epoch_draws(reshaped$distractor_2, epoch_len = epoch_len, cred_mass = credMass)

  out_1$Distractor <- "Distractor1"
  out_2$Distractor <- "Distractor2"

  dplyr::bind_rows(out_1, out_2)[, c("Distractor", "Epoch", "Mean", "Low", "High")]
}

get_PPC_epochsY <- function(fit, group_size = 576, epoch_len = 24, credMass = 0.95) {
  preds <- fit$draws(variables = "y_rep", format = "matrix")
  reshaped <- reshape_draws_by_distractor(preds, group_size = group_size)

  out_1 <- summarise_epoch_draws(reshaped$distractor_1, epoch_len = epoch_len, cred_mass = credMass)
  out_2 <- summarise_epoch_draws(reshaped$distractor_2, epoch_len = epoch_len, cred_mass = credMass)

  out_1$Distractor <- "Distractor1"
  out_2$Distractor <- "Distractor2"

  dplyr::bind_rows(out_1, out_2)[, c("Distractor", "Epoch", "Mean", "Low", "High")]
}

get_PPC_diff_epochs <- function(fit, group_size = 576, epoch_len = 24, credMass = 0.95) {
  preds <- fit$draws(variables = "y_rep", format = "matrix")
  reshaped <- reshape_draws_by_distractor(preds, group_size = group_size)

  avg_1 <- apply(reshaped$distractor_1, c(1, 2), mean)
  avg_2 <- apply(reshaped$distractor_2, c(1, 2), mean)
  diff_matrix <- avg_1 - avg_2

  n_draws <- nrow(diff_matrix)
  n_epochs <- ceiling(ncol(diff_matrix) / epoch_len)
  diff_epochs <- matrix(NA_real_, nrow = n_draws, ncol = n_epochs)

  for (epoch_idx in seq_len(n_epochs)) {
    trial_idx <- ((epoch_idx - 1) * epoch_len + 1):min(epoch_idx * epoch_len, ncol(diff_matrix))
    diff_epochs[, epoch_idx] <- rowMeans(diff_matrix[, trial_idx, drop = FALSE])
  }

  interval <- compute_hdi_matrix(diff_epochs, cred_mass = credMass)

  data.frame(
    Epoch = seq_len(n_epochs),
    MeanDiff = colMeans(diff_epochs),
    LowDiff = interval[, 1],
    HighDiff = interval[, 2]
  )
}

## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(dplyr, ggplot2, tidyr, tibble, cowplot, ggthemes, posterior, ggdist, ragg)

source("scripts/model_load_data.R")

if (!dir.exists("output/figures")) {
  dir.create("output/figures", recursive = TRUE)
}

results_path <- "output/models/model_comparison_results.rds"

if (!file.exists(results_path)) {
  stop("Could not find 'output/models/model_comparison_results.rds'. Run scripts/run_model_comparison.R first.")
}

model_results <- readRDS(results_path)
posterior_y_df <- model_results$posterior_predictions

if (is.null(posterior_y_df)) {
  stop("The results file does not contain posterior predictions. Run scripts/run_model_comparison.R first.")
}

posterior_y_df <- posterior_y_df %>%
  mutate(Epoch = as.numeric(as.character(Epoch)))

group_metadata <- tibble::tribble(
  ~group_code, ~group_id, ~group_label,
  "RevFB", "reversal_with_feedback", "Reversal with feedback",
  "Rev", "nominal_extinction", "Nominal extinction"
)

model_metadata <- tibble::tribble(
  ~model_code, ~model_id, ~model_label,
  "model_based", "model_based_revaluation", "Model-based",
  "model_free", "model_free_learning", "Model-free"
)

ppc_levels <- list(
  high = c("#C9333E", "#F29CA5"),
  low = c("#228BCB", "#9DDCF2")
)

# Both models are drawn in the same panel (model-based: solid line with its 95% credible
# interval; model-free: dashed line), so that their predictions at the onset of the
# reversal stage can be compared directly
prediction_panels <- list()

for (group_name in group_metadata$group_code) {
  group_row <- group_metadata[group_metadata$group_code == group_name, ]

  y_epochs <- posterior_y_df %>%
    filter(group == group_name) %>%
    mutate(
      Phase = factor(Phase, levels = c("Acquisition", "Reversal")),
      Model = factor(dplyr::recode(model, model_based = "Model-based", model_free = "Model-free"),
                     levels = c("Model-based", "Model-free"))
    )

  raw_plot <- dPlots[[group_name]] %>%
    mutate(
      Epoch = as.numeric(as.character(Epoch)),
      Phase = factor(Phase, levels = c("Acquisition", "Reversal"))
    )

  prediction_panels[[group_name]] <- ggplot(y_epochs, aes(x = Epoch, y = Mean * 100, colour = Singleton)) +
    geom_ribbon(
      data = filter(y_epochs, Model == "Model-based"),
      aes(ymin = Low * 100, ymax = High * 100, fill = Singleton),
      alpha = 0.3, colour = NA
    ) +
    geom_line(aes(linetype = Model), linewidth = 1.1) +
    geom_errorbar(
      data = raw_plot,
      aes(y = P * 100, ymin = (P - se) * 100, ymax = (P + se) * 100),
      width = 0.01, linewidth = 0.5
    ) +
    geom_point(data = raw_plot, aes(y = P * 100, fill = Singleton), shape = 21, size = 2.6, stroke = 0.8) +
    facet_grid(. ~ Phase, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = seq(1, 12, 1)) +
    scale_y_continuous(breaks = seq(10, 30, 10), labels = function(x) paste0(x, "%")) +
    scale_colour_manual(values = c("High" = ppc_levels$high[1], "Low" = ppc_levels$low[1])) +
    scale_fill_manual(values = c("High" = ppc_levels$high[2], "Low" = ppc_levels$low[2]), guide = "none") +
    scale_linetype_manual(values = c("Model-based" = "solid", "Model-free" = "22")) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(shape = 21, linetype = "solid",
                                                           fill = c(ppc_levels$high[2], ppc_levels$low[2]))),
      linetype = guide_legend(order = 2, override.aes = list(colour = "grey20"))
    ) +
    labs(title = group_row$group_label, x = "Epoch", y = "% of omission trials",
         colour = "Singleton", linetype = "Predictions") +
    theme_Publication(text_size = 10, base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 6)),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 15),
      legend.key.width = unit(1.4, "cm"),
      strip.text.x = element_text(size = 18, face = "bold"),
      strip.background = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      plot.margin = margin(5, 6, 6, 10)
    )
}

# Instructed reversal parameter (rho) of the model-based model ----
# Population mean of rho on its 0-1 scale, computed from mu and sigma as mu_beta[5]
# in the Stan model. With rho = .5, the instructions make both distractors equally valuable
rho_draws <- bind_rows(lapply(group_metadata$group_code, function(group_name) {
  draws <- as_draws_df(model_results$param_draws[[group_name]]$model_based)
  tibble(group = group_name, rho = pnorm(draws[["mu[5]"]] / sqrt(1 + draws[["sigma[5]"]]^2)))
}))

# The two groups were fitted separately, so their draws can be paired to compare them
rho_difference <- rho_draws$rho[rho_draws$group == "RevFB"] - rho_draws$rho[rho_draws$group == "Rev"]

rho_summary <- rho_draws %>%
  dplyr::summarise(median = median(rho), low = quantile(rho, 0.025), high = quantile(rho, 0.975),
                   p_above_half = mean(rho > 0.5), .by = group) %>%
  bind_rows(tibble(group = "RevFB_minus_Rev", median = median(rho_difference),
                   low = quantile(rho_difference, 0.025), high = quantile(rho_difference, 0.975),
                   p_above_half = mean(rho_difference > 0)))

write.csv(rho_summary, "output/models/rho_summary.csv", row.names = FALSE)

# Posterior of each group, with a gradient that follows the density (darker where the
# posterior is denser), and the median and the 66% and 95% intervals below
rho_panel <- rho_draws %>%
  dplyr::mutate(group = factor(group, levels = rev(group_metadata$group_code), labels = rev(group_metadata$group_label))) %>%
  ggplot(aes(x = rho, y = group)) +
  stat_slab(aes(fill_ramp = after_stat(pdf)), fill = "grey25", fill_type = "gradient",
            height = 0.9, justification = -0.05) +
  stat_pointinterval(.width = c(0.66, 0.95), point_interval = median_qi,
                     point_size = 3, interval_size_range = c(1.2, 2.6)) +
  scale_fill_ramp_continuous(from = "white", guide = "none") +
  # Reference line drawn over the distributions
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.8) +
  # Group names inside the panel, where there is no density, so that long axis labels do
  # not push panels A and B to the right when the panels are aligned
  annotate("text", x = 0.99, y = rev(group_metadata$group_label), label = rev(group_metadata$group_label),
           hjust = 1, vjust = -0.6, size = 5.5) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = expression(bold("Instructed reversal (" * rho * ")")), y = NULL) +
  theme_Publication(text_size = 10, base_size = 12) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(22, 22, 4, 10)
  )

# Shared legend above the whole figure, outside the panels
prediction_legend <- get_plot_component(
  prediction_panels[[1]] + theme(legend.position = "top", legend.box.margin = margin(0, 0, 0, 0)),
  "guide-box-top", return_all = TRUE
)
prediction_panels <- lapply(prediction_panels, function(x) x + theme(legend.position = "none"))

prediction_figure <- plot_grid(
  prediction_legend,
  plot_grid(
    plotlist = c(prediction_panels, list(rho_panel)),
    ncol = 1,
    labels = c("A", "B", "C"),
    label_size = 20,
    rel_heights = c(1, 1, 0.7),
    align = "v",
    axis = "lr"
  ),
  ncol = 1,
  rel_heights = c(0.05, 1)
)

# ragg draws the gradient of panel C as a true gradient (the default device draws stripes)
ggsave(
  filename = "output/figures/fig-6.png",
  plot = prediction_figure,
  width = 25,
  height = 30,
  units = "cm",
  dpi = 1200,
  device = ragg::agg_png,
  bg = "white"
)

print(prediction_figure)

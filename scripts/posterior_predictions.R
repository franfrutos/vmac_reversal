## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(dplyr, ggplot2, tidyr, tibble, cowplot, ggthemes)

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
  ~group_code, ~group_id, ~group_label, ~group_title, ~line_type, ~point_shape,
  "RevFB", "reversal_with_feedback", "Reversal with feedback", "Instructed reversal with feedback", "solid", 21,
  "Rev", "nominal_extinction", "Nominal extinction", "Instructed reversal in nominal extinction", "dashed", 22
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

prediction_panels <- list()
legend_plot <- NULL

for (group_name in group_metadata$group_code) {
  group_row <- group_metadata[group_metadata$group_code == group_name, ]

  y_epochs <- posterior_y_df %>%
    filter(group == group_name) %>%
    mutate(
      Phase = factor(Phase, levels = c("Acquisition", "Reversal")),
      model_label = dplyr::recode(
        model,
        model_based = "Model-based",
        model_free = "Model-free"
      ),
      model_label = factor(model_label, levels = c("Model-based", "Model-free"))
    )

  raw_plot <- dPlots[[group_name]] %>%
    mutate(
      Epoch = as.numeric(as.character(Epoch)),
      Phase = factor(Phase, levels = c("Acquisition", "Reversal"))
    ) %>%
    tidyr::crossing(model_label = factor(c("Model-based", "Model-free"), levels = c("Model-based", "Model-free")))

  pred_plot <- ggplot(
    y_epochs,
    aes(x = Epoch, y = Mean * 100, color = Singleton, fill = Singleton)
  ) +
    geom_line(linewidth = 1.1, linetype = group_row$line_type) +
    geom_ribbon(aes(ymin = Low * 100, ymax = High * 100), alpha = 0.3, color = NA) +
    geom_errorbar(
      data = raw_plot,
      aes(
        y = P * 100,
        ymin = (P - se) * 100,
        ymax = (P + se) * 100
      ),
      width = 0.01,
      linewidth = 0.5,
      position = position_dodge(width = 0.5)
    ) +
    geom_point(
      data = raw_plot,
      aes(y = P * 100),
      shape = group_row$point_shape,
      size = 2.05,
      stroke = 0.8,
      position = position_dodge(width = 0.5)
    ) +
    facet_grid(model_label ~ Phase, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = seq(1, 12, 1)) +
    scale_y_continuous(breaks = seq(10, 30, 10), labels = function(x) paste0(x, "%")) +
    scale_color_manual(values = c("High" = ppc_levels$high[1], "Low" = ppc_levels$low[1])) +
    scale_fill_manual(values = c("High" = ppc_levels$high[2], "Low" = ppc_levels$low[2]), guide = "none") +
    guides(
      color = guide_legend(
        override.aes = list(
          linetype = "solid",
          linewidth = 1.1,
          shape = 21,
          size = 3,
          fill = c(ppc_levels$high[2], ppc_levels$low[2]),
          alpha = 1
        )
      )
    ) +
    labs(
      title = group_row$group_title,
      x = "Epoch",
      y = "% of omission trials",
      color = "Singleton"
    ) +
    theme_Publication(text_size = 10, base_size = 12) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      legend.position = "bottom",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1.2, "lines"),
      legend.spacing.x = unit(0.2, "cm"),
      strip.text.x = element_text(size = 18, face = "bold"),
      strip.text.y.right = element_text(angle = 270, size = 18, face = "bold"),
      strip.background = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      plot.margin = margin(5, 2, 0, 2)
    )

  if (is.null(legend_plot)) {
    legend_plot <- cowplot::get_legend(
      pred_plot +
        theme(
          legend.position = "bottom",
          legend.justification = "center",
          legend.box.margin = margin(0, 0, 0, 0),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(0, 0, 0, 0)
        )
    )
  }

  prediction_panels[[group_name]] <- pred_plot + theme(legend.position = "none")
}

prediction_figure <- plot_grid(
  plotlist = prediction_panels,
  ncol = 1,
  labels = c("A", "B"),
  label_size = 20,
  align = "v",
  axis = "lr"
)

prediction_figure <- plot_grid(
  prediction_figure,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.055)
)

ggsave(
  filename = "output/figures/fig-6.png",
  plot = prediction_figure,
  width = 25,
  height = 25,
  units = "cm",
  dpi = 1200
)

print(prediction_figure)

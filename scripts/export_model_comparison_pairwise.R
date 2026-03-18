## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

results_path <- "output/models/model_comparison_results.rds"
output_path <- "output/models/model_comparison_summary.csv"

if (!file.exists(results_path)) {
  stop("Could not find 'output/models/model_comparison_results.rds'.")
}

model_results <- readRDS(results_path)

if (!is.null(model_results$comparison_summary)) {
  comparison_summary <- as.data.frame(model_results$comparison_summary)
  write.csv(comparison_summary, output_path, row.names = FALSE)
  print(comparison_summary)
  invisible(comparison_summary)
} else {
  comparisons <- model_results$comparisons

  if (is.null(comparisons) || length(comparisons) == 0) {
    stop("The results file does not contain model comparison objects.")
  }

  group_metadata <- data.frame(
    group = c("Rev", "RevFB"),
    group_id = c("nominal_extinction", "reversal_with_feedback"),
    group_label = c("Nominal extinction", "Reversal with feedback"),
    stringsAsFactors = FALSE
  )

  comparison_rows <- lapply(names(comparisons), function(group_name) {
    comparison_df <- as.data.frame(comparisons[[group_name]])
    comparison_df$comparison_row <- rownames(comparison_df)
    rownames(comparison_df) <- NULL
    group_row <- group_metadata[group_metadata$group == group_name, , drop = FALSE]

    if (nrow(comparison_df) < 2) {
      stop(sprintf("Model comparison for group '%s' does not contain the expected two rows.", group_name))
    }

    free_row <- comparison_df[2, , drop = FALSE]
    delta_elpd <- -free_row$elpd_diff[[1]]
    se_diff <- free_row$se_diff[[1]]

    data.frame(
      group = group_name,
      group_id = group_row$group_id,
      group_label = group_row$group_label,
      reference_model = "model_based",
      comparison_model = "model_free",
      comparison_row = free_row$comparison_row[[1]],
      elpd_diff = delta_elpd,
      se_diff = se_diff,
      p_superiority = stats::pnorm(delta_elpd / se_diff),
      stringsAsFactors = FALSE
    )
  })

  comparison_summary <- do.call(rbind, comparison_rows)
  write.csv(comparison_summary, output_path, row.names = FALSE)

  print(comparison_summary)
}

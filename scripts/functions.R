## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

required_packages <- c(
  "dplyr", "Hmisc", "grid", "ggthemes", "ggplot2", "marginaleffects",
  "Rmisc", "tidyr", "tibble", "plyr", "afex", "emmeans", "effectsize"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    sprintf(
      "Missing required packages: %s",
      paste(missing_packages, collapse = ", ")
    )
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Function to filter data:
filter_data <- function(data,
                        f.absent = T,
                        acc = T,
                        sd_filter = NULL,
                        fixed = T,
                        two_t = F,
                        phase = NULL,
                        experiment = "e1"
) {
  
  filt_trial <- NA
  if (two_t) {
    filt_trial <- c(seq(1, 577, 24), seq(2, 577, 24))
  }
  
  out <- data %>%
    filter(!(trial_num %in% filt_trial))

  if (!is.null(phase)) {
    out <- out %>% filter(Phase %in% phase)
  }
  
  if (experiment == "e1") {
    out <- out[which(!out$ID %in% acc_ex,),]
  } else {
    out <- out[which(!out$ID %in% acc_ex2),]
  }

  if (acc)
    out <- out[which(out$correct == 1),]
  
  
  if (f.absent)
    out <- out[which(out$Singleton != "Absent"),]
  
  if (fixed)
    out <- out[which(out$rt > 150 & out$rt < 1800),]
  
  if (is.numeric(sd_filter) & ifelse(is.null(sd_filter), F, sd_filter %in% 1:3)) {
    out <- out %>%
      group_by(ID) %>%
      dplyr::mutate(high_rt = mean(rt, na.rm = T) + sd(rt, na.rm = T)*sd_filter,
                    low_rt = mean(rt, na.rm = T) - sd(rt, na.rm = T)*sd_filter) %>%
      ungroup()  %>%
      filter(rt > low_rt, rt < high_rt) %>%
      select(-c(high_rt, low_rt))
  }
  out <- out[which(out$Phase != "Practice"),]
  return(out)
}

# Function to collapse blocks in different epochs.
create_epochs <- function(blocks, epoch = 2) {
  vapply(blocks, function(x, e = epoch) {
    ceiling(x / e)
  }, FUN.VALUE = numeric(1))
}

# Theme used in the plots
theme_Publication <-
  function(base_size = 12,
           base_family = "sans",
           text_size = 11) {
    library(ggthemes)
    (
      theme_foundation(base_size = base_size, base_family = base_family)
      + theme(
        plot.title = element_text(
          face = "bold",
          size = rel(1.2),
          hjust = 0.5
        ),
        text = element_text(size = text_size),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.2, "cm"),
        #legend.margin = margin(0, "cm"),
        #legend.title = element_text(face="italic"),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = NA, fill = NA),
        strip.text = element_text(face = "bold")
      )
    )
    
  }

transform_epoch <- function(epochs, length_v) {
  n <- length(epochs)
  paso <- length_v / n
  ajuste <- sapply(1:n, function(i) {
    mean(seq((i - 1) * paso + 1, i * paso))
  })
  return(ajuste)
}

# Functions to plot

# This function gets predictions for both models presented in the main text to plot (Figure 2)
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_predictions <- function(d1, d2, m1, m2, epoch = T) {
  library(marginaleffects)
  library(dplyr)
  
  if (epoch) {
    d1$Block_num <- create_epochs(d1$Block_num)
    d2$Block_num <- create_epochs(d2$Block_num)
  }
  raw_data <- rbind(
    d1 %>% group_by(ID, Phase, Block_num, Singleton) %>% dplyr::summarise(estimate = mean(rt)) %>%
      Rmisc::summarySEwithin(
        data = .,
        measurevar = "estimate",
        withinvars = c("Phase", "Block_num", "Singleton"),
        idvar = "ID"
      ) %>% mutate(Phase = "Acquisition"),
    d2 %>% group_by(ID, Phase, Block_num, Singleton) %>% dplyr::summarise(estimate = mean(rt)) %>%
      Rmisc::summarySEwithin(
        data = .,
        measurevar = "estimate",
        withinvars = c("Phase", "Block_num", "Singleton"),
        idvar = "ID"
      ) %>% mutate(Phase = "Reversal")
  )
  raw_data$Singleton <-
    factor(raw_data$Singleton, levels = c("High", "Low", "Absent"))
  raw_data$Block_num <- as.numeric(raw_data$Block_num) * 2 - .5
  
  model_preds <- rbind(
    predictions(
      m1,
      newdata = datagrid(
        Singleton = unique,
        Block_num = seq(1, 12, .01),
        ID = NA
      ),
      re.form = NA,
      transform = \(x) exp(x + (sigma(m1) ^ 2) / 2),
    ) %>% mutate(Phase = "Acquisition"),
    predictions(
      m2,
      newdata = datagrid(
        Singleton = unique,
        Block_num = seq(13, 24, .01),
        ID = NA
      ),
      re.form = NA,
      transform = \(x) exp(x + (sigma(m2) ^ 2) / 2)
    ) %>% mutate(Phase = "Reversal")
  )
  
  model_preds$Singleton <-
    factor(model_preds$Singleton, levels = c("High", "Low", "Absent"))
  
  return(list(raw = raw_data,
              mod = model_preds))
}

# Intermediate function to get VMAC and AC effects in get_comparisons()
get_raw_effect <- function(d, epoch = T) {
  if (epoch)
    d$Block_num <- create_epochs(d$Block_num)
  effs <- d %>%
    group_by(ID, Phase, Block_num, Singleton) %>%
    dplyr::summarise(RT = mean(rt)) %>%
    spread(Singleton, RT) %>%
    mutate(VMAC = High-Low,
           AC = Low-Absent) %>%
    ungroup() %>%
    drop_na()
  
  return(
    rbind(
      Rmisc::summarySEwithin(
        data = effs,
        measurevar = "VMAC",
        withinvars = c("Phase", "Block_num"),
        idvar = "ID"
      ) %>%
        dplyr::rename("estimate" = "VMAC") %>% mutate(Block_num = as.numeric(Block_num),
                                                      Effect = "VMAC"),
      Rmisc::summarySEwithin(
        data = effs,
        measurevar = "AC",
        withinvars = c("Phase", "Block_num"),
        idvar = "ID"
      ) %>% dplyr::rename("estimate" = "AC") %>%
        mutate(Block_num = as.numeric(Block_num),
               Effect = "Attentional Capture")
    )
  )
  
}

# This function gets conditional effects for each contrasts in both models presented in the main text to plot (Figure 3)
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_comparisons <- function(d1, d2, m1, m2, epoch = T) {
  raw_data <- rbind(get_raw_effect(d1, epoch),
                    get_raw_effect(d2, epoch))
  raw_data$Effect <-
    factor(raw_data$Effect, levels = c("VMAC", "Attentional Capture"))
  raw_data$Block_num[raw_data$Phase == "Reversal"] <-
    raw_data$Block_num[raw_data$Phase == "Reversal"] + 6
  raw_data$Block_num <- as.numeric(raw_data$Block_num) * 2 - .5
  
  
  model_comps <- bind_rows(
    comparisons(
      m1,
      variables = list(Singleton = "revsequential"),
      newdata = datagrid(Block_num = seq(1, 12, .01),
                         ID = NA),
      re.form = NA,
      comparison = \(hi, lo) exp(hi + (sigma(m1) ^ 2) / 2) - exp(lo +
                                                                   (sigma(m1) ^ 2) / 2)
    ) %>% mutate(
      Effect = rep(
      c("VMAC", "Attentional Capture"), each = nrow(.) / 2
    ), Phase = "Acquisition"),
    comparisons(
      m2,
      variables = list(Singleton = c("revsequential")),
      newdata = datagrid(Block_num = seq(13, 24, .01),
                         ID = NA),
      re.form = NA,
      comparison = \(hi, lo) exp(hi + (sigma(m2) ^ 2) / 2) - exp(lo +
                                                                   (sigma(m2) ^ 2) / 2)
    ) %>% 
      mutate(
         Effect = rep(
           c("VMAC", "Attentional Capture"), each = nrow(.) / 2
         ),
        Phase = "Reversal")
  )
  
   model_comps$Effect <-
     factor(model_comps$Effect, levels = c("VMAC", "Attentional Capture"))
  
  return(list(raw = raw_data,
              mod = model_comps))
}

sim_lme4 <- function(dat, model){
  library(lme4)
  library(lmerTest)
  dat$log_RT <- doSim(model)
  f <- model@call[["formula"]]
  tmp <- lmer(formula = f, data = dat,
              control = lmerControl(optimizer = "bobyqa"))
  return(list(main=summary(tmp)$coefficients[,5][2] < .05,# Main VMAC
              inter = summary(tmp)$coefficients[,5][5] < .05, # VMAC Block interaction
              conv = check_convergence(tmp)))
}

# Report trial exclusions: 
report_exclusions <- function(d, acc) {
  # Only the trials of the task (no practice or awareness test)
  d <- d[which(d$Phase %in% c("Acquisition", "Reversal")),]
  incorrect <- (1-nrow(d[which(d$correct == 1 & !d$ID %in% acc),])/
     nrow(d[which(!d$ID %in% acc),]))*100
  outliers <- (1-nrow(d[which(d$correct == 1 & !d$ID %in% acc & (d$rt > 150 & d$rt < 1800)),])/
      nrow(d[which(d$correct == 1 & !d$ID %in% acc),]))*100 
  
  cat(paste0("Incorrect responses: ", round(incorrect, 2), "%",
             "\nOutliers RTs: ", round(outliers, 2), "%"))
  # Percentages reported in the manuscript
  invisible(list(incorrect = incorrect, outliers = outliers))
}

# Functions to build ROC curve form probit
get_roc <- function(d, c) {
  fas <- seq(0, 1, 0.001)
  hts <- pnorm(d + qnorm(fas))
  return(list(dp = tibble(hts, fas), c = tibble(fas=pnorm(c),
                                                hts=pnorm(d+c))))
}

get_roc_probit <- function(model) {
  # Criterion (signal = 0.5) and d' (-0.5 minus 0.5) by Experiment, selected by value
  # because the row order of predictions() depends on the marginaleffects version
  c_s <- predictions(model, datagrid(ID = unique, Experiment = unique, signal = unique),
                     type = "link", re.form = NA, by = c("Experiment", "signal"))
  c_s <- c_s[c_s$signal == 0.5,]

  d_s <- predictions(model, datagrid(ID = unique, Experiment = unique, signal = unique),
                     type = "link", re.form = NA, by = c("Experiment", "signal"),
                     hypothesis = difference ~ revreference | Experiment)
  

  roc0 <- get_roc(d_s$estimate[1], c_s$estimate[1])
  roc1 <- get_roc(d_s$estimate[2], c_s$estimate[2])
  
  # dprime
  roc0$dp$upp <- get_roc(d_s$conf.high[1], 0)$dp$hts
  roc0$dp$low <- get_roc(d_s$conf.low[1], 0)$dp$hts
  
  roc1$dp$upp <- get_roc(d_s$conf.high[2], 0)$dp$hts
  roc1$dp$low <- get_roc(d_s$conf.low[2], 0)$dp$hts
  
  rocs_d <- rbind(roc0$dp, roc1$dp)
  
  rocs_d$Experiment <- as.factor(rep(c("Implicit reversal", "Nominal extinction"), each = nrow(rocs_d)/2))
  
  # c 
  
  rocs_c <- rbind(roc0$c, roc1$c)
  rocs_c$Experiment <- as.factor(c("Implicit reversal", "Nominal extinction"))
  
  # Raw data:
  raw_aw <- d_aw %>% 
    dplyr::summarise(Prop = mean(response_aw), .by = c(ID, Experiment, signal)) %>% drop_na() %>%
    spread(signal, Prop) %>%
    mutate(hts = `-0.5`,
           fas = (`0.5`))
  
  raw_aw$Experiment <- ifelse(raw_aw$Experiment== "A", "Nominal extinction", "Implicit reversal")
  
  plot <- ggplot(rocs_d, aes(1-fas, 1-hts, color = Experiment, shape = Experiment)) +
    geom_point(data = raw_aw, alpha = 2/3, aes(fas, hts)) +
    geom_line(aes(fas, hts, linetype = Experiment), linewidth = .8)+
    geom_ribbon(aes(ymin = upp, ymax = low, fill = Experiment, x = fas, y = hts), color = NA,
                alpha = .1)+
    geom_point(data = rocs_c, aes(fas, hts, fill = Experiment),
               size = 2) +
    geom_line(data = rocs_d %>% filter(((fas < pnorm(c_s$conf.high[1]) & fas > pnorm(c_s$conf.low[1])) & Experiment == "Implicit reversal") |
                                         ((fas < pnorm(c_s$conf.high[2]) & fas > pnorm(c_s$conf.low[2])) & Experiment == "Nominal extinction")),
              linewidth = 3, alpha = .3, lineend = "round", aes(fas, hts)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
    scale_color_manual(values = rev(c("#3B4994", "#A73061"))) +
    scale_fill_manual(values = rev(c("#7E87D4", "#D86491"))) +
    scale_shape_manual(values = rev(c(22, 23))) +
    scale_linetype_manual(values = c(
      "Nominal extinction" = "dashed",
      "Implicit reversal" = "twodash"
    )) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_Publication() +
    labs(x = "False alarm rate", y = "Hit rate") +
    theme(legend.position = "none",
          plot.margin = margin(5, 2, 5, 2))
  
  return(plot)
}

# Split-half reliability of d-primes
# Function to obtain dpimes in rapidsplithalf
dprimesfunc <- function(x) {
  
  if (mean(x) == 0 | mean(x) == 1) { # "unseen" in all trials
    # Hautus (1995) correction
    cell_corrected <- sum(x)+.5 #adds .5 to all cells
    Ntrials_corrected <- length(x)+1
    
    return(qnorm(cell_corrected / Ntrials_corrected))
  } 
  
  return(qnorm(mean(x)))
}

# Functions to get average and SE for raw data in mixed designs. The source code is in: :
summarySEwithin2 <- function (data = NULL,
                              measurevar,
                              betweenvars = NULL,
                              withinvars = NULL,
                              idvar = NULL,
                              na.rm = FALSE,
                              conf.interval = 0.95,
                              .drop = TRUE)
{
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop = FALSE], FUN = is.factor, FUN.VALUE = logical(1))
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message(
      "Automatically converting the following non-factors to factors: ",
      paste(nonfactorvars, collapse = ", ")
    )
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  datac <- summarySE(
    data,
    measurevar,
    groupvars = c(betweenvars, withinvars),
    na.rm = na.rm,
    conf.interval = conf.interval,
    .drop = .drop
  )
  
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop = .drop)
  measurevar_n <- paste(measurevar, "_norm", sep = "")
  ndatac <- summarySE(
    ndata,
    measurevar_n,
    groupvars = c(betweenvars, withinvars),
    na.rm = na.rm,
    conf.interval = conf.interval,
    .drop = .drop
  )
  nWithinGroups <- prod(vapply(ndatac[, withinvars, drop = FALSE], FUN = nlevels, FUN.VALUE = numeric(1)))
  correctionFactor <- sqrt(nWithinGroups / (nWithinGroups - 1))
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  merge(datac, ndatac)
}
normDataWithin <- function(data = NULL,
                           idvar,
                           measurevar,
                           betweenvars = NULL,
                           na.rm = FALSE,
                           .drop = TRUE) {
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- plyr::ddply(
    data,
    c(idvar, betweenvars),
    .drop = .drop,
    .fun = function(xx, col, na.rm) {
      c(subjMean = mean(xx[, col], na.rm = na.rm))
    },
    measurevar,
    na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep = "")
  data[, measureNormedVar] <- data[, measurevar] - data[, "subjMean"] +
    mean(data[, measurevar], na.rm = na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

# Reporting helpers ----
format_p_latex <- function(p, digits = 3) {
  if (is.na(p)) {
    return("NA")
  }
  threshold <- 10^(-digits)
  if (p < threshold) {
    return(paste0("< .", paste(rep("0", digits - 1), collapse = ""), "1"))
  }
  out <- formatC(p, digits = digits, format = "f")
  paste0("= ", sub("^0", "", out))
}

format_num_latex <- function(x, digits = 2, drop_leading_zero = FALSE) {
  if (is.na(x)) {
    return("NA")
  }
  out <- formatC(x, digits = digits, format = "f")
  if (drop_leading_zero) {
    out <- sub("^0", "", out)
    out <- sub("^-0", "-.", out)
  }
  out
}

format_beta_label_latex <- function(label) {
  if (is.null(label) || is.na(label) || !nzchar(label)) {
    return("")
  }

  format_beta_component <- function(part) {
    if (grepl("\\\\|\\^|_|\\{|\\}", part)) {
      return(part)
    }

    paste0("\\text{", part, "}")
  }

  parts <- strsplit(label, "\\s+\\\\times\\s+|\\s+x\\s+", perl = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]

  if (length(parts) <= 1) {
    return(format_beta_component(label))
  }

  paste(vapply(parts, format_beta_component, character(1)), collapse = " \\times ")
}

extract_model_term <- function(model, term, conf_level = 0.95, conf_method = "Wald") {
  sm <- coef(summary(model))
  if (!term %in% rownames(sm)) {
    stop(sprintf("Term '%s' was not found in the model summary.", term))
  }
  
  out <- as.data.frame(sm[term, , drop = FALSE])
  out$term <- term
  out$estimate <- out[[1]]
  out$std_error <- out[[2]]
  
  # Columns selected by name, because lmerTest adds a df column before the t value
  stat_col <- grep("^(t|z) value$", names(out), value = TRUE)[1]
  p_col <- names(out)[ncol(out) - 3 + 1]
  df_col <- "df"

  if ("Pr(>|t|)" %in% names(out)) {
    p_col <- "Pr(>|t|)"
  } else if ("Pr(>|z|)" %in% names(out)) {
    p_col <- "Pr(>|z|)"
  }
  
  out$statistic <- out[[stat_col]]
  out$statistic_label <- if (grepl("z", stat_col, ignore.case = TRUE)) "z" else "t"
  out$df <- if (df_col %in% names(out) && grepl("df", df_col, ignore.case = TRUE)) out[[df_col]] else NA_real_
  out$p_value <- if (p_col %in% names(out)) out[[p_col]] else NA_real_
  
  ci <- suppressMessages(suppressWarnings(confint(model, parm = term, level = conf_level, method = conf_method)))
  if (is.matrix(ci)) {
    out$conf_low <- ci[1, 1]
    out$conf_high <- ci[1, 2]
  } else {
    out$conf_low <- NA_real_
    out$conf_high <- NA_real_
  }
  
  out[, c("term", "estimate", "std_error", "statistic", "statistic_label", "df", "p_value", "conf_low", "conf_high")]
}

latex_model_term <- function(model,
                             term,
                             beta_label = NULL,
                             conf_level = 0.95,
                             estimate_digits = 3,
                             stat_digits = 2,
                             df_digits = 2,
                             ci_digits = 2,
                             conf_method = "Wald",
                             include_ci = FALSE) {
  res <- extract_model_term(model, term, conf_level = conf_level, conf_method = conf_method)
  beta_label <- if (is.null(beta_label)) term else beta_label
  
  beta_part <- paste0(
    "$\\beta_{", format_beta_label_latex(beta_label), "} = ",
    format_num_latex(res$estimate, digits = estimate_digits),
    "$"
  )
  
  stat_part <- if (is.na(res$df)) {
    paste0(
      "$", res$statistic_label, " = ",
      format_num_latex(res$statistic, digits = stat_digits),
      "$"
    )
  } else {
    paste0(
      "$", res$statistic_label, "(",
      format_num_latex(res$df, digits = df_digits),
      ") = ",
      format_num_latex(res$statistic, digits = stat_digits),
      "$"
    )
  }
  
  p_part <- paste0("$p ", format_p_latex(res$p_value), "$")
  
  if (!include_ci || anyNA(c(res$conf_low, res$conf_high))) {
    return(paste(beta_part, stat_part, p_part, sep = ", "))
  }
  
  ci_pct <- round(conf_level * 100)
  ci_part <- paste0(
    ci_pct, "\\% CI [",
    format_num_latex(res$conf_low, digits = ci_digits), ", ",
    format_num_latex(res$conf_high, digits = ci_digits), "]"
  )
  
  paste(beta_part, stat_part, p_part, ci_part, sep = ", ")
}

latex_model_terms <- function(model,
                              terms,
                              beta_labels = NULL,
                              conf_level = 0.95,
                              conf_method = "Wald",
                              include_ci = FALSE) {
  if (is.null(beta_labels)) {
    beta_labels <- terms
  }
  stats::setNames(
    vapply(seq_along(terms), function(i) {
      latex_model_term(
        model = model,
        term = terms[i],
        beta_label = beta_labels[i],
        conf_level = conf_level,
        conf_method = conf_method,
        include_ci = include_ci
      )
    }, FUN.VALUE = character(1)),
    terms
  )
}

extract_afex_term <- function(model, term) {
  anova_table <- model$anova_table
  if (is.null(anova_table) || !term %in% rownames(anova_table)) {
    stop(sprintf("Term '%s' was not found in the afex ANOVA table.", term))
  }

  out <- as.data.frame(anova_table[term, , drop = FALSE])
  out$term <- term
  out$num_df <- out$`num Df`
  out$den_df <- out$`den Df`
  out$f_value <- out$F
  out$mse <- out$MSE
  out$pes <- out$pes
  out$p_value <- out$`Pr(>F)`

  out[, c("term", "num_df", "den_df", "f_value", "mse", "pes", "p_value")]
}

latex_afex_term <- function(model,
                            term,
                            effect_label = NULL,
                            include_label = FALSE,
                            f_digits = 2,
                            df_digits = 0,
                            mse_digits = 2,
                            pes_digits = 3) {
  res <- extract_afex_term(model, term)
  effect_label <- if (is.null(effect_label)) term else effect_label

  stats_text <- paste0(
    "$F(",
    format_num_latex(res$num_df, digits = df_digits), ", ",
    format_num_latex(res$den_df, digits = df_digits),
    ") = ",
    format_num_latex(res$f_value, digits = f_digits),
    "$, $MSE = ",
    format_num_latex(res$mse, digits = mse_digits),
    "$, $p ",
    format_p_latex(res$p_value),
    "$, $\\eta_p^2 = ",
    format_num_latex(res$pes, digits = pes_digits, drop_leading_zero = TRUE),
    "$"
  )

  if (include_label) {
    paste0(effect_label, " ", stats_text)
  } else {
    stats_text
  }
}

analysis_results_path <- function(name, dir = "output/analysis") {
  file.path(dir, paste0(name, ".rds"))
}

ensure_project_dirs <- function() {
  dirs <- c("output/analysis", "output/figures", "output/models", "output/manuscript")
  invisible(lapply(dirs, function(dir) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }))
}

save_analysis_objects <- function(name, ..., dir = "output/analysis") {
  objects <- list(...)
  if (is.null(names(objects)) || any(names(objects) == "")) {
    stop("All objects passed to save_analysis_objects() must be named.")
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  saveRDS(objects, analysis_results_path(name, dir = dir))
  invisible(objects)
}

load_analysis_objects <- function(name, dir = "output/analysis") {
  path <- analysis_results_path(name, dir = dir)
  if (!file.exists(path)) {
    stop(sprintf("Could not find saved analysis objects at '%s'.", path))
  }
  readRDS(path)
}

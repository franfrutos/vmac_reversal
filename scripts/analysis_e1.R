## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

# Loading packages----
source("scripts/functions.R")
source("scripts/load_data.R")

if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(dplyr,# Data wrangling
       lme4,# LMMs
       lmerTest, # P-values in summary
       marginaleffects,# Model predictions and conditional effects
       hypr,# Set the contrast matrix
       sjPlot,# Tables
       afex,# ANOVA
       Rmisc,# Function to averaged the data in within-subject designs
       effectsize, # Cohen's d
       ggplot2, # Plotting
       cowplot,
       emmeans
)

ensure_project_dirs()

# Contrast for Singleton distractor condition (VMAC and AC effects):
HcRep <- hypr(
  VMAC = High ~ Low,
  AC = Low ~ Absent,
  levels = c("High", "Low", "Absent")
)

HcRepRev <- hypr(
  VMAC = High ~ Low, # Negative VMAC effect when there is a reversal
  AC = Low ~ Absent, # New low is compared against absent
  levels = c("High", "Low", "Absent")
)

# RT analysis: ----
# Acquisition model: 
d_RT <- filter_data(raw, f.absent = F) %>% filter(Phase == "Acquisition")
d_RT$log_RT <- log(d_RT$rt)
d_RT$Singleton <- factor(d_RT$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT$Singleton) <- HcRep

f.l <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT
  )

summary(f.l) 

f.l2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT
  )

summary(f.l2) # Maximal model

f.l2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT
  )

summary(f.l2_log) # Log-maximal model

anova(f.l2, f.l2_log) # |AIC_Log - AIC_Linear| = 200.3

# Selected model predictions
predictions(
  f.l2_log,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(f.l2_log) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
avg_comparisons(
  f.l2_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f.l2_log) ^ 2) / 2) - exp(lo + (sigma(f.l2_log) ^ 2) / 2)
)

# Reversal model:
d_RT_rev <- filter_data(raw, f.absent = F) %>% filter(Phase == "Reversal")
d_RT_rev$log_RT <- log(d_RT_rev$rt)
d_RT_rev$Singleton <- factor(d_RT_rev$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT_rev$Singleton) <- HcRepRev

fRev.l <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT_rev
  )

summary(fRev.l) 

fRev.l2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT_rev
  )

summary(fRev.l2) 

fRev.l2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d_RT_rev
  )

summary(fRev.l2_log) # Maximal model 

anova(fRev.l2, fRev.l2_log) ## |AIC_Log - AIC_Linear| = 6

# Selected model predictions
predictions(
  fRev.l2_log,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(fRev.l2_log) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
avg_comparisons(
  fRev.l2_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(fRev.l2_log) ^ 2) / 2) - exp(lo + (sigma(fRev.l2_log) ^ 2) / 2)
)

# Difference between phases? ----
phase_aovE1 <- filter_data(raw, f.absent = T, phase = c("Acquisition", "Reversal")) %>%
   .[.$Block_num %in% 11:14, ] %>%
   dplyr::summarise(RT = mean(rt), .by = c("Phase", "ID", "Singleton")) %>%
   aov_ez(dv = "RT", within = c("Singleton", "Phase"), id = "ID", anova_table = list(es = "pes"))

phase_aovE1 # Significant interaction

# Check interaction: 
emmeans(phase_aovE1, ~ Singleton | Phase) %>% pairs()
emmeans(phase_aovE1, ~ Singleton | Phase) %>% pairs() %>% confint()

emm_options(opt.digits = F)
t_to_d(4.656, df_error = 90, paired = T) # cohen's d for Acquisition (0.49 95% CI[0.27, 0.71])
t_to_d(0.086, df_error = 90, paired = T) # cohen's d for Reversal (0.01 95% CI[-0.20, 0.22])

# ACC analysis----
# Acquisition model: 
d_acc <- filter_data(raw, f.absent = F, acc = F) %>% filter(Phase == "Acquisition") %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(ID, Block_num, Singleton))
d_acc$Singleton <- factor(d_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_acc$Singleton) <- HcRep

f.acc <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d_acc
  ) # Singular fit


f.acc2 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d_acc
  ) # Singular fit


f.acc3 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d_acc
  )

summary(f.acc3) # Maximal model

# Reversal:
d_acc_rev <- filter_data(raw, f.absent = F, acc = F) %>% filter(Phase == "Reversal") %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(ID, Block_num, Singleton))
d_acc_rev$Singleton <- factor(d_acc_rev$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_acc_rev$Singleton) <- HcRepRev

fRev.acc <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_acc_rev
  )

summary(fRev.acc) 

fRev.acc2 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_acc_rev
  )

summary(fRev.acc2)

fRev.acc3 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_acc_rev
  )

summary(fRev.acc3) # Maximal model

# Plots ----
# Model predictions:
preds <- get_predictions(d_RT, d_RT_rev, f.l2_log, fRev.l2_log)
comps <- get_comparisons(d_RT, d_RT_rev, f.l2_log, fRev.l2_log)

preds_palette_fill <- c("High" = "#F29CA5", "Low" = "#9DDCF2", "Absent" = "grey80")
preds_palette_line <- c("High" = "#C9333E", "Low" = "#228BCB", "Absent" = "grey50")


# Gráfico base con más espacio entre facets
plot1_predsRT <- ggplot(data = preds[["mod"]],
                        aes(y = estimate, x = Block_num, color = Singleton, fill = Singleton)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  facet_grid(. ~ Phase, scales = "free_x", space = "free_x")+
  geom_errorbar(linewidth = 0.5, data = preds[["raw"]],
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = preds[["raw"]], aes(x = as.numeric(Block_num)),
             position = position_dodge(0.5), shape = 21) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +coord_cartesian(ylim = c(600,900))+
  scale_color_manual(values = preds_palette_line) +
  scale_fill_manual(values = preds_palette_fill) +
  labs(y = "Reaction time (ms)", x = "Block") +
  theme_Publication(text_size = 10) +
  ggtitle("Instructions and trial-by-trial feedback") +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "lines"),
    legend.spacing.x = unit(0.2, 'cm'),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5
    ))

pdata_1mod_acq <- comps[["mod"]] %>% filter(Effect == "VMAC", Phase == "Acquisition")
pdata_1raw_acq <- comps[["raw"]] %>% filter(Effect == "VMAC", Phase == "Acquisition")

plot1_compACQ <- ggplot(data = pdata_1mod_acq,
                        aes(y = estimate, x = Block_num, color = Effect, fill = Effect)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  geom_errorbar(linewidth = 0.5, data = pdata_1raw_acq,
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = pdata_1raw_acq, aes(x = as.numeric(Block_num)),
             position = position_dodge(0.5), shape = 21) +
  scale_x_continuous(breaks = seq(2, 12, 2)) +
  scale_y_continuous(limits = c(-40, 45), breaks = seq(-40, 40, 20)) +
  coord_cartesian(ylim = c(-40, 41)) +
  scale_color_manual(values = "#762A83") +
  scale_fill_manual(values = "#A83BB3") +
  labs(y = "VMAC (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_Publication(text_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  )

pdata_1mod_rev <- comps[["mod"]] %>% filter(Effect == "VMAC", Phase == "Reversal")
pdata_1raw_rev <- comps[["raw"]] %>% filter(Effect == "VMAC", Phase == "Reversal")

plot1_compREV <- ggplot(data = pdata_1mod_rev,
                        aes(y = estimate, x = Block_num, color = Effect, fill = Effect)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  geom_errorbar(linewidth = 0.5, data = pdata_1raw_rev,
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = pdata_1raw_rev, aes(x = as.numeric(Block_num)),
             position = position_dodge(0.5), shape = 21) +
  scale_x_continuous(breaks = seq(2, 24, 2))+
  scale_y_continuous(limits = c(-40, 45), breaks = seq(-40, 40, 20)) +
  coord_cartesian(ylim = c(-40, 41)) +
  scale_color_manual(values = "#762A83") +
  scale_fill_manual(values = "#A83BB3") +
  labs(y = "VMAC (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_Publication(text_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 12),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0, unit = "pt")
  )


subplot_left <- ggplotGrob(plot1_compACQ + theme(legend.position = "none"))
subplot_right <- ggplotGrob(plot1_compREV + theme(legend.position = "none"))

final_plot_1 <- ggdraw(plot1_predsRT) +
  # Subplot izquierdo (ajusta x/y para posición)
  draw_plot(subplot_left,
            x = 0.3, y = .555,   # Posición (x: izquierda/derecha, y: arriba/abajo)
            width = 0.175, height = 0.25) +  # Tamaño relativo
  # Subplot derecho
  draw_plot(subplot_right,
            x = 0.565, y = .565,
            width = 0.165, height = 0.25)

final_plot_1

ggsave(
  "output/figures/fig-2.png",
  height = 15,
  width = 25,
  dpi = 1200,
  units = "cm"
)

# Point estimation:
points <- raw %>%
  filter(!is.na(contingency_rating)) %>%
  dplyr::summarise(points = min(contingency_rating),
            .by = ID) %>%
  pull(points)

# Is the proportion of points estimated with either distractor different from 50%?
t.test(points, mu = 50)

save_analysis_objects(
  "analysis_e1",
  acquisition_rt_model = f.l2_log,
  reversal_rt_model = fRev.l2_log,
  acquisition_acc_model = f.acc3,
  reversal_acc_model = fRev.acc3,
  phase_anova = phase_aovE1,
  predictions = preds,
  comparisons = comps,
  contingency_points = points,
  final_plot = final_plot_1
)

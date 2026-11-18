## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

source("scripts/functions.R")
source("scripts/load_data.R") # Data of Experiments 1, 2a, and 2b, for the within-block analysis

required_packages <- c(
  "dplyr", "lme4", "lmerTest", "marginaleffects", "hypr", "sjPlot",
  "afex", "Rmisc", "ggplot2", "papaja", "emmeans", "effectsize"
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

ensure_project_dirs()

eyeD <- read.csv("Input/df_le2.csv") %>%
  filter(fixationTimeout != 1, trialPropGoodSamples > .25) %>%
  select(-group) # Not used (groups come from revGroup), and marginaleffects does not allow this name

group_df <- read.csv("Input/summary_Le.csv") %>% mutate(ID = subNum)

eyeD <- left_join(eyeD, group_df[, c("ID", "revGroup")], by = "ID")



# Contrast for VMAC
eyeD$distractType <- ifelse(eyeD$distractType == 1, "High", "Low")
eyeD$Singleton <- factor(eyeD$distractType, levels = c("High", "Low"))
contrasts(eyeD$Singleton) <- contr.sum(2)/2
colnames(contrasts(eyeD$Singleton)) <- "VMAC"

# Contrast for Groups
eyeD<- eyeD %>% mutate(Group = case_when(
  revGroup == 1 ~ "Rev",
  revGroup == 2 ~ "Control",
  T~"RevFB"
))

# eyeD$Group <- factor(eyeD$Group, levels = c("Rev", "Control", "RevFB"))
# 
# Gcontrast <- hypr(
#   Control_Rev = Control ~ Rev,
#   Rev_RevFB = Rev ~ RevFB,
#   levels = c("Rev", "Control", "RevFB")
# )
# 
# contrasts(eyeD$Group) <- Gcontrast

eyeD[eyeD$phase == 2,] %>%
  dplyr::summarise(P = mean(omissionTrial),
                   .by = c(ID, distractLoc, Singleton)) %>%
  spread(Singleton, P) %>% mutate(VMAC = High - Low,
                                  Pos = ifelse(distractLoc %in% c(1, 2, 6), "Top", "Bottom")) %>%
  summarySEwithin2(
    measurevar = "VMAC",
    withinvars = c("Pos"),
    idvar = "ID"
  ) %>%
  ggplot(aes(Pos, VMAC)) +
  geom_point() +
  geom_errorbar(aes(ymin = VMAC - ci, ymax = VMAC + ci), width = 0)

# Feedback ----
#Acquisition ----
# Training:
fit_Le_FB <- glmer(
  omissionTrial~Singleton*scale(block) + (Singleton*scale(block)|ID),
  data = eyeD[eyeD$phase == 2 & eyeD$Group == "RevFB",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)


fit_Le_log_FB <- glmer(
  omissionTrial~Singleton*scale(log(block)) + (Singleton*scale(log(block))|ID),
  data = eyeD[eyeD$phase == 2 & eyeD$Group == "RevFB",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)
anova(fit_Le_FB, fit_Le_log_FB)
summary(fit_Le_log_FB)

# Reversal -----
fit_LeTest_FB <- glmer(
  omissionTrial~Singleton*scale(block) + (Singleton+scale(block)|ID),
  data = eyeD[eyeD$phase == 3 & eyeD$Group == "RevFB",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)
fit_LeTest_log_FB <- glmer(
  omissionTrial~Singleton*scale(log(block)) + (Singleton+scale(log(block))|ID),
  data = eyeD[eyeD$phase == 3 & eyeD$Group == "RevFB",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)
anova(fit_LeTest_FB, fit_LeTest_log_FB)

summary(fit_LeTest_FB)

# Nominal ----
#Acquisition ----
# Training:
fit_Le <- glmer(
  omissionTrial~Singleton*scale(block) + (Singleton*scale(block)|ID),
  data = eyeD[eyeD$phase == 2 & eyeD$Group == "Rev",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)


summary(fit_Le)
fit_Le_log <- glmer(
  omissionTrial~Singleton*scale(log(block)) + (Singleton*scale(log(block))|ID),
  data = eyeD[eyeD$phase == 2 & eyeD$Group == "Rev",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)
anova(fit_Le, fit_Le_log)
summary(fit_Le_log)

# Reversal -----
fit_LeTest <- glmer(
  omissionTrial~Singleton*scale(block) + (Singleton*scale(block)|ID),
  data = eyeD[eyeD$phase == 3 & eyeD$Group == "Rev",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)
fit_LeTest_log <- glmer(
  omissionTrial~Singleton*scale(log(block)) + (Singleton*scale(log(block))|ID),
  data = eyeD[eyeD$phase == 3 & eyeD$Group == "Rev",], control = glmerControl(optimizer = "bobyqa"),
  family=binomial()
)
anova(fit_LeTest, fit_LeTest_log)
summary(fit_LeTest)
tab_model(fit_LeTest, digits = 3, digits.re = 4)

#ANOVA -----
phase_aovLe_rev <- eyeD[eyeD$Group == "Rev",]%>%
  filter(
    (phase == 2 & block %in% 15:16) | 
      (phase == 3 & block %in% 1:2)) %>%
  dplyr::summarise(DistractorTrial = mean(omissionTrial), .by = c("phase", "ID", "distractType"))%>%
  aov_ez(dv = "DistractorTrial", within = c("distractType", "phase"), id = "ID", anova_table = list(es = "pes"))

phase_aovLe_rev

# Check interaction: 
emmeans(phase_aovLe_rev, ~ distractType | phase) %>% pairs()
emmeans(phase_aovLe_rev, ~ distractType | phase) %>% pairs() %>% confint()

emm_options(opt.digits = F)
t_to_d(-5.354, df_error = 28, paired = T, alternative = "two.sided") 
t_to_d(-1.146, df_error = 28, paired = T, alternative = "two.sided") 

phase_aovLe_revfb <- eyeD[eyeD$Group == "RevFB",]%>%
  filter(
    (phase == 2 & block %in% 15:16) | 
      (phase == 3 & block %in% 1:2)) %>%
  dplyr::summarise(DistractorTrial = mean(omissionTrial), .by = c("phase", "ID", "distractType"))%>%
  aov_ez(dv = "DistractorTrial", within = c("distractType", "phase"), id = "ID", anova_table = list(es = "pes"))

phase_aovLe_revfb

# Check interaction: 
emmeans(phase_aovLe_revfb, ~ distractType | phase) %>% pairs()
emmeans(phase_aovLe_revfb, ~ distractType | phase) %>% pairs() %>% confint()

emm_options(opt.digits = F)
t_to_d(-5.088, df_error = 28, paired = T, alternative = "two.sided") 
t_to_d(1.061, df_error = 28, paired = T, alternative = "two.sided") 


#Plot----
# Model predictions for plot:

TrainComp <- rbind(
  comparisons(
  fit_Le_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(block = unique,
                     ID = NA),
  re.form = NA,
  by = c("block")) %>% mutate(Phase = "Acquisition",
                              Group = "Rev"),
  comparisons(
  fit_Le_log_FB,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(block = unique,
                     ID = NA),
  re.form = NA,
  by = c("block")) %>% mutate(Phase = "Acquisition",
                              Group = "RevFB")
  )


TestComp_rev <- rbind(
  comparisons(
    fit_LeTest_log,
    variables = list(Singleton = "revsequential"),
    newdata = datagrid(block = unique,
                       ID = NA),
    re.form = NA,
    by = c("block")) %>% mutate(Phase = "Reversal",
                                Group = "Rev"),
  comparisons(
    fit_LeTest_FB,
    variables = list(Singleton = "revsequential"),
    newdata = datagrid(block = unique,
                       ID = NA),
    re.form = NA,
    by = c("block")) %>% mutate(Phase = "Reversal",
                                Group = "RevFB")
)

TestComp <- rbind(TrainComp, TestComp_rev)

# Raw data for plot:----

eyeD$Block <- create_epochs(eyeD$block)

raw_Comps <- eyeD %>% filter(phase != 1, Group %in% c("Rev", "RevFB")) %>%
  mutate(Phase = ifelse(phase == 2, "Acquisition", "Reversal")) %>%
  group_by(Group, ID, Phase, Block, Singleton) %>%
  dplyr::summarise(Prop = mean(omissionTrial)) %>%
  spread(Singleton, Prop) %>%
  mutate(estimate = High - Low) %>%
  summarySEwithin2(
    measurevar = "estimate",
    betweenvars = "Group",
    withinvars = c("Phase", "Block"),
    idvar = "ID"
  )
raw_Comps$Phase <- factor(raw_Comps$Phase, levels = c("Acquisition", "Reversal"))
raw_Comps$Block <- as.numeric(raw_Comps$Block) * 2 - .5
raw_Comps$block <- raw_Comps$Block

# Plot:
plotLe_preds <- ggplot(data = TestComp,
                       aes(y = estimate * 100, x = block,
                           color = Group, fill = Group, shape = Group,
                           linetype = Group)) +  
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low * 100, ymax = conf.high * 100), 
              alpha = 0.3, color = NA) +
  facet_grid(. ~ Phase, scales = "free_x", space = "free_x") +
  geom_errorbar(linewidth = 0.7, data = raw_Comps,
                aes(x = as.numeric(block), 
                    ymin = (estimate - se) * 100, 
                    ymax = (estimate + se) * 100),  
                position = position_dodge(0.5), width = 0.01,
                linetype = "solid") +
  geom_point(size = 2.05, data = raw_Comps, 
             aes(x = as.numeric(block), y = estimate * 100),  
             position = position_dodge(0.5)) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_color_manual(values = c(
    "RevFB" = "#762A83",  # violeta
    "Rev" = "#3B4994"  # violeta-azul
  ),
  labels = c("Nominal extinction", "Feedback")) +
  scale_fill_manual(values = c(
    "RevFB" = "#A83BB3",
    "Rev" = "#7E87D4"
  ), labels = c("Nominal extinction", "Feedback")) +
  scale_shape_manual(values = c(
    "RevFB" = 21,
    "Rev" = 22
  ),
  labels =  c("Nominal extinction", "Feedback")) +
  scale_linetype_manual(values = c(
    "RevFB" = "solid",
    "Rev" = "dashed"
  ), labels =  c("Nominal extinction", "Feedback")) +
  
  scale_y_continuous(
    breaks = seq(-30, 30, by = 10), 
    limits = c(-20, 30)
  ) +
  labs(
    y = expression(atop(bold("VMAC effect"), bold("(distractor difference %)"))),
    x = "Block"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_Publication(text_size = 14, base_size = 12) +
  ggtitle("Re-analysis from Le et al. (2024)") +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    panel.spacing.x = unit(0.5, "lines"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )
plotLe_preds

ggsave(
  "output/figures/fig-5.png",
  units = "cm",
  height = 15,
  width = 25
)

# Trial vs. block analysis:----
# Creating a trial variable (1:24) within each block of trials
combined_d <- bind_rows(
  raw %>%
    filter(Phase == "Reversal") %>% 
    group_by(ID, Block_num) %>%
    dplyr::mutate(trial = 1:n()) %>% # plyr (loaded with Rmisc) masks mutate()
    filter_data(., f.absent = T) %>%
    mutate(Experiment = "1") %>% select(ID, rt, Singleton, trial, Block_num, Experiment) %>%
    ungroup(),
  raw2 %>%
    filter(Phase == "Reversal") %>%
    group_by(ID, Block_num) %>%
    dplyr::mutate(trial = 1:n()) %>%
    filter_data(., f.absent = T, experiment = "e2") %>%
    select(ID, rt, Singleton, trial, Block_num, Experiment) %>%
    ungroup()
) %>%
  mutate(ID = paste0(ID, Experiment))


# Log-RTs and effect coding for Singleton
combined_d$log_RT <- log(combined_d$rt)
combined_d$Singleton <- factor(combined_d$Singleton, levels = c("High", "Low"))
contrasts(combined_d$Singleton) <- contr.sum(2)/2

# Fit the model: 
f_combined <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * Experiment + # Effect of Block
      Singleton * scale(log(trial)) * Experiment + # Effect of trial
      (Singleton*scale(log(Block_num)) + scale(trial) | ID), # Random part
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = combined_d
  )

summary(f_combined)

# Plot:
d_plot <- comparisons(
  f_combined,
  variables = "Singleton",
  newdata = datagrid(
    Singleton = c("High", "Low"),
    Experiment = c("1", "A", "B"),
    trial = seq(1, 24, 0.01),
    ID = NA
  ),
  by = c("Experiment", "trial"),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_combined) ^ 2) / 2) - exp(lo + (sigma(f_combined) ^ 2) / 2),
)

trial_data <- combined_d %>%
  dplyr::summarise(
    RT = mean(rt),
    .by = c(Experiment, ID, Singleton, trial)
  ) %>%
  tidyr::spread(
    Singleton, RT
  ) %>%
  mutate(VMAC = High - Low)

trial_data$trial <- create_epochs(trial_data$trial, 4)

group_data <- summarySEwithin2(
  trial_data,
  measurevar = "VMAC",
  withinvars = "trial",
  betweenvars = "Experiment",
  na.rm = T
)

group_data$trial <- (as.numeric(group_data$trial) - 1) * 4 + 2.5

trial_plot <- ggplot(d_plot, aes(trial, -estimate, color = Experiment, linetype = Experiment, shape = Experiment)) +
  geom_line(linewidth = 1) +
  geom_point(data=group_data, aes(y = VMAC, color = Experiment, x = as.numeric(trial))) +
  geom_errorbar(data=group_data, aes(y = VMAC, color = Experiment, x = as.numeric(trial),
                                     ymin = VMAC - se, ymax = VMAC + se),
                width = 0, linetype = "solid") +
  geom_point(data=group_data, aes(y = VMAC, color = Experiment, x = as.numeric(trial),
                                  fill = Experiment),
             size = 2.5) +
  geom_ribbon(aes(ymin = -conf.low, ymax = -conf.high, fill = Experiment), color = NA, alpha = .2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = seq(1,24, 1)) +
  scale_color_manual(values = c(
    "1" = "#762A83",  # violeta
    "A" = "#3B4994",  # violeta-azul
    "B" = "#A73061"   # magenta
  ),
  labels = c("1", "2a", "2b")) +
  scale_fill_manual(values = c(
    "1" = "#A83BB3",
    "A" = "#7E87D4",
    "B" = "#D86491"
  ),
  labels = c("1", "2a", "2b")) +
  scale_shape_manual(values = c(
    "1" = 21,
    "A" = 22,
    "B" = 23
  ),
  labels = c("1", "2a", "2b")) +
  scale_linetype_manual(values = c(
    "1" = "solid",
    "A" = "dashed",
    "B" = "twodash"
  ), labels = c("1", "2a", "2b")) +
  labs(x = "Trial",
       y = "VMAC effect (ms)") +
  theme_Publication(text_size = 14, base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    panel.spacing.x = unit(0.5, "lines"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(5, 2, 0, 2)
  )

fig_7 <- trial_plot

ggsave(
  "output/figures/fig-7.png",
  fig_7,
  units = "cm",
  height = 16,
  width = 25
)

save_analysis_objects(
  "re_analysis_le2024",
  nominal_acquisition_model = fit_Le_log,
  nominal_reversal_model = fit_LeTest_log,
  feedback_acquisition_model = fit_Le_log_FB,
  feedback_reversal_model = fit_LeTest_FB,
  phase_anova_nominal = phase_aovLe_rev,
  phase_anova_feedback = phase_aovLe_revfb,
  within_block_model = f_combined,
  figure_5 = plotLe_preds,
  figure_7 = fig_7
)

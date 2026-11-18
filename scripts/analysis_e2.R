## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

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
       ggpubr,
       cowplot,
       patchwork,
       rapidsplithalf,
       emmeans
)

ensure_project_dirs()

#Sociodemographic data ----
descriptive_stats_e2 <- descriptive_e2 %>%
  dplyr::group_by(Experiment) %>%
  dplyr::summarise(
    n_participants = dplyr::n(),
    age = mean(Edad),
    age_sd = sd(Edad),
    age_max = max(Edad),
    females = sum(Sexo == 2),
    university = sum(University == 0),
    .groups = 'drop'
  )

descriptive_stats_e2

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


# Get data for analysis: 
d2 <- filter_data(raw2, f.absent = F, experiment = "e2")
d2$Epoch <- create_epochs(d2$Block_num, 2)
d2$log_RT <- log(d2$rt)
d2$Singleton <- factor(d2$Singleton, levels = c("High", "Low", "Absent"))

# Experiment 2A ----
contrasts(d2$Singleton) <- HcRep
## RTs----
### Acquisition-----
f_2A_l <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton*scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ]
  ) # Singular fit

f_2A_l2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton+scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ]
  ) #Maximal model

f_2A_l2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ]
  )
aic_2a_acq <- anova(f_2A_l2,f_2A_l2_log)
aic_2a_acq
summary(f_2A_l2_log)

# Selected model predictions
pred_2a_acq <- predictions(
  f_2A_l2_log,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(f_2A_l2_log) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
comp_2a_acq <- avg_comparisons(
  f_2A_l2_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_2A_l2_log) ^ 2) / 2) - exp(lo + (sigma(f_2A_l2_log) ^ 2) / 2)
)

### Reversal -----
contrasts(d2$Singleton) <- HcRepRev

f_2A_rev <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton*scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "A", ]
  ) # Singlular fit

f_2A_rev2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton+scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "A", ]
  ) # Maximal model

f_2A_rev2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "A", ]
  )
aic_2a_rev <- anova(f_2A_rev2, f_2A_rev2_log)
aic_2a_rev
summary(f_2A_rev2_log)

# Selected model predictions
pred_2a_rev <- predictions(
  f_2A_rev2_log,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(f_2A_rev2_log) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
comp_2a_rev <- avg_comparisons(
  f_2A_rev2_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_2A_rev2_log) ^ 2) / 2) - exp(lo + (sigma(f_2A_rev2_log) ^ 2) / 2)
)

# VMAC in the last two blocks of the reversal stage:
comp_2a_rev_last <- avg_comparisons(
  f_2A_rev2_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = 23:24,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_2A_rev2_log) ^ 2) / 2) - exp(lo + (sigma(f_2A_rev2_log) ^ 2) / 2)
)
comp_2a_rev_last

# Accuracy----
## Acquisition----
d2_acc <- filter_data(raw2, f.absent = F, acc = F) %>% filter(Phase == "Acquisition", Experiment == "A") %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(ID, Block_num, Singleton))
d2_acc$Singleton <- factor(d2_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d2_acc$Singleton) <- HcRep

f.acc_2A_l <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  ) # Singular fit

f.acc_2A_l2 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  ) # Singular fit


f.acc_2A_l3 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  )

summary(f.acc_2A_l3) # Maximal model
#p-values > .144

tab_model(f_2A_l2_log, f.acc_2A_l3, digits = 3, digits.re = 4)


avg_predictions(
  f.acc_2A_l3,
  newdata = datagrid(
    Singleton = unique, 
    ID = NA
  ),
  re.form = NA         
)

## Reversal----
d2_acc <- filter_data(raw2, f.absent = F, acc = F) %>% filter(Phase == "Reversal", Experiment == "A") %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(ID, Block_num, Singleton))
d2_acc$Singleton <- factor(d2_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d2_acc$Singleton) <- HcRep

fRev_2A_acc <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  ) # Singular fit

fRev_2A_acc2 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  ) #Singular fit


fRev_2A_acc3 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  ) # Singular fit

fRev_2A_acc4 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (scale(Block_num)| ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  )

fRev_2A_acc5 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (1| ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  )

summary(fRev_2A_acc5) # Maximal model

tab_model(f_2A_rev2_log,fRev_2A_acc5, digits = 3, digits.re = 4)

avg_predictions(
  fRev_2A_acc5,
  newdata = datagrid(
    Singleton = unique, 
    ID = NA
  ),
  re.form = NA
)
avg_comparisons(
  fRev_2A_acc5,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA)

# ANOVA ----
phase_aovE2A <- filter_data(raw2[raw2$Experiment == "A",], f.absent = T) %>%
  filter(Block_num %in% 11:14) %>%
  dplyr::summarise(RT = mean(rt), .by = c("Phase", "ID", "Singleton")) %>%
  aov_ez(dv = "RT", within = c("Singleton", "Phase"), id = "ID", anova_table = list(es = "pes"))

phase_aovE2A

# Check interaction: 
emmeans(phase_aovE2A, ~ Singleton | Phase) %>% pairs()
emmeans(phase_aovE2A, ~ Singleton | Phase) %>% pairs() %>% confint()

emm_options(opt.digits = F)
t_to_d(3.703, df_error = 62, paired = T) # cohen's d for Acquisition (0.49 95% CI[0.27, 0.71])
t_to_d(0.344, df_error = 62, paired = T) # cohen's d for Reversal (0.01 95% CI[-0.20, 0.22])


# Experiment 2B ----
## RTs----
### Acquisition-----
f_2B_l <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton*scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ]
  ) # Singular fit

f_2B_l2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton+scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ]
  ) #Maximal model

f_2B_l2_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ]
  )
aic_2b_acq <- anova(f_2B_l2,f_2B_l2_log)
aic_2b_acq
summary(f_2B_l2_log)

# Selected model predictions
pred_2b_acq <- predictions(
  f_2B_l2_log,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(f_2B_l2_log) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
comp_2b_acq <- avg_comparisons(
  f_2B_l2_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_2B_l2_log) ^ 2) / 2) - exp(lo + (sigma(f_2B_l2_log) ^ 2) / 2)
)

### Reversal -----
contrasts(d2$Singleton) <- HcRepRev

f_2B_rev <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton*scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "B", ]
  ) # Singlular fit

f_2B_rev2 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton+scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "B", ]
  )

f_2B_rev3 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton| ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "B", ]
  ) # Maximal model

f_2B_rev3_log <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton| ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "B", ]
  )
aic_2b_rev <- anova(f_2B_rev3, f_2B_rev3_log)
aic_2b_rev
summary(f_2B_rev3_log)

# Selected model predictions
pred_2b_rev <- predictions(
  f_2B_rev3_log,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(f_2B_rev3_log) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
comp_2b_rev <- avg_comparisons(
  f_2B_rev3_log,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_2B_rev3_log) ^ 2) / 2) - exp(lo + (sigma(f_2B_rev3_log) ^ 2) / 2)
)


#plotting
comps <- get_comparisons(d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "A", ],
                         f_2A_l2_log,
                         f_2A_rev2_log)

ggplot(data = comps[["mod"]],
       aes(
         y = estimate,
         x = Block_num,
         color = Effect,
         fill = Effect
       )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ Phase, scales = "free_x") +
  geom_point(
    data = comps[["raw"]],
    aes(x = as.numeric(Block_num)),
    position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]],
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  #scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

## Accuracy----
### Acquisition----
d2_acc <- filter_data(raw2, f.absent = F, acc = F) %>% filter(Phase == "Acquisition", Experiment == "B") %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(ID, Block_num, Singleton))
d2_acc$Singleton <- factor(d2_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d2_acc$Singleton) <- HcRep

f.acc_2B_l <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  ) # Singular fit

f.acc_2B_l2 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  ) # Singular fit


f.acc_2B_l3 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  ) # Singular fit

f.acc_2B_l4 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Block_num | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  ) #Maximal model

summary(f.acc_2B_l4)
tab_model(f_2B_l2_log, f.acc_2B_l4, digits = 3, digits.re = 4)
avg_predictions(
  f.acc_2B_l4,
  newdata = datagrid(
    Singleton = unique, 
    ID = NA
  ),
  re.form = NA,
  by = "Singleton"
)

avg_comparisons(
  f.acc_2B_l4,  
  variables = list(Singleton = "revsequential"), 
  newdata = datagrid(
    Singleton = unique, 
    Block_num = unique, 
    ID = NA            
  ),
  re.form = NA,
  comparison = \(hi, lo) lo - hi
)

#Exploring interactive effects.
d2_acc <- filter_data(raw2, f.absent = F, acc = F) %>% filter(Phase == "Acquisition", Experiment == "B") %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(ID, Block_num, Singleton))
d2_acc$Singleton <- factor(d2_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d2_acc$Singleton) <- contr.sum(3)
contr.sum(levels(d2_acc$Singleton))

f2.acc_2B_l4 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Block_num | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n, # Magic trick
    data = d2_acc
  )

summary(f2.acc_2B_l4)
tab_model(f2.acc_2B_l4, digits = 3, digits.re = 4)
### Reversal----
d2_acc <- filter_data(raw2, f.absent = F, acc = F) %>% filter(Phase == "Reversal", Experiment == "B") %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(ID, Block_num, Singleton))
d2_acc$Singleton <- factor(d2_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d2_acc$Singleton) <- HcRep

fRev_2B_acc <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  ) # Singular fit

fRev_2B_acc2 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  ) #Singular fit


fRev_2B_acc3 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  ) # Singular fit

fRev_2B_acc4 <-
  glmer(
    ACC ~ Singleton * scale(Block_num) + (scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d2_acc
  )

summary(fRev_2B_acc4) # Maximal model
tab_model(f_2B_rev3_log, fRev_2B_acc4, digits = 3, digits.re = 4)

avg_predictions(
  fRev_2B_acc4,
  newdata = datagrid(
    Singleton = unique, 
    ID = NA
  ),
  re.form = NA
)

avg_comparisons(
  fRev_2B_acc4,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA)
# ANOVA ----
phase_aovE2B <- filter_data(raw2[raw2$Experiment == "B",], f.absent = T) %>%
  filter(Block_num %in% 11:14) %>%
  dplyr::summarise(RT = mean(rt), .by = c("Phase", "ID", "Singleton")) %>%
  mutate(Singleton = factor(Singleton, levels = c("High", "Low"))) %>%
  aov_ez(dv = "RT", within = c("Singleton", "Phase"), id = "ID", anova_table = list(es = "pes"))

phase_aovE2B

# Check interaction: 
emmeans(phase_aovE2B, ~ Singleton | Phase) %>% pairs()
emmeans(phase_aovE2B, ~ Singleton | Phase) %>% pairs() %>% confint()

emm_options(opt.digits = F)
t_to_d(3.976, df_error = 61, paired = T) # cohen's d for Acquisition (0.51 95% CI[0.24, 0.77])
t_to_d(2.234, df_error = 61, paired = T) # cohen's d for Reversal (0.29 95% CI[-0.03, 0.54])

# Plots ----
### Experiment 2A ----
# Get data
preds <- get_predictions(
  d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ],
  d2[d2$Phase == "Reversal" & d2$Experiment == "A", ],
  f_2A_l2_log,
  f_2A_rev2_log
)

comps <- get_comparisons(
  d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ],
  d2[d2$Phase == "Reversal" & d2$Experiment == "A", ],
  f_2A_l2_log,
  f_2A_rev2_log
)

# Color palette for predictions
preds_palette_fill <- c("High" = "#F29CA5", "Low" = "#9DDCF2", "Absent" = "grey80")
preds_palette_line <- c("High" = "#C9333E", "Low" = "#228BCB", "Absent" = "grey50")

# Prediction plot
plot2A_predsRT <- ggplot(data = preds[["mod"]],
                         aes(y = estimate, x = Block_num, color = Singleton, fill = Singleton)) +
  geom_line(linewidth = 1.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  facet_grid(. ~ Phase, scales = "free_x", space = "free_x") +
  geom_errorbar(linewidth = 0.5, data = preds[["raw"]],
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = preds[["raw"]],
             aes(x = as.numeric(Block_num)), position = position_dodge(0.5), shape = 22) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  coord_cartesian(ylim = c(600, 900)) +
  scale_color_manual(values = preds_palette_line) +
  scale_fill_manual(values = preds_palette_fill) +
  labs(y = "Reaction time (ms)", x = "Block") +
  theme_Publication(text_size = 10) +
  ggtitle("Nominal extinction") +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.position = "none",
    #panel.spacing.x = unit(2, "lines"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# Conditional VMAC effect data
pdata_2Amod_acq <- comps[["mod"]] %>% filter(Effect == "VMAC", Phase == "Acquisition")
pdata_2Araw_acq <- comps[["raw"]] %>% filter(Effect == "VMAC", Phase == "Acquisition")
pdata_2Amod_rev <- comps[["mod"]] %>% filter(Effect == "VMAC", Phase == "Reversal")
pdata_2Araw_rev <- comps[["raw"]] %>% filter(Effect == "VMAC", Phase == "Reversal")

# Plot for conditional VMAC (acquisition)
plot2A_compACQ <- ggplot(data = pdata_2Amod_acq,
                         aes(y = estimate, x = Block_num,
                             color = Effect, fill = Effect)) +
  geom_line(linewidth = 1.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  geom_errorbar(linewidth = 0.5, data = pdata_2Araw_acq,
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                color = "#3B4994", position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = pdata_2Araw_acq,
             aes(x = as.numeric(Block_num)), color = "#3B4994", fill = "#7E87D4",
             position = position_dodge(0.5), shape = 22) +
  scale_x_continuous(breaks = seq(2, 12, 2)) +
  scale_y_continuous(limits = c(-40, 45), breaks = seq(-40, 40, 20)) +
  coord_cartesian(ylim = c(-40, 41)) +
  scale_color_manual(values = "#3B4994") +
  scale_fill_manual(values = "#7E87D4") +
  labs(y = "VMAC (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_Publication(text_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Plot for conditional VMAC (reversal)
plot2A_compREV <- ggplot(data = pdata_2Amod_rev,
                         aes(y = estimate, x = Block_num,
                             color = Effect, fill = Effect)) +
  geom_line(linewidth = 1.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  geom_errorbar(linewidth = 0.5, data = pdata_2Araw_rev,
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                color = "#3B4994", position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = pdata_2Araw_rev,
             aes(x = as.numeric(Block_num)), color = "#3B4994", fill = "#7E87D4",
             position = position_dodge(0.5), shape = 22) +
  scale_x_continuous(breaks = seq(2, 24, 2)) +
  scale_y_continuous(limits = c(-40, 45), breaks = seq(-40, 40, 20)) +
  coord_cartesian(ylim = c(-40, 41)) +
  scale_color_manual(values = "#3B4994") +
  scale_fill_manual(values = "#7E87D4") +
  labs(y = "VMAC (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_Publication(text_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Writting conditional effect (aligned, so both insets have the same panel size)
insets <- cowplot::align_plots(plot2A_compACQ, plot2A_compREV, align = "hv", axis = "tblr")
subplot_left <- insets[[1]]
subplot_right <- insets[[2]]

final_plot_2A <- ggdraw(plot2A_predsRT) +
  draw_plot(subplot_left, x = 0.3, y = .525, width = 0.19, height = 0.275) +
  draw_plot(subplot_right, x = 0.565, y = .525, width = 0.19, height = 0.275)


###Experiment 2b----
# Predictions:
preds <- get_predictions(d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "B", ],
                         f_2B_l2_log,
                         f_2B_rev3_log)
comps <- get_comparisons(d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "B", ],
                         f_2B_l2_log,
                         f_2B_rev3_log)
# Prediction plot
plot2B_predsRT <- ggplot(data = preds[["mod"]],
                         aes(y = estimate, x = Block_num, color = Singleton, fill = Singleton)) +
  geom_line(linewidth = 1.1, linetype = "twodash") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  facet_grid(. ~ Phase, scales = "free_x", space = "free_x") +
  geom_errorbar(linewidth = 0.5, data = preds[["raw"]],
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = preds[["raw"]],
             aes(x = as.numeric(Block_num)), position = position_dodge(0.5), shape = 23) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  coord_cartesian(ylim = c(600, 900)) +
  scale_color_manual(values = preds_palette_line) +
  scale_fill_manual(values = preds_palette_fill) +
  labs(y = "Reaction time (ms)", x = "Block") +
  theme_Publication(text_size = 10) +
  ggtitle("Omission of instructions") +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1.2, "lines"),
    legend.spacing.x = unit(0.2, 'cm')
  )

# Conditional VMAC effect data
pdata_2Bmod_acq <- comps[["mod"]] %>% filter(Effect == "VMAC", Phase == "Acquisition")
pdata_2Braw_acq <- comps[["raw"]] %>% filter(Effect == "VMAC", Phase == "Acquisition")
pdata_2Bmod_rev <- comps[["mod"]] %>% filter(Effect == "VMAC", Phase == "Reversal")
pdata_2Braw_rev <- comps[["raw"]] %>% filter(Effect == "VMAC", Phase == "Reversal")

# Plot for conditional VMAC (acquisition)
plot2B_compACQ <- ggplot(data = pdata_2Bmod_acq,
                         aes(y = estimate, x = Block_num,
                             color = Effect, fill = Effect)) +
  geom_line(linewidth = 1.1, linetype = "twodash") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  geom_errorbar(linewidth = 0.5, data = pdata_2Braw_acq,
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                color = "#A73061", position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = pdata_2Braw_acq,
             aes(x = as.numeric(Block_num)), color = "#A73061", fill = "#D86491",
             position = position_dodge(0.5), shape = 23) +
  scale_x_continuous(breaks = seq(2, 12, 2)) +
  scale_y_continuous(limits = c(-40, 45), breaks = seq(-40, 40, 20)) +
  coord_cartesian(ylim = c(-40, 41)) +
  scale_color_manual(values = "#A73061") +
  scale_fill_manual(values = "#D86491") +
  labs(y = "VMAC (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_Publication(text_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.position = "none"
  )

# Plot for conditional VMAC (reversal)
plot2B_compREV <- ggplot(data = pdata_2Bmod_rev,
                         aes(y = estimate, x = Block_num,
                             color = Effect, fill = Effect)) +
  geom_line(linewidth = 1.1, linetype = "twodash") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3, color = NA) +
  geom_errorbar(linewidth = 0.5, data = pdata_2Braw_rev,
                aes(x = as.numeric(Block_num), ymin = estimate - se, ymax = estimate + se),
                color = "#A73061", position = position_dodge(0.5), width = 0.01) +
  geom_point(size = 2.05, data = pdata_2Braw_rev,
             aes(x = as.numeric(Block_num)), color = "#A73061", fill = "#D86491",
             position = position_dodge(0.5), shape = 23) +
  scale_x_continuous(breaks = seq(2, 24, 2)) +
  scale_y_continuous(limits = c(-40, 45), breaks = seq(-40, 40, 20)) +
  coord_cartesian(ylim = c(-40, 41)) +
  scale_color_manual(values = "#A73061") +
  scale_fill_manual(values = "#D86491") +
  labs(y = "VMAC (ms)") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  theme_Publication(text_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 12),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.position = "none"
  )

# Writting conditional effect (aligned, so both insets have the same panel size)
insets <- cowplot::align_plots(plot2B_compACQ, plot2B_compREV, align = "hv", axis = "tblr")
subplot_left <- insets[[1]]
subplot_right <- insets[[2]]

final_plot_2B <- ggdraw(plot2B_predsRT) +
  draw_plot(subplot_left, x = 0.3, y = .57, width = 0.19, height = 0.25) +
  draw_plot(subplot_right, x = 0.565, y = .57, width = 0.19, height = 0.25)

ggarrange(
  final_plot_2A, 
  final_plot_2B,
  ncol = 1,
  nrow = 2,
  heights = c(9, 10),  
  common.legend = TRUE,
  labels = c("A","B"),
  font.label = list(size = 20, face = "bold", color = "black"),
  legend = "bottom"
)

ggsave(
  "output/figures/fig-3.png",
  height = 30,
  width = 25,
  dpi = 1200,
  units = "cm"
)

# Awareness test ----
# Data
d_aw <- raw2[!is.na(raw2$Phase) & raw2$Phase == "Awareness" & !raw2$ID %in% acc_ex2 & raw2$Singleton != "Absent",]

# Setting contrasts
d_aw$response_aw <- ifelse(d_aw$response == "c", 1, 0)
d_aw$signal <- ifelse(d_aw$Singleton == "High", .5, -.5)
d_aw$Experiment <- factor(d_aw$Experiment, levels = unique(d_aw$Experiment))
contrasts(d_aw$Experiment) <- contr.sum(2)/2

# Hierarchical Probit SDT model
fit_awA <- glmer(response_aw~signal+(signal|ID),data=d_aw[d_aw$Experiment == "A",],
                 family=binomial(link="probit"), control = glmerControl(optimizer = "bobyqa"))

# Model results
summary(fit_awA)

# CI for response criteria
-avg_predictions(fit_awA,
           type = "link",
           re.form = NA)

# CI for d'
avg_slopes(fit_awA,
       variables = "signal",
       type = "link",
       re.form = NA)

fit_awB <- glmer(response_aw~signal+(signal|ID),data=d_aw[d_aw$Experiment == "B",],
                 family=binomial(link="probit"), control = glmerControl(optimizer = "bobyqa"))

summary(fit_awB)

# CI for response criteria
-avg_predictions(fit_awB,
                 type = "link",
                 re.form = NA)

# CI for d'
avg_slopes(fit_awB,
           variables = "signal",
           type = "link",
           re.form = NA)


fit_aw <- glmer(response_aw~signal*Experiment+(signal|ID),data=d_aw,family=binomial(link="probit"),
                control = glmerControl(optimizer = "bobyqa"))

summary(fit_aw)

# Plot
roc_plot <- get_roc_probit(fit_aw)

# Individual differences in awareness ----
d2 <- filter_data(raw2, f.absent = F, experiment = "e2") %>% filter(Experiment == "B", Singleton != "Absent")
d2$Epoch <- create_epochs(d2$Block_num, 2)
d2$log_RT <- log(d2$rt)
d2$Singleton <- factor(d2$Singleton, levels = c("High", "Low"))
#d2$Singleton <- factor(d2$Singleton, levels = c("High", "Low", "Absent"))
d2$ID <- as.factor(d2$ID)
#contrasts(d2$Singleton) <- HcRepRev
contrasts(d2$Singleton) <- contr.sum(2)/2

ranef_df <- coef(fit_awB)$ID %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID")%>%
  mutate(ID = as.factor(ID)) %>%
  left_join(d2, ., by = "ID")

f_d <-
  lmer(
    log_RT ~ Singleton*signal*scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = ranef_df[ranef_df$Phase == "Reversal" & ranef_df$Experiment == "B", ]
  )

summary(f_d)

# VMAC effect when d' = 0
avg_comparisons(
  f_d,
  variables = "Singleton",
  datagrid(
    signal = 0,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_d) ^ 2) / 2) - exp(lo + (sigma(f_d) ^ 2) / 2),
)

# Reliability:
set.seed(2025) # Random splits
# d-prime
(splithalf_dprimes <- rapidsplit(data = d_aw[d_aw$Experiment == "B",], 
                                subjvar = "ID", 
                                aggvar = "response_aw",
                                diffvars = "signal", # 1: present (signal); 0: absent (noise)
                                aggfunc = dprimesfunc,
                                splits = 5000))

# VMAC effect: 
(rel <- rapidsplit(data = ranef_df[ranef_df$Phase == "Reversal" & ranef_df$Experiment == "B", ],
                   subjvar = "ID",
                   aggvar = "rt",
                   diffvars = "Singleton", # High minus low
                   aggfunc = "means",
                   splits = 5000))


# Expected attenuation
expected_attenuation <- (1 - sqrt(rel$r * splithalf_dprimes$r)) * 100
expected_attenuation

# Maximal expected attenuation
expected_attenuation_low <- (1 - sqrt(rel$ci[1] * splithalf_dprimes$ci[1])) * 100
expected_attenuation_low

# Maximal expected attenuation
expected_attenuation_high <- (1 - sqrt(rel$ci[2] * splithalf_dprimes$ci[2])) * 100
expected_attenuation_high

# Plot:
new_data <- expand.grid(
  Singleton = unique(ranef_df$Singleton),
  Block_num =  13:24,
  signal = unique(ranef_df$signal),
  ID = NA# marginaleffects needs a valid value for ID, letter we will set re.form 
)

# Step 3: Generate Predictions
# Generate predictions for the specific levels of 'signal' for each participant
preds <- avg_comparisons(f_d, newdata = new_data, re.form = NA, variables = "Singleton",
                         by = "signal", 
                         comparison  = \(hi, lo) exp(hi + (sigma(f_d) ^ 2) / 2) - exp(lo + (sigma(f_d) ^ 2) / 2))

dprimes <- coef(fit_awB)$ID %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID")%>%
  mutate(ID = as.factor(ID)) 

raw_effects <- ranef_df %>% filter(Phase == "Reversal") %>%
  dplyr::summarise(rt = mean(rt), .by = c(ID, Singleton)) %>%
  spread(Singleton, rt) %>% mutate(VMAC = Low - High) %>%
  left_join(., dprimes, by = "ID")

group_effects <- cbind(
  summarySE(
    data=raw_effects,
    measurevar = "VMAC"
  ) %>%
    mutate(
      CI_vmac = ci
    ) %>%
    select(c(VMAC, CI_vmac)),
  summarySE(
    data=raw_effects,
    measurevar = "signal"
  ) %>%
    mutate(
      CI_signal = ci
    ) %>%
    select(c(signal, CI_signal))
)

d_plot <- ggplot(preds, aes(x = signal, y = -estimate)) +
  geom_line(color = "#A73061", linewidth = 1, linetype = "twodash") +
  geom_point(data=raw_effects, aes(y = -VMAC, x = signal), alpha = 2/3, color = "#A73061", shape = 23) +
  geom_errorbar(data=group_effects, aes(y = -VMAC, ymax = -VMAC + CI_vmac, ymin = -VMAC - CI_vmac),
                color = "#A73061", size = 0.5, width = 0) +
  geom_errorbarh(data=group_effects, aes(y = -VMAC, xmax = signal + CI_signal, xmin = signal - CI_signal),
                color = "#A73061", size = 0.5, height = 0) +
  geom_point(data=group_effects, aes(y = -VMAC, x = signal), color = "#A73061", fill = "#D86491",
             size = 3, shape = 23) +
  geom_ribbon(aes(ymin = -conf.low, ymax = -conf.high), color = NA, alpha = .2,
              fill = "#A73061") +
  scale_y_continuous(breaks = seq(-60, 90, 30)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = expression(bolditalic(d) * "'"),
       y = "VMAC effect (ms)") +
  theme_Publication()


# Between experiment comparison ----
combined_d <- bind_rows(
  filter_data(raw, f.absent = T) %>% 
    filter(Phase == "Reversal") %>% 
    mutate(Experiment = "1") %>% select(ID, rt, Singleton, Block_num, Experiment),
  filter_data(raw2, f.absent = T, experiment = "e2") %>%
    filter(Phase == "Reversal") %>%
    select(ID, rt, Singleton, Block_num, Experiment),
) %>%
  mutate(ID = paste0(ID, Experiment))

combined_d$log_RT <- log(combined_d$rt)
combined_d$Singleton <- factor(combined_d$Singleton, levels = c("High", "Low"))
contrasts(combined_d$Singleton) <- contr.sum(2)/2

f_combined <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * Experiment + (Singleton*scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = combined_d
  )

summary(f_combined)
f_combined_rt <- f_combined

pred_combined <- predictions(
  f_combined,
  newdata = datagrid(
    Experiment = unique,
    ID = NA,
    Block_num = 13:24
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(f_combined) ^ 2) / 2),
  by = "Experiment"
)

# Plot for comparison: 
comps <- comparisons(
  f_combined,
  variables = "Singleton",
  newdata = datagrid(
    Experiment = unique,
    ID = NA,
    Block_num = 13:24
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(f_combined) ^ 2) / 2) - exp(lo + (sigma(f_combined) ^ 2) / 2),
  by = c("Experiment", "Block_num")
)

combined_d$Epoch <- create_epochs(combined_d$Block_num, 2)
raw_d <- combined_d %>%
  dplyr::summarise(RT = mean(rt),
            .by = c(Experiment, ID, Epoch, Singleton)) %>%
  spread(Singleton, RT) %>%
  mutate(VMAC = High - Low)

group_data <- summarySEwithin2(
  data = raw_d,
  measurevar = "VMAC",
  withinvars = "Epoch",
  betweenvars = "Experiment"
)

group_data$Block_num <- (as.numeric(group_data$Epoch) + 6) *2 -0.5

comp_plot <- ggplot(data=comps, aes(y=-estimate, x=Block_num, color = Experiment, linetype = Experiment, shape = Experiment)) +
  geom_line(linewidth = 1) +
  geom_errorbar(data = group_data, aes(y = VMAC, ymin = VMAC- se, ymax=VMAC+se),
                width = 0, position = position_dodge(.5), linetype = "solid") +
  geom_point(data = group_data, aes(y = VMAC, fill = Experiment), size = 2.5,
             position = position_dodge(.5)) +
  geom_ribbon(aes(ymin = -conf.low, ymax = -conf.high, fill = Experiment), color = NA,
              alpha = .2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = 13:24) +
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
  labs(y = "VMAC effect (ms)", x = "Block") +
  theme_Publication()

fig_4 <- ggarrange(
  ggarrange(
    roc_plot, d_plot,
    labels = c("A", "B"),
    ncol = 2,
    widths = c(1, 1),
    align = "hv"
  ),
  comp_plot,
  labels = c("", "C"),
  ncol = 1,
  heights = c(1, 1.3)
)

ggsave(
  "output/figures/fig-4.png",
  fig_4,
  units = "cm",
  height = 20,
  width = 18
)

# Accuracy:
# Data:
combined_d <- bind_rows(
  filter_data(raw, f.absent = T, acc = F) %>% 
    filter(Phase == "Reversal") %>% 
    mutate(Experiment = "1") %>% select(ID, correct, Singleton, Block_num, Experiment),
  filter_data(raw2, f.absent = T, experiment = "e2", acc = F) %>%
    filter(Phase == "Reversal") %>%
    select(ID, correct, Singleton, Block_num, Experiment),
) %>%
  mutate(ID = paste0(ID, Experiment))

# Set Contrasts:
combined_d <- combined_d %>%
  dplyr::summarise(ACC = mean(correct), n = n(), .by = c(Experiment, ID, Block_num, Singleton))
combined_d$Singleton <- factor(combined_d$Singleton, levels = c("High", "Low"))
contrasts(combined_d$Singleton) <- contr.sum(2)/2

# Fit model: 
f_combined <-
  glmer(
    ACC ~ Singleton * scale(Block_num) * Experiment + (Singleton * scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = combined_d
  )

summary(f_combined)
f_combined_acc <- f_combined

save_analysis_objects(
  "analysis_e2",
  descriptive_stats = descriptive_stats_e2,
  e2a_acquisition_rt_model = f_2A_l2_log,
  e2a_reversal_rt_model = f_2A_rev2_log,
  e2b_acquisition_rt_model = f_2B_l2_log,
  e2b_reversal_rt_model = f_2B_rev3_log,
  e2a_acquisition_acc_model = f.acc_2A_l3,
  e2a_reversal_acc_model = fRev_2A_acc5,
  e2b_acquisition_acc_model = f.acc_2B_l4,
  e2b_reversal_acc_model = fRev_2B_acc4,
  e2a_phase_anova = phase_aovE2A,
  e2b_phase_anova = phase_aovE2B,
  awareness_model_a = fit_awA,
  awareness_model_b = fit_awB,
  awareness_model_combined = fit_aw,
  awareness_vmac_model = f_d,
  awareness_reliability = list(
    dprime_split_half = splithalf_dprimes,
    vmac_split_half = rel,
    expected_attenuation = expected_attenuation,
    expected_attenuation_low = expected_attenuation_low,
    expected_attenuation_high = expected_attenuation_high
  ),
  comparison_rt_model = f_combined_rt,
  comparison_acc_model = f_combined_acc,
  figure_3a = final_plot_2A,
  figure_3b = final_plot_2B,
  figure_4 = fig_4,
  roc_plot = roc_plot,
  awareness_plot = d_plot,
  comparison_plot = comp_plot,
  # Values reported in the text
  reported = list(
    demographics = descriptive_e2 %>%
      dplyr::group_by(Experiment) %>%
      dplyr::summarise(
        n = dplyr::n(),
        females = sum(Sexo == 2),
        age_min = min(Edad),
        age_max = max(Edad),
        age_mean = mean(Edad),
        age_sd = sd(Edad),
        malaga = sum(University == 0),
        .groups = 'drop'
      ) %>%
      as.data.frame(),
    exclusions = list(e2a = exclusions_e2a, e2b = exclusions_e2b),
    aic_diff = c(
      e2a_acquisition = abs(diff(aic_2a_acq$AIC)),
      e2a_reversal = abs(diff(aic_2a_rev$AIC)),
      e2b_acquisition = abs(diff(aic_2b_acq$AIC)),
      e2b_reversal = abs(diff(aic_2b_rev$AIC))
    ),
    predictions = list(
      e2a_acquisition = as.data.frame(pred_2a_acq),
      e2a_reversal = as.data.frame(pred_2a_rev),
      e2b_acquisition = as.data.frame(pred_2b_acq),
      e2b_reversal = as.data.frame(pred_2b_rev),
      combined_by_experiment = as.data.frame(pred_combined)
    ),
    comparisons = list(
      e2a_acquisition = as.data.frame(comp_2a_acq),
      e2a_reversal = as.data.frame(comp_2a_rev),
      e2a_reversal_last_two_blocks = as.data.frame(comp_2a_rev_last),
      e2b_acquisition = as.data.frame(comp_2b_acq),
      e2b_reversal = as.data.frame(comp_2b_rev)
    )
  )
)


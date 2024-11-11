source("scripts/functions.R") # This will load functions and data into the workspace

if (!require(pacman)) {
  install.packages(pacman)
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
       ggplot2 # Plotting
)


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

# RT ----
# Get data for analysis: 
d2 <- filter_data(raw2, f.absent = F, experiment = "e2")
d2$Epoch <- create_epochs(d2$Block_num, 2)
d2$log_RT <- log(d2$rt)
d2$Singleton <- factor(d2$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d2$Singleton) <- HcRep

# Experiment 2A
f_r <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ]
  )

summary(f_r)

contrasts(d2$Singleton) <- HcRepRev

f_rev <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "A", ]
  )

summary(f_rev)

contrasts(d2$Singleton) <- HcRep

# Experiment 2B
f_rb <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "B" & d2$Block_num > 6, ]
  )

f_rb <-
  lmer(
    log_RT ~ Singleton + (Singleton | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Acquisition" & d2$Experiment == "B" & d2$Block_num > 6, ]
  )

summary(f_rb)

contrasts(d2$Singleton) <- HcRepRev

f_revb <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = d2[d2$Phase == "Reversal" & d2$Experiment == "B", ]
  )

summary(f_revb)

# Accuracy (do properly)----
d_acc <- filter_data(raw2, f.absent = F, acc = F, experiment = "e2") %>% filter(Phase != "Awareness")
d_acc$Singleton <- factor(d_acc$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_acc$Singleton) <- HcRep

# Trick to improve performance (evaluating likelihood as binomial(p, n), which gave same results)
d_opt <- d_acc %>%
  dplyr::summarise(correct = mean(correct),
                   n = n(),
                   .by = c(Experiment, ID, Phase, Block_num, Singleton))
  
f_r_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (Singleton * scale(log(Block_num)) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_opt[which(d_opt$Phase == "Acquisition" & d_opt$Experiment == "A"), ]
  ) # Maximal model

summary(f_r_acc) # Singular fit

# Dropping interaction from RE
f_r_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_opt[which(d_opt$Phase == "Acquisition" & d_opt$Experiment == "A"), ]
  ) 

summary(f_r_acc) # Singular fit

# Dropping Block from RE
f_r_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    weights = n,
    data = d_opt[which(d_opt$Phase == "Acquisition" & d_opt$Experiment == "A"), ]
  ) 

summary(f_r_acc) # Maximal feasible model

contrasts(d_acc$Singleton) <- HcRepRev

f_rev_acc <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (1 | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc[which(d_acc$Phase == "Reversal" & d_acc$Experiment == "A"), ]
  )

summary(f_rev_acc)

contrasts(d_acc$Singleton) <- HcRep

f_r_accb <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (scale(log(Block_num)) | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc[which(d_acc$Phase == "Acquisition" & d_acc$Experiment == "B"), ]
  )

summary(f_r_accb)

contrasts(d_acc$Singleton) <- HcRepRev

f_rev_accb <-
  glmer(
    correct ~ Singleton * scale(log(Block_num)) + (1 | ID),
    control = glmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    family = binomial(),
    data = d_acc[which(d_acc$Phase == "Reversal" & d_acc$Experiment == "B"), ]
  )

summary(f_rev_accb)

plot_predictions(
  f_rev_acc,
  newdata = insight::get_data(f_r_acc),
  condition = c("Singleton"),
  re.form = NA
)

plot_comparisons(
  f_r_accb,
  variables = list(Singleton = "revsequential"),
  condition = c("Block_num"),
  re.form = NA
)


# Is there persistence at the beginning of the reversal phase compared to the end of acquisition? ----

d2[d2$Block_num %in% 11:14, ] %>%
  mutate(Experiment = ifelse(Experiment == "A", "Experiment 2A", "Experiment 2B"),
         Phase = ifelse(Phase == "Rewarded", "Acquisition", Phase),
         Singleton = case_when(
           T~Singleton
         )) %>%
  dplyr::summarise(RT = mean(rt), .by = c(University, Experiment, ID, Phase, Singleton)) %>%
  mutate(RT =RT)-> d_aov

d_aov$Experiment <- factor(d_aov$Experiment, levels = c("Experiment 2A", "Experiment 2B"))

afex::aov_ez(data = d_aov %>% filter(Singleton != "Absent", Experiment == "Experiment 2A"), id = "ID", dv = "RT",
             within = c("Phase", "Singleton"), anova_table = list(es = "pes")) # Significant interaction

afex::aov_ez(data = d_aov %>% filter(Singleton != "Absent", Experiment == "Experiment 2B"), id = "ID", dv = "RT",
             within = c("Phase", "Singleton"), anova_table = list(es = "pes")) # Non-significant interaction


# Comparison of experiments plot
d_aov %>%
  mutate(ID = as.factor(paste(ID, Experiment))) %>%
  filter(Singleton != "Absent") %>%
  bind_rows(.,
    filter_data(raw, f.absent = F) %>% mutate(Experiment = "Experiment 1") %>%
      filter(Block_num %in% 11:14) %>%
      mutate(
             Phase = ifelse(Phase == "Rewarded", "Acquisition", Phase),
             Singleton = case_when(
               Phase == "Reversal" & Singleton == "Low"~"High",
               Phase == "Reversal" & Singleton == "High"~"Low",
               T~Singleton
             )) %>%
      dplyr::summarise(RT = mean(rt), .by = c(Experiment, ID, Phase, Singleton)) %>%
  filter(Singleton != "Absent") %>%
  mutate(
    ID = as.factor(paste(ID)),
    Singleton = as.factor(Singleton),
    Phase = as.factor(Phase),
    #Experiment = factor(Experiment, levels = c("Experiment 1", "Experiment 2A","Experiment 2B")),
    RT = RT
  )) %>%
  spread(Singleton, RT) %>%
  mutate(VMAC = High-Low)-> merge_data

merge_data$Experiment <- factor(merge_data$Experiment, levels = c("Experiment 1", "Experiment 2A","Experiment 2B"))
#merge_data$Singleton <- factor(merge_data$Singleton, levels = c("High", "Low"))


plot <- afex::aov_ez(data = merge_data, id = "ID", dv = "VMAC",
               within = c("Phase"), between = "Experiment", anova_table = list(es = "pes")) %>%
  afex_plot(., x="Experiment", trace = "Phase", error = "CMO", mapping = c("color", "shape"),
            dodge = .1,
            data_alpha = 0.25,
            data_arg = list(
              position = 
                ggplot2::position_jitterdodge(
                  jitter.width = .05, 
                  jitter.height = 0, 
                  dodge.width = 0.1  ## needs to be same as dodge
                ))) +
  theme_Publication() +
  scale_y_continuous(breaks = seq(-200, 200, 50)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "VMAC effect (ms)") +
  scale_color_brewer(palette = "Set1")

plot[["layers"]][[3]] <- NULL
plot[["layers"]][[5]] <- NULL
plot
ggsave(
  "plots/reversal_exps.png",
  height = 15,
  width = 17,
  dpi = 1200,
  units = "cm"
) 

# Plots ----
# Experiment 2A
# Predictions:
preds <- get_predictions(d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "A", ],
                         f_r, f_rev)
ggplot(data = preds[["mod"]],
       aes(
         y = estimate,
         x = Block_num,
         color = Singleton,
         fill = Singleton
       )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ Phase, scales = "free_x") +
  geom_point(data = preds[["raw"]],
             aes(x = as.numeric(Block_num)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = preds[["raw"]],
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(600, 900, 50)) +
  scale_color_brewer(palette = "Set1", name = "Singleton type") +
  scale_fill_brewer(palette = "Set1", name = "Singleton type") +
  labs(y = "Response time (ms)", x = "Block") +
  theme_Publication(text_size = 10) + 
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/preds2A_RT.png",
  height = 15,
  width = 20,
  dpi = 1200,
  units = "cm"
) 

# Conditional effects: 
comps <- get_comparisons(d2[d2$Phase == "Acquisition" & d2$Experiment == "A", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "A", ],
                         f_r, f_rev)

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
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(x = as.numeric(Block_num)),
    position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/cond2A_RT.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
) 

# Experiment 2B
# Predictions
predsB <- get_predictions(d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "B", ],
                         f_rb, f_revb)

ggplot(data = predsB[["mod"]],
       aes(
         y = estimate,
         x = Block_num,
         color = Singleton,
         fill = Singleton
       )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ Phase, scales = "free_x") +
  geom_point(data = predsB[["raw"]],
             aes(x = as.numeric(Block_num)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = predsB[["raw"]],
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(600, 900, 50)) +
  scale_color_brewer(palette = "Set1", name = "Singleton type") +
  scale_fill_brewer(palette = "Set1", name = "Singleton type") +
  labs(y = "Response time (ms)", x = "Block") +
  theme_Publication(text_size = 10) + 
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/preds2B_RT.png",
  height = 15,
  width = 20,
  dpi = 1200,
  units = "cm"
) 

comps <- get_comparisons(d2[d2$Phase == "Acquisition" & d2$Experiment == "B", ],
                         d2[d2$Phase == "Reversal" & d2$Experiment == "B", ],
                         f_rb, f_revb)

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
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(x = as.numeric(Block_num)),
    position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]] %>% mutate(Phase = ifelse(Phase == "Acquisition", "Acquisition",
                                                    "Reversal")),
    aes(
      x = as.numeric(Block_num),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_brewer(palette = "Set1", name = "Contrast") +
  scale_fill_brewer(palette = "Set1", name = "Contrast") +
  labs(y = "Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "plots/cond2B_RT.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
) 

# Awareness ----
d_aw <- raw2[raw2$Phase == "Awareness" & !raw2$ID %in% acc_ex2 & raw2$Singleton != "Absent",]
d_aw$response_aw <- ifelse(d_aw$response == "c", 1, 0)
d_aw$signal <- ifelse(d_aw$Singleton == "High", .5, -.5)
d_aw$Experiment <- factor(d_aw$Experiment, levels = unique(d_aw$Experiment))
contrasts(d_aw$Experiment) <- contr.sum(2)/2

fit_aw <- glmer(response_aw~signal*Experiment+(signal|ID),data=d_aw,family=binomial(link="probit"),
                control = glmerControl(optimizer = "bobyqa"))

summary(fit_aw)

get_roc_probit(fit_aw) 

ggsave("plots/plots_roc.png", height = 15, width = 15, dpi = 900,
       units = "cm")

# Individual differences in awareness ----
# Try with and without random effects
fit_awb <- glmer(response_aw~signal+(signal|ID),data=d_aw[d_aw$Experiment == "A",],family=binomial(link="probit"),
                control = glmerControl(optimizer = "bobyqa"))


HcRepRev <- hypr(
  VMAC = High ~ Low, # Negative VMAC effect when there is a reversal
  AC = Low ~ Absent, # New low is compared against absent
  levels = c("High", "Low", "Absent")
)

d2 <- filter_data(raw2, f.absent = F, experiment = "e2") %>% filter(Experiment == "A", Singleton != "Absent")
d2$Epoch <- create_epochs(d2$Block_num, 2)
d2$log_RT <- log(d2$rt)
d2$Singleton <- factor(d2$Singleton, levels = c("High", "Low"))
#d2$Singleton <- factor(d2$Singleton, levels = c("High", "Low", "Absent"))
d2$ID <- as.factor(d2$ID)
#contrasts(d2$Singleton) <- HcRepRev
contrasts(d2$Singleton) <- contr.sum(2)/2

ranef_df <- coef(fit_awb)$ID %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID")%>%
  mutate(ID = as.factor(ID)) %>%
  left_join(d2, ., by = "ID")

#ranef_df$Group<- ifelse(ranef_df$signal >= quantile(ranef_df$signal, probs = .5), "High", "Low")

f_d <-
  lmer(
    log_RT ~ Singleton*signal*scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = ranef_df[ranef_df$Phase == "Reversal" & ranef_df$Experiment == "A", ]
  )

summary(f_d)

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

coef(fit_awb)$ID %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ID")%>%
  mutate(ID = as.factor(ID)) -> a
a<- ranef_df %>% filter(Phase == "Reversal") %>%
  dplyr::summarise(rt = mean(rt), .by = c(ID, Singleton)) %>%
  spread(Singleton, rt) %>% mutate(VMAC = Low - High) %>%
  left_join(., a, by = "ID")

ggplot(preds, aes(x = signal, y = estimate)) +
  geom_line() +
  geom_point(data=a, aes(y = VMAC, x = signal)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), color = NA, alpha = .2) +
  scale_y_continuous(breaks = seq(-60, 90, 30)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = expression(bolditalic(d) * "'"),
       y = "VMAC effect (ms)",
       color = "Singleton") +
  theme_Publication()

# Paso 2: Crear las Columnas de Hits y Falsas Alarmas en los Datos
# Suponiendo que `response_aw` == 1 indica una respuesta correcta y `signal` == 1 indica la presencia de la señal
d_awb <- d_aw[d_aw$Experiment == "B" & d_aw$Singleton != "Absent",] %>%
  mutate(
    hits = ifelse(response_aw == 1 & signal == 0.5, 1, 0),           # Acierto cuando hay señal (signal == 1)
    fas = ifelse(response_aw == 1 & signal == -0.5, 1, 0)    # Falsa alarma cuando no hay señal (signal == 0)
  )

# Paso 3: Calcular las Tasas de Hits y Falsas Alarmas por Participante
rates_df <- d_awb %>%
  group_by(ID) %>%
  dplyr::summarise(
    hit_rate = mean(hits),
    fa_rate = mean(fas),
    n = n()
  )

calculate_d_prime <- function(hit_rate, fa_rate) {
  # Aplicar correcciones de continuidad según el número de ensayos
  hit_rate <- ifelse(hit_rate == 1, 1 - 1 / (2 * 40), 
                     ifelse(hit_rate == 0, 1 / (2 * 40), hit_rate))
  
  fa_rate <- ifelse(fa_rate == 1, 1 - 1 / (2 * 40), 
                    ifelse(fa_rate == 0, 1 / (2 * 40), fa_rate))
  
  # Calcular d' (d-prime)
  z_hit_rate <- qnorm(hit_rate)
  z_fa_rate <- qnorm(fa_rate)
  d_prime <- z_hit_rate - z_fa_rate
  
  return(d_prime)
}

# Paso 4: Calcular d' para Cada Participante
rates_df <- rates_df %>%
  mutate(d_prime = calculate_d_prime(hit_rate, fa_rate))

# Ver los resultados
print(rates_df)

sample_df <- rates_df %>% 
  mutate(ID = as.factor(ID)) %>%
  left_join(d2, ., by = "ID") %>%
  drop_na()

f_d2 <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num))*scale(d_prime) + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = sample_df[sample_df$Phase == "Reversal" & sample_df$Experiment == "B", ]
  )

summary(f_d2)

ggplot(rates_df, aes(d_prime)) +
  geom_histogram(aes(y=after_stat(density))) +
  geom_density()


# Between experiment comparison ----
combined_d <- bind_rows(
  filter_data(raw, f.absent = F) %>% 
    filter(Phase == "Reversal") %>% 
    mutate(Experiment = "1") %>% select(ID, rt, Singleton, Block_num, Experiment),
  filter_data(raw2, f.absent = F, experiment = "e2") %>%
    filter(Phase == "Reversal") %>%
    select(ID, rt, Singleton, Block_num, Experiment),
)

combined_d$Epoch <- create_epochs(combined_d$Block_num, 2)
combined_d$log_RT <- log(combined_d$rt)
combined_d$Singleton <- factor(combined_d$Singleton, levels = c("High", "Low"))
contrasts(d2$Singleton) <- contr.sum(2)/2

f_combined <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) * Experiment + (Singleton*scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa', optCtrl=list(maxfun=1e5)),
    data = combined_d
  )

summary(f_combined)


## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

# Load packages ----
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(
  dplyr,
  Rmisc
)

source("scripts/model_functions.R")

library(dplyr)
library(Rmisc)  # para summarySEwithin
# Suponemos que create_epochs() ya está definido en tu entorno

# ———————— Cargar datos y asignar grupos ————————
d <- read.csv("Input/modeling/le2023_e2.csv")
group_df <- read.csv("Input/modeling/summary_Le.csv") %>%
  mutate(ID = subNum)

d <- left_join(d, group_df[, c("ID", "revGroup")], by = "ID") %>%
  mutate(Group = case_when(
    revGroup == 1 ~ "Rev",
    revGroup == 2 ~ "Control",
    revGroup == 3 ~ "RevFB"
  )) %>%
  filter(phase %in% 2:3) %>%
  mutate(
    reward       = case_when(
      phase == 2 ~trialPay / 500,
      phase == 3 & Group == "RevFB" & distractType == 1~10/500,
      phase == 3 & Group == "RevFB" & distractType == 2~500/500,
      phase == 3 & (Group =="Rev" | Group == "Control")~0
      ),
    change_phase = ifelse(phase == 3 & block == 1 & trial == 1, 1, 0)
  )%>%
  group_by(
    ID
  ) %>%
  dplyr::mutate(trial = 1:n()) %>%
  ungroup()


# ———————— Preparar listas vacías ————————
dStan  <- list()
dPlots <- list()
group_metadata <- tibble::tribble(
  ~group_code, ~group_id, ~group_label,
  "Rev", "nominal_extinction", "Nominal extinction",
  "RevFB", "reversal_with_feedback", "Reversal with feedback"
)

# ———————— Iterar sobre cada grupo ————————
groups_to_process <- c("RevFB", "Rev")

for (grp in groups_to_process) {
  # 1) Filtrar y reenumerar IDs
  d_grp <- d %>%
    filter(Group == grp) %>%
    mutate(ID = as.integer(as.factor(ID))) %>%
    arrange(ID, phase, trial) 
  
  
  dat <- list(
    I         = table(d_grp$ID)[1],
    J         = length(unique(d_grp$ID)),
    N         = nrow(d_grp),
    idx_start = which(!duplicated(d_grp$ID)),
    singleton = d_grp$distractType,
    reward    = d_grp$reward,
    valid     = ifelse(d_grp$timeout != 1 & d_grp$trialPropGoodSamples > 0.25, 1, 0),
    change    = d_grp$change_phase,
    phase     = d_grp$phase - 1,
    Omission  = d_grp$omissionTrial
  )

  # # IDs consecutivos 1..J
  # subj_idx <- as.integer(factor(d_grp$ID))
  # 
  # # ensayos válidos (tu misma lógica)
  # valid_vec <- ifelse(d_grp$timeout != 1 & d_grp$trialPropGoodSamples > 0.25, 1L, 0L)
  # 
  # # filtro para Omission (solo válidos, en orden)
  # Omission_valid <- as.integer(d_grp$omissionTrial[valid_vec == 1L])
  # 
  # # Singleton 1/2 y con nombre que espera Stan
  # Singleton_vec <- as.integer(d_grp$distractType)  # asegúrate que es 1/2
  # 
  # # phase en 1/2 (¡no restes 1!)
  # phase_vec <- as.integer(d_grp$phase) -1            # 1 = adquisición, 2 = reversión
  # 
  # # trial como real (centrado fuera si quieres)
  # trial = (d_grp$trial - mean(d_grp$trial)) / (max(d$trial) - min(d$trial))
  # 
  # dat <- list(
  #   N         = nrow(d_grp),
  #   N_valid   = sum(valid_vec),
  #   J         = length(unique(subj_idx)),
  #   subj      = subj_idx,
  #   trial     = trial,
  #   Singleton = Singleton_vec,
  #   reward    = as.numeric(d_grp$reward),
  #   valid     = valid_vec,
  #   change    = as.integer(d_grp$change_phase),
  #   phase     = phase_vec,
  #   Omission  = Omission_valid
  # )

  
  dStan[[grp]] <- dat
  
  # 3) Construir resumen para plots
  d_plot_grp <- d_grp %>%
    mutate(
      trial = rep(1:(384 + 192), length.out = nrow(.)),
      Epoch = create_epochs(trial, 48)
    ) %>%
    filter(timeout != 1 & trialPropGoodSamples > 0.25) %>%
    dplyr::summarise(
      P = mean(omissionTrial),
      .by = c(ID, phase, Epoch, distractType)
    ) %>%
    dplyr::rename("Singleton" = distractType) %>%
    mutate(
      Singleton = ifelse(Singleton == 1, "High", "Low"),
      Data  = grp,
      group_id = group_metadata$group_id[group_metadata$group_code == grp],
      group_label = group_metadata$group_label[group_metadata$group_code == grp]
    ) %>%
    summarySEwithin(
      measurevar  = "P",
      withinvars  = c("Singleton", "phase", "Epoch"),
      idvar       = "ID"
    )  %>%
    mutate(Phase = ifelse(phase == 2, "Acquisition", "Reversal"))
  
  dPlots[[grp]] <- d_plot_grp
}



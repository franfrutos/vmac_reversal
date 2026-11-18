## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

# Installing (if needed) and loading the packages ----
if (!require(pacman)) {
  install.packages("pacman")
  library(pacman)
}

p_load(here, dplyr, tidyr, afex)

setwd(here::here())

if (!exists("report_exclusions")) {
  report_exclusions <- function(d, acc) {
    # Only the trials of the task (no practice or awareness test)
    d <- d[which(d$Phase %in% c("Acquisition", "Reversal")), ]
    incorrect <- (1 - nrow(d[which(d$correct == 1 & !d$ID %in% acc), ]) /
      nrow(d[which(!d$ID %in% acc), ])) * 100
    outliers <- (1 - nrow(d[which(d$correct == 1 & !d$ID %in% acc &
      (d$rt > 150 & d$rt < 1800)), ]) /
      nrow(d[which(d$correct == 1 & !d$ID %in% acc), ])) * 100

    cat(
      paste0(
        "Incorrect responses: ", round(incorrect, 2), "%",
        "\nOutliers RTs: ", round(outliers, 2), "%"
      )
    )
    invisible(list(incorrect = incorrect, outliers = outliers))
  }
}

# Experiment 1----
raw <- read.csv("Input/data_e1.csv") %>%
  mutate(Phase = ifelse(Phase == "Reward", "Acquisition", Phase))

# Check 
# There are participants with less observations than the max number of observations?
check <- raw %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID, Phase) %>% filter(Phase %in% c("Acquisition", "Reversal")) %>%
  filter(n < 288 | is.na(n))

length(unique(check$ID)) # N = 0 participants with missing data

length(unique(raw$ID)) # N = 92 participants

# Participants with less than .7 of accuracy?

acc_ex <- raw %>%
  filter(Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(correct)) %>% filter(mean_ACC < .7) %>% pull(ID)

length(unique(acc_ex)) # 1 participants with less than .7 of accuracy that will be excluded from the analysis

# How many participants?
raw %>%
  filter(!ID %in% acc_ex, Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  group_by(Phase) %>%
  dplyr::count() # 91 participants with complete observations and ACC > .7

# Report trial exclusions:
print("Experiment 1:")
exclusions_e1 <- report_exclusions(raw, acc_ex)

# Experiments 2a and 2b ----
raw2 <- read.csv("Input/data_e2.csv") %>%
  mutate(Phase = ifelse(Phase == "Reward", "Acquisition", Phase),
         Singleton = case_when(
           Singleton == "High" & Phase == "Reversal"~"Low",
           Singleton == "Low" & Phase == "Reversal"~"High",
           T~Singleton
         ))

check <- raw2 %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID, Phase) %>% filter(Phase %in% c("Acquisition", "Reversal")) %>%
  filter(n < 288 | is.na(n))

length(unique(check$ID)) # N = 0 participants with missing data

# How many participants?
acc_ex2 <- raw2 %>%
  filter(Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(correct)) %>%
  filter(mean_ACC < .7) %>%
  pull(ID) # All participants above cut-off


raw2 %>%
  filter(!ID %in% acc_ex2, Phase %in% c("Acquisition", "Reversal")) %>%
  group_by(Experiment, ID, Phase) %>%
  dplyr::count() %>%
  group_by(Experiment, Phase) %>%
  dplyr::count() 

# Report trial exclusions for experiment 2A and 2B:
print("Experiment 2A:")
exclusions_e2a <- report_exclusions(raw2[which(raw2$Experiment == "A"),], acc_ex2)

print("Experiment 2B:")
exclusions_e2b <- report_exclusions(raw2[which(raw2$Experiment == "B"),], acc_ex2)


# Merged dataset ----
raw_all <- bind_rows(
  raw %>%
    mutate(Experiment = "Experiment 1") %>% 
    select(ID, Block_num, trial_num, Phase, Singleton, rt, correct, Experiment),
  raw2 %>% mutate(Experiment = ifelse(Experiment == "A", "Experiment 2A", "Experiment 2B")) %>%
    select(ID, Block_num, trial_num, Phase, Singleton, rt, correct, Experiment)
)

descriptive_e1 <- read.csv("Input/descriptive_e1.csv")
descriptive_e2 <- read.csv("Input/descriptive_e2.csv")

descriptive_e2 <- descriptive_e2 %>%
  mutate(ID = factor(ID)) %>%
  inner_join(raw2 %>% mutate(ID = factor(ID)) %>%
               select(ID, Experiment, University) %>% distinct(), by = "ID")

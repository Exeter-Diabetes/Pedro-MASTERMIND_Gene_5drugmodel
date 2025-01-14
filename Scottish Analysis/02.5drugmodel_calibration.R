##############################################################################
#
# This file tests the 5 drug model in the T&F dataset 
#
##############################################################################

######
# load libraries
library(tidyverse)
library(rms)
library(MatchIt)

######

##
## Change this line of code to load the Scottish dataset
##

# load dataset
original.dataset <- readRDS("../Synthetic Data/full.synth.dataset.rds")

# load the model
load("../fivedrugmodel_5knot_share_20230823.Rdata")  # name: m1.5.final

# load functions
source("functions.R")


######


######
# Intercept vs Intercept + Slope vs Nothing test

## SGLT2
closed_loop_test_results_SGLT2 <- closedtest_continuous_function(
  cohort = "SGLT2 subcohort",
  dataset = original.dataset %>% filter(drugclass == "SGLT2"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

### save object
saveRDS(closed_loop_test_results_SGLT2, "02.closed_loop_test_results_SGLT2.rds")

## GLP1
closed_loop_test_results_GLP1 <- closedtest_continuous_function(
  cohort = "GLP1 subcohort",
  dataset = original.dataset %>% filter(drugclass == "GLP1"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

### save object
saveRDS(closed_loop_test_results_GLP1, "02.closed_loop_test_results_GLP1.rds")

## DPP4
closed_loop_test_results_DPP4 <- closedtest_continuous_function(
  cohort = "DPP4 subcohort",
  dataset = original.dataset %>% filter(drugclass == "DPP4"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

### save object
saveRDS(closed_loop_test_results_DPP4, "02.closed_loop_test_results_DPP4.rds")

## SU
closed_loop_test_results_SU <- closedtest_continuous_function(
  cohort = "SU subcohort",
  dataset = original.dataset %>% filter(drugclass == "SU"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

### save object
saveRDS(closed_loop_test_results_SU, "02.closed_loop_test_results_SU.rds")

## TZD
closed_loop_test_results_TZD <- closedtest_continuous_function(
  cohort = "TZD subcohort",
  dataset = original.dataset %>% filter(drugclass == "TZD"),
  original_model = m1.5.final,
  outcome_name = "posthba1cfinal",
  p_value = 0.05
)

### save object
saveRDS(closed_loop_test_results_TZD, "02.closed_loop_test_results_TZD.rds")


# Make predictions for each treatment for all patients
interim.dataset <- original.dataset %>%
  cbind(
    pred.SGLT2 = predict_with_modelchoice_function(closed_loop_test_results_SGLT2, original.dataset %>% mutate(drugclass = "SGLT2")),
    pred.GLP1 = predict_with_modelchoice_function(closed_loop_test_results_GLP1, original.dataset %>% mutate(drugclass = "GLP1")),
    pred.DPP4 = predict_with_modelchoice_function(closed_loop_test_results_DPP4, original.dataset %>% mutate(drugclass = "DPP4")),
    pred.SU = predict_with_modelchoice_function(closed_loop_test_results_SU, original.dataset %>% mutate(drugclass = "SU")),
    pred.TZD = predict_with_modelchoice_function(closed_loop_test_results_TZD, original.dataset %>% mutate(drugclass = "TZD"))
  ) %>%
  mutate(
    pred.current = ifelse(drugclass == "SGLT2", pred.SGLT2,
                          ifelse(drugclass == "GLP1", pred.GLP1,
                                 ifelse(drugclass == "DPP4", pred.DPP4,
                                        ifelse(drugclass == "SU", pred.SU, pred.TZD))))
  )



######
# Concordant vs discordant calibration

### 
# Overall calibration

# select best treatments
for (i in 1:nrow(interim.dataset)) {
  
  # Best drug
  ## add value/name to interim.dataset
  interim.dataset$first_best_drug_value[i] <- sort(interim.dataset %>% select(pred.SGLT2, pred.GLP1, pred.DPP4, pred.SU, pred.TZD) %>% slice(i) %>% unlist())[1]
  interim.dataset$first_best_drug_name[i] <- gsub("pred.", "", names(sort(interim.dataset %>% select(pred.SGLT2, pred.GLP1, pred.DPP4, pred.SU, pred.TZD) %>% slice(i) %>% unlist())[1]))
  
  # Second best drug
  ## add value/name to interim.dataset
  interim.dataset$second_best_drug_value[i] <- sort(interim.dataset %>% select(pred.SGLT2, pred.GLP1, pred.DPP4, pred.SU, pred.TZD) %>% slice(i) %>% unlist())[2]
  interim.dataset$second_best_drug_name[i] <- gsub("pred.", "", names(sort(interim.dataset %>% select(pred.SGLT2, pred.GLP1, pred.DPP4, pred.SU, pred.TZD) %>% slice(i) %>% unlist())[2]))
  
}

# add concordant vs discordant label
interim.dataset <- interim.dataset %>%
  mutate(
    conc_disc_label = ifelse(drugclass == first_best_drug_name, "Concordant", "Discordant"),
    conc_disc_label_numeric = ifelse(drugclass == first_best_drug_name, 1, 0),
    hba1c_group_percentile = cut(prehba1c, quantile(prehba1c, prob = 0:20 / 20, names = FALSE), include = TRUE)
  )


# calculate benefit as overall measure 
## Take patient A (concordant) and patient B (discordant)
## To calculate predicted benefit (x-axis)
### Take the difference of predicted benefit (for patient A) on drug taken by patient A 
###   and predicted benefit (for patient A) on drug taken by patient B.
## To calculate observed benefit* (y-axis)
### Take the difference of observed response for Patient A and observed response for Patient B.


# check the formula for matching (only keep categorical variables with two or more unique values)
formula <- paste0("conc_disc_label_numeric ~ t2dmduration + prebmi + agetx + prealt + preegfr + pretotalcholesterol + prehdl")

cat_vars <- c("smoke", "imd5", "ncurrtx", "drugline")
for (var in cat_vars) {
  
  if (length(unique(interim.dataset %>% select(all_of(var)) %>% unlist())) > 1) {
    formula <- paste0(formula, "+", var)
  }
  
}

## Matching model
matching_model <- MatchIt::matchit(
  formula = as.formula(formula),
  data = interim.dataset,
  method = "nearest",
  distance = "mahalanobis", # try this next
  replace = TRUE,
  
  # caliper = 0.05,
  # antiexact = c("drugclass"),
  exact = c("sex", "hba1c_group_percentile", "first_best_drug_name")
)


# Number of groupings (edit this number depending on the number of concordant/discordant pairs you have (usually about 100 per group, more is better))
group_num = 5



# matched dataset
matched_interim.dataset <- MatchIt::get_matches(matching_model, data = interim.dataset)

matched_interim.dataset <- matched_interim.dataset %>%
  group_by(subclass) %>%
  mutate(
    discordant_drugclass = paste(drugclass[2], collapse = ","),
    calibration_obs = diff(posthba1cfinal)
  ) %>%
  ungroup() %>%
  distinct(subclass, .keep_all = TRUE) %>%
  mutate(
    calibration_interim = ifelse(discordant_drugclass == "SGLT2", pred.SGLT2, 
                                 ifelse(discordant_drugclass == "DPP4", pred.DPP4, 
                                        ifelse(discordant_drugclass == "GLP1", pred.GLP1, 
                                               ifelse(discordant_drugclass == "SU", pred.SU, pred.TZD)))),
    calibration_pred = calibration_interim - pred.current
  ) %>%
  select(calibration_pred, calibration_obs) %>%
  mutate(
    grouping = as.numeric(cut(calibration_pred, unique(quantile(calibration_pred, prob = 0:group_num / group_num, names = FALSE)), include = TRUE))
  )


# initiate vectors for values
coef <- rep(0, group_num)
coef_low <- rep(0, group_num)
coef_high <- rep(0, group_num)
mean <- rep(0, group_num)
n_vector <- rep(0, group_num)
benefit_col <- "calibration_pred"
# iterate through each group
for (i in 1:group_num) {
  
  # select patients in this group
  group.data <- matched_interim.dataset %>%
    filter(grouping == i)
  
  # calculate the mean benefit for this group
  mean[i] <- mean(group.data %>% select(all_of(benefit_col)) %>% unlist(), na.rm = TRUE)
  
  lm <- lm(as.formula("calibration_obs ~ calibration_pred"), group.data)
  
  
  # predictions
  predictions_vector = predict(lm, newdata = data.frame(calibration_pred = mean(group.data %>% select(all_of(benefit_col)) %>% unlist(), na.rm = TRUE)), interval = "confidence")
  
  # add coefficients
  coef[i] <- predictions_vector[1]
  coef_low[i] <- predictions_vector[2]
  coef_high[i] <- predictions_vector[3]
  n_vector[i] <- nrow(group.data)
  
}

overall_calibration_table <- data.frame(mean, coef, coef_low, coef_high, 
                                        n = n_vector, 
                                        total_conc = interim.dataset %>% filter(conc_disc_label == "Concordant") %>% nrow(),
                                        total_disconc = interim.dataset %>% filter(conc_disc_label == "Discordant") %>% nrow())






## treatments being compared
drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4")

## Iterating through combinations
# SGLT2 vs GLP1
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "GLP1")) %>%
  mutate(drugclass = factor(drugclass, levels = c("GLP1", "SGLT2"))) %>%
  mutate(
    benefit = pred.SGLT2 - pred.GLP1
  )

SGLT2_GLP1_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
SGLT2_GLP1_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# SGLT2 vs TZD
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "TZD")) %>%
  mutate(drugclass = factor(drugclass, levels = c("TZD", "SGLT2"))) %>%
  mutate(
    benefit = pred.SGLT2 - pred.TZD
  )

SGLT2_TZD_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
SGLT2_TZD_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# SGLT2 vs SU
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "SU")) %>%
  mutate(drugclass = factor(drugclass, levels = c("SU", "SGLT2"))) %>%
  mutate(
    benefit = pred.SGLT2 - pred.SU
  )

SGLT2_SU_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
SGLT2_SU_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# SGLT2 vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "DPP4")) %>%
  mutate(drugclass = factor(drugclass, levels = c("DPP4", "SGLT2"))) %>%
  mutate(
    benefit = pred.SGLT2 - pred.DPP4
  )

SGLT2_DPP4_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
SGLT2_DPP4_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# GLP1 vs TZD
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("GLP1", "TZD")) %>%
  mutate(drugclass = factor(drugclass, levels = c("TZD", "GLP1"))) %>%
  mutate(
    benefit = pred.GLP1 - pred.TZD
  )

GLP1_TZD_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
GLP1_TZD_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# GLP1 vs SU
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("GLP1", "SU")) %>%
  mutate(drugclass = factor(drugclass, levels = c("SU", "GLP1"))) %>%
  mutate(
    benefit = pred.GLP1 - pred.SU
  )

GLP1_SU_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
GLP1_SU_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# GLP1 vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("GLP1", "DPP4")) %>%
  mutate(drugclass = factor(drugclass, levels = c("DPP4", "GLP1"))) %>%
  mutate(
    benefit = pred.GLP1 - pred.DPP4
  )

GLP1_DPP4_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
GLP1_DPP4_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# TZD vs SU
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("TZD", "SU")) %>%
  mutate(drugclass = factor(drugclass, levels = c("SU", "TZD"))) %>%
  mutate(
    benefit = pred.TZD - pred.SU
  )

TZD_SU_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
TZD_SU_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# TZD vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("TZD", "DPP4")) %>%
  mutate(drugclass = factor(drugclass, levels = c("DPP4", "TZD"))) %>%
  mutate(
    benefit = pred.TZD - pred.DPP4
  )

TZD_DPP4_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
TZD_DPP4_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")

# SU vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SU", "DPP4")) %>%
  mutate(drugclass = factor(drugclass, levels = c("DPP4", "SU"))) %>%
  mutate(
    benefit = pred.SU - pred.DPP4
  )

SU_DPP4_5_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 5, "benefit")
SU_DPP4_3_conc_disc_object <- conc_disc_validation_function(combination.interim , "drugclass", 3, "benefit")


# Plot

# plot_calibration_concordant_discordant <- SGLT2_GLP1_5_conc_disc_object %>%
#   ggplot() +
#   geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red") +
#   geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "red") +
#   geom_abline(aes(intercept = 0, slope = 1)) +
#   geom_errorbar(aes(x = mean, y = coef, ymin = coef_low, ymax = coef_high)) +
#   labs(x = "CATE", y = "ATE", title = "SGLT2 - GLP1 concordant vs discordant") +
#   theme_bw()

# Output files
output_table <- SGLT2_GLP1_5_conc_disc_object %>%
  mutate(grouping = 5) %>%
  rbind(
    SGLT2_GLP1_3_conc_disc_object %>%
      mutate(grouping = 3)
  ) %>%
  mutate(
    drug1 = "SGLT2",
    drug2 = "GLP1"
  ) %>%
  rbind(
    
    rbind(
      SGLT2_TZD_5_conc_disc_object %>%
        mutate(grouping = 5),
      SGLT2_TZD_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "SGLT2",
        drug2 = "TZD"
      ),
    rbind(
      SGLT2_SU_5_conc_disc_object %>%
        mutate(grouping = 5),
      SGLT2_SU_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "SGLT2",
        drug2 = "SU"
      ),
    rbind(
      SGLT2_DPP4_5_conc_disc_object %>%
        mutate(grouping = 5),
      SGLT2_DPP4_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "SGLT2",
        drug2 = "DPP4"
      ),
    rbind(
      GLP1_TZD_5_conc_disc_object %>%
        mutate(grouping = 5),
      GLP1_TZD_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "GLP1",
        drug2 = "TZD"
      ),
    rbind(
      GLP1_SU_5_conc_disc_object %>%
        mutate(grouping = 5),
      GLP1_SU_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "GLP1",
        drug2 = "SU"
      ),
    rbind(
      GLP1_DPP4_5_conc_disc_object %>%
        mutate(grouping = 5),
      GLP1_DPP4_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "GLP1",
        drug2 = "DPP4"
      ),
    rbind(
      TZD_SU_5_conc_disc_object %>%
        mutate(grouping = 5),
      TZD_SU_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "TZD",
        drug2 = "SU"
      ),
    rbind(
      TZD_DPP4_5_conc_disc_object %>%
        mutate(grouping = 5),
      TZD_DPP4_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "TZD",
        drug2 = "DPP4"
      ),
    rbind(
      SU_DPP4_5_conc_disc_object %>%
        mutate(grouping = 5),
      SU_DPP4_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "SU",
        drug2 = "DPP4"
      )
  )

# save output table
saveRDS(overall_calibration_table, "02.5drugmodel_overall_calibration.rds")
saveRDS(output_table, "02.5drugmodel_calibration_conc_disc.rds")

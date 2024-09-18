##############################################################################
#
# This file tests the 5 drug model in the T&F dataset 
#
##############################################################################

######
# load libraries
library(tidyverse)
library(rms)

######

##
## Change this line of code to load the Scottish dataset
##

# load dataset
original.dataset <- readRDS("../Synthetic Data/full.synth.dataset.rds")

# load the model
load("../fivedrugmodel_5knot_share_20230823.Rdata")  # name: m1.5.final

######
# Function for doing the concordant vs discordant validation
conc_disc_validation_function <- function(data, drug_col, group_num, benefit_col) {
  
  # split benefit into groups
  interim.data <- data %>%
    mutate(
      quantile = ntile(data %>% select(all_of(benefit_col)), group_num)
    )
  
  # initiate vectors for values
  coef <- rep(0, group_num)
  coef_low <- rep(0, group_num)
  coef_high <- rep(0, group_num)
  mean <- rep(0, group_num)
  
  # iterate through each group
  for (i in 1:group_num) {
    
    # select patients in this group
    group.data <- interim.data %>%
      filter(quantile == i)
    
    # calculate the mean benefit for this group
    mean[i] <- mean(group.data %>% select(all_of(benefit_col)) %>% unlist(), na.rm = TRUE)
    
    # fit the regression
    ########## Need to check whether the categorical variables are elements, otherwise remove
    formula <- as.formula(paste0("posthba1cfinal ~ ", drug_col ," + sex + t2dmduration + prebmi + prehba1c + agetx + 
                prealt + preegfr + pretotalcholesterol + prehdl + ethnicity + smoke + imd5 + 
                hba1cmonth + ncurrtx + drugline"))
    
    lm <- glm(formula, data = group.data)
    
    # add coefficients
    coef[i] <- coef(lm)[2]
    coef_low[i] <- confint(lm)[2,1]
    coef_high[i] <- confint(lm)[2,2]
    
  }
  
  # return output
  return(data.frame(mean, coef, coef_low, coef_high))
  
}


######


# Make predictions for each treatment for all patients
interim.dataset <- original.dataset %>%
  cbind(
    pred.SGLT2 = predict(m1.5.final, original.dataset %>% mutate(drugclass = "SGLT2")),
    pred.GLP1 = predict(m1.5.final, original.dataset %>% mutate(drugclass = "GLP1")),
    pred.DPP4 = predict(m1.5.final, original.dataset %>% mutate(drugclass = "DPP4")),
    pred.SU = predict(m1.5.final, original.dataset %>% mutate(drugclass = "SU")),
    pred.TZD = predict(m1.5.final, original.dataset %>% mutate(drugclass = "TZD"))
  )


######
# Intercept vs Intercept + Slope vs Nothing test




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
    conc_disc_label = ifelse(drugclass == first_best_drug_name, "Concordant", "Discordant")
  )

# calculate benefit
## If concordant:
##  - Patient took best drug, therefore compare best 1 vs best 2
## If discordant:
##  - Patient took something else, therefore compare drug taken vs best 1

for (i in 1:nrow(interim.dataset)) {
  
  if (interim.dataset$conc_disc_label[i] == "Concordant") {
    # If patient is concordant
    ## Add predicted benefit
    interim.dataset$benefit[i] <- interim.dataset$second_best_drug_value[i] - interim.dataset$first_best_drug_value[i]
    
  } else {
    # If patient is discordant
    ## drug taken
    drug_taken = interim.dataset$drugclass[i]
    ## column needed for benefit calculation
    column_name = paste0("pred.", drug_taken)
    ## Add predicted benefit
    interim.dataset$benefit[i] <- interim.dataset %>% select(all_of(column_name)) %>% slice(i) %>% unlist() - interim.dataset$first_best_drug_value[i]
    
  }
  
}

# add grouping variable = this can be changed to include more or less groups
interim.dataset <- interim.dataset %>%
  mutate(
    quantile = ntile(benefit, 5)
  )

# Running function for overall tests
overall_10_conc_disc_object <- conc_disc_validation_function(interim.dataset, "conc_disc_label", 10, "benefit")
overall_5_conc_disc_object <- conc_disc_validation_function(interim.dataset, "conc_disc_label", 5, "benefit")
overall_3_conc_disc_object <- conc_disc_validation_function(interim.dataset, "conc_disc_label", 3, "benefit")


## treatments being compared
drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4")

## Iterating through combinations
# SGLT2 vs GLP1
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "GLP1")) %>%
  mutate(
    benefit = pred.SGLT2 - pred.GLP1
  )

SGLT2_GLP1_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
SGLT2_GLP1_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# SGLT2 vs TZD
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "TZD")) %>%
  mutate(
    benefit = pred.SGLT2 - pred.TZD
  )

SGLT2_TZD_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
SGLT2_TZD_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# SGLT2 vs SU
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "SU")) %>%
  mutate(
    benefit = pred.SGLT2 - pred.SU
  )

SGLT2_SU_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
SGLT2_SU_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# SGLT2 vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SGLT2", "DPP4")) %>%
  mutate(
    benefit = pred.SGLT2 - pred.DPP4
  )

SGLT2_DPP4_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
SGLT2_DPP4_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# GLP1 vs TZD
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("GLP1", "TZD")) %>%
  mutate(
    benefit = pred.GLP1 - pred.TZD
  )

GLP1_TZD_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
GLP1_TZD_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# GLP1 vs SU
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("GLP1", "SU")) %>%
  mutate(
    benefit = pred.GLP1 - pred.SU
  )

GLP1_SU_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
GLP1_SU_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# GLP1 vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("GLP1", "DPP4")) %>%
  mutate(
    benefit = pred.GLP1 - pred.DPP4
  )

GLP1_DPP4_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
GLP1_DPP4_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# TZD vs SU
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("TZD", "SU")) %>%
  mutate(
    benefit = pred.TZD - pred.SU
  )

TZD_SU_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
TZD_SU_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# TZD vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("TZD", "DPP4")) %>%
  mutate(
    benefit = pred.TZD - pred.DPP4
  )

TZD_DPP4_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
TZD_DPP4_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")

# SU vs DPP4
combination.interim <- interim.dataset %>%
  filter(drugclass %in% c("SU", "DPP4")) %>%
  mutate(
    benefit = pred.SU - pred.DPP4
  )

SU_DPP4_5_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 5, "benefit")
SU_DPP4_3_conc_disc_object <- conc_disc_validation_function(interim.dataset , "drugclass", 3, "benefit")


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
output_table <- overall_10_conc_disc_object %>%
  mutate(grouping = 10) %>%
  rbind(
    overall_5_conc_disc_object %>%
      mutate(grouping = 5),
    overall_3_conc_disc_object %>%
      mutate(grouping = 3)
  ) %>%
  mutate(
    drug1 = "Concordant",
    drug2 = "Discordant"
  ) %>%
  rbind(
    rbind(
      SGLT2_GLP1_5_conc_disc_object %>%
        mutate(grouping = 5),
      SGLT2_GLP1_3_conc_disc_object %>%
        mutate(grouping = 3)
    ) %>%
      mutate(
        drug1 = "SGLT2",
        drug2 = "GLP1"
      ),
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
saveRDS(output_table, "02.5drugmodel_calibration_conc_disc.rds")

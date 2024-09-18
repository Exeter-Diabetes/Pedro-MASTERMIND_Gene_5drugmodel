##############################################################################
#
# This file confirms the SNP effect and saves a table with effect
#
##############################################################################

######
# load libraries
library(tidyverse)

######

##
## Change this line of code to load the Scottish dataset
##

# load dataset
original.dataset <- readRDS("../Synthetic Data/full.synth.dataset.rds")


######

#
# The following code is designed to test one singular SNP at a time.
#
# Currently the code is set-up to test CYP2C9*2 or CYP2C9*3
#   In order to test difference SNPs, you can change the column being renamed 
#   to the pre-specified name 'testing_SNP' and the code should run smoothly.
#   You can also change the filtering so that it only testing the correct strata.
#   At the end of the analysis, a table is created containing the betas.
#   Don't forget to change the name of the output table so that they don't
#   overwrite each other.
#

# rename SNP variable
interim.dataset <- original.dataset %>%
  
  # Filter to only analysis patients with SU
  filter(drugclass == "SU") %>%
  
  # The code below create the variable since the synthetic dataset does not have it
  mutate(testing_SNP = rep(c(0, 1), each = n()/2))

  # # The code below renames the existing column, replace '_____' with the actual column
  # rename("testing_SNP" = "_________")


# running models testing effect size
## unadjusted model
unadjusted_model <- lm(posthba1cfinal ~ testing_SNP, data = interim.dataset)

## adjusted model
adjusted_model <- lm(posthba1cfinal ~ testing_SNP + prehba1c + agetx + sex, data = interim.dataset)


# save table with effect sizes
output_table <- as.data.frame(
  # combine coefficients with confidence intervals
  unadjusted_model$coefficients %>%
    cbind(
      confint(unadjusted_model, levels = 0.95)
    )
) %>%
  mutate(Model = "unadjusted") %>%
  # turn row names into a column
  rownames_to_column() %>%
  # add adjusted models
  rbind(
    as.data.frame(
      # combine coefficients with confidence intervals
      adjusted_model$coefficients %>%
        cbind(
          confint(adjusted_model, levels = 0.95)
        )
    ) %>%
      mutate(Model = "adjusted") %>%
      # turn row names into a column
      rownames_to_column()
  ) %>%
  # rename columns to understandable names
  rename(
    "Coefficient" = "rowname",
    "mean" = "."
  )

# save the file - change name appropriately
saveRDS(output_table, "01.snp_effect_SU_CYP2C92.rds")




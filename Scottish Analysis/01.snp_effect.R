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
  
  # The code below create the variable since the synthetic dataset does not have it 
  #   (remove after renaming the column with line below)
  mutate(testing_SNP = rep(c(0, 1), each = n()/2))

  # # The code below renames the existing column, replace '_____' with the actual column
  # KEEP THIS VARIABLE CONTINUOUS using 0 or 1. This will make sure the code doesn't break.
  # rename("testing_SNP" = "_________")

# Output file - going to continuously add elements to this table
output_table <- NULL



# running models testing effect size
## unadjusted model
unadjusted_model_SGLT2 <- glm(posthba1cfinal ~ testing_SNP, data = interim.dataset %>% filter(drugclass == "SGLT2"))
if(length(unique(interim.dataset %>% filter(drugclass == "SGLT2") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- as.data.frame(summary(unadjusted_model_SGLT2)$coefficients[2,] %>% t() %>% cbind(t(confint(unadjusted_model_SGLT2)[2,])) %>% as.data.frame() %>%
                       mutate(Model = "unadjusted", drug = "SGLT2", interaction = NA))
} else {
  output_table <- as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "unadjusted", drug = "SGLT2", interaction = NA))
}

unadjusted_model_GLP1 <- glm(posthba1cfinal ~ testing_SNP, data = interim.dataset %>% filter(drugclass == "GLP1"))
if(length(unique(interim.dataset %>% filter(drugclass == "GLP1") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(unadjusted_model_GLP1)$coefficients[2,] %>% t() %>% cbind(t(confint(unadjusted_model_GLP1)[2,])) %>% as.data.frame() %>%
                                  mutate(Model = "unadjusted", drug = "GLP1", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "unadjusted", drug = "GLP1", interaction = NA)))
}

unadjusted_model_TZD <- glm(posthba1cfinal ~ testing_SNP, data = interim.dataset %>% filter(drugclass == "TZD"))
if(length(unique(interim.dataset %>% filter(drugclass == "TZD") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(unadjusted_model_TZD)$coefficients[2,] %>% t() %>% cbind(t(confint(unadjusted_model_TZD)[2,])) %>% as.data.frame() %>%
                                  mutate(Model = "unadjusted", drug = "TZD", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "unadjusted", drug = "TZD", interaction = NA)))
}

unadjusted_model_SU <- glm(posthba1cfinal ~ testing_SNP, data = interim.dataset %>% filter(drugclass == "SU"))
if(length(unique(interim.dataset %>% filter(drugclass == "SU") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(unadjusted_model_SU)$coefficients[2,] %>% t() %>% cbind(t(confint(unadjusted_model_SU)[2,])) %>% as.data.frame() %>%
                                  mutate(Model = "unadjusted", drug = "SU", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "unadjusted", drug = "SU", interaction = NA)))
}

unadjusted_model_DPP4 <- glm(posthba1cfinal ~ testing_SNP, data = interim.dataset %>% filter(drugclass == "DPP4"))
if(length(unique(interim.dataset %>% filter(drugclass == "DPP4") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(unadjusted_model_DPP4)$coefficients[2,] %>% t() %>% cbind(t(confint(unadjusted_model_DPP4)[2,])) %>% as.data.frame() %>%
                                  mutate(Model = "unadjusted", drug = "DPP4", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "unadjusted", drug = "DPP4", interaction = NA)))
}

## adjusted simple model
adjusted_simple_model_SGLT2 <- lm(posthba1cfinal ~ testing_SNP + prehba1c, data = interim.dataset %>% filter(drugclass == "SGLT2"))
if(length(unique(interim.dataset %>% filter(drugclass == "SGLT2") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_model_SGLT2)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_model_SGLT2)[2,])) %>% as.data.frame() %>%
                                  mutate(Model = "adjusted_simple", drug = "SGLT2", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple", drug = "SGLT2", interaction = NA)))
}

adjusted_simple_model_GLP1 <- lm(posthba1cfinal ~ testing_SNP + prehba1c, data = interim.dataset %>% filter(drugclass == "GLP1"))
if(length(unique(interim.dataset %>% filter(drugclass == "GLP1") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_model_GLP1)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_model_GLP1)[2,])) %>% as.data.frame() %>%
                                  mutate(Model = "adjusted_simple", drug = "GLP1", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple", drug = "GLP1", interaction = NA)))
}

adjusted_simple_model_TZD <- lm(posthba1cfinal ~ testing_SNP + prehba1c, data = interim.dataset %>% filter(drugclass == "TZD"))
if(length(unique(interim.dataset %>% filter(drugclass == "TZD") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_model_TZD)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_model_TZD)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple", drug = "TZD", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple", drug = "TZD", interaction = NA)))
}

adjusted_simple_model_SU <- lm(posthba1cfinal ~ testing_SNP + prehba1c, data = interim.dataset %>% filter(drugclass == "SU"))
if(length(unique(interim.dataset %>% filter(drugclass == "SU") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_model_SU)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_model_SU)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple", drug = "SU", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple", drug = "SU", interaction = NA)))
}

adjusted_simple_model_DPP4 <- lm(posthba1cfinal ~ testing_SNP + prehba1c, data = interim.dataset %>% filter(drugclass == "DPP4"))
if(length(unique(interim.dataset %>% filter(drugclass == "DPP4") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_model_DPP4)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_model_DPP4)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple", drug = "DPP4", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple", drug = "DPP4", interaction = NA)))
}

## adjusted simple interation model
adjusted_simple_interaction_model_SGLT2 <- lm(posthba1cfinal ~ testing_SNP * prehba1c, data = interim.dataset %>% filter(drugclass == "SGLT2"))
if(length(unique(interim.dataset %>% filter(drugclass == "SGLT2") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_interaction_model_SGLT2)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_SGLT2)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple_interaction", drug = "SGLT2", interaction = "testing_SNP")),
                        as.data.frame(summary(adjusted_simple_interaction_model_SGLT2)$coefficients[4,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_SGLT2)[4,])) %>% as.data.frame() %>%
                                        mutate(Model = "adjusted_simple_interaction", drug = "SGLT2", interaction = "testing_SNP:prehba1c")))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "SGLT2", interaction = "testing_SNP")),
                        as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "SGLT2", interaction = "testing_SNP:prehba1c")))
}

adjusted_simple_interaction_model_GLP1 <- lm(posthba1cfinal ~ testing_SNP * prehba1c, data = interim.dataset %>% filter(drugclass == "GLP1"))
if(length(unique(interim.dataset %>% filter(drugclass == "GLP1") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_interaction_model_GLP1)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_GLP1)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple_interaction", drug = "GLP1", interaction = "testing_SNP")),
                        as.data.frame(summary(adjusted_simple_interaction_model_GLP1)$coefficients[4,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_GLP1)[4,])) %>% as.data.frame() %>%
                                        mutate(Model = "adjusted_simple_interaction", drug = "GLP1", interaction = "testing_SNP:prehba1c")))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "GLP1", interaction = "testing_SNP")),
                        as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "GLP1", interaction = "testing_SNP:prehba1c")))
}

adjusted_simple_interaction_model_TZD <- lm(posthba1cfinal ~ testing_SNP * prehba1c, data = interim.dataset %>% filter(drugclass == "TZD"))
if(length(unique(interim.dataset %>% filter(drugclass == "TZD") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_interaction_model_TZD)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_TZD)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple_interaction", drug = "TZD", interaction = "testing_SNP")),
                        as.data.frame(summary(adjusted_simple_interaction_model_TZD)$coefficients[4,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_TZD)[4,])) %>% as.data.frame() %>%
                                        mutate(Model = "adjusted_simple_interaction", drug = "TZD", interaction = "testing_SNP:prehba1c")))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "TZD", interaction = "testing_SNP")),
                        as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "TZD", interaction = "testing_SNP:prehba1c")))
}

adjusted_simple_interaction_model_SU <- lm(posthba1cfinal ~ testing_SNP * prehba1c, data = interim.dataset %>% filter(drugclass == "SU"))
if(length(unique(interim.dataset %>% filter(drugclass == "SU") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_interaction_model_SU)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_SU)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple_interaction", drug = "SU", interaction = "testing_SNP")),
                        as.data.frame(summary(adjusted_simple_interaction_model_SU)$coefficients[4,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_SU)[4,])) %>% as.data.frame() %>%
                                        mutate(Model = "adjusted_simple_interaction", drug = "SU", interaction = "testing_SNP:prehba1c")))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "SU", interaction = "testing_SNP")),
                        as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "SU", interaction = "testing_SNP:prehba1c")))
}

adjusted_simple_interaction_model_DPP4 <- lm(posthba1cfinal ~ testing_SNP * prehba1c, data = interim.dataset %>% filter(drugclass == "DPP4"))
if(length(unique(interim.dataset %>% filter(drugclass == "DPP4") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_simple_interaction_model_DPP4)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_DPP4)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_simple_interaction", drug = "DPP4", interaction = "testing_SNP")),
                        as.data.frame(summary(adjusted_simple_interaction_model_DPP4)$coefficients[4,] %>% t() %>% cbind(t(confint(adjusted_simple_interaction_model_DPP4)[4,])) %>% as.data.frame() %>%
                                        mutate(Model = "adjusted_simple_interaction", drug = "DPP4", interaction = "testing_SNP:prehba1c")))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "DPP4", interaction = "testing_SNP")),
                        as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_simple_interaction", drug = "DPP4", interaction = "testing_SNP:prehba1c")))
}

## adjusted complex model
adjusted_complex_model_SGLT2 <- lm(posthba1cfinal ~ testing_SNP + prehba1c + agetx + sex, data = interim.dataset %>% filter(drugclass == "SGLT2"))
if(length(unique(interim.dataset %>% filter(drugclass == "SGLT2") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_complex_model_SGLT2)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_complex_model_SGLT2)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_complex", drug = "SGLT2", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_complex", drug = "SGLT2", interaction = NA)))
}

adjusted_complex_model_GLP1 <- lm(posthba1cfinal ~ testing_SNP + prehba1c + agetx + sex, data = interim.dataset %>% filter(drugclass == "GLP1"))
if(length(unique(interim.dataset %>% filter(drugclass == "GLP1") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_complex_model_GLP1)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_complex_model_GLP1)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_complex", drug = "GLP1", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_complex", drug = "GLP1", interaction = NA)))
}

adjusted_complex_model_TZD <- lm(posthba1cfinal ~ testing_SNP + prehba1c + agetx + sex, data = interim.dataset %>% filter(drugclass == "TZD"))
if(length(unique(interim.dataset %>% filter(drugclass == "TZD") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_complex_model_TZD)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_complex_model_TZD)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_complex", drug = "TZD", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_complex", drug = "TZD", interaction = NA)))
}

adjusted_complex_model_SU <- lm(posthba1cfinal ~ testing_SNP + prehba1c + agetx + sex, data = interim.dataset %>% filter(drugclass == "SU"))
if(length(unique(interim.dataset %>% filter(drugclass == "SU") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_complex_model_SU)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_complex_model_SU)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_complex", drug = "SU", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_complex", drug = "SU", interaction = NA)))
}

adjusted_complex_model_DPP4 <- lm(posthba1cfinal ~ testing_SNP + prehba1c + agetx + sex, data = interim.dataset %>% filter(drugclass == "DPP4"))
if(length(unique(interim.dataset %>% filter(drugclass == "DPP4") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_complex_model_DPP4)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_complex_model_DPP4)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_complex", drug = "DPP4", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_complex", drug = "DPP4", interaction = NA)))
}

## adjusted full model (all variables in the model)
adjusted_full_model_SGLT2 <- lm(posthba1cfinal ~ testing_SNP + sex + t2dmduration + prebmi + prehba1c + agetx + 
                             prealt + preegfr + pretotalcholesterol + prehdl + ethnicity + smoke + imd5 + 
                             hba1cmonth + ncurrtx + drugline, data = interim.dataset %>% filter(drugclass == "SGLT2"))
if(length(unique(interim.dataset %>% filter(drugclass == "SGLT2") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_full_model_SGLT2)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_full_model_SGLT2)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_full", drug = "SGLT2", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_full", drug = "SGLT2", interaction = NA)))
}

adjusted_full_model_GLP1 <- lm(posthba1cfinal ~ testing_SNP + sex + t2dmduration + prebmi + prehba1c + agetx + 
                             prealt + preegfr + pretotalcholesterol + prehdl + ethnicity + smoke + imd5 + 
                             hba1cmonth + ncurrtx + drugline, data = interim.dataset %>% filter(drugclass == "GLP1"))
if(length(unique(interim.dataset %>% filter(drugclass == "GLP1") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_full_model_GLP1)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_full_model_GLP1)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_full", drug = "GLP1", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_full", drug = "GLP1", interaction = NA)))
}

adjusted_full_model_TZD <- lm(posthba1cfinal ~ testing_SNP + sex + t2dmduration + prebmi + prehba1c + agetx + 
                             prealt + preegfr + pretotalcholesterol + prehdl + ethnicity + smoke + imd5 + 
                             hba1cmonth + ncurrtx + drugline, data = interim.dataset %>% filter(drugclass == "TZD"))
if(length(unique(interim.dataset %>% filter(drugclass == "TZD") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_full_model_TZD)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_full_model_TZD)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_full", drug = "TZD", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_full", drug = "TZD", interaction = NA)))
}

adjusted_full_model_SU <- lm(posthba1cfinal ~ testing_SNP + sex + t2dmduration + prebmi + prehba1c + agetx + 
                             prealt + preegfr + pretotalcholesterol + prehdl + ethnicity + smoke + imd5 + 
                             hba1cmonth + ncurrtx + drugline, data = interim.dataset %>% filter(drugclass == "SU"))
if(length(unique(interim.dataset %>% filter(drugclass == "SU") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_full_model_SU)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_full_model_SU)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_full", drug = "SU", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_full", drug = "SU", interaction = NA)))
}

adjusted_full_model_DPP4 <- lm(posthba1cfinal ~ testing_SNP + sex + t2dmduration + prebmi + prehba1c + agetx + 
                             prealt + preegfr + pretotalcholesterol + prehdl + ethnicity + smoke + imd5 + 
                             hba1cmonth + ncurrtx + drugline, data = interim.dataset %>% filter(drugclass == "DPP4"))
if(length(unique(interim.dataset %>% filter(drugclass == "DPP4") %>% select(testing_SNP) %>% unlist())) == 2) {
  output_table <- rbind(output_table, as.data.frame(summary(adjusted_full_model_DPP4)$coefficients[2,] %>% t() %>% cbind(t(confint(adjusted_full_model_DPP4)[2,])) %>% as.data.frame() %>%
                                                      mutate(Model = "adjusted_full", drug = "DPP4", interaction = NA)))
} else {
  output_table <- rbind(output_table, as.data.frame(cbind(`Estimate` = NA, `Std. Error` = NA, `t value` = NA, `Pr(>|t|)` = NA, `2.5 %` = NA, `97.5 %` = NA, Model = "adjusted_full", drug = "DPP4", interaction = NA)))
}



# save the file - change name appropriately
## Original format is 
##  - '01.snp_effect' helps distinguish the results tables
##  - '_SU' helps understand if there is any filtering of patients (replace as needed)
##  - '_CYP2C92' helps understand the SNP being tested (replace as needed)
##  - '.rds' is the file format
saveRDS(output_table, "01.snp_effect_SU_CYP2C92.rds")

######


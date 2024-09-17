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

# Make predictions for each treatment for all patients
interim.dataset <- original.dataset %>%
  cbind(
    pred.SGLT2 = predict(m1.5.final, original.dataset %>% mutate(drugclass = "SGLT2")),
    pred.GLP1 = predict(m1.5.final, original.dataset %>% mutate(drugclass = "GLP1")),
    pred.DPP4 = predict(m1.5.final, original.dataset %>% mutate(drugclass = "DPP4")),
    pred.SU = predict(m1.5.final, original.dataset %>% mutate(drugclass = "SU")),
    pred.TZD = predict(m1.5.final, original.dataset %>% mutate(drugclass = "TZD"))
  )

# Concordant vs discordant calibration

# ??????????????

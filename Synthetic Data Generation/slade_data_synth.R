############################################################################
#
# This file contains the code for generating synth dataset for the 5-drug model
#
############################################################################

# load libraries
library(tidyverse)
library(rms)

# load objects needed
## CPRD dataset
load("/slade/CPRD_data/mastermind_2022/John/data/match.t.Rdata")
## Five drug model
load("/slade/CPRD_data/mastermind_2022/John/data/fivedrugmodel_5knot.Rdata")

# wrangle the dataset
interim.dataset <- match.t %>%
  # select the correct columns
  select(
    sex, t2dmduration, prebmi, prehba1c, agetx, prealt, preegfr, 
    pretotalcholesterol, prehdl, ethnicity, smoke, imd5, hba1cmonth, 
    ncurrtx, drugline, drugclass, yrdrugstart, posthba1cfinal
  ) %>%
  # drop any missingness
  drop_na() %>%
  # sample n patients
  sample_n(10000)

# test whether the dataset can be used for predicting from the model
## this makes sure the dataset is the right one
predict(m1.5, interim.dataset, conf.int = 0.95)


############################################################################
# Create synthetic dataset

# load libraries
library(synthpop)

# generate synth data
synth.dataset <- syn(interim.dataset, seed = 123)$syn

# save full synth dataset
saveRDS(synth.dataset, "full.synth.dataset.rds")

# save shorten synth dataset
saveRDS(synth.dataset %>% slice(1:100), "short.synth.dataset.rds")

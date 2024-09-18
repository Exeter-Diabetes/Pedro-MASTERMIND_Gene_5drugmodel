##############################################################################
#
# This file loads the dataset and creates/saves descriptive tables
#
##############################################################################

######
# load libraries
library(tidyverse)
library(tableone)

######

##
## Change this line of code to load the Scottish dataset
##

# load dataset
original.dataset <- readRDS("../Synthetic Data/full.synth.dataset.rds")


######

# Generate a column for HbA1c change
interim.dataset <- original.dataset %>%
  mutate(hba1cresp = prehba1c - posthba1cfinal)

# continuous variables for the table
vars <- c(
  # General descriptives
  "agetx", "sex", "t2dmduration", "ethnicity", "imd5", "smoke",
  # Laboratory measurements
  "prebmi", "prehba1c", "preegfr", "pretotalcholesterol", "prehdl", "prealt",
  # Anti-hyperglycaemic treatment
  "drugline", "ncurrtx", "yrdrugstart", "hba1cmonth", "posthba1cfinal", "hba1cresp"
  )

# categorical variables for the table
vars_cat <- c(
  # General descriptives
  "sex", "ethnicity", "imd5", "smoke",
  # Anti-hyperglycaemic treatment
  "drugline", "ncurrtx"
)

# Generate table by drug taken
table_characteristics <- tableone::CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("drugclass"),
  data = interim.dataset,
  test = FALSE
)

# turn the table into a printable format
table_characteristics_print <- print(
  table_characteristics,
  exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1
)

# save the table into an excel file. Currently, it will save in the working directory.
write.csv(table_characteristics_print, file = "00.table_characteristics.csv")



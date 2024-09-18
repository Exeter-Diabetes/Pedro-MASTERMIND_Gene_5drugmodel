library(tidyverse)
library(rms)

load("/slade/CPRD_data/mastermind_2022/John/data/match.t.Rdata")
load("/slade/CPRD_data/mastermind_2022/John/data/fivedrugmodel_5knot.Rdata")


interim.dataset <- match.t %>%
  select(sex, t2dmduration, prebmi, prehba1c, agetx, prealt, preegfr, pretotalcholesterol, prehdl, ethnicity, smoke, imd5, hba1cmonth, ncurrtx, drugline, drugclass, yrdrugstart, posthba1cfinal) %>%
  drop_na() %>%
  sample_n(10000)


predict(m1.5, interim.dataset, conf.int = 0.95)



## Create synthetic dataset
library(synthpop)

synth.dataset <- syn(interim.dataset, seed = 123)$syn

saveRDS(synth.dataset, "full.synth.dataset.rds")



saveRDS(synth.dataset %>% slice(1:100), "short.synth.dataset.rds")
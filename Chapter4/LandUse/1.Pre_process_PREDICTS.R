## Pre-process PREDICTS database

library(dplyr)

Predicts <- readRDS("../Data/PredictsVertebrates.rds") %>%
  filter(Best_guess_binomial !="")

Predicts <- predictsFunctions::MergeSites(Predicts)
Predicts <- Predicts %>% droplevels(Predicts$Predominant_land_use)
Predicts <- Predicts %>% droplevels(Predicts$Use_intensity)
Predicts <- Predicts %>% droplevels(Predicts$Class)

saveRDS(Predicts, "../Results/Predicts_merged_sites.rds")


## Adding phylogenetic eigenvectors to trait datasets (using phylogenetic eigenvectors already extracted in the FD analyses)

X <- c("dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ape", "treeio", "ngram","phylobase")
lapply(X, library, character.only=TRUE); rm(X)
source("Functions.R")

## Load data (with taxonomic correction)
AmphibiansPhyEi <- read.csv("../../../1.2.Trait_imputations_FD_using_Global_Gaps_data/Results/Datasets_FD_imputations/Amphibians.csv")
BirdsPhyEi <- read.csv("../../../1.2.Trait_imputations_FD_using_Global_Gaps_data/Results/Datasets_FD_imputations/Birds.csv")
MammalsPhyEi <- read.csv("../../../1.2.Trait_imputations_FD_using_Global_Gaps_data/Results/Datasets_FD_imputations/Mammals.csv")
ReptilesPhyEi <- read.csv("../../../1.2.Trait_imputations_FD_using_Global_Gaps_data/Results/Datasets_FD_imputations/Reptiles.csv")

## Load trait data with diet
Amphibians <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Amphibians.csv")
Birds <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Birds.csv")
Mammals <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Mammals.csv")
Reptiles <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Reptiles.csv")

## Data with phylogenetic eigenvectors

Amphibians <- cbind(Amphibians, AmphibiansPhyEi[, c("EV_1", "EV_2", "EV_3", "EV_4", "EV_5", "EV_6", "EV_7", "EV_8", "EV_9", "EV_10")])
Birds <- cbind(Birds, BirdsPhyEi[, c("EV_1", "EV_2", "EV_3", "EV_4", "EV_5", "EV_6", "EV_7", "EV_8", "EV_9", "EV_10")])
Mammals <- cbind(Mammals, MammalsPhyEi[, c("EV_1", "EV_2", "EV_3", "EV_4", "EV_5", "EV_6", "EV_7", "EV_8", "EV_9", "EV_10")])
Reptiles <- cbind(Reptiles, ReptilesPhyEi[, c("EV_1", "EV_2", "EV_3", "EV_4", "EV_5", "EV_6", "EV_7", "EV_8", "EV_9", "EV_10")])


# quick check
(Amphibians$Best_guess_binomial == AmphibiansPhyEi$Best_guess_binomial) %>%  unique
(Birds$Best_guess_binomial == BirdsPhyEi$Best_guess_binomial) %>%  unique
(Mammals$Best_guess_binomial == MammalsPhyEi$Best_guess_binomial) %>%  unique
(Reptiles$Best_guess_binomial == ReptilesPhyEi$Best_guess_binomial) %>%  unique

write.csv(Amphibians, "../../Results/Traits_with_phy_eigenvectors/Amphibians.csv", row.names = FALSE)
write.csv(Birds, "../../Results/Traits_with_phy_eigenvectors/Birds.csv", row.names = FALSE)
write.csv(Mammals, "../../Results/Traits_with_phy_eigenvectors/Mammals.csv", row.names = FALSE)
write.csv(Reptiles, "../../Results/Traits_with_phy_eigenvectors/Reptiles.csv", row.names = FALSE)

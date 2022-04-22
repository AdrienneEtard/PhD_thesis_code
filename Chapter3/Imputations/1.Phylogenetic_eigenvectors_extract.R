## Adding phylogenetic eigenvectors to trait datasets for imputations of missing trait values

X <- c("dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ape", "treeio", "ngram","phylobase")
lapply(X, library, character.only=TRUE); rm(X)
source("Functions.R")

## Load trait data (with taxonomic correction)

Amphibians <- read.csv("../Data/GlobalGaps_traitdata/Amphibians.csv")
Birds <- read.csv("../Data/GlobalGaps_traitdata/Birds.csv")
Mammals <- read.csv("../Data/GlobalGaps_traitdata/Mammals.csv")
Reptiles <- read.csv("../Data/GlobalGaps_traitdata/Reptiles.csv")

All_species <- c(as.character(Amphibians$Best_guess_binomial),
                 as.character(Birds$Best_guess_binomial),
                 as.character(Mammals$Best_guess_binomial),
                 as.character(Reptiles$Best_guess_binomial))

## Load PREDICTS (with taxonomic corrections)
Predicts <- readRDS("../Data/PredictsVertebrates.rds")
Predicts_species <- unique(Predicts$Best_guess_binomial)

## intersection: verify that all PREDICTS species found in trait datasets.
setdiff(Predicts_species, All_species)

## Load phylogenies, the ones used in Global Gaps article (consensus trees).
Phylo_Mammals <- read.nexus("../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Mammals_complete_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Birds <- read.nexus("../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Birds_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Amphibians <- read.nexus("../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Amphibians_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Reptiles <- read.nexus("../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Reptiles_TreeAnnotator.nex")  %>% .Format_tiplabels()

## species representation in the phylogenies
RepMammals <- intersect(as.character(Mammals$Best_guess_binomial), Phylo_Mammals$tip.label)
length(RepMammals)/nrow(Mammals) * 100

RepBids <- intersect(as.character(Birds$Best_guess_binomial), Phylo_Birds$tip.label)
length(RepBids)/nrow(Birds) * 100

RepAmphibians <- intersect(as.character(Amphibians$Best_guess_binomial), Phylo_Amphibians$tip.label)
length(RepAmphibians)/nrow(Amphibians) * 100

RepReptiles <- intersect(as.character(Reptiles$Best_guess_binomial), Phylo_Reptiles$tip.label)
length(RepReptiles)/nrow(Reptiles) * 100

## Phylogenetic eigenvectors extraction, when possible (when the species is present in the phylogeny).
## I add the first 10 phylogenetic eigenvectors to the trait datasets, following Penone et al. 2014
## (enough to maximise the accuracy of further imputations.)

EV_Amphibians <- Extract_eigenvectors(Amphibians, Phylo_Amphibians, 10)
EV_Reptiles <- Extract_eigenvectors(Reptiles, Phylo_Reptiles, 10)
EV_Birds <- Extract_eigenvectors(Birds, Phylo_Birds, 10)
EV_Mammals <- Extract_eigenvectors(Mammals, Phylo_Mammals, 10)

## Add eigenvectors to datasets and select columns for imputations. Here, we impute without dietary trait (primary diet and diet breadth).

# Select traits for imputations - for FD analysis
Habitats <- colnames(Amphibians)[c(22:34)]

# exclude svl_length.
Amphibians %>%  filter(!is.na(Svl_length_mm)) %>% filter(is.na(Body_length_mm)) %>%  nrow()

Amphibians <- Amphibians %>%
  dplyr::select(Order, Family, Genus, Best_guess_binomial,
                Body_mass_g, Body_length_mm,
                Maturity_d, Max_longevity_d,
                Litter_size,
                Trophic_level, Diel_activity,
                Habitat_breadth_IUCN,
                Specialisation,
                all_of(Habitats))

Birds <- Birds %>%
  dplyr::select(Order, Family, Genus, Best_guess_binomial,
                                  Body_mass_g,
                                  Generation_length_d,
                                  Litter_size,
                                  Trophic_level, Diel_activity,
                                  Habitat_breadth_IUCN,
                                  Specialisation,
                                  all_of(Habitats))

Mammals <- Mammals %>%
  dplyr::select(Order, Family, Genus, Best_guess_binomial,
                Body_mass_g,
                Generation_length_d,
                Litter_size,
                Trophic_level, Diel_activity,
                Habitat_breadth_IUCN,
                Specialisation,
                all_of(Habitats))


# exclude adult svl for reptiles
Reptiles %>%  filter(!is.na(Adult_svl_cm)) %>%  filter(is.na(Body_mass_g)) %>%  nrow()

Reptiles <- Reptiles %>%
  dplyr::select(Order, Family, Genus, Best_guess_binomial,
                Body_mass_g,
                Maturity_d, Max_longevity_d, Longevity_d,
                Litter_size,
                Trophic_level, Diel_activity,
                Habitat_breadth_IUCN,
                Specialisation,
                all_of(Habitats))

# Adding eigenvectors to corrected trait datasets
Amphibians <- merge(Amphibians, EV_Amphibians, by="Best_guess_binomial", all = TRUE)
Birds <- merge(Birds, EV_Birds, by="Best_guess_binomial", all = TRUE)
Reptiles <- merge(Reptiles, EV_Reptiles, by="Best_guess_binomial", all = TRUE)
Mammals <- merge(Mammals, EV_Mammals, by="Best_guess_binomial", all = TRUE)

nrow(Amphibians[!is.na(Amphibians$Maturity_d),])
nrow(Amphibians[!is.na(Amphibians$Max_longevity_d),])

nrow(Reptiles[!is.na(Reptiles$Maturity_d),])
nrow(Reptiles[!is.na(Reptiles$Max_longevity_d),])
nrow(Reptiles[!is.na(Reptiles$Longevity_d),])

write.csv(Amphibians, "../Results/Datasets_FD_imputations/Amphibians.csv", row.names = FALSE)
write.csv(Birds, "../Results/Datasets_FD_imputations/Birds.csv", row.names = FALSE)
write.csv(Reptiles, "../Results/Datasets_FD_imputations/Reptiles.csv", row.names = FALSE)
write.csv(Mammals, "../Results/Datasets_FD_imputations/Mammals.csv", row.names = FALSE)

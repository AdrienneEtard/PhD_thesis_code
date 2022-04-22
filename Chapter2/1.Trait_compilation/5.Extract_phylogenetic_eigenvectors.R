## Adding phylogenetic eigenvectors to trait datasets, using both corrected and uncorrected phylogenies -
# phylogenies that have been treated for pseudoreplication in (4.)

# # # # P r e a m b l e 
X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ape", "treeio", "ngram","phylobase")
lapply(X, library, character.only=TRUE); rm(X)
source("Functions_for_phylogenies.R")

## Load data

# No taxonomic correction
UN_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Amphibians.csv")
UN_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Birds.csv")
UN_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Mammals.csv")
UN_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Reptiles.csv")

# With taxonomic correction
C_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Amphibians.csv")
C_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Birds.csv")
C_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Mammals.csv")
C_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Reptiles.csv")

## Load PREDICTS with and without taxonomic correction
UN_Predicts <- readRDS("../../Data/PREDICTS_database.rds") %>%
  filter(Class %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"))
C_Predicts <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")


## Load phylogenies

# Corrected
C_Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk")  %>% .Format_tiplabels()
C_Phylo_Birds <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk")  %>% .Format_tiplabels()
C_Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk")  %>% .Format_tiplabels()
C_Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk")  %>% .Format_tiplabels()

# Uncorrected and without having attached any species its genus
UN_Phylo_Mammals <- read.newick("../../Data/Mammals/Phylogenies/TTOL_mammals_smoothed_interpolated.nwk")  %>% .Format_tiplabels()
UN_Phylo_Birds <-  read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk")  %>% .Format_tiplabels()
UN_Phylo_Amphibians <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk")  %>% .Format_tiplabels()
UN_Phylo_Reptiles <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk")  %>% .Format_tiplabels()


## Phylogenetic eigenvectors extraction, when possible (when the species is present in the phylogeny).
## I add the first 10 eigenvectors (enough to maximise the accuracy of further imputations.)

# For corrected datasets
C.EV_Amphibians <- Extract_eigenvectors(C_Amphibians, C_Phylo_Amphibians, 10)
C.EV_Reptiles <- Extract_eigenvectors(C_Reptiles, C_Phylo_Reptiles, 10)
C.EV_Birds <- Extract_eigenvectors(C_Birds, C_Phylo_Birds, 10)
C.EV_Mammals <- Extract_eigenvectors(C_Mammals, C_Phylo_Mammals, 10)

# For uncorrected datasets
UN.EV_Amphibians <- Extract_eigenvectors(UN_Amphibians, UN_Phylo_Amphibians, 10)
UN.EV_Mammals <- Extract_eigenvectors(UN_Mammals, UN_Phylo_Mammals, 10)
UN.EV_Reptiles <- Extract_eigenvectors(UN_Reptiles, UN_Phylo_Reptiles, 10)
UN.EV_Birds <- Extract_eigenvectors(UN_Birds, UN_Phylo_Birds, 10)

# For corrected datasets but uncorrected phylogenies
C.EV_Amphibians_unphylo <- Extract_eigenvectors(C_Amphibians, UN_Phylo_Amphibians, 10)
C.EV_Mammals_unphylo <- Extract_eigenvectors(C_Mammals, UN_Phylo_Mammals, 10)
C.EV_Reptiles_unphylo <- Extract_eigenvectors(C_Reptiles, UN_Phylo_Reptiles, 10)
C.EV_Birds_unphylo <- Extract_eigenvectors(C_Birds, UN_Phylo_Birds, 10)


## Save eigenvectors
write.csv(C.EV_Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Amphibians.csv", row.names=FALSE)
write.csv(C.EV_Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Reptiles.csv", row.names=FALSE)
write.csv(C.EV_Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Birds.csv", row.names=FALSE)
write.csv(C.EV_Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Mammals.csv", row.names=FALSE)

write.csv(UN.EV_Amphibians, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Amphibians.csv", row.names=FALSE)
write.csv(UN.EV_Mammals, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Mammals.csv", row.names=FALSE)
write.csv(UN.EV_Reptiles, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Reptiles.csv", row.names=FALSE)
write.csv(UN.EV_Birds, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Birds.csv", row.names=FALSE)

write.csv(C.EV_Amphibians_unphylo, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/AmphibiansEV.csv", row.names=FALSE)
write.csv(C.EV_Reptiles_unphylo, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/ReptilesEV.csv", row.names=FALSE)
write.csv(C.EV_Birds_unphylo, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/BirdsEV.csv", row.names=FALSE)
write.csv(C.EV_Mammals_unphylo, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/MammalsEV.csv", row.names=FALSE)


## Load eigenvectors

# Corrected
C.EV_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Amphibians.csv")
C.EV_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Birds.csv")
C.EV_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Mammals.csv")
C.EV_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Reptiles.csv")

# Uncorrected
UN.EV_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Amphibians.csv")
UN.EV_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Birds.csv")
UN.EV_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Mammals.csv")
UN.EV_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Eigenvectors/Reptiles.csv")


# hybrid
C.EV_Amphibians_unphylo <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/AmphibiansEV.csv")
C.EV_Reptiles_unphylo <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/ReptilesEV.csv")
C.EV_Birds_unphylo <- read.csv ("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/BirdsEV.csv")
C.EV_Mammals_unphylo <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/MammalsEV.csv")


## Add eigenvectors to datasets and order columns

# Adding eigenvectors to corrected trait datasets
C_Amphibians <- Add_eigenvectors(C_Amphibians, C.EV_Amphibians)
C_Birds <- Add_eigenvectors(C_Birds, C.EV_Birds)
C_Reptiles <- Add_eigenvectors(C_Reptiles, C.EV_Reptiles)
C_Mammals <- Add_eigenvectors(C_Mammals, C.EV_Mammals)

# Adding eigenvectors to uncorrected trait datasets
UN_Amphibians <- Add_eigenvectors(UN_Amphibians, UN.EV_Amphibians)
UN_Birds <- Add_eigenvectors(UN_Birds, UN.EV_Birds)
UN_Reptiles <- Add_eigenvectors(UN_Reptiles, UN.EV_Reptiles)
UN_Mammals <- Add_eigenvectors(UN_Mammals, UN.EV_Mammals)

# Adding uncorrected eigenvectors to corrected dataset for further imputation comparisons
Cdf_UNphy_Amphibians <-  Add_eigenvectors(C_Amphibians, C.EV_Amphibians_unphylo)
Cdf_UNphy_Birds <-  Add_eigenvectors(C_Birds, C.EV_Birds_unphylo)
Cdf_UNphy_Mammals <-  Add_eigenvectors(C_Mammals, C.EV_Mammals_unphylo)
Cdf_UNphy_Reptiles <-  Add_eigenvectors(C_Reptiles, C.EV_Reptiles_unphylo)


## Rearranging column order and cleaning datasets

# Merging together maturity and generation length, and longevity and max longevity

# Reptiles 
C_Reptiles$Longevity_d <- apply(C_Reptiles[, c("Longevity_d","Max_longevity_d")], 1, median, na.rm=TRUE)
C_Reptiles <- C_Reptiles %>% select(-Max_longevity_d)

UN_Reptiles$Longevity_d <- apply(UN_Reptiles[, c("Longevity_d","Max_longevity_d")], 1, median, na.rm=TRUE)
UN_Reptiles <- UN_Reptiles %>% select(-Max_longevity_d)
Cdf_UNphy_Reptiles$Longevity_d<- apply(Cdf_UNphy_Reptiles[, c("Longevity_d","Max_longevity_d")], 1, median, na.rm=TRUE)
Cdf_UNphy_Reptiles <- Cdf_UNphy_Reptiles %>% select(-Max_longevity_d)

C_Reptiles$Primary_diet <- NA; C_Reptiles$Diet_breadth <- NA
UN_Reptiles$Primary_diet <- NA; UN_Reptiles$Diet_breadth <- NA
Cdf_UNphy_Reptiles$Diet_breadth <- NA; Cdf_UNphy_Reptiles$Primary_diet <- NA

C_Reptiles$Class <- "Reptilia"
UN_Reptiles$Class <- "Reptilia"
Cdf_UNphy_Reptiles$Class <- "Reptilia"

C_Reptiles$Maturity_d[C_Reptiles$Maturity_d==0] <- NA
UN_Reptiles$Maturity_d[UN_Reptiles$Maturity_d==0] <- NA
Cdf_UNphy_Reptiles$Maturity_d[Cdf_UNphy_Reptiles$Maturity_d==0] <- NA


# Amphibians
colnames(C_Amphibians)[9] <- "Longevity_d"
colnames(UN_Amphibians)[9] <- "Longevity_d"
colnames(Cdf_UNphy_Amphibians)[9] <- "Longevity_d"

C_Amphibians$Class <- "Amphibia"
UN_Amphibians$Class <- "Amphibia"
Cdf_UNphy_Amphibians$Class <- "Amphibia"


# Birds
C_Birds$Longevity_d <- apply(C_Birds[, c("Max_longevity_d","Longevity_d")], 1, median, na.rm=TRUE)
C_Birds <- C_Birds %>% select(-Max_longevity_d)
UN_Birds$Longevity_d <- apply(UN_Birds[, c("Max_longevity_d","Longevity_d")], 1, median, na.rm=TRUE)
UN_Birds <- UN_Birds %>% select(-Max_longevity_d)
Cdf_UNphy_Birds$Longevity_d <- apply(Cdf_UNphy_Birds[, c("Max_longevity_d","Longevity_d")], 1, median, na.rm=TRUE)
Cdf_UNphy_Birds <- Cdf_UNphy_Birds %>% select(-Max_longevity_d)

C_Birds$Class <- "Aves"
UN_Birds$Class <- "Aves"
Cdf_UNphy_Birds$Class <- "Aves"

# Mammals
C_Mammals$Longevity_d <- apply(C_Mammals[, c("Max_longevity_d","Longevity_d")], 1, median, na.rm=TRUE)
C_Mammals <- C_Mammals %>% select(-Max_longevity_d)
UN_Mammals$Longevity_d <- apply(UN_Mammals[, c("Max_longevity_d","Longevity_d")], 1, median, na.rm=TRUE)
UN_Mammals <- UN_Mammals %>% select(-Max_longevity_d)
Cdf_UNphy_Mammals$Longevity_d <-  apply(Cdf_UNphy_Mammals[, c("Max_longevity_d","Longevity_d")], 1, median, na.rm=TRUE)
Cdf_UNphy_Mammals <- Cdf_UNphy_Mammals %>% select(-Max_longevity_d)

C_Mammals$Class <- "Mammalia"
UN_Mammals$Class <- "Mammalia"
Cdf_UNphy_Mammals$Class <- "Mammalia"

## Rearrange columns
Ev <- vector()
for (i in 1:10) {Ev <- c(Ev, paste("EV", i, sep="_")) }
Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")
Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas",
             "Caves.and.subterranean","Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")
Taxinfo <- c("Class", "Order", "Family", "Genus", "Best_guess_binomial")


# Corrected datasets
C_Reptiles <- C_Reptiles[, c(Taxinfo,
                         "Body_mass_g", "Adult_svl_cm", "Maturity_d", "Longevity_d",
                         "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", "Specialisation",
                         "Habitat_breadth_IUCN", Habitat, Ev)]

C_Amphibians <- C_Amphibians[, c(Taxinfo,
                             "Body_mass_g", "Body_length_mm","Svl_length_mm", "Maturity_d", "Longevity_d",
                             "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet",  Diet,"Specialisation",
                             "Habitat_breadth_IUCN", Habitat, Ev)]

C_Birds <- C_Birds[, c(Taxinfo,
                   "Body_mass_g", "Adult_svl_cm", "Generation_length_d","Maturity_d", "Longevity_d",
                   "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth","Primary_diet",Diet,"Specialisation",
                   "Habitat_breadth_IUCN", Habitat, Ev)]


C_Mammals <- C_Mammals[, c(Taxinfo,
                       "Body_mass_g", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm",
                       "Generation_length_d","Maturity_d", "Longevity_d", "AFR_d",
                       "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", Diet, "Specialisation",
                       "Habitat_breadth_IUCN", Habitat, Ev)]

# Uncorrected datasets
UN_Reptiles <- UN_Reptiles[, c("Class", "Best_guess_binomial",
                             "Body_mass_g", "Adult_svl_cm", "Maturity_d", "Longevity_d",
                             "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", "Specialisation",
                             "Habitat_breadth_IUCN", Habitat, Ev)]

UN_Amphibians <- UN_Amphibians[, c("Class", "Best_guess_binomial",
                                 "Body_mass_g", "Body_length_mm","Svl_length_mm", "Maturity_d", "Longevity_d",
                                 "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet",  Diet,"Specialisation",
                                 "Habitat_breadth_IUCN", Habitat, Ev)]

UN_Birds <- UN_Birds[, c("Class","Best_guess_binomial",
                       "Body_mass_g", "Adult_svl_cm", "Generation_length_d","Maturity_d", "Longevity_d",
                       "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth","Primary_diet",Diet,"Specialisation",
                       "Habitat_breadth_IUCN", Habitat, Ev)]


UN_Mammals <- UN_Mammals[, c("Class","Best_guess_binomial",
                           "Body_mass_g", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm",
                           "Generation_length_d","Maturity_d", "Longevity_d", "AFR_d",
                           "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", Diet, "Specialisation",
                           "Habitat_breadth_IUCN", Habitat, Ev)]


# Hybrid datasets
Cdf_UNphy_Reptiles <- Cdf_UNphy_Reptiles[, c(Taxinfo,
                             "Body_mass_g", "Adult_svl_cm", "Maturity_d", "Longevity_d",
                             "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", "Specialisation",
                             "Habitat_breadth_IUCN", Habitat, Ev)]

Cdf_UNphy_Amphibians <- Cdf_UNphy_Amphibians[, c(Taxinfo,
                                 "Body_mass_g", "Body_length_mm","Svl_length_mm", "Maturity_d", "Longevity_d",
                                 "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet",  Diet,"Specialisation",
                                 "Habitat_breadth_IUCN", Habitat, Ev)]


Cdf_UNphy_Birds <- Cdf_UNphy_Birds[, c(Taxinfo,
                       "Body_mass_g", "Adult_svl_cm", "Generation_length_d","Maturity_d", "Longevity_d",
                       "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth","Primary_diet",Diet,"Specialisation",
                       "Habitat_breadth_IUCN", Habitat, Ev)]


Cdf_UNphy_Mammals <- Cdf_UNphy_Mammals[, c(Taxinfo,
                           "Body_mass_g", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm",
                           "Generation_length_d","Maturity_d", "Longevity_d", "AFR_d",
                           "Litter_size", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", Diet, "Specialisation",
                           "Habitat_breadth_IUCN", Habitat, Ev)]

## Save files
write.csv(C_Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv", row.names = FALSE)
write.csv(C_Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv", row.names = FALSE)
write.csv(C_Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv", row.names = FALSE)
write.csv(C_Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv", row.names = FALSE)

write.csv(UN_Amphibians, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv", row.names = FALSE)
write.csv(UN_Birds, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv", row.names = FALSE)
write.csv(UN_Mammals, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv", row.names = FALSE)
write.csv(UN_Reptiles, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv", row.names = FALSE)



# hybrid files
write.csv(Cdf_UNphy_Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Cor_datasets_uncor_phy/Amphibians.csv", row.names = FALSE)
write.csv(Cdf_UNphy_Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Cor_datasets_uncor_phy/Birds.csv", row.names = FALSE)
write.csv(Cdf_UNphy_Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Cor_datasets_uncor_phy/Mammals.csv", row.names = FALSE)
write.csv(Cdf_UNphy_Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Cor_datasets_uncor_phy/Reptiles.csv", row.names = FALSE)








## Compute distance matrices from selected trait dataset

source("Functions.R")
library(dplyr)
library(FD)
library(usedist)
library(beepr)

## Load PREDICTS data

# Predicts <- readRDS("../Data/PredictsVertebrates.rds") %>%
#   filter(Best_guess_binomial !="")
#
# Predicts <- CorrectSamplingEffort(Predicts)
# Predicts <- predictsFunctions::MergeSites(Predicts)
# length(unique(Predicts$SS))
#
# Predicts <- Predicts %>% droplevels(Predicts$Predominant_land_use)
# Predicts <- Predicts %>% droplevels(Predicts$Use_intensity)
# Predicts <- Predicts %>% droplevels(Predicts$Class)
#
# saveRDS(Predicts, "../Results/Predicts_merged_sites.rds")

Predicts <- readRDS("../../Results/Predicts_merged_sites.rds")
X <- unique(Predicts$Best_guess_binomial)

## Get Gower distances for each set of imputed dataset (subsetted for PREDICTS species)
Imputed8 <- readRDS("../../Results/Imputed_zscored_standardised/All_imputed_sets.rds")

gc()
memory.limit(size=50000)

Selected_traits <- c("log10_Body_mass_g",
                     "log10_Lifespan_proxy",
                     "log10_Litter_size",
                     "sqrt_Habitat_breadth_IUCN",
                     "Specialisation",
                     "Diel_activity",
                     "Trophic_level")

Gower_dist_Predicts_8 <- list()

for(i in 1:8){
  Gower_dist_i <- GowerDist(Imputed8[[i]], c("Best_guess_binomial", Selected_traits))
  Gower_dist_Predicts_8[[i]] <- dist_subset(Gower_dist_i, X)
  rm(Gower_dist_i)
  print(i)
}

## Save Gower distance matrices calculated across all ter ver
saveRDS(Gower_dist_Predicts_8, "../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets.rds")


## Obtaining distance matrices by class for the 8 imputed trait datasets
Class_imputed <- readRDS("../../Results/Imputed_zscored_standardised/All_imputed_sets_byClass.rds")

# distance matrices by class
Classes <- lapply(Class_imputed, function(x){
  y <- split(x, f=as.factor(x$Class))
})

Mammals <- lapply(Classes, function(x){y <- x[["Mammals"]]; return(y)})
Birds <- lapply(Classes, function(x){y <- x[["Birds"]]; return(y)})
Reptiles <- lapply(Classes, function(x){y <- x[["Reptiles"]]; return(y)})
Amphibians <- lapply(Classes, function(x){y <- x[["Amphibians"]]; return(y)})

XMammals <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Mammalia"])
XBirds <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Aves"])
XAmphibians <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Amphibia"])
XReptiles <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Reptilia"])

GetGowerList <- function(Data, species){
  ListGower <- list()
  for(i in 1:8){
   Gower_dist_i <- GowerDist(Data[[i]], c("Best_guess_binomial", Selected_traits))
   ListGower[[i]] <- dist_subset(Gower_dist_i, species)
   rm(Gower_dist_i)
   print(i)
  }
  return(ListGower)
}

GowerMammals8 <- GetGowerList(Mammals, XMammals)
GowerBirds8 <- GetGowerList(Birds, XBirds)
GowerAmphibians8 <- GetGowerList(Amphibians, XAmphibians)
GowerReptiles8 <- GetGowerList(Reptiles, XReptiles)

## save class-specific Gower distance matrices
saveRDS(GowerMammals8, "../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Mammals.rds")
saveRDS(GowerBirds8, "../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Birds.rds")
saveRDS(GowerAmphibians8, "../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Amphibians.rds")
saveRDS(GowerReptiles8, "../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Reptiles.rds")

## Gower distance matrix across all terrestrial vertebrates for the 8th imputed trait dataset - to use for ecoregional simulations
Gower_dist_terver <- GowerDist(Imputed8[[8]], c("Best_guess_binomial", Selected_traits))
saveRDS(Gower_dist_terver, "../../Results/Gower_distances/Gower_distances_terver8.rds")

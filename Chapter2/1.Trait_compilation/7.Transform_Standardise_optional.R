## Transform and standardise compiled datasets before imputations.
# standardise traits (log10 and sqrt transformations, and z-scoring)

# Only for corrected datasets as they will be the ones imputed on. 
# The imputations will be run on both standardised and non standardised traits in two different sets, to compare.

# # PREAMBLE

X <-c("data.table", "plyr", "dplyr", "tidyr", "magrittr", "reshape", "reshape2", "stringr", "stringi", "lazyeval", "rlang", "PerformanceAnalytics") 
lapply(X, library, character.only=TRUE); rm(X)

Transform_zscore <- function(TraitDF, Trait, Transf) {
  
  if (Transf=="log10"){
    TraitDF[,paste("log10", Trait, sep="_")] <- as.numeric(log10(TraitDF[,Trait]))
    TraitDF[, paste("log10", Trait, sep="_")] <- scale(TraitDF[, paste("log10", Trait, sep="_")], center=TRUE, scale=TRUE)
    TraitDF[, paste("log10", Trait, sep="_")] <- as.numeric(TraitDF[, paste("log10", Trait, sep="_")])
  }
  
  if(Transf=="sqrt") {
    TraitDF[,Trait]  <- as.numeric(sqrt(TraitDF[,Trait]))
    TraitDF[, paste("sqrt", Trait, sep="_")] <- scale(TraitDF[,Trait] , center=TRUE, scale=TRUE)
    TraitDF[, paste("sqrt", Trait, sep="_")] <- as.numeric(TraitDF[, paste("sqrt", Trait, sep="_")])
    
  }
  
  return(TraitDF)
}

## Load trait data 
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

## Predicts
PredictsVertebrates <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

## Standardising traits

# Mammals
Mammals <- Transform_zscore(Mammals, "Body_mass_g", "log10")
Mammals <- Transform_zscore(Mammals, "Adult_svl_cm", "log10")
Mammals <- Transform_zscore(Mammals, "Forearm_length_mm", "log10")
Mammals <- Transform_zscore(Mammals, "Head_length_mm", "log10")
Mammals <- Transform_zscore(Mammals, "Generation_length_d", "log10")
Mammals <- Transform_zscore(Mammals, "Longevity_d", "log10")
Mammals <- Transform_zscore(Mammals, "Maturity_d", "log10")
Mammals <- Transform_zscore(Mammals, "AFR_d", "log10")
Mammals <- Transform_zscore(Mammals, "Litter_size", "log10")
Mammals <- Transform_zscore(Mammals, "Habitat_breadth_IUCN", "sqrt")
Mammals <- Transform_zscore(Mammals, "Diet_breadth", "sqrt")

# Reptiles
Reptiles <- Transform_zscore(Reptiles, "Body_mass_g", "log10")
Reptiles <- Transform_zscore(Reptiles, "Adult_svl_cm", "log10")
Reptiles <- Transform_zscore(Reptiles, "Maturity_d", "log10")
Reptiles <- Transform_zscore(Reptiles, "Longevity_d", "log10")
Reptiles <- Transform_zscore(Reptiles, "Litter_size", "log10")
Reptiles <- Transform_zscore(Reptiles, "Habitat_breadth_IUCN", "sqrt")

# Birds
Birds <- Transform_zscore(Birds, "Body_mass_g", "log10")
Birds <- Transform_zscore(Birds, "Adult_svl_cm", "log10")
Birds <- Transform_zscore(Birds, "Maturity_d", "log10")
Birds <- Transform_zscore(Birds, "Longevity_d", "log10")
Birds <- Transform_zscore(Birds, "Generation_length_d", "log10")
Birds <- Transform_zscore(Birds, "Litter_size", "log10")
Birds <- Transform_zscore(Birds, "Habitat_breadth_IUCN", "sqrt")
Birds <- Transform_zscore(Birds, "Diet_breadth", "sqrt")

# Amphibians
Amphibians <- Transform_zscore(Amphibians, "Body_mass_g", "log10")
Amphibians <- Transform_zscore(Amphibians, "Body_length_mm", "log10")
Amphibians <- Transform_zscore(Amphibians, "Svl_length_mm", "log10")
Amphibians <- Transform_zscore(Amphibians, "Maturity_d", "log10")
Amphibians <- Transform_zscore(Amphibians, "Longevity_d", "log10")
Amphibians <- Transform_zscore(Amphibians, "Litter_size", "log10")
Amphibians <- Transform_zscore(Amphibians, "Habitat_breadth_IUCN", "sqrt")
Amphibians <- Transform_zscore(Amphibians, "Diet_breadth", "sqrt")


## Reorganising columns
Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")
Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas",
             "Caves.and.subterranean","Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")
Ev <- vector()
for (i in 1:10) {Ev <- c(Ev, paste("EV", i, sep="_")) }

Reptiles$sqrt_Diet_breadth <- NA
Reptiles <- Reptiles[, c("Class", "Order", "Family", "Genus", "Best_guess_binomial",
                         "log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Maturity_d", "log10_Longevity_d",
                         "log10_Litter_size", "Diel_activity", "Trophic_level", "sqrt_Diet_breadth", "Primary_diet", "Specialisation",
                         "sqrt_Habitat_breadth_IUCN", Habitat, Ev)]

Amphibians <- Amphibians[, c("Class", "Order", "Family", "Genus", "Best_guess_binomial",
                         "log10_Body_mass_g", "log10_Body_length_mm","log10_Svl_length_mm", "log10_Maturity_d", "log10_Longevity_d",
                         "log10_Litter_size", "Diel_activity", "Trophic_level", "sqrt_Diet_breadth", "Primary_diet",  Diet,"Specialisation",
                         "sqrt_Habitat_breadth_IUCN", Habitat, Ev)]

Birds <- Birds[, c("Class", "Order", "Family", "Genus", "Best_guess_binomial",
                             "log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Generation_length_d","log10_Maturity_d", "log10_Longevity_d",
                             "log10_Litter_size", "Diel_activity", "Trophic_level", "sqrt_Diet_breadth","Primary_diet",Diet,"Specialisation",
                             "sqrt_Habitat_breadth_IUCN", Habitat, Ev)]


Mammals <- Mammals[, c("Class", "Order", "Family", "Genus", "Best_guess_binomial",
                   "log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Forearm_length_mm","log10_Head_length_mm",
                   "log10_Generation_length_d","log10_Maturity_d", "log10_Longevity_d", "log10_AFR_d",
                   "log10_Litter_size", "Diel_activity", "Trophic_level", "sqrt_Diet_breadth", "Primary_diet", Diet, "Specialisation",
                   "sqrt_Habitat_breadth_IUCN", Habitat, Ev)]


# Saving results
write.csv(Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Reptiles.csv", row.names=FALSE)
write.csv(Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Mammals.csv", row.names=FALSE)
write.csv(Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Birds.csv", row.names=FALSE)
write.csv(Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Amphibians.csv", row.names=FALSE)


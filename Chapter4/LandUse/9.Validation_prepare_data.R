## Load trait data with phylogenetic imformation as eigenvectors
Mammals <- read.csv("../../Results/Traits_with_phy_eigenvectors/Mammals.csv") %>% 
  dplyr::select(-Trophic_level.Elton)
Birds <- read.csv("../../Results/Traits_with_phy_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/Traits_with_phy_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/Traits_with_phy_eigenvectors/Reptiles.csv")

## set lifespan proxy
colnames(Mammals)[10] <- "Lifespan_proxy"
colnames(Birds)[11] <- "Lifespan_proxy"

nrow(Amphibians[is.na(Amphibians$Max_longevity_d),])
nrow(Amphibians[is.na(Amphibians$Maturity_d),])
colnames(Amphibians)[10] <- "Lifespan_proxy"

nrow(Reptiles[is.na(Reptiles$Max_longevity_d),])
nrow(Reptiles[is.na(Reptiles$Longevity_d),])
nrow(Reptiles[is.na(Reptiles$Maturity_d),])
colnames(Reptiles)[8] <- "Lifespan_proxy"

## transform all traits
Transform <- function(TraitDF, Trait, Transf) {
  
  if (Transf=="log10"){
    TraitDF[,paste("log10", Trait, sep="_")] <- as.numeric(log10(TraitDF[,Trait]))
  }
  
  if(Transf=="sqrt") {
    TraitDF[,Trait] <- as.numeric(TraitDF[,Trait])
    TraitDF[, paste("sqrt", Trait, sep="_")]  <- as.numeric(sqrt(TraitDF[,Trait]))
  }
  
  return(TraitDF)
}

Candidate_traits <-  c("Body_mass_g",
                       "Lifespan_proxy",
                       "Litter_size",
                       "Habitat_breadth_IUCN",
                       "Diet_breadth")

transf <- c("log10", "log10", "log10", "sqrt", "sqrt")

for(t in 1:length(Candidate_traits)){
  Birds <- Transform(Birds, Trait=Candidate_traits[t], Transf = transf[t])
  Mammals <- Transform(Mammals, Trait=Candidate_traits[t], Transf = transf[t])
  Amphibians <- Transform(Amphibians, Trait=Candidate_traits[t], Transf = transf[t])
  Reptiles <- Transform(Reptiles, Trait=Candidate_traits[t], Transf = transf[t])
}

## combine traits and range sizes
Index <- read.csv("E:/PhD/PhD_R_projects/Range_maps_work/Results/7_GenerateSpeciesIndices/SpeciesIndices.csv")
RS <- read.csv("../../Data/Range_sizes.csv")
RS$Species <- paste0("sp", RS$Index)
RS <- RS %>% 
  filter(Range_area_sq_km!=0)
RS <- left_join(RS, Index)
colnames(RS)[5] <- "Best_guess_binomial"

Birds <- left_join(Birds, RS)
Birds$log10_Range_area <- log10(Birds$Range_area_sq_km)

Mammals <- left_join(Mammals, RS)
Mammals$log10_Range_area <- log10(Mammals$Range_area_sq_km)

Amphibians <- left_join(Amphibians, RS)
Amphibians$log10_Range_area <- log10(Amphibians$Range_area_sq_km)

Reptiles <- left_join(Reptiles, RS)
Reptiles$log10_Range_area <- log10(Reptiles$Range_area_sq_km)

## select all relevant traits
Traits_cat <- c("Diel_activity", "Specialisation", "Primary_diet")
Traits_cont <- c("log10_Body_mass_g","log10_Lifespan_proxy", "log10_Litter_size", "sqrt_Diet_breadth", "sqrt_Habitat_breadth_IUCN", "log10_Range_area")
Traits <- c(Traits_cat, Traits_cont, "Best_guess_binomial")

Birds <- Birds[, Traits]
Birds$Class <- "Birds"
Mammals <- Mammals[, Traits]
Mammals$Class <- "Mammals"
Amphibians <- Amphibians[, Traits]
Amphibians$Class <- "Amphibians"
Reptiles <- Reptiles[, Traits]
Reptiles$Class<- "Reptiles"

AllClasses <- rbind(Amphibians, Birds, Mammals, Reptiles)

## save complete trait data ready for use in validation
write.csv(AllClasses, "../../Results/Data_complete_traits/Vertebrate_complete.csv", row.names = FALSE)





## Investigating the effects of traits on species probability of occurrence in disturbed land uses
## taxonomic effects as random effects
## land use and use intensity in interaction with other traits


### VALIDATIONS ON COMPLETE TRAIT DATA SUBSET

## within-class models

library(dplyr)
library(StatisticalModels)
library(lme4)
library(MCMCglmm)
library(ggplot2)
library(ggpubr)
library(scales)
library(viridis)


#setwd("E:/3.Explanatory_traits/Code/")

## load model data
ModelData <- readRDS( "../../Results/Data_complete_traits/Model_data_within_class.rds")

## check and set levels
ModelData$LandUseGrouped %>% levels() %>% print()
ModelData$Use_intensity %>% levels() %>%  print()

any(is.na(ModelData$LandUseGrouped))
any(is.na(ModelData$Use_intensity))

## check primary diet
unique(ModelData$Primary_diet)
table(ModelData$Primary_diet)

# ## filter out specis with unknown range area
# ModelData <- ModelData %>% 
#   filter(!is.na(log10_Range_area))


#################################################################################################################
##### fitting models: Here, class-specific models with use intensity  #### GLMER version

FitGLMER_Mammals_Birds <- function(Data, VClass) {
  
  # subset for given class
  Data <- Data %>% 
    dplyr::filter(Class==VClass)
  print(length(unique(Data$Best_guess_binomial)))
  
  Start <- Sys.time()
  
    Model <- lme4::glmer(Occurrence ~ 
                           log10_Body_mass_g +
                           log10_Litter_size +
                           log10_Lifespan_proxy +
                           log10_Range_area +
                           sqrt_Habitat_breadth_IUCN +
                           sqrt_Diet_breadth +
                           
                           LandUseGrouped +
                           Use_intensity +
                           Specialisation +
                           Diel_activity +
                           Primary_diet +
                           
                           LandUseGrouped:log10_Body_mass_g +
                           LandUseGrouped:log10_Litter_size +
                           LandUseGrouped:log10_Lifespan_proxy +
                           LandUseGrouped:log10_Range_area +
                           LandUseGrouped:sqrt_Habitat_breadth_IUCN +
                           LandUseGrouped:sqrt_Diet_breadth +
                           LandUseGrouped:Use_intensity +
                           LandUseGrouped:Specialisation +
                           LandUseGrouped:Diel_activity +
                           LandUseGrouped:Primary_diet +
                           
                           Use_intensity:log10_Body_mass_g +
                           Use_intensity:log10_Litter_size +
                           Use_intensity:log10_Lifespan_proxy +
                           Use_intensity:log10_Range_area +
                           Use_intensity:sqrt_Habitat_breadth_IUCN +
                           Use_intensity:sqrt_Diet_breadth +
                           Use_intensity:Specialisation +
                           Use_intensity:Diel_activity +
                           Use_intensity:Primary_diet +
                           
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         
                         data = Data,
                         family = "binomial")

  End <- Sys.time()
  print(Start-End)
  return(Model)
}

FitGLMER_Amphibians <- function(Data, VClass) {
  
  # subset for given class
  Data <- Data %>% 
    dplyr::filter(Class==VClass)
  print(length(unique(Data$Best_guess_binomial)))
  
  Start <- Sys.time()
  
  Model <- lme4::glmer(Occurrence ~
                         #log10_Body_mass_g +
                         #log10_Litter_size +
                         #log10_Lifespan_proxy +
                         log10_Range_area +
                         sqrt_Habitat_breadth_IUCN +
                         sqrt_Diet_breadth +
                         
                         LandUseGrouped +
                         Use_intensity +
                         Specialisation +
                         Diel_activity +
                         
                         #LandUseGrouped:log10_Body_mass_g +
                         #LandUseGrouped:log10_Litter_size +
                         #LandUseGrouped:log10_Lifespan_proxy +
                         LandUseGrouped:log10_Range_area +
                         LandUseGrouped:sqrt_Habitat_breadth_IUCN +
                         LandUseGrouped:sqrt_Diet_breadth +
                         LandUseGrouped:Use_intensity +
                         LandUseGrouped:Specialisation +
                         LandUseGrouped:Diel_activity +
                         
                         #Use_intensity:log10_Body_mass_g +
                         #Use_intensity:log10_Litter_size +
                         #Use_intensity:log10_Lifespan_proxy +
                         Use_intensity:log10_Range_area +
                         Use_intensity:sqrt_Habitat_breadth_IUCN +
                         Use_intensity:sqrt_Diet_breadth +
                         Use_intensity:Specialisation +
                         Use_intensity:Diel_activity +
                         
                         (1|SS) +
                         (1|SSBS) +
                         (1|Best_guess_binomial),
                       
                       data = Data,
                       family = "binomial")
  
  End <- Sys.time()
  print(Start-End)
  return(Model)
}

FitGLMER_Reptiles <- function(Data, VClass) {
  
  # subset for given class
  Data <- Data %>% 
    dplyr::filter(Class==VClass)
  print(length(unique(Data$Best_guess_binomial)))
  
  Start <- Sys.time()
  
  Model <- lme4::glmer(Occurrence ~
                         #log10_Body_mass_g +
                         log10_Litter_size +
                         log10_Lifespan_proxy +
                         log10_Range_area +
                         sqrt_Habitat_breadth_IUCN +
                         sqrt_Diet_breadth +
                         
                         LandUseGrouped +
                         Use_intensity +
                         Specialisation +
                         Diel_activity +
                         
                         #LandUseGrouped:log10_Body_mass_g +
                         LandUseGrouped:log10_Litter_size +
                         LandUseGrouped:log10_Lifespan_proxy +
                         LandUseGrouped:log10_Range_area +
                         LandUseGrouped:sqrt_Habitat_breadth_IUCN +
                         LandUseGrouped:sqrt_Diet_breadth +
                         LandUseGrouped:Use_intensity +
                         LandUseGrouped:Specialisation +
                         LandUseGrouped:Diel_activity +
                         
                         #Use_intensity:log10_Body_mass_g +
                         Use_intensity:log10_Litter_size +
                         Use_intensity:log10_Lifespan_proxy +
                         Use_intensity:log10_Range_area +
                         Use_intensity:sqrt_Habitat_breadth_IUCN +
                         Use_intensity:sqrt_Diet_breadth +
                         Use_intensity:Specialisation +
                         Use_intensity:Diel_activity +
                         
                         (1|SS) +
                         (1|SSBS) +
                         (1|Best_guess_binomial),
                       
                       data = Data,
                       family = "binomial")
  
  End <- Sys.time()
  print(Start-End)
  return(Model)
}

## Reptiles ~ 5 minutes
Reptiles <- FitGLMER_Reptiles(ModelData, "Reptilia")
saveRDS(Reptiles, "../Results/Validations/Reptiles_Diet.rds")

RepData <- ModelData[ModelData$Class=="Reptilia",]
nrow(RepData[is.na(RepData$log10_Lifespan_proxy),]) /nrow(RepData) * 100
nrow(RepData[is.na(RepData$log10_Body_mass_g),]) /nrow(RepData) * 100
nrow(RepData[is.na(RepData$log10_Litter_size),]) /nrow(RepData) * 100
nrow(RepData[is.na(RepData$Diel_activity),]) /nrow(RepData) * 100
nrow(RepData[is.na(RepData$sqrt_Diet_breadth),]) /nrow(RepData) * 100
nrow(RepData[is.na(RepData$sqrt_Habitat_breadth_IUCN),]) /nrow(RepData) * 100
nrow(RepData[is.na(RepData$Specialisation),]) /nrow(RepData) * 100
nrow(RepData[is.na(RepData$log10_Range_area),]) /nrow(RepData) * 100


## mammals ~ 1.3 hours
Mammals <- FitGLMER_Mammals_Birds(ModelData, "Mammalia")
saveRDS(Mammals, "../Results/Validations/Mammals_Diet.rds")

## Amphibians ~ 25 minutes 
AmphData <- ModelData[ModelData$Class=="Amphibia",]
nrow(AmphData[is.na(AmphData$log10_Lifespan_proxy),]) /nrow(AmphData) * 100
nrow(AmphData[is.na(AmphData$log10_Body_mass_g),]) /nrow(AmphData) * 100
nrow(AmphData[is.na(AmphData$log10_Litter_size),]) /nrow(AmphData) * 100
nrow(AmphData[is.na(AmphData$Diel_activity),]) /nrow(AmphData) * 100
nrow(AmphData[is.na(AmphData$sqrt_Diet_breadth),]) /nrow(AmphData) * 100
nrow(AmphData[is.na(AmphData$sqrt_Habitat_breadth_IUCN),]) /nrow(AmphData) * 100
nrow(AmphData[is.na(AmphData$Specialisation),]) /nrow(AmphData) * 100
nrow(AmphData[is.na(AmphData$log10_Range_area),]) /nrow(AmphData) * 100

Amphibians <- FitGLMER_Amphibians(ModelData, "Amphibia")
saveRDS(Amphibians, "../Results/Validations/Amphibians_Diet.rds")

## Birds
Birds <- FitGLMER_Mammals_Birds(ModelData, "Aves")
saveRDS(Birds, "../Results/Validations/Birds_Diet.rds")

#########################################################################################################


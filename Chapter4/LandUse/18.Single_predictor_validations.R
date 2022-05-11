## Investigating the effects of traits on species probability of occurrence in disturbed land uses
## taxonomic effects as random effects
## land use and use intensity in interaction with other traits

## VALIDATIONS ON COMPLETE TRAIT DATA

## within-class models

library(dplyr)
library(StatisticalModels)
library(lme4)
library(MCMCglmm)

#setwd("E:/3.Explanatory_traits/Code/")

## load model data
ModelData <- readRDS( "../Results/Data_complete_traits/Model_data_within_class.rds")

## check and set levels
ModelData$LandUseGrouped %>% levels() %>% print()
ModelData$Use_intensity %>% levels() %>%  print()

any(is.na(ModelData$LandUseGrouped))
any(is.na(ModelData$Use_intensity))

## check primary diet
unique(ModelData$Primary_diet)
table(ModelData$Primary_diet[ModelData$Class=="Amphibia"])
table(ModelData$Primary_diet[ModelData$Class=="Reptilia"])
table(ModelData$Primary_diet[ModelData$Class=="Aves"])
table(ModelData$Primary_diet[ModelData$Class=="Mammalia"])

table(is.na(ModelData$Primary_diet[ModelData$Class=="Amphibia"]))
table(is.na(ModelData$Primary_diet[ModelData$Class=="Reptilia"]))
table(is.na(ModelData$Primary_diet[ModelData$Class=="Mammalia"]))
table(is.na(ModelData$Primary_diet[ModelData$Class=="Aves"]))

## filter out specis with unknown range area
ModelData <- ModelData %>% 
  filter(!is.na(log10_Range_area))

#################################################################################################################
##### fitting models: Here, class-specific models with use intensity  #### GLMER version with single predictor 

## CATEGORICAL PREDICTORS

FitGLMER_single_cat_pred <- function(Data, VClass, Cat) {
  
  # subset for given class
  Data <- Data %>% 
    dplyr::filter(Class==VClass)
  print(length(unique(Data$Best_guess_binomial)))
  
  Start <- Sys.time()
  
  if(Cat=="Primary_diet"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           Primary_diet +
                           LandUseGrouped:Primary_diet +
                           Use_intensity:Primary_diet +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  if(Cat=="Diel_activity"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           Diel_activity +
                           LandUseGrouped:Diel_activity +
                           Use_intensity:Diel_activity +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  if(Cat=="Specialisation"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           Specialisation +
                           LandUseGrouped:Specialisation +
                           Use_intensity:Specialisation +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  End <- Sys.time()
  print(Start-End)
  return(Model)
}

## Amphibians
Amphibians1 <- FitGLMER_single_cat_pred(ModelData, "Amphibia", Cat = "Specialisation")
Amphibians2 <- FitGLMER_single_cat_pred(ModelData, "Amphibia", Cat = "Diel_activity")
Amphibians3 <- FitGLMER_single_cat_pred(ModelData, "Amphibia", Cat = "Primary_diet")
saveRDS(Amphibians1, "../Results/Single_predictor_models/Validations/Amphibians_specialisation.rds")
saveRDS(Amphibians2, "../Results/Single_predictor_models/Validations/Amphibians_diel_activity.rds")
saveRDS(Amphibians3, "../Results/Single_predictor_models/Validations/Amphibians_diet.rds")

## Reptiles 
Reptiles1 <- FitGLMER_single_cat_pred(ModelData, "Reptilia", Cat = "Specialisation")
Reptiles2 <- FitGLMER_single_cat_pred(ModelData, "Reptilia", Cat = "Diel_activity")
Reptiles3 <- FitGLMER_single_cat_pred(ModelData, "Reptilia", Cat = "Primary_diet")
saveRDS(Reptiles1, "../Results/Single_predictor_models/Validations/Reptiles_specialisation.rds")
saveRDS(Reptiles2, "../Results/Single_predictor_models/Validations/Reptiles_diel_activity.rds")
saveRDS(Reptiles3, "../Results/Single_predictor_models/Validations/Reptiles_diet.rds")

## Mammals
Mammals1 <- FitGLMER_single_cat_pred(ModelData, "Mammalia", Cat = "Specialisation")
Mammals2 <- FitGLMER_single_cat_pred(ModelData, "Mammalia", Cat = "Diel_activity")
Mammals3 <- FitGLMER_single_cat_pred(ModelData, "Mammalia", Cat = "Primary_diet")
saveRDS(Mammals1, "../Results/Single_predictor_models/Validations/Mammals_specialisation.rds")
saveRDS(Mammals2, "../Results/Single_predictor_models/Validations/Mammals_diel_activity.rds")
saveRDS(Mammals3, "../Results/Single_predictor_models/Validations/Mammals_diet.rds")

## Birds 
Birds1 <- FitGLMER_single_cat_pred(ModelData, "Aves", Cat = "Specialisation")
Birds2 <- FitGLMER_single_cat_pred(ModelData, "Aves", Cat = "Diel_activity")
Birds3 <- FitGLMER_single_cat_pred(ModelData, "Aves", Cat = "Primary_diet")
saveRDS(Birds1, "../Results/Single_predictor_models/Validations/Birds_specialisation.rds")
saveRDS(Birds2, "../Results/Single_predictor_models/Validations/Birds_diel_activity.rds")
saveRDS(Birds3, "../Results/Single_predictor_models/Validations/Birds_diet.rds")


## CONTINUOUS PREDICTORS (RA, BM, HB, DB, LCS, LP)

FitGLMER_single_cont_pred <- function(Data, VClass, Cont) {
  
  # subset for given class
  Data <- Data %>% 
    dplyr::filter(Class==VClass)
  print(length(unique(Data$Best_guess_binomial)))
  
  Start <- Sys.time()
  
  if(Cont=="log10_Range_area"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           log10_Range_area +
                           LandUseGrouped:log10_Range_area +
                           Use_intensity:log10_Range_area +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  if(Cont=="log10_Body_mass_g"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           log10_Body_mass_g +
                           LandUseGrouped:log10_Body_mass_g +
                           Use_intensity:log10_Body_mass_g +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  if(Cont=="log10_Litter_size"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           log10_Litter_size +
                           LandUseGrouped:log10_Litter_size +
                           Use_intensity:log10_Litter_size +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  if(Cont=="sqrt_Habitat_breadth_IUCN"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           sqrt_Habitat_breadth_IUCN +
                           LandUseGrouped:sqrt_Habitat_breadth_IUCN +
                           Use_intensity:sqrt_Habitat_breadth_IUCN +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  if(Cont=="sqrt_Diet_breadth"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           sqrt_Diet_breadth +
                           LandUseGrouped:sqrt_Diet_breadth +
                           Use_intensity:sqrt_Diet_breadth +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  if(Cont=="log10_Lifespan_proxy"){
    Model <- lme4::glmer(Occurrence ~ 
                           LandUseGrouped +
                           Use_intensity +
                           log10_Lifespan_proxy +
                           LandUseGrouped:log10_Lifespan_proxy +
                           Use_intensity:log10_Lifespan_proxy +
                           (1|SS) +
                           (1|SSBS) + 
                           (1|Best_guess_binomial),
                         data = Data,
                         family = "binomial")
  }
  
  
  End <- Sys.time()
  print(Start-End)
  return(Model)
}

Run_on_all <- function(Class) {
  
  Traits <- c("log10_Range_area", "log10_Body_mass_g", "log10_Litter_size", "log10_Lifespan_proxy", "sqrt_Habitat_breadth_IUCN", "sqrt_Diet_breadth")
  
  for(t in Traits){
    cat("Fitting for", t)
    Model <- FitGLMER_single_cont_pred(ModelData, VClass = Class, Cont = t)
    Path <- paste0("../Results/Single_predictor_models/Validations/", Class, "_", t, ".rds")
    saveRDS(Model, Path)
  }
}

## Amphibians 

Run_on_all("Amphibia")
Run_on_all("Reptilia")
Run_on_all("Mammalia")
Run_on_all("Aves")


# Amphibians1 <- FitGLMER_single_cont_pred(ModelData, "Amphibia", Cont = "log10_Body_mass_g")
# saveRDS(Amphibians1, "../Results/Single_predictor_models/Amphibians_Body_mass.rds")
# 
# Mammals1 <- FitGLMER_single_cont_pred(ModelData, "Mammalia", Cont = "log10_Body_mass_g")
# saveRDS(Mammals1, "../Results/Single_predictor_models/Mammals_Body_mass.rds")
# 
# Reptiles1 <- FitGLMER_single_cont_pred(ModelData, "Reptilia", Cont = "log10_Body_mass_g")
# saveRDS(Reptiles1, "../Results/Single_predictor_models/Reptiles_Body_mass.rds")
# 
# Birds1 <- FitGLMER_single_cont_pred(ModelData, "Aves", Cont = "log10_Body_mass_g")
# saveRDS(Birds1, "../Results/Single_predictor_models/Birds_Body_mass.rds")



## PGLS species-level models ## validation using complete trait data
library(ape)
library(geiger)
library(phytools)
library(dplyr)
library(caper)
library(car)
library(ggplot2)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  return(Phylogeny)
}

stepwise.vif <- function (dataset,
                          metrics,
                          vif.threshold = 5,
                          verbose = F){
  

  
  dataset$dummy <- rnorm(nrow(dataset)) # generates a dummy normal distribution
  output <- metrics # those are all the candidate variables
  step.count <- 1
  output.results <- list() 
  
  repeat {
    
    # VIF scores for the linear model: dummy ~ all variables; returns VIF or GVIF, depending on Degrees of freedom (if any variables has more than one then GVIF is computed rather than VIF)
    vif.scores <- vif(lm(as.formula(paste0(
      "dummy~", paste0(output,
                       collapse = "+")
    )), data = dataset))
    na.coefficients <- Reduce('|', is.nan(vif.scores)) # for scores that are NA -- if so stop the function
    if (na.coefficients) {
      stop("NA coefficient in a regression model.")
    }
    
    #browser()
    # Select VIF scores
    #vif.scores <- vif.scores[,1] 
    
    ## output.results stores VIF calculated at each step and output stores variable names that are selected 
    output.results[[step.count]] <-
      sort(vif.scores, decreasing = F) # sort VIF scores
    
    vif.scores <- vif.scores[vif.scores >= vif.threshold]
    
    ## If all VIF scores are under the threshold values then keep all variables and stop the function. Else drop the variable
    # that has a VIF above the threshold and proceed iteratively
    if (length(vif.scores) == 0)
      break
    
    drop.var <-
      names(vif.scores[vif.scores == max(vif.scores)]) # select variable names that have VIF score more than the threshold
    if (verbose) {
      print(paste0(
        "Step ",
        step.count,
        " - Exclude ",
        drop.var,
        " (VIF = ",
        max(vif.scores),
        ")"
      ))
    }
    step.count <- step.count + 1
    output <- output[!output %in% drop.var]
  }
  
  names(output.results) <- paste0("Iteration ", 1:step.count)
  names(output.results)[length(output.results)] <- "Final"
  
  if(length(output.results)==1) {output.results <- as.data.frame(output.results)}
  
  return(list(Selected_vars=output, VIF=output.results))
}

GetVIF_scores <- function(Data, Class, Diet) {
  Data_Class <- Data[Data$Class==Class,]
  
  if(Diet){
    Metrics <-  c("log10_Body_mass_g",
                  "log10_Lifespan_proxy",
                  "log10_Litter_size",
                  "log10_Range_area",
                  "sqrt_Habitat_breadth_IUCN",
                  "sqrt_Diet_breadth",
                  "Specialisation", 
                  "Diel_activity", 
                  "Primary_diet")
  } else{
    Metrics <-  c("log10_Body_mass_g",
                  "log10_Lifespan_proxy",
                  "log10_Litter_size",
                  "log10_Range_area",
                  "sqrt_Habitat_breadth_IUCN",
                  "Specialisation", 
                  "Diel_activity",
                  "Trophic_level")   
  }
  
  print(nrow(Data_Class))
  
  return(stepwise.vif(dataset = Data_Class, vif.threshold = 5, verbose = TRUE, metrics = Metrics)
  )
}


#######################################################################################################

## load complete trait data
Traits <- read.csv("../Data/Vertebrate_complete_traits_validation.csv")

## load CENFA results
Index <- read.csv("../../Range_maps_work/Results/7_GenerateSpeciesIndices/SpeciesIndices.csv")
Index$Species <- paste0("sp",Index$Index)
Sensitivity_5 <- read.csv("../Results/2.CENFA_summary_dataframes/Vertebrates_5km.csv")
Sensitivity_5 <- left_join(Sensitivity_5, Index[,c(2,4)], by="Species")
colnames(Sensitivity_5)[5] <- "Best_guess_binomial"

## match
Data <- left_join(Traits, Sensitivity_5[,c(1,5)], by="Best_guess_binomial")
ggplot(Data, aes(log10_Range_area, log10(sensitivity))) + geom_point() + facet_wrap(~Class)

## filter for species with RA>100km2 (log10(100)=2)
Data_filtered <- Data %>% 
  dplyr::filter(log10_Range_area>2)
ggplot(Data_filtered, aes(log10_Range_area, log10(sensitivity))) + geom_point() + facet_wrap(~Class)


## checking the data
levels(as.factor(Data_filtered$Primary_diet))
Data_filtered$Primary_diet <- as.character(Data_filtered$Primary_diet)
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="SE"] <- "PL|SE"
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="PL"] <- "PL|SE"
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="NE"] <- "FR|NE"
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="FR"] <- "FR|NE"
levels(as.factor(Data_filtered$Primary_diet))


## checking VIF score within each class with stepwise.vif
GetVIF_scores(Data_filtered, "Mammals", TRUE) 
GetVIF_scores(Data_filtered, "Birds", TRUE)  

Data_filtered$Primary_diet[Data_filtered$Class=="Reptiles"] %>%  table()
GetVIF_scores(Data_filtered, "Reptiles", TRUE)  # reptiles: exclude primary diet

Data_filtered$Primary_diet[Data_filtered$Class=="Amphibians"] %>%  table()
Data_filtered$Diel_activity[Data_filtered$Class=="Amphibians"] %>%  table()
Data_filtered$Specialisation[Data_filtered$Class=="Amphibians"] %>%  table()
Data_filtered$sqrt_Diet_breadth[Data_filtered$Class=="Amphibians"] %>%  table()
Data_filtered$log10_Litter_size[Data_filtered$Class=="Amphibians" & !is.na(Data_filtered$log10_Litter_size)] %>%  length()
Data_filtered$log10_Lifespan_proxy[Data_filtered$Class=="Amphibians" & !is.na(Data_filtered$log10_Lifespan_proxy)] %>%  length()
Data_filtered$log10_Body_mass_g[Data_filtered$Class=="Amphibians" & !is.na(Data_filtered$log10_Body_mass_g)] %>%  length()
Data_filtered$sqrt_Habitat_breadth_IUCN[Data_filtered$Class=="Amphibians" & !is.na(Data_filtered$sqrt_Habitat_breadth_IUCN)] %>%  length()

GetVIF_scores(Data_filtered, "Amphibians", TRUE) # amphibians: exclude primary diet


#######################################################################################################

## fitting PGLS models


## ! for amphibians / diet breadth and primary diet need to be considered in separate models

# phylogenies
Phylo_mammals <- read.newick("../Data/phylogenies/consensus/Mammals.nwk") %>% .Format_tiplabels
Phylo_amphibians <- read.newick("../Data/phylogenies/consensus/Amphibians.nwk") %>% .Format_tiplabels
Phylo_birds <- read.newick("../Data/phylogenies/consensus/Birds.nwk") %>% .Format_tiplabels
Phylo_reptiles <- read.newick("../Data/phylogenies/consensus/Reptiles.nwk") %>% .Format_tiplabels

## assemble data for PGLS models 

Assemble_data <- function(Data, Class, Phylo){
  Data <- Data[Data$Class==Class,]
  
  # dropping tips from tree (for species that are not in the dataset but that are represented in the tree)
  to_drop <- setdiff(Phylo$tip.label, Data$Best_guess_binomial)
  Phylo <- drop.tip(Phylo, to_drop)
  
  return(comparative.data(phy = Phylo, data = Data,
                          names.col = "Best_guess_binomial",
                          vcv = TRUE, vcv.dim = 2,
                          na.omit = FALSE,
                          warn.dropped = TRUE))
}

Traits %>%  group_by(Class) %>%  summarise(C=n())

## for species with RS < 100km square
gc()
memory.limit(size = 500000000)
Mammals_data_filtered <- Assemble_data(Data_filtered, "Mammals", Phylo_mammals)
saveRDS(Mammals_data_filtered, "../Results/9.Data_PGLS_validations/Mammals_filtered.rds")
gc()
memory.limit(size = 100000000)
Amphibians_data_filtered <- Assemble_data(Data_filtered, "Amphibians", Phylo_amphibians)
saveRDS(Amphibians_data_filtered, "../Results/9.Data_PGLS_validations/Amphibians_filtered.rds")
gc()
memory.limit(size = 100000000)
Reptiles_data_filtered <- Assemble_data(Data_filtered, "Reptiles", Phylo_reptiles)
saveRDS(Reptiles_data_filtered, "../Results/9.Data_PGLS_validations/Reptiles_filtered.rds")
gc()
memory.limit(size = 500000000)
Birds_data_filtered <- Assemble_data(Data_filtered, "Birds", Phylo_birds)
saveRDS(Birds_data_filtered, "../Results/9.Data_PGLS_validations/Birds_filtered.rds")



#######################################################################################################

## running the class-specific PGLS models with class-specific phylogenies
## for species whose range sizes is >100km2

memory.limit(size = 7000000)

## WITH DIET INFORMATION

FitPGLS <- function(AssembledData, Diet) {
  
  if(Diet){
    
    Model <-  pgls(log10(sensitivity) ~ 
                     log10_Body_mass_g +
                     log10_Lifespan_proxy +
                     log10_Litter_size +
                     log10_Range_area +
                     sqrt_Habitat_breadth_IUCN +
                     sqrt_Diet_breadth +
                     Specialisation +
                     Diel_activity +
                     Primary_diet,
                   data = AssembledData,
                   lambda = "ML")
  }else{
    Model <-  pgls(log10(sensitivity) ~ 
                     log10_Body_mass_g +
                     log10_Lifespan_proxy +
                     log10_Litter_size +
                     log10_Range_area +
                     sqrt_Habitat_breadth_IUCN +
                     sqrt_Diet_breadth +
                     Specialisation +
                     Diel_activity,
                   
                   data = AssembledData,
                   lambda = "ML")
  }
  return(Model)
}

# mammals (~26 minutes)
Mammals_filtered <- readRDS("../Results/9.Data_PGLS_validations/Mammals_filtered.rds") # 4714 species
Mammals_filtered$data %>%  nrow()
Mammals_filtered$phy $tip.label %>%  length
gc()
Start <- Sys.time() 
PGLS_Mammals <- FitPGLS(Mammals_filtered, Diet= TRUE)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Mammals, "../Results/10.PGLS_complete_models_results/Mammals.rds")
summary(PGLS_Mammals)

# reptiles
Reptiles_filtered <- readRDS("../Results/9.Data_PGLS_validations/Reptiles_filtered.rds") # 7333 species
Reptiles_filtered$data %>%  nrow()
Reptiles_filtered$phy $tip.label %>%  length
gc()
Start <- Sys.time() 
PGLS_Reptiles <- FitPGLS(Reptiles_filtered, Diet = FALSE)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Reptiles, "../Results/10.PGLS_complete_models_results/Reptiles.rds")

# birds 
Birds_filtered <- readRDS("../Results/9.Data_PGLS_validations/Birds_filtered.rds") # 10198 species
Birds_filtered$data %>%  nrow()
Birds_filtered$phy $tip.label %>%  length
gc()
Start <- Sys.time() 
PGLS_Birds <- FitPGLS(Birds_filtered, Diet = TRUE)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Birds, "../Results/10.PGLS_complete_models_results/Birds.rds")

# amphibians
Amphibians_filtered <- readRDS("../Results/9.Data_PGLS_validations/Amphibians_filtered.rds") # 4537 species
Amphibians_filtered$data %>% nrow()
Amphibians_filtered$phy 
gc()
Start <- Sys.time() 
PGLS_Amphibians <- FitPGLS(Amphibians_filtered, Diet=FALSE)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Amphibians, "../Results/10.PGLS_complete_models_results/Amphibians.rds")

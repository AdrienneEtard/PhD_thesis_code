## PGLS species-level models
library(ape)
library(geiger)
library(phytools)
library(dplyr)
library(caper)
library(car)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  return(Phylogeny)
}

stepwise.vif <- function (dataset,
                          metrics,
                          vif.threshold = 5,
                          verbose = F){
  
  #browser()
  
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
    
    # Select VIF scores
    vif.scores <- vif.scores[,1] 
    
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

#setwd("D:/PhD/PhD_R_work/5.Climatic_niche_space/Code/")

#######################################################################################################

# trait + CENFA data
Data_all <- read.csv("../Results/Traits_CENFA_RS_all.csv") ## with all species
levels(as.factor(Data_all$Primary_diet))
Data_all$Primary_diet <- as.character(Data_all$Primary_diet)
Data_all$Primary_diet[Data_all$Primary_diet=="SE"] <- "PL|SE"
Data_all$Primary_diet[Data_all$Primary_diet=="PL"] <- "PL|SE"
Data_all$Primary_diet[Data_all$Primary_diet=="NE"] <- "FR|NE"
Data_all$Primary_diet[Data_all$Primary_diet=="FR"] <- "FR|NE"

#######################################################################################################
## checking diet -- particularly for herptiles

## amphibians
Data_all[Data_all$Class.x=="Amphibians",] %>%  group_by(Primary_diet) %>%  summarise(C=n())
Data_all[Data_all$Class.x=="Amphibians",] %>%  group_by(Trophic_level) %>%  summarise(C=n())

X <- Data_all[Data_all$Class.x=="Amphibians",] %>%  group_by(Primary_diet) 
X$Best_guess_binomial[X$Primary_diet=="PL|SE"]
X$Best_guess_binomial[X$Primary_diet=="VE"]

# for amphibians, need to drop: Plant/seed eaters; vertebrate eaters (only 5 species identified as such) and also, herbivores

## reptiles
Data_all[Data_all$Class.x=="Reptiles",] %>%  group_by(Primary_diet) %>%  summarise(C=n())
Data_all[Data_all$Class.x=="Reptiles",] %>%  group_by(Trophic_level) %>%  summarise(C=n())
Data_all[Data_all$Class.x=="Reptiles",] %>%  group_by(Trophic_level, Primary_diet) %>%  summarise(C=n())

## birds
Data_all[Data_all$Class.x=="Birds",] %>%  group_by(Primary_diet, Trophic_level) %>%  summarise(C=n()) %>% arrange(desc(C))

## mammals
Data_all[Data_all$Class.x=="Mammals",] %>%  group_by(Primary_diet, Trophic_level) %>%  summarise(C=n()) %>% arrange(desc(C))

#######################################################################################################


Data_filtered <- read.csv("../Results/Traits_CENFA_RS_filtered.csv") ## with species that have range size > 100 square kilometers
levels(as.factor(Data_filtered$Primary_diet))
Data_filtered$Primary_diet <- as.character(Data_filtered$Primary_diet)
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="SE"] <- "PL|SE"
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="PL"] <- "PL|SE"
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="NE"] <- "FR|NE"
Data_filtered$Primary_diet[Data_filtered$Primary_diet=="FR"] <- "FR|NE"

Data_all$log10_Range_area <- log10(Data_all$Range_area_sq_km)
Data_filtered$log10_Range_area <- log10(Data_filtered$Range_area_sq_km)

Data_all <- Data_all %>%  dplyr::select(-Class.y)
Data_filtered <- Data_filtered %>%  dplyr::select(-Class.y)

colnames(Data_all)[4] <- "Class"
colnames(Data_filtered)[4] <- "Class"

# checking VIF score within each class with stepwise.vif
GetVIF_scores(Data_all, "Mammals", TRUE) 
GetVIF_scores(Data_all, "Mammals", FALSE) 

GetVIF_scores(Data_all, "Birds", TRUE)  
GetVIF_scores(Data_all, "Birds", FALSE)  

GetVIF_scores(Data_all, "Amphibians", TRUE) 
GetVIF_scores(Data_all, "Amphibians", FALSE) 

GetVIF_scores(Data_all, "Reptiles", TRUE)  
GetVIF_scores(Data_all, "Reptiles", FALSE) 

GetVIF_scores(Data_filtered, "Mammals", TRUE) 
GetVIF_scores(Data_filtered, "Mammals", FALSE) 

GetVIF_scores(Data_filtered, "Birds", TRUE)  
GetVIF_scores(Data_filtered, "Birds", FALSE)  

GetVIF_scores(Data_filtered, "Amphibians", TRUE) 
GetVIF_scores(Data_filtered, "Amphibians", FALSE) 

GetVIF_scores(Data_filtered, "Reptiles", TRUE)  
GetVIF_scores(Data_filtered, "Reptiles", FALSE) 

## need to consider trophic level and primary diet in different models. Otherwise no other multicollinearity problem.
## use of trophic levels only (as PD was collected initially for PREDICTS species?)

## for amphibians, diet breadth and primary diet need to be considered in different models

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


## for species with RS < 100km square
gc()
memory.limit(size = 500000000)
Mammals_data_filtered <- Assemble_data(Data_filtered, "Mammals", Phylo_mammals)
saveRDS(Mammals_data_filtered, "../Results/5.Data_for_PGLS_models/Mammals_filtered.rds")
gc()
memory.limit(size = 100000000)
Amphibians_data_filtered <- Assemble_data(Data_filtered, "Amphibians", Phylo_amphibians)
saveRDS(Amphibians_data_filtered, "../Results/5.Data_for_PGLS_models/Amphibians_filtered.rds")
gc()
memory.limit(size = 100000000)
Reptiles_data_filtered <- Assemble_data(Data_filtered, "Reptiles", Phylo_reptiles)
saveRDS(Reptiles_data_filtered, "../Results/5.Data_for_PGLS_models/Reptiles_filtered.rds")
gc()
memory.limit(size = 500000000)
Birds_data_filtered <- Assemble_data(Data_filtered, "Birds", Phylo_birds)
saveRDS(Birds_data_filtered, "../Results/5.Data_for_PGLS_models/Birds_filtered.rds")


## for all species (robusteness test)
gc()
memory.limit(size = 500000000)
Mammals_data_all <- Assemble_data(Data_all, "Mammals", Phylo_mammals)
saveRDS(Mammals_data_all, "../Results/5.Data_for_PGLS_models/Mammals_all.rds")
gc()
memory.limit(size = 100000000)
Amphibians_data_all <- Assemble_data(Data_all, "Amphibians", Phylo_amphibians)
saveRDS(Amphibians_data_all, "../Results/5.Data_for_PGLS_models/Amphibians_all.rds")
gc()
memory.limit(size = 100000000)
Reptiles_data_all <- Assemble_data(Data_all, "Reptiles", Phylo_reptiles)
saveRDS(Reptiles_data_all, "../Results/5.Data_for_PGLS_models/Reptiles_all.rds")
gc()
memory.limit(size = 500000000)
Birds_data_all <- Assemble_data(Data_all, "Birds", Phylo_birds)
saveRDS(Birds_data_all, "../Results/5.Data_for_PGLS_models/Birds_all.rds")


#######################################################################################################

## running the class-specific PGLS models with class-specific phylogenies
## for species whose range sizes is >100km2

memory.limit(size = 7000000)

## NO DIET INFORMATION  --  TROPHIC LEVELS INSTEAD

FitPGLS_no_diet <- function(AssembledData, Degree) {
  
  Model_no_diet <- pgls(log10(sensitivity) ~ 
                                  poly(log10_Body_mass_g,Degree) +
                                  poly(log10_Lifespan_proxy,Degree) +
                                  poly(log10_Litter_size,Degree) +
                                  poly(log10_Range_area,Degree) +
                                  poly(sqrt_Habitat_breadth_IUCN,Degree) +
                                  Specialisation +
                                  Diel_activity +
                                  Trophic_level,
                                data = AssembledData,
                                lambda = "ML")
  
  return(Model_no_diet)
}

## WITH DIET INFORMATION

FitPGLS_diet <- function(AssembledData, Diet_breadth, Degree) {
  
  if(Diet_breadth){
  Model_diet <-  pgls(log10(sensitivity) ~ 
                        poly(log10_Body_mass_g,Degree) +
                        poly(log10_Lifespan_proxy,Degree) +
                        poly(log10_Litter_size,Degree) +
                        poly(log10_Range_area,Degree) +
                        poly(sqrt_Habitat_breadth_IUCN,Degree) +
                        poly(sqrt_Diet_breadth,Degree)+
                        Specialisation +
                        Diel_activity +
                        Primary_diet,
                      data = AssembledData,
                      lambda = "ML")
  }else{
    Model_diet <-  pgls(log10(sensitivity) ~ 
                          poly(log10_Body_mass_g,Degree) +
                          poly(log10_Lifespan_proxy,Degree) +
                          poly(log10_Litter_size,Degree) +
                          poly(log10_Range_area,Degree) +
                          poly(sqrt_Habitat_breadth_IUCN,Degree) +
                          #poly(sqrt_Diet_breadth,Degree)+
                          Specialisation +
                          Diel_activity +
                          Primary_diet,
                        data = AssembledData,
                        lambda = "ML")
  }
  return(Model_diet)
}

# mammals: 2 hours each model
Mammals_filtered <- readRDS("../Results/5.Data_for_PGLS_models/Mammals_filtered.rds") # 4714 species
Mammals_filtered$data %>%  nrow()
Mammals_filtered$phy $tip.label %>%  length

# gc()
# Start <- Sys.time() 
# PGLS_Mammals_filtered3 <- FitPGLS_no_diet(Mammals_filtered, Degree = 3)
# End <- Sys.time()
# print(Start-End)
# saveRDS(PGLS_Mammals_filtered3, "../Results/5.PGLS_models_results/Mammals_no_diet_v2.rds")

## degree 3
gc()
Start <- Sys.time() 
PGLS_Mammals_filtered_diet3 <- FitPGLS_diet(Mammals_filtered, Diet_breadth = TRUE, Degree=3)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Mammals_filtered_diet3, "../Results/5.PGLS_models_results/Mammals_diet_v2.rds")

# ## degree 2
# gc()
# Start <- Sys.time() 
# PGLS_Mammals_filtered2 <- FitPGLS_diet(Mammals_filtered,  Diet_breadth = TRUE, Degree = 2)
# End <- Sys.time()
# print(Start-End)
# saveRDS(PGLS_Mammals_filtered2, "../Results/5.PGLS_models_results/Other_degrees/Mammals_diet_degree2.rds")


# reptiles
Reptiles_filtered <- readRDS("../Results/5.Data_for_PGLS_models/Reptiles_filtered.rds") # 7333 species
Reptiles_filtered$data %>%  nrow()
Reptiles_filtered$phy $tip.label %>%  length
# gc()
# Start <- Sys.time() 
# PGLS_Reptiles_filtered <- FitPGLS_no_diet(Reptiles_filtered)
# End <- Sys.time()
# print(Start-End)
# saveRDS(PGLS_Reptiles_filtered, "../Results/5.PGLS_models_results/Reptiles_no_diet_v2.rds")

gc()
Start <- Sys.time() 
PGLS_Reptiles_filtered_diet <- FitPGLS_diet(Reptiles_filtered, Diet_breadth = TRUE)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Reptiles_filtered_diet, "../Results/5.PGLS_models_results/Reptiles_diet_v2.rds")

# birds 
Birds_filtered <- readRDS("../Results/5.Data_for_PGLS_models/Birds_filtered.rds") # 10198 species
Birds_filtered$data %>%  nrow()
Birds_filtered$phy $tip.label %>%  length
# gc()
# Start <- Sys.time() 
# PGLS_Birds_filtered <- FitPGLS_no_diet(Birds_filtered)
# End <- Sys.time()
# print(Start-End)
# saveRDS(PGLS_Birds_filtered, "../Results/5.PGLS_models_results/Birds_no_diet_v2.rds")

gc()
Start <- Sys.time() 
PGLS_Birds_filtered_diet <- FitPGLS_diet(Birds_filtered, Diet_breadth = TRUE)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Birds_filtered_diet, "../Results/5.PGLS_models_results/Birds_diet_v2.rds")


# amphibians
Amphibians_filtered <- readRDS("../Results/5.Data_for_PGLS_models/Amphibians_filtered.rds") # 4537 species
Amphibians_filtered$data %>% nrow()
Amphibians_filtered$phy 

# gc()
# Start <- Sys.time()
# print(Start)
# PGLS_Amphibians_filtered <- FitPGLS_no_diet(Amphibians_filtered)
# End <- Sys.time()
# print(Start-End)
# saveRDS(PGLS_Amphibians_filtered, "../Results/5.PGLS_models_results/Amphibians_no_diet_v2.rds")

gc()
Start <- Sys.time() 
PGLS_Amphibians_filtered_diet <- FitPGLS_diet(Amphibians_filtered, Diet_breadth = FALSE) # need to remove diet breadth here because of multicollinearity problem
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Amphibians_filtered_diet, "../Results/5.PGLS_models_results/Amphibians_diet_without_DB_v2.rds")


#######################################################################################################

## running the class-specific PGLS models with class-specific phylogenies
## for all species 



memory.limit(size = 7000000)


## mammals
Mammals_all <- readRDS("../Results/5.Data_for_PGLS_models/Mammals_All.rds") # 4844 species
Mammals_all$data %>%  nrow()
Mammals_all$phy $tip.label %>%  length
gc()
Start <- Sys.time() 
PGLS_Mammals_all_diet3 <- FitPGLS_diet(Mammals_all, Diet_breadth = TRUE, Degree=3)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Mammals_all_diet3, "../Results/5.PGLS_models_results/All_species/Mammals_diet_v2.rds")


## birds (~20hours)
Birds_all <- readRDS("../Results/5.Data_for_PGLS_models/Birds_all.rds") # xxxxx species
Birds_all$data %>%  nrow()
Birds_all$phy $tip.label %>%  length
gc()
Start <- Sys.time() 
PGLS_Birds_all_diet <- FitPGLS_diet(Birds_all, Diet_breadth = TRUE, Degree=3)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Birds_all_diet, "../Results/5.PGLS_models_results/All_species/Birds_diet_v2.rds")

## reptiles
Reptiles_all <- readRDS("../Results/5.Data_for_PGLS_models/Reptiles_all.rds") # 10340 species
Reptiles_all$data %>%  nrow()
Reptiles_all$phy $tip.label %>%  length
gc()
Start <- Sys.time() 
PGLS_Reptiles_all_diet <- FitPGLS_diet(Reptiles_all, Diet_breadth = TRUE, Degree=3)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Reptiles_all_diet, "../Results/5.PGLS_models_results/All_species/Reptiles_diet_v2.rds")

## amphibians
Amphibians_all <- readRDS("../Results/5.Data_for_PGLS_models/Amphibians_all.rds") # xxxxx species
Amphibians_all$data %>%  nrow()
Amphibians_all$phy $tip.label %>%  length
gc()
Start <- Sys.time() 
PGLS_Amphibians_all_diet <- FitPGLS_diet(Amphibians_all, Diet_breadth = FALSE, Degree=3)
End <- Sys.time()
print(Start-End)
saveRDS(PGLS_Amphibians_all_diet, "../Results/5.PGLS_models_results/All_species/Amphibians_diet_v2.rds")


#######################################################################################################

## Looking at the models -- filtered for species whose RS > threshold (100km2)

## mammals, with diet
MammalsModel <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Mammals_diet_v2.rds")
MammalsModel$call
MammalsModel$fitted %>%  nrow() ## sample size: 4212
as.data.frame(summary(MammalsModel)$coefficients)
MammalsModel$aic
MammalsModel$aicc

## reptiles, with diet
ReptilesModel <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Reptiles_diet_v2.rds")
ReptilesModel$call
ReptilesModel$data$data$Primary_diet %>% table
ReptilesModel$data$data$sqrt_Diet_breadth %>% table
ReptilesModel$fitted %>%  nrow() ## sample size: 7330

## birds, with diet
BirdsModel <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Birds_diet_v2.rds")
BirdsModel$call
BirdsModel$fitted %>%  nrow() ## sample size: 10198

## amphibians, without diet breadth (multicollinearity issues with DB)
AmphibiansModel <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Amphibians_diet_without_DB_v2.rds")
AmphibiansModel$call
AmphibiansModel$data$data$Primary_diet %>% table
AmphibiansModel$data$data$sqrt_Diet_breadth %>% table
AmphibiansModel$fitted %>%  nrow() ## sample size: 4537

########## ANOVAS on the models

## anova on the models
AnovaMammals <- anova(MammalsModel)
AnovaAmphibians <- anova(AmphibiansModel)
saveRDS(AnovaMammals, "../Results/5.PGLS_models_results/ANOVA_diet_models/Mammals_ANOVA.rds")
saveRDS(AnovaAmphibians, "../Results/5.PGLS_models_results/ANOVA_diet_models/Amphibians_ANOVA.rds")

gc()
memory.limit(size = 500000000)
AnovaBirds <- anova(BirdsModel)
saveRDS(AnovaBirds, "../Results/5.PGLS_models_results/ANOVA_diet_models/Birds_ANOVA.rds")

gc()
AnovaReptiles <- anova(ReptilesModel)
saveRDS(AnovaReptiles, "../Results/5.PGLS_models_results/ANOVA_diet_models/Reptiles_ANOVA.rds")


#######################################################################################################
########## Effects of model with RS>100km versus all species

MammalsF <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Mammals_diet_v2.rds")
AmphibiansF <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Amphibians_diet_without_DB_v2.rds")
ReptilesF <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Reptiles_diet_v2.rds")
BirdsF <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Birds_diet_v2.rds")

ReptilesA <- readRDS("../Results/5.PGLS_models_results/All_species/Reptiles_diet_v2.rds")
BirdsA <- readRDS("../Results/5.PGLS_models_results/All_species/Birds_diet_v2.rds")
AmphibiansA <- readRDS("../Results/5.PGLS_models_results/All_species/Amphibians_diet_v2.rds")
MammalsA <- readRDS("../Results/5.PGLS_models_results/All_species/Mammals_diet_v2.rds")

## difference in species number
nrow(MammalsA$fitted) - nrow(MammalsF$fitted)
nrow(BirdsA$fitted) - nrow(BirdsF$fitted)
nrow(AmphibiansA$fitted) - nrow(AmphibiansF$fitted)
nrow(ReptilesA$fitted) - nrow(ReptilesF$fitted)

Plot_estimated_effects <- function(DFALIST, DFFLIST, ClassList) {
  
  par(family='serif', tcl=0.2, cex.lab=1.3, mfrow=c(4,3), oma=c(0.5,0.5,2,0.5))
  
  for(i in 1:4){
    
    DFA <- DFALIST[[i]]
    DFF <- DFFLIST[[i]]
    Class <- ClassList[[i]]
  
  DFA <- as.data.frame(summary(DFA)$coefficients)
  DFA$Which <- "All species"
  colnames(DFA)[c(1,2)] <- c("Estimate_all", "StdError_all")
  
  DFF <- as.data.frame(summary(DFF)$coefficients)
  DFF$Which <- "Species with range area>100km2"
  colnames(DFF)[c(1,2)] <- c("Estimate_filtered", "StdError_filtered")
  
  S <- cbind(DFA[,c(1,2)], DFF[,c(1,2)])
  
  r1 <- which(rownames(S)=="poly(log10_Range_area, Degree)1")
  r2 <- which(rownames(S)=="poly(log10_Range_area, Degree)2")
  r3 <- which(rownames(S)=="poly(log10_Range_area, Degree)3")
  
  plot(S$Estimate_all[-c(r1,r2,r3)]~S$Estimate_filtered[-c(r1,r2,r3)], pch=19, 
       xlab="Excl. species", ylab="All species",
       main=paste(Class, "\nAll effects (except range area)"))
  abline(a=0, b=1, lty="dashed", col="blue")
  
  plot(S$Estimate_all[c(r1,r2,r3)]~S$Estimate_filtered[c(r1,r2,r3)], pch=19, col="darkred", 
       xlab="Excl.species", ylab="All species",
       main=paste(Class, "\nEffects for range area"))
  abline(a=0, b=1, lty="dashed", col="blue")

  plot(log(S$StdError_all)~log(S$StdError_filtered), pch=19, 
       xlab="Excl. species", ylab="All species",
       main=paste(Class, "\nEstimated standard errors"))
  abline(a=0, b=1, lty="dashed", col="blue")
  
  # mtext(Class,                   
  #       side = 3,
  #       line = 0,
  #       outer = TRUE, 
  #       cex = 1.2,
  #       adj=0.1)
  
  }
}

pdf("G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/ClimateChangeModelsEstimates.pdf",
    width=6, height=8, family="Times", pointsize=11)
Plot_estimated_effects(list(AmphibiansA, BirdsA, MammalsA, ReptilesA),
                       list(AmphibiansF, BirdsF, MammalsF, ReptilesF),
                       Class=list("Amphibians", "Birds", "Mammals", "Reptiles"))
dev.off()


#######################################################################################################
## summaries for SI

library(stargazer)
stargazer(summary(AmphibiansF)$coefficients, summary=FALSE, digits=2)
stargazer(summary(BirdsF)$coefficients, summary=FALSE, digits=2)
stargazer(summary(MammalsF)$coefficients, summary=FALSE, digits=2)
stargazer(summary(ReptilesF)$coefficients, summary=FALSE, digits=2)

stargazer(summary(AmphibiansA)$coefficients, summary=FALSE, digits=2)
stargazer(summary(BirdsA)$coefficients, summary=FALSE, digits=2)
stargazer(summary(MammalsA)$coefficients, summary=FALSE, digits=2)
stargazer(summary(ReptilesA)$coefficients, summary=FALSE, digits=2)


ReptilesA$fitted %>%  nrow()

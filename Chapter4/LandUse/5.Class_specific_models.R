## Investigating the effects of traits on species probability of occurrence in disturbed land uses
## taxonomic effects as random effects
## land use and use intensity in interaction with other traits

## within-class models

library(dplyr)
library(StatisticalModels)
library(lme4)
library(MCMCglmm)

setwd("E:/3.Explanatory_traits/Code/")

## load model data
ModelData <- readRDS("../Results/Model_data_within_class.rds")

## check and set levels
ModelData$LandUseGrouped %>% levels() %>% print()
ModelData$Use_intensity %>% levels() %>%  print()

any(is.na(ModelData$LandUseGrouped))
any(is.na(ModelData$Use_intensity))

## check primary diet
unique(ModelData$Primary_diet)
table(ModelData$Primary_diet)

## filter out specis with unknown ramge area
ModelData <- ModelData %>% 
  filter(!is.na(log10_Range_area))

saveRDS(ModelData, "../Data/Data_for_class_specific_models.rds")

#################################################################################################################
##### fitting models: Here, class-specific models with use intensity  #### GLMER version

FitGLMER <- function(Data, VClass, Diet) {
 
  # subset for given class
  Data <- Data %>% 
    dplyr::filter(Class==VClass)
  print(length(unique(Data$Best_guess_binomial)))
  
  Start <- Sys.time()
  
  if(Diet){
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
  }

  if(!Diet){
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

                           LandUseGrouped:log10_Body_mass_g +
                           LandUseGrouped:log10_Litter_size +
                           LandUseGrouped:log10_Lifespan_proxy +
                           LandUseGrouped:log10_Range_area +
                           LandUseGrouped:sqrt_Habitat_breadth_IUCN +
                           LandUseGrouped:sqrt_Diet_breadth +
                           LandUseGrouped:Use_intensity +
                           LandUseGrouped:Specialisation +
                           LandUseGrouped:Diel_activity +

                           Use_intensity:log10_Body_mass_g +
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
  }

  End <- Sys.time()
  print(Start-End)
  return(Model)
}


## Amphibians sample size 307 ~25 minutes
Amphibians <- FitGLMER(ModelData, "Amphibia", Diet = FALSE)
saveRDS(Amphibians, "../Results/GLMER_Models/Class_specific/Use_intensity/Amphibians_Diet.rds")

## Reptiles sample size 305 ~ 50 minutes
Reptiles <- FitGLMER(ModelData, "Reptilia", Diet = FALSE)
saveRDS(Reptiles, "../Results/GLMER_Models/Class_specific/Use_intensity/Reptiles_Diet.rds")

## mammals sample size 532 ~ 2 hours
Mammals <- FitGLMER(ModelData, "Mammalia", Diet = TRUE)
saveRDS(Mammals, "../Results/GLMER_Models/Class_specific/Use_intensity/Mammals_Diet.rds")

## Birds sample size ~ 2963
Birds <- FitGLMER(ModelData, "Aves", Diet = TRUE)
saveRDS(Birds, "../Results/GLMER_Models/Class_specific/Use_intensity/Birds_Diet.rds")





#################################################################################################################
##### fitting models: Here, class-specific models with use intensity  #### mcmcglmm version

FitBayes <- function(Data, VClass, Npar) {
  
  Priors <- list(B = list(mu = rep(0,Npar), V = diag(Npar) *(1 + pi^2/3)), # n_parameters is the number of estimated parameters, need prior information for each of these
                 R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                          G3 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000))) # weakly informative
  
  Data <- Data %>%
    dplyr::filter(Class==VClass) %>%
    as.data.frame()
  
  print(length(unique(Data$Best_guess_binomial)))
  print(levels(Data$LandUseGrouped))
  print(unique(Data$LandUseGrouped))
  
  gc()
  memory.limit(size=90000000)
  
  if(VClass=="Mammalia"|VClass=="Aves"){
    
    cat("Starting mcmcmGLMM for", VClass, "...")
    
    Start <- Sys.time()
    print(Start)
    Model <- MCMCglmm::MCMCglmm(fixed=Occurrence ~ 
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
                                  Use_intensity:Primary_diet,
                                
                                random= ~ SS + SSBS + Best_guess_binomial,
                                data=Data,
                                family="categorical",
                                prior=Priors,
                                nitt=100000,
                                thin=20,
                                burnin=5000,
                                pl=TRUE,  # saving posterior distribution for latent variables
                                pr=TRUE,
                                singular.ok = TRUE)
    End <- Sys.time()
  }
  
  if(VClass=="Amphibia"|VClass=="Reptilia"){
    
    cat("Starting mcmcmGLMM for", VClass, "...")
    
    Start <- Sys.time()
    print(Start)
    Model <- MCMCglmm::MCMCglmm(fixed=Occurrence ~
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
                                  
                                  LandUseGrouped:log10_Body_mass_g +
                                  LandUseGrouped:log10_Litter_size +
                                  LandUseGrouped:log10_Lifespan_proxy +
                                  LandUseGrouped:log10_Range_area +
                                  LandUseGrouped:sqrt_Habitat_breadth_IUCN +
                                  LandUseGrouped:sqrt_Diet_breadth +
                                  LandUseGrouped:Use_intensity +
                                  LandUseGrouped:Specialisation +
                                  LandUseGrouped:Diel_activity +
                                  
                                  Use_intensity:log10_Body_mass_g +
                                  Use_intensity:log10_Litter_size +
                                  Use_intensity:log10_Lifespan_proxy +
                                  Use_intensity:log10_Range_area +
                                  Use_intensity:sqrt_Habitat_breadth_IUCN +
                                  Use_intensity:sqrt_Diet_breadth +
                                  Use_intensity:Specialisation +
                                  Use_intensity:Diel_activity,
                                
                                random= ~ SS + SSBS + Best_guess_binomial,
                                data=Data,
                                family="categorical",
                                prior=Priors,
                                nitt=100000,
                                thin=20,
                                burnin=5000,
                                pl=TRUE,  # saving posterior distribution for latent variables
                                pr=TRUE,
                                singular.ok = TRUE)
    End <- Sys.time()
  }
  
  print(Start-End)
  return(Model)
}

Amphibians <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Amphibians_Diet.rds")
Birds <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Birds_Diet.rds")
Mammals <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Mammals_Diet.rds")
Reptiles <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Reptiles_Diet.rds")

## npar = n estimated par from lme4 + n dropped coefficients
NPar_Amphibians <- nrow(summary(Amphibians)$coefficients)+1
NPar_Birds <- nrow(summary(Birds)$coefficients)
NPar_Reptiles <- nrow(as.data.frame(summary(BayesAmphibians)$solutions))
NPar_Mammals <- nrow(summary(Mammals)$coefficients)+2

## checking predictors
any(is.na(ModelData$LandUseGrouped))
any(is.na(ModelData$Use_intensity))
any(is.na(ModelData$log10_Body_mass_g))
any(is.na(ModelData$log10_Lifespan_proxy))
any(is.na(ModelData$log10_Litter_size))
any(is.na(ModelData$sqrt_Diet_breadth))
any(is.na(ModelData$sqrt_Habitat_breadth_IUCN))
any(is.na(ModelData$Diel_activity))
any(is.na(ModelData$Specialisation))

any(is.na(ModelData$log10_Range_area))
ModelData <- ModelData[!is.na(ModelData$log10_Range_area),]

## fitting Bayesian model for amphibians ## 45 minutes
gc()
BayesAmphibians <- FitBayes(Data = ModelData, VClass = "Amphibia", Npar = NPar_Amphibians)
saveRDS(BayesAmphibians, "../Results/MCMC_glmm_models/Class_specific/Amphibians_Diet.rds")

## fitting Bayesian model for mammals ## 
gc()
BayesMammals <- FitBayes(Data = ModelData, VClass = "Mammalia", Npar = NPar_Mammals)
saveRDS(BayesMammals, "../Results/MCMC_glmm_models/Class_specific/Mammals_Diet.rds")

## fitting Bayesian model for birds ## 20 hours 
gc()
BayesBirds <- FitBayes(Data = ModelData, VClass = "Aves", Npar = NPar_Birds)
saveRDS(BayesBirds, "../Results/MCMC_glmm_models/Class_specific/Birds_Diet.rds")

## reptiles ~ 2 hours
gc()
BayesReptiles <- FitBayes(Data = ModelData, VClass = "Reptilia", Npar = 71)
saveRDS(BayesReptiles, "../Results/MCMC_glmm_models/Class_specific/Reptiles_Diet.rds")


####

## plotting estimates GLMER against mcmcglmm

Plot_glmm_vs_mcmc <- function(Model_GLMER, Model_Bayesian) {
  
  SummaryGLMER <- as.data.frame(summary(Model_GLMER)$coefficients)
  SummaryBayesian <- as.data.frame(summary(Model_Bayesian)$solutions)
  
  ## checking order of estimate coefficients
  print(unique(rownames(SummaryGLMER)==rownames(SummaryBayesian)))
  
  if(nrow(SummaryGLMER)!=nrow(SummaryBayesian)){
    to_drop <- setdiff(rownames(SummaryBayesian),rownames(SummaryGLMER))
    print(to_drop)
    x <- which(rownames(SummaryBayesian) %in% to_drop)
    SummaryBayesian <- SummaryBayesian[-x,]
  }
  
  ## plotting
  SummaryGLMER$Which <- "GLMER"
  SummaryBayesian$Which <- "MCMCGLMM"
  
  SummaryGLMER$CI_up <- SummaryGLMER$Estimate + 1.96*SummaryGLMER$`Std. Error`
  SummaryGLMER$CI_low <- SummaryGLMER$Estimate - 1.96*SummaryGLMER$`Std. Error`
  
  colnames(SummaryGLMER)
  colnames(SummaryBayesian)
  
  colnames(SummaryBayesian)[c(1,2,3)] <- c("Estimate", "CI_low", "CI_up")
  
  #browser()
  
  plot(SummaryBayesian$pMCMC ~ SummaryGLMER$`Pr(>|z|)`, pch=19)
  
  SummaryBayesian <- SummaryBayesian %>% 
    dplyr::select(Estimate, CI_up, CI_low, Which)
  
  SummaryGLMER <- SummaryGLMER %>% 
    dplyr::select(Estimate, CI_up, CI_low, Which)
  
  Summaries <- rbind(SummaryGLMER, SummaryBayesian)
  
  ## plotting estimates 
  #dev.off()
  
  #pdf(file = "../Results/Figures/8.Occurrence_models/SI_Coefficients_Bayesian_against_GLMER.pdf", width=8, height=7, pointsize = 13)
  par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,1), oma=c(2,1,2,1))
  par(mfrow=c(2,2))
  
  plot(log(Summaries$Estimate[Summaries$Which=="MCMCGLMM"]),
       log(Summaries$Estimate[Summaries$Which=="GLMER"]),
       pch=19, xlab="mcmcglmm", ylab="GLMER")
  abline(a=0, b=1, lty="dashed", col="red")
  title("(a) Coefficient estimates", adj = 0, line = 1)
  
  plot(Summaries$CI_up[Summaries$Which=="MCMCGLMM"],
       Summaries$CI_up[Summaries$Which=="GLMER"],
       pch=19, xlab="mcmcglmm", ylab="GLMER")
  abline(a=0, b=1, lty="dashed", col="red")
  title("(b) Confidence/credible interval - upper bound", adj = 0, line = 1)
  
  plot(Summaries$CI_low[Summaries$Which=="MCMCGLMM"],
       Summaries$CI_low[Summaries$Which=="GLMER"],
       pch=19, xlab="mcmcglmm", ylab="GLMER")
  title("(c) Confidence/credible interval - lower bound", adj = 0, line = 1)
  abline(a=0, b=1, lty="dashed", col="red")
  
  #dev.off() 
}

Plot_glmm_vs_mcmc(Birds, BayesBirds)
Plot_glmm_vs_mcmc(Mammals, BayesMammals)
Plot_glmm_vs_mcmc(Amphibians, BayesAmphibians)
Plot_glmm_vs_mcmc(Reptiles, BayesReptiles)


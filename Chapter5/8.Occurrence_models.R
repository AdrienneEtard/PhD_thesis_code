## Occurrence model: do species with lower residual RMR persist better in disturbed land uses?

## testing for effect of residual RMR on occurrence probability.

###############################################################################################################

library(dplyr)
library(lme4)
library(ggplot2)
library(StatisticalModels)
library(performance)
library(ggmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)
library(raster)
library(viridis)
library(car)
library(MCMCglmm)


## Function to check variance inflation factors (stepwise selection)
# https://rdrr.io/github/software-analytics/Rnalytica/src/R/stepwise.vif.R
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

## for plotting

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

Predict_effects <- function(Newdata, Model, rescale, Cont_or_Cat, seMultiplier, LU_n) {
  
  # function for backtransforming
  logit2prob <- function(logit){
    return(exp(logit)/(1+exp(logit)))
  }
  
  # which coefs are not estimated?
  coefs <- mvrnorm(n=1, mu=fixef(Model), Sigma=vcov(Model))
  mm <- model.matrix(terms(Model), Newdata)
  print(setdiff(names(coefs), colnames(mm)))
  
  # predictions, for categories
  if(Cont_or_Cat=="categorical"){
    
    preds <- sapply(X=1:1000, FUN=function(i){
      
      coefs <- mvrnorm(n=1, mu=fixef(Model), Sigma=vcov(Model))
      mm <- model.matrix(terms(Model), Newdata)
      
      y <- mm %*% coefs
      
      # backtranforming
      y <- logit2prob(y)
      
      # rescaling
      if(rescale){
        
        # inittialisation
        seq <- 1:LU_n
        y[seq] <- y[seq]/y[seq[1]]*100
        
        # for loop to rescale all values
        for(i in 1:(nrow(Newdata)/LU_n-1)){
          seq <- seq + LU_n
          y[seq] <- y[seq]/y[seq[1]]*100
        }
      }
      
      return(y)
    })
    
    preds <- data.frame(Median=apply(X=preds, MARGIN=1, FUN=median),
                        Upper=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.975),
                        Lower=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.025))
    
    preds <- cbind(preds, Newdata)
    return(preds)
  }
  
  if(Cont_or_Cat=="continuous"){
    
    mm <- model.matrix(terms(Model), Newdata)
    
    # predictions
    preds <- mm %*% fixef(Model)
    
    # backtransform
    Newdata$Estimate <- logit2prob(preds)
    
    # getting the errors around the predicted values
    # variance covariance matrix and matrix multiplication
    VarCoVar <- as.matrix(vcov(Model))
    VarCoVar <- base::tcrossprod(VarCoVar, mm)
    
    # matrix multiplication to get estimates and take diagonal matrix
    pvar <- diag(mm %*% VarCoVar)
    
    # estimate error of predictions and backtransform
    Lower <- preds - seMultiplier*sqrt(pvar)
    Upper <- preds + seMultiplier*sqrt(pvar)
    
    Newdata$Lower <- logit2prob(Lower)
    Newdata$Upper <- logit2prob(Upper)
    
    return(Newdata)
    
  }
}

Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")


Labels=c("PV", "SV", "PF", "PA", "CR", "UR")

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols)


################################################################################################################################################

## loading and checking PREDICTS data
PredictsTraits <- "../Results/PredictsTraits.rds" %>% 
  readRDS() 
levels(PredictsTraits$LandUse)
levels(PredictsTraits$Use_intensity)

## addind in Diet data
DietData <- readRDS("E:/3.Explanatory_traits/Results/Model_data_allvertebrates.rds")
DietData <- unique(DietData[, c("Best_guess_binomial", "sqrt_Diet_breadth", "Primary_diet")])
DietData$Primary_diet %>%  table()
PredictsTraits <- left_join(PredictsTraits, DietData)
rm(DietData)

PredictsTraits$Primary_diet  %>%  table()
PredictsTraits$Primary_diet  <- as.character(PredictsTraits$Primary_diet)
PredictsTraits$Primary_diet[PredictsTraits$Primary_diet=="SE"] <- "PL|SE"
PredictsTraits$Primary_diet[PredictsTraits$Primary_diet=="PL"] <- "PL|SE"
PredictsTraits$Primary_diet[PredictsTraits$Primary_diet=="FR"] <- "FR|NE"
PredictsTraits$Primary_diet[PredictsTraits$Primary_diet=="NE"] <- "FR|NE"
PredictsTraits$Primary_diet  %>%  table()


## checking VIF scores for predictors - is there a multicollinearity problem? -- NO

VIFVal <- stepwise.vif(dataset=PredictsTraits, 
                       metrics = c("Residual_BMR_log_log", 
                                   "Annual_mean_temperature",
                                   "LandUse", 
                                   "Use_intensity",
                                   "Trophic_level"), 
                       vif.threshold = 5, 
                       verbose = TRUE)$VIF

print(VIFVal)

VIFVal <- stepwise.vif(dataset=PredictsTraits, 
                       metrics = c("Residual_BMR_log_log", 
                                   "Annual_mean_temperature",
                                   "LandUse", 
                                   "Use_intensity",
                                   "Primary_diet"), 
                       vif.threshold = 5, 
                       verbose = TRUE)$VIF

print(VIFVal)


################################################################################################################################################
#################################################################################################################################################

## WITH TROPHIC LEVEL


## FIRST: MODEL SELECTION TO FIND BEST-FITTING MODEL FOR OCCURRENCE PROBABILITY
## there are 4 candidate main effects and 6 possible interaction terms.

Start <- Sys.time() 
SelectModel <- GLMERSelect(modelData = PredictsTraits, 
                           responseVar = "Occurrence",   
                           fitFamily = "binomial",
                           fixedFactors = c("LandUse", "Use_intensity", "Trophic_level"),
                           fixedTerms = list(Residual_BMR_log_log=1),
                           fitInteractions=TRUE,
                           randomStruct = "(1|SS)+(1|SSBS)+(1|Best_guess_binomial)",
                           verbose = TRUE)
End <- Sys.time()
print(Start-End) ## take a bit more than 2 days

SelectModel$final.call

"Occurrence~LandUse+Use_intensity+Trophic_level+poly(Residual_BMR_log_log,1)+
LandUse:Use_intensity+
LandUse:Trophic_level+
LandUse:poly(Residual_BMR_log_log,1)+
Use_intensity:Trophic_level+
Use_intensity:poly(Residual_BMR_log_log,1)+
(1|SS)+(1|SSBS)+(1|Best_guess_binomial)"

saveRDS(SelectModel, "../Results/8.Occurrence_models/ModelSelect.rds")

SelectModel <- readRDS("../Results/8.Occurrence_models/ModelSelect.rds")
BestModel <- SelectModel$model
summary(BestModel)

#################################################################################################################################################
## check full model against:
## model with the 3-way interaction between land use, trophic level and residual RMR.
## model with the 3-way interaction between land use intensity, trophic level and residual RMR.

## full model (4 hours):
Start <- Sys.time()
print(Start)
Model_Full <- lme4::glmer(Occurrence ~
                             LandUse + 
                             Use_intensity +
                             Trophic_level +
                             Residual_BMR_log_log +
                             LandUse:Use_intensity +
                             LandUse:Trophic_level +
                             LandUse:Residual_BMR_log_log +
                             Use_intensity:Trophic_level +
                             Use_intensity:Residual_BMR_log_log +
                             Trophic_level:Residual_BMR_log_log +
                             (1|SS) +
                             (1|SSBS) + 
                             (1|Best_guess_binomial),
                           data=PredictsTraits,
                           family="binomial")
End <- Sys.time()
print(Start-End)

saveRDS(Model_Full, "../Results/8.Occurrence_models/Model_occurrence_full.rds")

## model with three-way interaction LU:TL:ResRMR (6 hours):
Start <- Sys.time()
print(Start)
Model_3_way <- lme4::glmer(Occurrence ~
                            LandUse + 
                            Use_intensity +
                            Trophic_level +
                            Residual_BMR_log_log +
                            LandUse:Use_intensity +
                            LandUse:Trophic_level +
                            LandUse:Residual_BMR_log_log +
                            Use_intensity:Trophic_level +
                            Use_intensity:Residual_BMR_log_log +
                            Trophic_level:Residual_BMR_log_log +
                            LandUse:Trophic_level:Residual_BMR_log_log +
                            (1|SS) +
                            (1|SSBS) + 
                            (1|Best_guess_binomial),
                          data=PredictsTraits,
                          family="binomial")
End <- Sys.time()
print(Start-End)

saveRDS(Model_3_way, "../Results/8.Occurrence_models/Model_occurrence_3_way_LU_TL_Res.rds")

summary(Model_3_way)
anova(Model_3_way)

## which of the full model versus the model with the 3-way interaction fits better?
## test.

anova(Model_Full, Model_3_way)

## model 3-way fits better

## adding in 3-way interaction between Use intensity, trophic level and residual RMR

Start <- Sys.time()
print(Start)
Model_3_way_UI <- lme4::glmer(Occurrence ~
                             LandUse + 
                             Use_intensity +
                             Trophic_level +
                             Residual_BMR_log_log +
                             LandUse:Use_intensity +
                             LandUse:Trophic_level +
                             LandUse:Residual_BMR_log_log +
                             Use_intensity:Trophic_level +
                             Use_intensity:Residual_BMR_log_log +
                             Trophic_level:Residual_BMR_log_log +
                             LandUse:Trophic_level:Residual_BMR_log_log +
                             Use_intensity:Trophic_level:Residual_BMR_log_log +
                             (1|SS) +
                             (1|SSBS) + 
                             (1|Best_guess_binomial),
                           data=PredictsTraits,
                           family="binomial")
End <- Sys.time()
print(Start-End)

saveRDS(Model_3_way_UI, "../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")


## which of the 3-way models fits better?

anova(Model_3_way, Model_3_way_UI)

## the model with both 3-way interaction fits better than the model with only one 3-way interaction. 

summary(Model_3_way_UI)

# #################################################################################################################################################
# ## model ~6 hours with residual BMR (log - log) + TEMPERATURE as an additional variable
# Start <- Sys.time()
# print(Start)
# Model_temp <- lme4::glmer(Occurrence ~
#                        LandUse + 
#                        Use_intensity +
#                        Trophic_level +
#                        Residual_BMR_log_log +
#                        Annual_mean_temperature +
#                        LandUse:Use_intensity +
#                        LandUse:Trophic_level +
#                        LandUse:Residual_BMR_log_log +
#                        LandUse:Annual_mean_temperature +
#                        Use_intensity:Trophic_level +
#                        Use_intensity:Residual_BMR_log_log +
#                        Use_intensity:Annual_mean_temperature +
#                        Trophic_level:Annual_mean_temperature +
#                        Residual_BMR_log_log:Annual_mean_temperature +
#                        (1|SS) +
#                        (1|SSBS) + 
#                        (1|Best_guess_binomial),
#                      data=PredictsTraits,
#                      family="binomial")
# End <- Sys.time()
# print(Start-End)
# summary(Model_temp)
# saveRDS(Model_temp, "../Results/8.Occurrence_models/ModelTemperature.rds")

##############################################################################################################################

## model ~11 hours with residual BMR (log - log) WITH A BAYESIAN FRAMEWORK to check for robustness of estimated coefficients

Model3way <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")

## number of parameters to estimate (determined from the model summary lme4)
nB <- nrow(summary(Model3way)$coefficients)

## priors
Priors <- list(B = list(mu = rep(0, nB), V = diag(nB) *(1 + pi^2/3)), # nB is the number of estimated parameters, need prior information for each of these
                R = list(V = 1, fix = 1),
                G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                         G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                         G3 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000))) # weakly informative

# PRIORS: 
# residual variance is fixed at 1 (R = list(V = 1, fix = 1))
# for B: using a prior that is relatively flat on the probability scale (cf Course Notes, mcmcglmm package)
# G: parameter expansion to improve mixing of the chains. G structure = random effects structure.

any(is.na(PredictsTraits$LandUse))
any(is.na(PredictsTraits$Use_intensity))

any(is.na(PredictsTraits$Trophic_level))
any(is.na(PredictsTraits$Use_intensity))
any(is.na(PredictsTraits$Residual_BMR_log_log))

PredictsTraits <- PredictsTraits %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity))

gc()
memory.limit(size=90000000)
Start <- Sys.time()
print(Start)
Model <- MCMCglmm::MCMCglmm(fixed=Occurrence ~
                              LandUse + 
                              Use_intensity +
                              Trophic_level +
                              Residual_BMR_log_log +
                              LandUse:Use_intensity +
                              LandUse:Trophic_level +
                              LandUse:Residual_BMR_log_log +
                              Use_intensity:Trophic_level +
                              Use_intensity:Residual_BMR_log_log +
                              Trophic_level:Residual_BMR_log_log +
                              LandUse:Trophic_level:Residual_BMR_log_log +
                              Use_intensity:Trophic_level:Residual_BMR_log_log,
                            random= ~ SS + SSBS + Best_guess_binomial,
                            data=PredictsTraits,
                            family="categorical",
                            prior=Priors,
                            nitt=70000,
                            thin=20,
                            burnin=5000,
                            pl=TRUE,  # saving posterior distribution for latent variables
                            pr=TRUE,
                            singular.ok = TRUE) # saving posterior distribution for random effects   

End <- Sys.time() ##12 hours
print(Start-End)
summary(Model)
saveRDS(Model, "../Results/8.Occurrence_models/ModelOccurrenceBayesian.rds")


## checking coefficients - Bayesian against lme4
## plotting coefficients estimated from GLMER against coefficients estimated from Bayesian model

Model_GLMER <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")
Model_Bayesian <- readRDS("../Results/8.Occurrence_models/ModelOccurrenceBayesian.rds")

## checking robustness of Bayesian model (trace plots)
plot(Model_Bayesian)

gc()
SummaryGLMER <- as.data.frame(summary(Model_GLMER)$coefficients)
SummaryBayesian <- as.data.frame(summary(Model_Bayesian)$solutions)

## checking order of estimate coefficients
rownames(SummaryGLMER)==rownames(SummaryBayesian)

## plotting
SummaryGLMER$Which <- "GLMER"
SummaryBayesian$Which <- "MCMCGLMM"

SummaryGLMER$CI_up <- SummaryGLMER$Estimate + 1.96*SummaryGLMER$`Std. Error`
SummaryGLMER$CI_low <- SummaryGLMER$Estimate - 1.96*SummaryGLMER$`Std. Error`

colnames(SummaryGLMER)
colnames(SummaryBayesian)

SummaryGLMER <- SummaryGLMER %>% 
  dplyr::select(Estimate, CI_up, CI_low, Which)

colnames(SummaryBayesian)[c(1,2,3)] <- c("Estimate", "CI_low", "CI_up")

SummaryBayesian <- SummaryBayesian %>% 
  dplyr::select(Estimate, CI_up, CI_low, Which)

Summaries <- rbind(SummaryGLMER, SummaryBayesian)

## plotting estimates 
dev.off()

pdf(file = "../Results/Figures/8.Occurrence_models/SI_Coefficients_Bayesian_against_GLMER.pdf", width=8, height=7, pointsize = 13)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,1), oma=c(2,1,2,1))
par(mfrow=c(2,2))

plot(Summaries$Estimate[Summaries$Which=="MCMCGLMM"],
     Summaries$Estimate[Summaries$Which=="GLMER"],
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

dev.off()

Diffs <- SummaryGLMER$Estimate-SummaryBayesian$Estimate
hist(Diffs)
max(Diffs)
min(Diffs)



################################################################################################################################################
#################################################################################################################################################

## WITH DIET


## FIRST: MODEL SELECTION TO FIND BEST-FITTING MODEL FOR OCCURRENCE PROBABILITY
## there are 4 candidate main effects and 6 possible interaction terms.

## not run
# Start <- Sys.time() 
# SelectModel <- GLMERSelect(modelData = PredictsTraits, 
#                            responseVar = "Occurrence",   
#                            fitFamily = "binomial",
#                            fixedFactors = c("LandUse", "Use_intensity", "Primary_diet"),
#                            fixedTerms = list(Residual_BMR_log_log=1),
#                            fitInteractions=TRUE,
#                            randomStruct = "(1|SS)+(1|SSBS)+(1|Best_guess_binomial)",
#                            verbose = TRUE)
# End <- Sys.time()
# print(Start-End) ## take a bit more than 2 days
# 
# SelectModel$final.call
# 
# "Occurrence~LandUse+Use_intensity+Trophic_level+poly(Residual_BMR_log_log,1)+
# LandUse:Use_intensity+
# LandUse:Trophic_level+
# LandUse:poly(Residual_BMR_log_log,1)+
# Use_intensity:Trophic_level+
# Use_intensity:poly(Residual_BMR_log_log,1)+
# (1|SS)+(1|SSBS)+(1|Best_guess_binomial)"
# 
# saveRDS(SelectModel, "../Results/8.Occurrence_models/ModelSelect.rds")
# 
# SelectModel <- readRDS("../Results/8.Occurrence_models/ModelSelect.rds")
# BestModel <- SelectModel$model
# summary(BestModel)

#################################################################################################################################################
## check full model against:
## model with the 3-way interaction between land use, trophic level and residual RMR.
## model with the 3-way interaction between land use intensity, trophic level and residual RMR.

## full model (8 hours):
Start <- Sys.time()
print(Start)
Model_Full_PD <- lme4::glmer(Occurrence ~
                            LandUse + 
                            Use_intensity +
                            Primary_diet +
                            Residual_BMR_log_log +
                            LandUse:Use_intensity +
                            LandUse:Primary_diet +
                            LandUse:Residual_BMR_log_log +
                            Use_intensity:Primary_diet +
                            Use_intensity:Residual_BMR_log_log +
                            Primary_diet:Residual_BMR_log_log +
                            (1|SS) +
                            (1|SSBS) + 
                            (1|Best_guess_binomial),
                          data=PredictsTraits,
                          family="binomial")
End <- Sys.time()
print(Start-End)

saveRDS(Model_Full_PD, "../Results/8.Occurrence_models/Model_occurrence_full_PD.rds")

## model with three-way interaction LU:PD:ResRMR (14 hours):
Start <- Sys.time()
print(Start)
Model_3_way_PD <- lme4::glmer(Occurrence ~
                             LandUse + 
                             Use_intensity +
                             Primary_diet +
                             Residual_BMR_log_log +
                             LandUse:Use_intensity +
                             LandUse:Primary_diet +
                             LandUse:Residual_BMR_log_log +
                             Use_intensity:Primary_diet +
                             Use_intensity:Residual_BMR_log_log +
                             Primary_diet:Residual_BMR_log_log +
                             LandUse:Primary_diet:Residual_BMR_log_log +
                             (1|SS) +
                             (1|SSBS) + 
                             (1|Best_guess_binomial),
                           data=PredictsTraits,
                           family="binomial")
End <- Sys.time()
print(Start-End)

saveRDS(Model_3_way_PD, "../Results/8.Occurrence_models/Model_occurrence_3_way_LU_TL_Res_PD.rds")


## which of the full model versus the model with the 3-way interaction fits better?
## test.

Model_3_way_PD <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_TL_Res_PD.rds")
Model_Full_PD <- readRDS("../Results/8.Occurrence_models/Model_occurrence_full_PD.rds")
  
anova(Model_Full_PD, Model_3_way_PD)
## model 3-way fits better

## adding in 3-way interaction between Use intensity, primary diet and residual RMR ## 15 hours

Start <- Sys.time()
print(Start)
Model_3_way_UI_PD <- lme4::glmer(Occurrence ~
                                LandUse + 
                                Use_intensity +
                                Primary_diet +
                                Residual_BMR_log_log +
                                LandUse:Use_intensity +
                                LandUse:Primary_diet +
                                LandUse:Residual_BMR_log_log +
                                Use_intensity:Primary_diet +
                                Use_intensity:Residual_BMR_log_log +
                                Primary_diet:Residual_BMR_log_log +
                                LandUse:Primary_diet:Residual_BMR_log_log +
                                Use_intensity:Primary_diet:Residual_BMR_log_log +
                                (1|SS) +
                                (1|SSBS) + 
                                (1|Best_guess_binomial),
                              data=PredictsTraits,
                              family="binomial")
End <- Sys.time()
print(Start-End)

saveRDS(Model_3_way_UI_PD, "../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res_PD.rds")


Model_3_way_PD <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_TL_Res_PD.rds")
Model_3_way_UI_PD <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res_PD.rds")

## which of the 3-way models fits better?
anova(Model_3_way_PD, Model_3_way_UI_PD)

## the model with both 3-way interaction fits better than the model with only one 3-way interaction. 

## not run (Bayesian with PD)

##############################################################################################################################

## model ~11 hours with residual BMR (log - log) WITH A BAYESIAN FRAMEWORK to check for robustness of estimated coefficients

Model3way <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")

## number of parameters to estimate (determined from the model summary lme4)
nB <- nrow(summary(Model3way)$coefficients)

## priors
Priors <- list(B = list(mu = rep(0, nB), V = diag(nB) *(1 + pi^2/3)), # nB is the number of estimated parameters, need prior information for each of these
               R = list(V = 1, fix = 1),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                        G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000),
                        G3 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*1000))) # weakly informative

# PRIORS: 
# residual variance is fixed at 1 (R = list(V = 1, fix = 1))
# for B: using a prior that is relatively flat on the probability scale (cf Course Notes, mcmcglmm package)
# G: parameter expansion to improve mixing of the chains. G structure = random effects structure.

any(is.na(PredictsTraits$LandUse))
any(is.na(PredictsTraits$Use_intensity))

any(is.na(PredictsTraits$Trophic_level))
any(is.na(PredictsTraits$Use_intensity))
any(is.na(PredictsTraits$Residual_BMR_log_log))

PredictsTraits <- PredictsTraits %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity))

gc()
memory.limit(size=90000000)
Start <- Sys.time()
print(Start)
Model <- MCMCglmm::MCMCglmm(fixed=Occurrence ~
                              LandUse + 
                              Use_intensity +
                              Trophic_level +
                              Residual_BMR_log_log +
                              LandUse:Use_intensity +
                              LandUse:Trophic_level +
                              LandUse:Residual_BMR_log_log +
                              Use_intensity:Trophic_level +
                              Use_intensity:Residual_BMR_log_log +
                              Trophic_level:Residual_BMR_log_log +
                              LandUse:Trophic_level:Residual_BMR_log_log +
                              Use_intensity:Trophic_level:Residual_BMR_log_log,
                            random= ~ SS + SSBS + Best_guess_binomial,
                            data=PredictsTraits,
                            family="categorical",
                            prior=Priors,
                            nitt=70000,
                            thin=20,
                            burnin=5000,
                            pl=TRUE,  # saving posterior distribution for latent variables
                            pr=TRUE,
                            singular.ok = TRUE) # saving posterior distribution for random effects   

End <- Sys.time() ##12 hours
print(Start-End)
summary(Model)
saveRDS(Model, "../Results/8.Occurrence_models/ModelOccurrenceBayesian.rds")

##############################################################################################################################

## checking coefficients - Bayesian against lme4
## plotting coefficients estimated from GLMER against coefficients estimated from Bayesian model

gc()
memory.limit(size=90000000000)

Model_GLMER <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")
Model_Bayesian <- readRDS("../Results/8.Occurrence_models/ModelOccurrenceBayesian.rds")

## checking robustness of Bayesian model (trace plots)
# plot(Model_Bayesian)

gc()
SummaryGLMER <- as.data.frame(summary(Model_GLMER)$coefficients)
SummaryBayesian <- as.data.frame(summary(Model_Bayesian)$solutions)

## checking order of estimate coefficients
rownames(SummaryGLMER)==rownames(SummaryBayesian)

## plotting
SummaryGLMER$Which <- "GLMER"
SummaryBayesian$Which <- "MCMCGLMM"

SummaryGLMER$CI_up <- SummaryGLMER$Estimate + 1.96*SummaryGLMER$`Std. Error`
SummaryGLMER$CI_low <- SummaryGLMER$Estimate - 1.96*SummaryGLMER$`Std. Error`

colnames(SummaryGLMER)
colnames(SummaryBayesian)

SummaryGLMER <- SummaryGLMER %>% 
  dplyr::select(Estimate, CI_up, CI_low, Which)

colnames(SummaryBayesian)[c(1,2,3)] <- c("Estimate", "CI_low", "CI_up")

SummaryBayesian <- SummaryBayesian %>% 
  dplyr::select(Estimate, CI_up, CI_low, Which)

Summaries <- rbind(SummaryGLMER, SummaryBayesian)

## plotting estimates 
dev.off()

pdf(file = "../Results/Figures/8.Occurrence_models/SI_Coefficients_Bayesian_against_GLMER.pdf", width=8, height=7, pointsize = 13)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,1), oma=c(2,1,2,1))
par(mfrow=c(2,2))

plot(Summaries$Estimate[Summaries$Which=="MCMCGLMM"],
     Summaries$Estimate[Summaries$Which=="GLMER"],
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

dev.off()

Diffs <- SummaryGLMER$Estimate-SummaryBayesian$Estimate
hist(Diffs)
max(Diffs)
min(Diffs)



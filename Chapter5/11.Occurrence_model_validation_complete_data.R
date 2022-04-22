## Occurrence model: validation of results using complete data subset

library(dplyr)
library(lme4)
library(ggplot2)
library(StatisticalModels)
library(performance)
library(ggpubr)
library(phytools)
library(PVR)
library(viridis)
library(ggplot2)
library(ggpubr)
library(picante)


GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

# function for backtransforming from log-odds to probablities
logit2prob <- function(logit){
  return(exp(logit)/(1+exp(logit)))
}

# function for predicting effects
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
      
      # # drop coefs that couldn't be estimated
      print(setdiff(colnames(mm), names(coefs)))
      to_drop <- setdiff(colnames(mm), names(coefs))
      
      if(length(to_drop)!=0){
        mm <- as.data.frame(mm)
        mm <- mm[, -which(colnames(mm) %in% to_drop)]
        mm <- as.matrix(mm)
      }
      
      y <- mm %*% coefs
      
      # backtranforming
      #y <- logit2prob(y)
      
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



######################################################################################################################

## RESIDUAL EXTRACTION - FOR KNOWN RMR VALUES ONLY

BMR_data <- read.csv("../Results/2.BMR_data_to_impute.csv") %>% 
  filter(!is.na(Family))
colnames(BMR_data)[1] <- "Best_guess_binomial"

## first, get residual variation in RMR, for the known values only (no imputed values included).
BMR_data$log_BMR <- log(BMR_data$BMR_ml_O2_per_h)

# model for extraction of residuals
model1 <- lmer(log_BMR~log(Body_mass_g) + (1|Class/Order/Family), data=BMR_data)
Residual_BMR_log_log <- residuals(model1)

# insert residual values in dataframe
BMR_data$log_BMR[!is.na(BMR_data$log_BMR)] %>%  length()
length(Residual_BMR_log_log)
BMR_data$Residual_BMR_log_log[!is.na(BMR_data$log_BMR)] <- Residual_BMR_log_log


######################################################################################################################

## Load PREDICTS data anc check
## loading and checking PREDICTS data
PredictsTraits <- "../Results/PredictsTraits.rds" %>% 
  readRDS() 

levels(PredictsTraits$LandUse)
levels(PredictsTraits$Use_intensity)

## remove previous estimates of residual RMRs, etc
PredictsTraits_complete <- PredictsTraits %>% 
  dplyr::filter(!is.na(LandUse), !is.na(Use_intensity)) %>% 
  dplyr::select(-Residual_BMR_log_log, -log_BMR, -Body_mass_g, -Thermoregulation)

## join with new RMR data
BMR_data <- BMR_data %>% 
  dplyr::select(-Family, -Order, -Class)
PredictsTraits_complete <- left_join(PredictsTraits_complete, BMR_data, by="Best_guess_binomial")
colnames(PredictsTraits_complete)

## checking
any(is.na(PredictsTraits_complete$Body_mass_g))
any(is.na(PredictsTraits_complete$Trophic_level))

any(is.na(PredictsTraits_complete$log_BMR))
any(is.na(PredictsTraits_complete$Residual_BMR_log_log))

######################################################################################################################

# sample sizes (number of species for which RMR is known)
PredictsTraits_complete$Best_guess_binomial[!is.na(PredictsTraits_complete$Residual_BMR_log_log)] %>%  unique %>%  length() # 489 species
PredictsTraits_complete$Best_guess_binomial[!is.na(PredictsTraits_complete$Residual_BMR_log_log) & PredictsTraits_complete$Class=="Mammalia"] %>%  unique %>%  length() # 144 species
PredictsTraits_complete$Best_guess_binomial[!is.na(PredictsTraits_complete$Residual_BMR_log_log) & PredictsTraits_complete$Class=="Aves"] %>%  unique %>%  length() # 307 species
PredictsTraits_complete$Best_guess_binomial[!is.na(PredictsTraits_complete$Residual_BMR_log_log) & PredictsTraits_complete$Class=="Reptilia"] %>%  unique %>%  length() # 22 species
PredictsTraits_complete$Best_guess_binomial[!is.na(PredictsTraits_complete$Residual_BMR_log_log) & PredictsTraits_complete$Class=="Amphibia"] %>%  unique %>%  length() # 16 species

# distribution of residual RMR
PredictsTraits_complete$Residual_BMR_log_log %>%  hist()


######################################################################################################################
## Running occurrence model -- with secondary vegetation grouped

PredictsTraits_complete$Use_intensity %>%  levels()
PredictsTraits_complete$LandUse %>%  levels

## fitting the model with the three-way interactions

Start <- Sys.time()
print(Start)
Model <- lme4::glmer(Occurrence ~
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
                     data=PredictsTraits_complete,
                     family="binomial")

End <- Sys.time()
print(Start-End)

saveRDS(Model, "../Results/10.Occurrence_model_validation/Occurrence_model_complete_RMR.rds")
summary(Model)

##################################################################################

Model_complete <- readRDS( "../Results/10.Occurrence_model_validation/Occurrence_model_complete_RMR.rds")
Model_imputed <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")

Summary_complete <- as.data.frame(summary(Model_complete)$coefficients)
Summary_imputed <- as.data.frame(summary(Model_imputed)$coefficients)

## checking order of estimate coefficients
rownames(Summary_complete)==rownames(Summary_imputed)

## plotting
# Summary_imputed$Which <- "Imputed"
# Summary_complete$Which <- "Complete"
# Summaries <- rbind(Summary_imputed, Summary_complete)

# Summaries$CI_up <- Summaries$Estimate + 1.96*Summaries$`Std. Error`
# Summaries$CI_low <- Summaries$Estimate - 1.96*Summaries$`Std. Error`

## plotting estimates 
dev.off()

pdf(file = "../Results/Figures/10.Occurrence_model_validation/SI_Coefficients_Complete_against_Imputed.pdf", width=8, height=7, pointsize = 13)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,1), oma=c(2,1,2,1))
par(mfrow=c(2,1))

plot(Summary_complete$Estimate,
     Summary_imputed$Estimate,
     pch=19, main="Coefficient estimates", xlab="Estimate: imputed", ylab="Estimate: complete")
abline(a=0, b=1, lty="dashed", col="red")

plot(Summary_complete$`Std. Error`,
     Summary_imputed$`Std. Error`,
     pch=19, main="Coefficient estimates", xlab="Estimate: imputed", ylab="Estimate: complete")
abline(a=0, b=1, lty="dashed", col="red")


dev.off()

Diffs <- Summary_imputed$Estimate-Summary_complete$Estimate
hist(Diffs)
max(Diffs)
min(Diffs)



# ########################### PREDCITIONS FOR COMPLETE AGAINST PREDICTIONS FOR IMPUTED + COLLECTED









##################################################################################

# ## variance explained
# anova(Model)
# af <- anova(Model)
# afss <- af$"Sum Sq"
# afss <- cbind(af,PctExp=afss/sum(afss)*100)
# print(afss[order(afss$PctExp),])

# #######################################################################################################################
# # plotting slope of the relationship between occurrence probability and residual RMR - for each land-use intensity
# # since residual RMR doesn't interact with trophic level, the slopes are similar for different trohic groups
# 
# ModelSummary <- as.data.frame(summary(Model)$coefficients)
# 
# # subset for slopes
# ModelSummary_slopes <- ModelSummary[grepl("BMR", rownames(ModelSummary)),]
# 
# # slopes for reference level of use intensity 
# ModelSummary_slopes_Minimal <- ModelSummary_slopes[!grepl("Use_intensity", rownames(ModelSummary_slopes)),]
# ModelSummary_slopes_Minimal$Use_intensity <- "Minimal use"
# 
# # get slope estimates for minimally-used land uses
# ModelSummary_slopes_Minimal$Estimate[2:nrow(ModelSummary_slopes_Minimal)] <- 
#   ModelSummary_slopes_Minimal$Estimate[2:nrow(ModelSummary_slopes_Minimal)] +
#   ModelSummary_slopes_Minimal$Estimate[1]
# 
# # get slope estimates for lightly-used land uses
# ModelSummary_slopes_Light <- ModelSummary_slopes_Minimal
# ModelSummary_slopes_Light$Estimate <- ModelSummary_slopes_Light$Estimate + ModelSummary_slopes$Estimate[grepl("Light use", rownames(ModelSummary_slopes))]
# 
# ModelSummary_slopes_Intense <- ModelSummary_slopes_Minimal
# ModelSummary_slopes_Intense$Estimate <- ModelSummary_slopes_Intense$Estimate + ModelSummary_slopes$Estimate[grepl("Intense use", rownames(ModelSummary_slopes))]
# 
# ModelSummary_slopes_Light$Use_intensity <- "Light use"
# ModelSummary_slopes_Intense$Use_intensity <- "Intense use"
# 
# # binding slopes together
# Slopes <- rbind(ModelSummary_slopes_Minimal, ModelSummary_slopes_Light, ModelSummary_slopes_Intense)
# Slopes$Use_intensity <- factor(Slopes$Use_intensity, levels=c("Minimal use", "Light use", "Intense use"))
# 
# # add land uses
# FactorLU <- factor(c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban"), 
#                    levels = c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban"))
# Slopes$LandUse <- rep(FactorLU, 3)
# Slopes$LandUse
# 
# # add confidence intervals
# Slopes$Upper <- Slopes$Estimate + Slopes$`Std. Error`*1.96
# Slopes$Lower <- Slopes$Estimate - Slopes$`Std. Error`*1.96
# 
# Slopes$Lower[Slopes$LandUse=="Primary" & Slopes$Use_intensity=="Minimal use"] <- NA
# Slopes$Upper[Slopes$LandUse=="Primary" & Slopes$Use_intensity=="Minimal use"] <- NA
# 
# ## plotting
# scales::show_col(viridis(option = "plasma", n=8))
# cols <- viridis(option = "plasma", n=8)[c(2:5,7)]
# cbPalette <- c("#000000", cols)
# scales::show_col(cbPalette)
# 
# 
# # for adding horizontal line
# Slopes$h_line <- rep(Slopes$Estimate[Slopes$LandUse=="Primary"], each=6)
# 
# 
# Slopes_plot <-
#   ggplot(Slopes,
#          aes(LandUse, Estimate, ymin = Lower, ymax = Upper, col=LandUse)) +
#   geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   geom_rect(xmin=8.5, xmax=9.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   ylab("") + xlab("") +
#   geom_hline(aes(yintercept=h_line), col="darkgrey", linetype="dashed") +
#   geom_hline(yintercept=Slopes$Estimate[Slopes$LandUse=="Primary" & Slopes$Use_intensity=="Minimal use"],linetype="dashed") +
#   #geom_hline(yintercept=0,linetype="dotted", col="blue") +
#   geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
#   geom_point(size=2, position=position_dodge(width = 0.7)) +
#   GGPoptions +
#   ggtitle("") +
#   facet_grid(~Use_intensity, scales="free") +
#   theme(panel.spacing = unit(0, "lines")) +
#   theme( strip.text.x = element_text(size = 12, face = "bold"),
#          strip.text.y = element_text(size = 12, face = "bold")) +
#   scale_colour_manual(values=cbPalette, name="Use intensity") +
#   theme(legend.position = "right") +
#   scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
#   guides(colour=FALSE) +
#   theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
#   ylab("Slope estimate (+/- 95% CI)")
# 
# 
# #ggsave(Slopes_plot, filename="../Results/Figures/8.Occurrence_models/Best_model_Slopes.pdf", width=9, height = 4.5)
# 
# #######################################################################################################################
# 
# ## plotting continuous effetcs
# Newdata <- expand.grid(LandUse=levels(Model@frame$LandUse), 
#                        Use_intensity=levels(Model@frame$Use_intensity), 
#                        Trophic_level=levels(Model@frame$Trophic_level),  
#                        Residual_BMR_log_log=seq(from=min(Model@frame$Residual_BMR_log_log),
#                                                 to=max(Model@frame$Residual_BMR_log_log),
#                                                 length.out=100),  
#                        Occurrence=1)
# 
# 
# ## predictions
# Predictions <- Predict_effects(Newdata = Newdata, Model = Model, rescale = FALSE, Cont_or_Cat = "continuous", seMultiplier = 1.96, LU_n = 6)
# 
# scales::show_col(viridis(option = "plasma", n=8))
# cols <- viridis(option = "plasma", n=8)
# cbPalette_1 <- c("#000000", cols[c(4,7)])
# scales::show_col(cbPalette_1)
# cbPalette_2 <- cols[c(2,5,8)]
# scales::show_col(cbPalette_2)
# 
# 
# ## plotting for PV, SV, PF
# Predictions_PV_SV_PF <- subset(Predictions, LandUse %in% c("Primary vegetation", "Secondary vegetation", "Plantation forest"))
# ggplot(Predictions_PV_SV_PF, aes(Residual_BMR_log_log, Estimate, col=LandUse, fill=LandUse, ymin=Lower, ymax=Upper)) + 
#   geom_line() +
#   scale_colour_manual(values=cbPalette_1, name="Land use") +
#   scale_fill_manual(values=cbPalette_1, name="Land use") +
#   geom_ribbon(alpha=0.2, col=NA) + 
#   facet_grid(Use_intensity~Trophic_level) +
#   xlab("Residual RMR") + ylab("Occurrence probability") +
#   GGPoptions
# 
# Predictions_PA_CR_UR <- subset(Predictions, LandUse %in% c("Pasture", "Cropland", "Urban"))
# ggplot(Predictions_PA_CR_UR, aes(Residual_BMR_log_log, Estimate, col=LandUse, fill=LandUse, ymin=Lower, ymax=Upper)) + 
#   geom_line() + 
#   scale_fill_manual(values=cbPalette_2, name="Land use") +
#   geom_ribbon(alpha=0.3, col=NA) + 
#   facet_grid(Use_intensity~Trophic_level)+
#   xlab("Residual RMR") + ylab("Occurrence probability") +
#   GGPoptions
# 
# 
# ## predictions are very similar within trophic levels - so just plot for carnivores (reference level)
# Predictions_PV_SV_PF_Car <- subset(Predictions_PV_SV_PF, Trophic_level=="Carnivore")
# Predictions_PA_CR_UR_Car <- subset(Predictions_PA_CR_UR, Trophic_level=="Carnivore")
# 
# Predictions_PV_SV_PF_Car$LandUse <- factor(Predictions_PV_SV_PF_Car$LandUse,
#                                            levels = c("Primary vegetation",
#                                                       "Secondary vegetation",
#                                                       "Plantation forest"),
#                                            labels=c("Primary", "Secondary", "Plantation"))
# 
# Predictions_PA_CR_UR_Car$LandUse <- factor(Predictions_PA_CR_UR_Car$LandUse,
#                                            levels = c("Pasture",
#                                                       "Cropland",
#                                                       "Urban"),
#                                            labels=c("Pasture", "Cropland", "Urban"))
# 
# ## plotting for PV, SV, PF
# p1 <- 
#   ggplot(Predictions_PV_SV_PF_Car, aes(Residual_BMR_log_log, 
#                                        Estimate, col=LandUse,
#                                        fill=LandUse, 
#                                        ymin=Lower, ymax=Upper)) + 
#   geom_line(size=1) +
#   scale_colour_manual(values=cbPalette_1, name="Land use") +
#   scale_fill_manual(values=cbPalette_1, name="Land use") +
#   geom_ribbon(alpha=0.1, col=NA) + 
#   facet_grid(~Use_intensity) +
#   xlab("Residual RMR") + ylab("Occurrence probability") +
#   GGPoptions +
#   theme(panel.spacing = unit(0, "lines")) +
#   theme( strip.text.x = element_text(size = 12, face = "bold"),
#          strip.text.y = element_text(size = 12, face = "bold")) +
#   ggtitle("(A) Primary, secondary, Plantation")
# 
# p2 <- 
#   ggplot(Predictions_PA_CR_UR_Car, aes(Residual_BMR_log_log, Estimate, col=LandUse, fill=LandUse, ymin=Lower, ymax=Upper)) + 
#   geom_line(size=1) + 
#   scale_fill_manual(values=cbPalette_2, name="Land use") +
#   scale_colour_manual(values=cbPalette_2, name="Land use") +
#   geom_ribbon(alpha=0.1, col=NA) + 
#   facet_grid(~Use_intensity)+
#   xlab("Residual RMR") + ylab("Occurrence probability") +
#   GGPoptions +
#   theme(panel.spacing = unit(0, "lines")) +
#   theme( strip.text.x = element_text(size = 12, face = "bold"),
#          strip.text.y = element_text(size = 12, face = "bold")) +
#   ggtitle("(B) Pasture, Cropland, Urban")
# 
# ggarrange(p1,p2,  nrow=2)
# 
# ## trying a different facetting
# Predictions_sub <- subset(Predictions, Trophic_level=="Carnivore")
# Predictions_sub$Facet[Predictions_sub$LandUse %in% c("Primary vegetation", "Secondary vegetation", "Plantation forest")] <- "Primary, secondary, plantation"
# Predictions_sub$Facet[Predictions_sub$LandUse %in% c("Pasture", "Cropland", "Urban")] <- "Pasture, cropland, urban"
# Predictions_sub$Facet <- factor(Predictions_sub$Facet, levels = c("Primary, secondary, plantation", "Pasture, cropland, urban"))
# 
# Predictions_sub$LandUse <- factor(Predictions_sub$LandUse,
#                                   labels=c("Primary", "Secondary", "Plantation", 
#                                            "Pasture", "Cropland", "Urban"))
# 
# scales::show_col(viridis(option = "plasma", n=8))
# cols <- viridis(option = "plasma", n=8)[c(2:5,7)]
# cbPalette <- c("#000000", cols)
# scales::show_col(cbPalette)
# 
# p_continuous_effects <- 
#   ggplot(Predictions_sub, aes(Residual_BMR_log_log, Estimate, col=LandUse, fill=LandUse, ymin=Lower, ymax=Upper)) + 
#   geom_line(size=1) + 
#   geom_ribbon(alpha=0.1, col=NA) + 
#   facet_grid(Facet~Use_intensity)+
#   xlab("Residual RMR") + ylab("Occurrence probability") +
#   GGPoptions +
#   scale_colour_manual(values=cbPalette, name="Land use") +
#   scale_fill_manual(values=cbPalette, name="Land use") +
#   theme(panel.spacing = unit(0, "lines")) +
#   theme( strip.text.x = element_text(size = 12, face = "bold"),
#          strip.text.y = element_text(size = 12, face = "bold")) +
#   labs(col='Land use', fill="Land use") +
#   ggtitle("(a) Effect of residual RMR on occurrence probability") +
#   theme(strip.text.y = element_blank())
# # scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.9) +
# # scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.9)
# 
# 
# 
# #################### Arranging plots for manuscript
# 
# p_manuscript <- 
#   ggarrange(p_continuous_effects,
#             Slopes_plot +
#               ggtitle("(b) Slope of the relationship"), 
#             nrow=2,
#             common.legend = TRUE,
#             legend = "right",
#             heights = c(0.55, 0.45))
# 
# 
# ggsave(p_manuscript, 
#        filename = "../Results/Figures/10.Occurrence_model_validation/Figure_predictions_slopes.pdf", height=7, width=8)
# ggsave(p_manuscript, 
#        filename = "../Results/Figures/10.Occurrence_model_validation/Figure_predictions_slopes.png", height=7, width=8)
# 
# 
# # # plotting using Tim's function
# # PlotGLMERContinuous(model = SelectModel$model,
# #                     data = SelectModel$data,
# #                     effects = "Residual_BMR_log_log",
# #                     otherFactors = list(LandUse="Primary Vegetation", Use_intensity="Intense use", Trophic_level="Herbivore"),
# #                     xlab = "Residual RMR",
# #                     ylab = "Occurrence probability",
# #                     byFactor = "LandUse",
# #                     logLink = "b",
# #                     line.cols=c("#E6AB02","#D95F02","#7570B3","#66A61E","#1B9E77","#E7298A"))
# # 
# # PlotGLMERContinuous(model = SelectModel$model,
# #                     data = SelectModel$data,
# #                     effects = "Residual_BMR_log_log",
# #                     otherFactors = list(LandUse="Primary Vegetation", Use_intensity="Intense use", Trophic_level="Carnivore"),
# #                     xlab = "Residual RMR",
# #                     ylab = "Occurrence probability",
# #                     byFactor = "LandUse",
# #                     logLink = "b",
# #                     line.cols=c("#E6AB02","#D95F02","#7570B3","#66A61E","#1B9E77","#E7298A"))
# # 
# # 
# # legend(1.5,0.35,c("Primary","Secondary","Plantation","Pasture","Cropland","Urban"),
# #        col=c("#E6AB02","#D95F02","#7570B3","#66A61E","#1B9E77","#E7298A"),bty="n",lty=1)
# # 
# # 

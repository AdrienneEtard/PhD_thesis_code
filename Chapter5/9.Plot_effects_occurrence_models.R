## plotting effects for occurrence models 

library(StatisticalModels)
library(ggplot2)
library(viridis)
library(scales)
library(ggpubr)
library(performance)
library(DHARMa)


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
Predict_effects <- function(Newdata, Model, rescale, Cont_or_Cat, seMultiplier, LU_n, BT) {
  
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
      if(BT){
      y <- logit2prob(y)
      }
      
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


###########################################################################################################################################
Model <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")
Coefs <- summary(Model)$coefficients
Model@call

## checking residual distributions of the model
## diagnostic plots with DHarma
simulationOutput <- simulateResiduals(fittedModel = Model, plot = TRUE, asFactor=TRUE, n = 500)
testOutliers(simulationOutput, type="bootstrap")
plot(simulationOutput)
testDispersion(simulationOutput)

SampleSize <- Model@frame[!is.na(Model@frame$LandUse) &!is.na(Model@frame$Use_intensity),]
SampleSize$Best_guess_binomial %>%  unique %>%  length()
SampleSize$SS %>%  unique %>%  length()
SampleSize$SSBS %>%  unique %>%  length()


####################################################################################################################################
# plotting slope of the relationship between occurrence probability and residual RMR - for each land-use intensity and trophic level

## getting slopes by resampling in the models coefficients

Newdata <- expand.grid(LandUse=levels(Model@frame$LandUse), 
                       Use_intensity=levels(Model@frame$Use_intensity), 
                       Trophic_level=levels(Model@frame$Trophic_level),  
                       Residual_BMR_log_log=1,  
                       Occurrence=1)


## getting slope of the relationship between occurrence probability and residual RMR for each row of Newdata

# subsetting model's coefficients and variance-covariance matrix
Fixed_effects <- fixef(Model) 
Fixed_effects_slopes <- Fixed_effects[grepl("BMR", names(Fixed_effects))]

V_Cov <- vcov(Model)
# subset relevant rows
V_Cov_sub <- V_Cov[grepl("BMR", rownames(V_Cov)),] 
# subset relevant columns
Col_retain <- which(grepl("BMR", colnames(V_Cov_sub)))
length(Col_retain)==nrow(V_Cov_sub)
V_Cov_sub <- V_Cov_sub[, Col_retain]
nrow(V_Cov_sub)==ncol(V_Cov_sub)

## getting slopes and 95% CI by re sampling
preds <- sapply(X=1:5000, FUN=function(i){
  coefs <- mvrnorm(n=1, mu=Fixed_effects_slopes, Sigma=V_Cov_sub)
  mm <- model.matrix(terms(Model), Newdata)
  ## drop coefs that shouldn't be estimated (the ones that correspond to intercepts)
  #print(setdiff(colnames(mm), names(coefs)))
  to_drop <- setdiff(colnames(mm), names(coefs))
  if(length(to_drop)!=0){
    mm <- as.data.frame(mm)
    mm <- mm[, -which(colnames(mm) %in% to_drop)]
    mm <- as.matrix(mm)
  }
  ## get slope estimates
  y <- mm %*% coefs
  return(y)
})

preds_Slopes <- data.frame(Median=apply(X=preds, MARGIN=1, FUN=median),
                           Upper=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.975),
                           Lower=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.025))

preds_Slopes <- cbind(preds_Slopes, Newdata)

## plotting slopes
scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[c(2:5,7)]
cbPalette <- c("#000000", cols)
scales::show_col(cbPalette)

# setting levels for land use
preds_Slopes$LandUse <- factor(preds_Slopes$LandUse, 
                               labels = c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban"))

# for adding horizontal line
preds_Slopes$h_line <- rep(preds_Slopes$Median[preds_Slopes$LandUse=="Primary" &
                                                 preds_Slopes$Use_intensity=="Minimal use"], each=18)

# # remove CIs for Primary vegetation, minimal use
preds_Slopes$Lower[preds_Slopes$LandUse=="Primary"&
                     preds_Slopes$Use_intensity=="Minimal use"] <- NA

preds_Slopes$Upper[preds_Slopes$LandUse=="Primary"&
                     preds_Slopes$Use_intensity=="Minimal use"] <- NA

# adding in whether the slope differ within land use intensities for each land use
Signif_UI <- function(preds_Slopes, TL) {
  
  sub <- subset(preds_Slopes, Trophic_level==TL)
  
  # split Sub by land use
  SubSplit <- split(x=sub, f=sub$LandUse)
  
  # apply function to each element
  FunToApply <- function(X) {
    ValueMinimal <- X$Median[X$Use_intensity=="Minimal use"]
    X$Significance <- NA
    X$Significance[2] <- !dplyr::between(ValueMinimal, X$Lower[2], X$Upper[2])      
    X$Significance[3] <- !dplyr::between(ValueMinimal, X$Lower[3], X$Upper[3])      
    return(X)
  }
  
  SubSplitRes <- lapply(SubSplit, FunToApply)
  SubSplitRes <- data.table::rbindlist(SubSplitRes)
  
  return(SubSplitRes)
}

preds_SlopesC <- Signif_UI(preds_Slopes, "Carnivore")
preds_SlopesO <- Signif_UI(preds_Slopes, "Omnivore")
preds_SlopesH <- Signif_UI(preds_Slopes, "Herbivore")

preds_Slopes <- rbind(preds_SlopesC, preds_SlopesO, preds_SlopesH)
preds_Slopes$Significance[is.na(preds_Slopes$Significance)] <- TRUE

preds_Slopes$Significance <- factor(preds_Slopes$Significance, levels=c("TRUE", "FALSE"))

Slopes_plot <-
  ggplot(preds_Slopes,
         aes(LandUse, Median, ymin = Lower, ymax = Upper, group=Use_intensity, shape=Use_intensity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=8.5, xmax=9.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  ylab("") + xlab("") +
  geom_hline(aes(yintercept=h_line), col="black") +
  geom_hline(yintercept = 0, col="#808080", linetype="dashed") +
  #geom_hline(yintercept=0, col="darkgrey") +
  geom_hline(yintercept=preds_Slopes$Estimate[preds_Slopes$LandUse=="Primary" & preds_Slopes$Use_intensity=="Minimal use"],linetype="dashed") +
  #geom_hline(yintercept=0,linetype="dotted", col="blue") +
  geom_errorbar(width=.2, size=0.6, position=position_dodge(width = 0.7), stat="identity") + #aes(linetype=Significance)
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=c("black", "#808080", "#A0A0A0"), name="Use intensity")+
                     # name="Within land uses, is slope \nsignificantly different from \nminimal use-intensity level") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  #guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Slope estimate (+/- 95% CI)") +
  theme(legend.position = "bottom") 
  #scale_linetype_manual(values=c(1,8), name="Within land uses, is slope \nsignificantly different from \nminimal use-intensity level")
  
Slopes_plot
ggsave(Slopes_plot, filename="../Results/Figures/8.Occurrence_models/Model_3_ways_Slopes.pdf", width=9, height = 4)

## checking which relationships are significant for minimal uses
preds2 <-  data.frame(Median=apply(X=preds, MARGIN=1, FUN=median),
                      Upper=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.975),
                      Lower=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.025))
preds2 <- cbind(preds2, Newdata)
preds2 %>%  
  dplyr::filter(Use_intensity=="Minimal use", LandUse=="Primary vegetation") %>% 
  dplyr::filter(Upper<0|Lower>0)


#######################################################################################################################

# OLDER CODE ## reference level is carnivores, minimal use
# 
# ModelSummary <- as.data.frame(summary(Model)$coefficients)
# write.csv(ModelSummary, "../Results/8.Occurrence_models/Occurrence_model_3_ways_summary.csv", row.names = TRUE)
# 
# # subset for slopes
# ModelSummary_slopes <- ModelSummary[grepl("BMR", rownames(ModelSummary)),]
# 
# # slopes for reference level of use intensity (minimal use), for carnivores
# ModelSummary_slopes_Minimal <- ModelSummary_slopes[!grepl("Use_intensity", rownames(ModelSummary_slopes)),]
# ModelSummary_slopes_Minimal$Use_intensity <- "Minimal use"
# 
# ModelSummary_slopes_Minimal_carnivores <- ModelSummary_slopes_Minimal[!grepl("Trophic_level", rownames(ModelSummary_slopes_Minimal)),]
# # get slope estimates for minimally-used land uses for Carnivores
# ModelSummary_slopes_Minimal_carnivores$Estimate[2:nrow(ModelSummary_slopes_Minimal_carnivores)] <- 
#   ModelSummary_slopes_Minimal_carnivores$Estimate[2:nrow(ModelSummary_slopes_Minimal_carnivores)] +
#   ModelSummary_slopes_Minimal_carnivores$Estimate[1]
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


#######################################################################################################################

rm(Newdata)

## plotting continuous effects
Newdata <- expand.grid(LandUse=levels(Model@frame$LandUse), 
                       Use_intensity=levels(Model@frame$Use_intensity), 
                       Trophic_level=levels(Model@frame$Trophic_level),  
                       Residual_BMR_log_log=seq(from=min(Model@frame$Residual_BMR_log_log),
                                                to=max(Model@frame$Residual_BMR_log_log),
                                                length.out=100),  
                       Occurrence=1)


## predictions
Predictions <- Predict_effects(Newdata = Newdata,
                               Model = Model,
                               rescale = FALSE,
                               Cont_or_Cat = "continuous",
                               seMultiplier = 1.96,
                               LU_n = 6)

# setting levels for land use
Predictions$LandUse <- factor(Predictions$LandUse, 
                              labels = c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban"))

# setting 3 colors
scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)
cbPalette_1 <- c("#000000", cols[c(4,7)])
scales::show_col(cbPalette_1)

# # plot for minimal uses for carnivores and omnivores, and for minimal and light uses for herbivores (in Primary, secondary, plantation)
# Predictions_MU <- subset(Predictions, Use_intensity %in% c("Minimal use"))
# Predictions_LI <- subset(Predictions, Use_intensity =="Light use" &
#                            LandUse %in% c("Primary", "Secondary", "Plantation") &
#                            Trophic_level=="Herbivore")
# 
# PredsToplot <- rbind(Predictions_MU, Predictions_LI)
# 
# PredsToplot$Significance_from_PV <- NA
# PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Carnivore" & PredsToplot$LandUse=="Urban"] <- TRUE
# PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Omnivore" & PredsToplot$LandUse %in% c("Secondary", "Cropland", "Urban")] <- TRUE
# PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Herbivore" & PredsToplot$LandUse %in% c("Cropland", "Urban")] <- TRUE
# PredsToplot$Significance_from_PV[is.na(PredsToplot$Significance_from_PV)] <- FALSE
# PredsToplot$Significance_from_PV[PredsToplot$LandUse=="Primary"] <- TRUE

#PredsToplot$Significance_from_PV <- factor(PredsToplot$Significance_from_PV, levels=c("TRUE", "FALSE"))

Predictions_minimal <- subset(Predictions, Use_intensity=="Minimal use")

## add significance of the relationship for the minimal uses 
# herbivore only for PV, cropland and urban 
# urban carnivore
#omn: secondary, cropland, urban

Predictions_minimal$Significant <- FALSE
Predictions_minimal$Significant[Predictions_minimal$LandUse %in% c("Primary", "Cropland", "Urban") &
                                  Predictions_minimal$Trophic_level=="Herbivore"] <- TRUE
Predictions_minimal$Significant[Predictions_minimal$LandUse %in% c("Urban") &
                                  Predictions_minimal$Trophic_level=="Carnivore"] <- TRUE
Predictions_minimal$Significant[Predictions_minimal$LandUse %in% c("Secondary","Cropland", "Urban") &
                                  Predictions_minimal$Trophic_level=="Omnivore"] <- TRUE

p_continuous_effects <- 
ggplot(Predictions_minimal, aes(Residual_BMR_log_log, Estimate, ymin=Lower, ymax=Upper))+#, linetype=Significance_from_PV)) + 
  geom_line(aes(linetype=Significant), size=0.8) +
  scale_linetype_manual(values=c("dashed", "solid"))+
  scale_fill_manual(values=c("black", "black"), name="Use intensity") +
  scale_colour_manual(values=c("black", "darkgrey"), name="Use intensity") +
  geom_ribbon(alpha=0.05, col=NA) + 
  facet_grid(LandUse~Trophic_level)+
  xlab("Residual RMR") + ylab("Occurrence probability") +
  GGPoptions +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.25)) +
  theme(legend.position = "none")

p_continuous_effects

#################### Arranging plots for manuscript

Slopes_plot2 <- Slopes_plot +
  theme(legend.margin=margin(0,0,45,0),
legend.box.margin=margin(-5,-5,-5,-5))

p_manuscript <- 
  ggarrange(  Slopes_plot2 +
                ggtitle("(a) Slope of the relationships between \noccurrence probability (log-odds) and residual RMR"), 
              p_continuous_effects +
              ggtitle("(b) Predicted effects of residual RMR on \noccurrence probability (minimal use intensity)"),
          
            nrow=2,
            #common.legend = TRUE,
            heights = c(0.40, 0.55))
p_manuscript

ggsave(p_manuscript, 
       filename = "../Results/Figures/8.Occurrence_models/Figure_predictions_slopes.pdf", height=11, width=9)

ggsave(p_manuscript, 
       filename = "../Results/Figures/8.Occurrence_models/Figure_predictions_slopes.png", height=11, width=9)


ggsave(p_manuscript, 
       filename = "C:/Users/adrie/OneDrive/Desktop/Thesis/figures/Chapter5/Figure4.pdf", height=12, width=7)





################################################# Discretizing effects of residual RMR for effects on occurrence probability - rescaled compared to minimal PV

ResRMR <- unique(Model@frame[, c("Best_guess_binomial", "Residual_BMR_log_log")])
Low <- quantile(ResRMR$Residual_BMR_log_log)["25%"]
Median <- quantile(ResRMR$Residual_BMR_log_log)["50%"]
High <- quantile(ResRMR$Residual_BMR_log_log)["75%"]

Newdata <- expand.grid(LandUse=levels(Model@frame$LandUse), 
                       Residual_BMR_log_log=c(Low, Median, High),  
                       Use_intensity=levels(Model@frame$Use_intensity), 
                       Trophic_level=levels(Model@frame$Trophic_level),  
                       Occurrence=1)

Effects_discrete <- Predict_effects(Newdata = Newdata,
                                    Model = Model,
                                    rescale = TRUE, 
                                    Cont_or_Cat = "categorical", 
                                    seMultiplier = 1.96, 
                                    LU_n = 18, 
                                    BT=TRUE)


Effects_discrete$Median <- Effects_discrete$Median - 100
Effects_discrete$Lower <- Effects_discrete$Lower - 100
Effects_discrete$Upper <- Effects_discrete$Upper - 100

# subset for minimal uses
Effects_discreteSub <- Effects_discrete[Effects_discrete$Use_intensity=="Minimal use" & 
                                          Effects_discrete$LandUse %in% c("Primary vegetation", "Cropland", "Urban"),]
Effects_discreteSub$Residual_BMR_log_log <- factor(Effects_discreteSub$Residual_BMR_log_log, labels=c("Low (25% quantile)",
                                                                                                      "Median", 
                                                                                                      "High (75% quantile)"))


Effects_discreteSub$Signif <- TRUE
Effects_discreteSub$Signif[Effects_discreteSub$Trophic_level=="Carnivore" & Effects_discreteSub$LandUse=="Cropland"] <- FALSE
Effects_discreteSub$Signif[Effects_discreteSub$Trophic_level=="Carnivore" & Effects_discreteSub$LandUse=="Primary vegetation"] <- FALSE
Effects_discreteSub$Signif[Effects_discreteSub$Trophic_level=="Omnivore" & Effects_discreteSub$LandUse=="Primary vegetation"] <- FALSE


pOc <- 
ggplot(Effects_discreteSub,
       aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Residual_BMR_log_log, group=Residual_BMR_log_log, col=Signif)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(labels=c("Primary", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_shape_manual(values=c(16,3,4), name="Residual RMR") +
  theme(legend.position = "right") +
  scale_colour_manual(values=c("grey48", "black"), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Occurrence probability\n(% difference from primary)")
pOc

ggsave(pOc, 
       filename = "../Results/Figures/8.Occurrence_models/Figure_Occurrence_probability.pdf", height=3.5, width=8)

ggsave(pOc, 
       filename = "../Results/Figures/8.Occurrence_models/Figure_Occurrence_probability.png", height=3.5, width=8)


ggsave(pOc, 
       filename = "C:/Users/adrie/OneDrive/Desktop/Thesis/figures/Chapter5/Figure5.pdf", height=3.5, width=8)




# ## plotting for PV, SV, P
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
#                                      Estimate, col=LandUse,
#                                      fill=LandUse, 
#                                      ymin=Lower, ymax=Upper)) + 
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
#   # scale_color_viridis(discrete = TRUE, option = "viridis", end = 0.9) +
#   # scale_fill_viridis(discrete = TRUE, option = "viridis", end = 0.9)
# 



# # plotting using Tim's function
# PlotGLMERContinuous(model = SelectModel$model,
#                     data = SelectModel$data,
#                     effects = "Residual_BMR_log_log",
#                     otherFactors = list(LandUse="Primary Vegetation", Use_intensity="Intense use", Trophic_level="Herbivore"),
#                     xlab = "Residual RMR",
#                     ylab = "Occurrence probability",
#                     byFactor = "LandUse",
#                     logLink = "b",
#                     line.cols=c("#E6AB02","#D95F02","#7570B3","#66A61E","#1B9E77","#E7298A"))
# 
# PlotGLMERContinuous(model = SelectModel$model,
#                     data = SelectModel$data,
#                     effects = "Residual_BMR_log_log",
#                     otherFactors = list(LandUse="Primary Vegetation", Use_intensity="Intense use", Trophic_level="Carnivore"),
#                     xlab = "Residual RMR",
#                     ylab = "Occurrence probability",
#                     byFactor = "LandUse",
#                     logLink = "b",
#                     line.cols=c("#E6AB02","#D95F02","#7570B3","#66A61E","#1B9E77","#E7298A"))
# 
# 
# legend(1.5,0.35,c("Primary","Secondary","Plantation","Pasture","Cropland","Urban"),
#        col=c("#E6AB02","#D95F02","#7570B3","#66A61E","#1B9E77","#E7298A"),bty="n",lty=1)
# 
# 



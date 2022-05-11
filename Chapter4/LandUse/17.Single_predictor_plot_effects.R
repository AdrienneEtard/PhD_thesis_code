## effects of categorical traits -- single predictor models

library(dplyr)
library(ggplot2)
library(ggarrange)
library(viridis)
library(StatisticalModels)

## functions to plot effects

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

# function for predicting effects -- for categorical predictors
Predict_effects <- function(Newdata, Model, rescale, Cont_or_Cat, seMultiplier, LU_n, BT) {
  
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
      #print(setdiff(colnames(mm), names(coefs)))
      to_drop <- setdiff(colnames(mm), names(coefs))
      
      if(length(to_drop)!=0){
        mm <- as.data.frame(mm)
        mm <- mm[, -which(colnames(mm) %in% to_drop)]
        mm <- as.matrix(mm)
      }
      
      y <- mm %*% coefs
      
      # backtransforming
      if(BT){
        y <- logit2prob(y)
      }
      
      # rescaling
      if(rescale){
        
        # initialization
        seq <- 1:LU_n
        y[seq] <- y[seq]/y[seq[1]]*100
        
        
        if(LU_n!=nrow(Newdata)){
        # for loop to rescale all values
        N <- (nrow(Newdata)/LU_n-1)
        for(i in 1:N){
          seq <- seq + LU_n
          y[seq] <- y[seq]/y[seq[1]]*100
          }
        }
        
        if(LU_n==nrow(Newdata)){
          y <- y/y[1]*100
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
  
}


##########################################################################################################################################
##### PREDICTING FOR SPECIALISATION
Spe_A <- readRDS( "../Results/Single_predictor_models/Amphibians_specialisation.rds")
Spe_B <- readRDS( "../Results/Single_predictor_models/Birds_specialisation.rds")
Spe_M <- readRDS( "../Results/Single_predictor_models/Mammals_specialisation.rds")
Spe_R <- readRDS( "../Results/Single_predictor_models/Reptiles_specialisation.rds")

Newdata <- expand.grid(LandUseGrouped=levels(Spe_A@frame$LandUseGrouped), 
                       Use_intensity=levels(Spe_A@frame$Use_intensity), 
                       Specialisation =levels(Spe_A@frame$Specialisation),  
                       Occurrence=1)

Effects_A <- Predict_effects(Newdata, Spe_A, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_A$Class <- "Amphibians"
Effects_B <- Predict_effects(Newdata, Spe_B, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_B$Class <- "Birds"
Effects_M <- Predict_effects(Newdata, Spe_M, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_M$Class <- "Mammals"
Effects_R <- Predict_effects(Newdata, Spe_R, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_R$Class <- "Reptiles"
Effects_R$Median[Effects_R$LandUseGrouped=="Urban"] <- NA
Effects_R$Lower[Effects_R$LandUseGrouped=="Urban"] <- NA
Effects_R$Upper[Effects_R$LandUseGrouped=="Urban"] <- NA

Effects_Spe <- rbind(Effects_A, Effects_B, Effects_M, Effects_R)
rm(Effects_A, Effects_B, Effects_M, Effects_R)
Effects_Spe$Median <- Effects_Spe$Median - 100
Effects_Spe$Lower <- Effects_Spe$Lower - 100
Effects_Spe$Upper <- Effects_Spe$Upper - 100

## plotting

Effects_Spe$Class <- factor(Effects_Spe$Class, labels = c("(a) Amphibians", "(b) Birds", "(c) Mammals", "(d) Reptiles"))
Effects_Spe$Specialisation <- factor(Effects_Spe$Specialisation, labels=c("Artificial habitats user", "Natural habitats specialist"))

Plot_Spe <- 
ggplot(Effects_Spe,
       aes(LandUseGrouped, Median, ymin = Lower, ymax = Upper, group=interaction(Use_intensity,Specialisation), shape=Use_intensity, col=Specialisation)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  xlab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary", 
                            "Secondary",
                            "Plantation", 
                            "Agricultural",
                            "Urban")) +
  GGPoptions +
  ylab("Occurrence probability, % difference from primary (+/- 95% confidence interval)") + 
  facet_wrap(~Class, scales = "free_y") +
  scale_colour_manual(values = c("black", "#D6604D"), name="Use of artificial habitats") +
  theme(panel.spacing.y =  unit(0, "lines")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.05, hjust=1)) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 12, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")))
  
ggsave(Plot_Spe, filename="../Results/Single_predictor_models/Figures/Validations/Specialisation.pdf", width = 10, height=7)
  
##########################################################################################################################################
##### PREDICTING FOR DIEL ACTIVITY
DA_A <- readRDS( "../Results/Single_predictor_models/Amphibians_diel_activity.rds")
DA_B <- readRDS( "../Results/Single_predictor_models/Birds_diel_activity.rds")
DA_M <- readRDS( "../Results/Single_predictor_models/Mammals_diel_activity.rds")
DA_R <- readRDS( "../Results/Single_predictor_models/Reptiles_diel_activity.rds")

Urban_Mammals <- DA_M@frame %>%  
  filter(LandUseGrouped=="Urban")
Urban_Mammals <- unique(Urban_Mammals[,c("Best_guess_binomial", "Diel_activity")])

Newdata <- expand.grid(LandUseGrouped=levels(DA_A@frame$LandUseGrouped), 
                       Use_intensity=levels(DA_A@frame$Use_intensity), 
                       Diel_activity =levels(DA_A@frame$Diel_activity),  
                       Occurrence=1)

Effects_A <- Predict_effects(Newdata, DA_A, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_A$Class <- "Amphibians"
Effects_B <- Predict_effects(Newdata, DA_B, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_B$Class <- "Birds"
Effects_M <- Predict_effects(Newdata, DA_M, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_M$Class <- "Mammals"
Effects_R <- Predict_effects(Newdata, DA_R, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_R$Class <- "Reptiles"
Effects_R$Median[Effects_R$LandUseGrouped=="Urban"] <- NA
Effects_R$Lower[Effects_R$LandUseGrouped=="Urban"] <- NA
Effects_R$Upper[Effects_R$LandUseGrouped=="Urban"] <- NA

Effects_DA <- rbind(Effects_A, Effects_B, Effects_M, Effects_R)
rm(Effects_A, Effects_B, Effects_M, Effects_R)
Effects_DA$Median <- Effects_DA$Median - 100
Effects_DA$Lower <- Effects_DA$Lower - 100
Effects_DA$Upper <- Effects_DA$Upper - 100

# ## removing effects of mammals in urban for non-nocturnal species (null effects, error bars too large)
# Effects_DA$Median <- Effects_DA$Median - 100
# Effects_DA$Lower <- Effects_DA$Lower - 100
# Effects_DA$Upper <- Effects_DA$Upper - 100


## plotting

Effects_DA$Class <- factor(Effects_DA$Class, labels = c("(a) Amphibians", "(b) Birds", "(c) Mammals", "(d) Reptiles"))
Effects_DA$Diel_activity <- factor(Effects_DA$Diel_activity, labels=c("Nocturnal", "Non-nocturnal"))

Plot_DA <-
ggplot(Effects_DA,
       aes(LandUseGrouped, Median, ymin = Lower, ymax = Upper, group=interaction(Use_intensity,Diel_activity), shape=Use_intensity, col=Diel_activity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  xlab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary", 
                            "Secondary",
                            "Plantation", 
                            "Agricultural",
                            "Urban")) +
  GGPoptions +
  ylab("Occurrence probability, % difference from primary (+/- 95% confidence interval)") + 
  facet_wrap(~Class, scales = "free_y") +
  scale_colour_manual(values = c("black", "#2166AC"), name="Diel activity") +
  theme(panel.spacing.y =  unit(0, "lines")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.05, hjust=1)) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
theme( strip.text.x = element_text(size = 12, face = "bold"),
       strip.text.y = element_text(size = 12, face = "bold")) +
  theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 12, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")))

ggsave(Plot_DA, filename="../Results/Single_predictor_models/Figures/Validations/Diel_activity.pdf", width = 10, height=7)


##########################################################################################################################################
##### PREDICTING FOR PRIMARY DIET
Diet_B <- readRDS( "../Results/Single_predictor_models/Birds_diet.rds")
Diet_M <- readRDS( "../Results/Single_predictor_models/Mammals_diet.rds")

Urban_Mammals <- Diet_M@frame %>%  
  filter(LandUseGrouped=="Urban")
Urban_Mammals <- unique(Urban_Mammals[,c("Best_guess_binomial", "Primary_diet")])
table(Urban_Mammals$Primary_diet)

Diet_A <- readRDS("../Results/Single_predictor_models/Amphibians_diet.rds")
Urban_A <- Diet_A@frame %>%  
  filter(LandUseGrouped=="Urban")
Urban_A <- unique(Urban_A[,c("Best_guess_binomial", "Primary_diet")])
table(Urban_A$Primary_diet)

Diet_R <- readRDS("../Results/Single_predictor_models/Reptiles_diet.rds")

## birds and mammals
Newdata <- expand.grid(LandUseGrouped=levels(Diet_B@frame$LandUseGrouped), 
                       Use_intensity=levels(Diet_B@frame$Use_intensity), 
                       Primary_diet =levels(Diet_B@frame$Primary_diet),  
                       Occurrence=1)

Effects_B <- Predict_effects(Newdata, Diet_B, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_B$Class <- "Birds"
Effects_M <- Predict_effects(Newdata, Diet_M, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_M$Class <- "Mammals"

Effects_M$Median[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="PL|SE"] <- NA
Effects_M$Lower[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="PL|SE"] <- NA
Effects_M$Upper[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="PL|SE"] <- NA

Effects_M$Median[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="VE"] <- NA
Effects_M$Lower[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="VE"] <- NA
Effects_M$Upper[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="VE"] <- NA

Effects_M$Median[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="FR|NE"] <- NA
Effects_M$Lower[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="FR|NE"] <- NA
Effects_M$Upper[Effects_M$LandUseGrouped=="Urban" & Effects_M$Primary_diet=="FR|NE"] <- NA



## amphibians
Newdata <- expand.grid(LandUseGrouped=levels(Diet_A@frame$LandUseGrouped), 
                       Use_intensity=levels(Diet_A@frame$Use_intensity), 
                       Primary_diet =levels(Diet_A@frame$Primary_diet),  
                       Occurrence=1)

Effects_A <- Predict_effects(Newdata, Diet_A, rescale = TRUE, "categorical", NULL, LU_n = 15, BT=TRUE)
Effects_A$Class <- "Amphibians"

## reptiles
Newdata <- expand.grid(LandUseGrouped=levels(Diet_R@frame$LandUseGrouped), 
                       Use_intensity=levels(Diet_R@frame$Use_intensity), 
                       Primary_diet =levels(Diet_R@frame$Primary_diet),  
                       Occurrence=1)

Effects_R <- Predict_effects(Newdata, Diet_R, rescale = TRUE, "categorical", NULL, LU_n = 12, BT=TRUE)
Effects_R$Class <- "Reptiles"


## all classes together

Effects_Diet <- rbind(Effects_A, Effects_B, Effects_M, Effects_R)
#rm(Effects_B, Effects_M)
Effects_Diet$Median <- Effects_Diet$Median - 100
Effects_Diet$Lower <- Effects_Diet$Lower - 100
Effects_Diet$Upper <- Effects_Diet$Upper - 100

## creating empty plots for amphibians (FR/NE, PL/SE, VE)
Effects_Sub <- subset(Effects_B, Primary_diet %in% c("FR|NE", "PL|SE", "VE"))
Effects_Sub$Class <- "Amphibians"
Effects_Sub$Median <- NA
Effects_Sub$Lower <- NA
Effects_Sub$Upper <- NA

## creating empty plots for reptiles (FR/NE, PL/SE)
Effects_Sub2 <- subset(Effects_B, Primary_diet %in% c("FR|NE", "PL|SE"))
Effects_Sub2$Class <- "Reptiles"
Effects_Sub2$Median <- NA
Effects_Sub2$Lower <- NA
Effects_Sub2$Upper <- NA

Effects_Diet <- rbind(Effects_Diet, Effects_Sub, Effects_Sub2)


## plotting

# Effects_Diet$Primary_diet <- factor(Effects_Diet$Primary_diet, 
#                                     levels=c("FR|NE", "PL|SE", "IN", "VE", "OM"),
#                                     labels=c("(a) Fruit/nectar",  "(b) Plants/seeds", "(c) Invertebrates", "(d) Vertebrates","(e) Omnivore"))

Effects_Diet$Primary_diet <- factor(Effects_Diet$Primary_diet,
                                    levels=c("FR|NE", "PL|SE", "IN", "VE", "OM"),
                                    labels=c("Fruit/nectar",  "Plants/seeds", "Invertebrates", "Vertebrates","Omnivore"))


Effects_Diet$Which <- paste0(Effects_Diet$Class, ", ",tolower(Effects_Diet$Primary_diet))
Effects_Diet$Which2 <- paste0(Effects_Diet$Primary_diet, ", ",tolower(Effects_Diet$Class))

P_diet <- 
ggplot(Effects_Diet,
       aes(LandUseGrouped, Median, ymin = Lower, ymax = Upper, group=interaction(Use_intensity,Primary_diet), 
           shape=Use_intensity, col=Primary_diet)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  xlab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=1, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary", 
                            "Secondary",
                            "Plantation", 
                            "Agricultural",
                            "Urban")) +
  GGPoptions +
  ylab("Occurrence probability, % difference from primary (+/- 95% CI)") + 
  #facet_grid(Primary_diet~Class, scales="free_y") +
  facet_wrap(~Which2, scales="free_y", ncol=4) +
  theme(panel.spacing =  unit(0, "lines")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.05, hjust=1)) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 11, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
  theme(legend.position = "right") +
  viridis::scale_color_viridis(discrete = TRUE,
                               option = "B", 
                               end = 0.7, 
                               name="Primary diet",
                               labels=c("Fruit/nectar", 
                                        "Plants/seeds", 
                                        "Invertebrates", 
                                        "Vertebrates", 
                                        "Omnivore"))


ggsave(P_diet, filename="../Results/Single_predictor_models/Figures/Diet.pdf", width = 11, height=10)

##########################################################################################################################################
##########################################################################################################################################

## DISCRETIZING CONTINUOUS PREDICTORS AND PLOTTING EFFECTS

# Traits <- c("log10_Range_area", 
#             "log10_Body_mass_g", 
#             "log10_Litter_size", 
#             "log10_Lifespan_proxy", 
#             "sqrt_Habitat_breadth_IUCN", 
#             "sqrt_Diet_breadth")

GetNewdata <- function(Model, Trait, Class){

  traitvals <- unique(Model@frame[, c("Best_guess_binomial", Trait)])
  traitvals <- traitvals[,Trait]
  Low <- quantile(traitvals)["25%"]
  Median <- quantile(traitvals)["50%"]
  High <- quantile(traitvals)["75%"]
  
  Newdata <- expand.grid(LandUseGrouped=levels(Model@frame$LandUseGrouped), 
                         Use_intensity=levels(Model@frame$Use_intensity), 
                         Trait=c(Low, Median, High),  
                         Occurrence=1)
  
  colnames(Newdata)[3] <- Trait
  Newdata$Class <- Class
  return(Newdata)  
}

Function_Plot_effects_continuous_discretised <- function(Trait){
  
  PathA <- paste0("../Results/Single_predictor_models/","Amphibia_", Trait, ".rds")
  PathB <- paste0("../Results/Single_predictor_models/","Aves_", Trait, ".rds")
  PathM <- paste0("../Results/Single_predictor_models/","Mammalia_", Trait, ".rds")
  PathR <- paste0("../Results/Single_predictor_models/","Reptilia_", Trait, ".rds")
  
  Model_A <- readRDS(PathA)
  Model_B <- readRDS(PathB)
  Model_M <- readRDS(PathM)
  Model_R <- readRDS(PathR)
  
  print(levels(Model_R@frame$LandUseGrouped))
  
  Newdata_A <- GetNewdata(Model_A, Trait = Trait, "(a) Amphibians")
  Newdata_B <- GetNewdata(Model_B, Trait = Trait, "(b) Birds")
  Newdata_M <- GetNewdata(Model_M, Trait = Trait, "(c) Mammals")
  Newdata_R <- GetNewdata(Model_R, Trait = Trait, "(d) Reptiles")
  
  Effects_A <- Predict_effects(Newdata_A, Model_A, rescale = TRUE, "categorical", NULL, LU_n = 15*3, BT=TRUE)
  Effects_B <- Predict_effects(Newdata_B, Model_B, rescale = TRUE, "categorical", NULL, LU_n = 15*3, BT=TRUE)
  Effects_M <- Predict_effects(Newdata_M, Model_M, rescale = TRUE, "categorical", NULL, LU_n = 15*3, BT=TRUE)
  Effects_R <- Predict_effects(Newdata_R, Model_R, rescale = TRUE, "categorical", NULL, LU_n = 12*3, BT=TRUE)
  
  SetFac <- function(Effects, Trait) {
    Effects[, Trait] <- factor(Effects[, Trait], labels=c("Low (25% quantile)", "Median", "High (75% quantile)"))
    return(Effects)
  }
  Effects_A <- SetFac(Effects_A, Trait)
  Effects_B <- SetFac(Effects_B, Trait)
  Effects_M <- SetFac(Effects_M, Trait)
  Effects_R <- SetFac(Effects_R, Trait)
  Effects <- rbind(Effects_A, Effects_B, Effects_M, Effects_R)
  rm(Effects_A, Effects_B, Effects_M, Effects_R)
  Effects$Median <- Effects$Median - 100
  Effects$Lower <- Effects$Lower - 100
  Effects$Upper <- Effects$Upper - 100
  
  colnames(Effects)[colnames(Effects)==Trait] <- "The_trait"
  
  if(Trait=="log10_Body_mass_g"){ 
    Name <- "Body mass"}
  if(Trait=="log10_Range_area"){ 
    Name <- "Range area"}
  if(Trait=="log10_Litter_size"){ 
    Name <- "Litter/clutch size"}
  if(Trait=="log10_Lifespan_proxy"){ 
    Name <- "Lifespan proxy"}
  if(Trait=="sqrt_Habitat_breadth_IUCN"){ 
    Name <- "Habitat breadth"}
  if(Trait=="sqrt_Diet_breadth"){ 
    Name <- "Diet breadth"}
  
  ggplot(Effects,
         aes(LandUseGrouped, Median, ymin = Lower, ymax = Upper, group=interaction(Use_intensity,The_trait), shape=The_trait)) +
    geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
    geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
    geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
    geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
    xlab("") +
    geom_hline(yintercept = 0, col="black", linetype="dashed") +
    geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
    geom_point(size=1, position=position_dodge(width = 0.6)) +
    scale_x_discrete(labels=c("Primary", 
                              "Secondary",
                              "Plantation", 
                              "Agricultural",
                              "Urban")) +
    GGPoptions +
    ylab("Occurrence probability, % difference from primary (+/- 95% confidence interval)") + 
    facet_grid(Use_intensity~Class, scales="free_y") +
    theme(panel.spacing =  unit(0, "lines")) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1.05, hjust=1)) +
    scale_shape_manual(values=c(19, 17, 8), name=Name) +
    theme( strip.text.x = element_text(size = 12, face = "bold"),
           strip.text.y = element_text(size = 12, face = "bold")) +
    theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
          axis.text.x = element_text(size = 11, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
    theme(legend.position = "top")
  
  
}

## FOR BODY MASS
Function_Plot_effects_continuous_discretised("log10_Body_mass_g")

## FOR LITTER/CLUTCH SIZE
Function_Plot_effects_continuous_discretised("log10_Litter_size")

## FOR LIFESPAN
Function_Plot_effects_continuous_discretised("log10_Lifespan_proxy")

## FOR HABITAT BREADTH
Function_Plot_effects_continuous_discretised("sqrt_Habitat_breadth_IUCN")

## FOR DIET BREADTH
Function_Plot_effects_continuous_discretised("sqrt_Diet_breadth")

## FOR RANGE AREA
Function_Plot_effects_continuous_discretised("log10_Range_area")

##

# summary(readRDS("../Results/Single_predictor_models/Mammalia_log10_Range_area.rds"))
# summary(readRDS("../Results/Single_predictor_models/Aves_log10_Range_area.rds"))
# summary(readRDS("../Results/Single_predictor_models/Amphibia_log10_Range_area.rds"))
# summary(readRDS("../Results/Single_predictor_models/Mammalia_log10_Range_area.rds"))
# summary(readRDS("../Results/Single_predictor_models/Mammalia_sqrt_Habitat_breadth_IUCN.rds"))


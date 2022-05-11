library(ggplot2)
library(dplyr)
library(ggpubr)
library(StatisticalModels)
library(scales)
library(viridis)
library(dplyr)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=11, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


Get_effects_resampling <- function(Model, PD, Trait, Class){
  
  
  ## create new data 
  if(PD) {
    
    Newdata <- expand.grid(LandUseGrouped=levels(Model@frame$LandUseGrouped ), 
                           Use_intensity=levels(Model@frame$Use_intensity), 
                           Specialisation =levels(Model@frame$Specialisation),  
                           Diel_activity =levels(Model@frame$Diel_activity),  
                           Primary_diet =levels(Model@frame$Primary_diet),  
                           log10_Body_mass_g=mean(Model@frame$log10_Body_mass_g),  
                           log10_Litter_size=mean(Model@frame$log10_Litter_size),  
                           log10_Lifespan_proxy=mean(Model@frame$ log10_Lifespan_proxy),  
                           log10_Range_area=mean(Model@frame$log10_Range_area),  
                           sqrt_Habitat_breadth_IUCN=mean(Model@frame$sqrt_Habitat_breadth_IUCN),
                           sqrt_Diet_breadth =mean(Model@frame$ sqrt_Diet_breadth),
                           Occurrence=1)
  } else {
    
    Newdata <- expand.grid(LandUseGrouped=levels(Model@frame$LandUseGrouped ), 
                           Use_intensity=levels(Model@frame$Use_intensity), 
                           Specialisation =levels(Model@frame$Specialisation),  
                           Diel_activity =levels(Model@frame$Diel_activity),  
                           #Primary_diet =levels(Model@frame$Primary_diet),  
                           log10_Body_mass_g=mean(Model@frame$log10_Body_mass_g),  
                           log10_Litter_size=mean(Model@frame$log10_Litter_size),  
                           log10_Lifespan_proxy=mean(Model@frame$ log10_Lifespan_proxy),  
                           log10_Range_area=mean(Model@frame$log10_Range_area),  
                           sqrt_Habitat_breadth_IUCN=mean(Model@frame$sqrt_Habitat_breadth_IUCN),
                           sqrt_Diet_breadth =mean(Model@frame$ sqrt_Diet_breadth),
                           Occurrence=1)
  }
  
  
  # subsetting model's coefficients and variance-covariance matrix for the trait
  Fixed_effects <- fixef(Model) 
  Fixed_effects_Trait <- Fixed_effects[grepl(Trait, names(Fixed_effects))]
  
  V_Cov <- vcov(Model)
  # subset relevant rows
  V_Cov_sub <- V_Cov[grepl(Trait, rownames(V_Cov)),] 
  # subset relevant columns
  Col_retain <- which(grepl(Trait, colnames(V_Cov_sub)))
  
  # sanity check
  if(length(Col_retain)!=nrow(V_Cov_sub)){
    print("Problem!")
    stop()
  } else{
    print("Sanity check ok")
  }
  
  V_Cov_sub <- V_Cov_sub[, Col_retain]
  
  # sanity check n2
  if(nrow(V_Cov_sub)!=ncol(V_Cov_sub)){
    print("Problem!")
    stop()
  } else{
    print("Sanity check number 2 ok")
  }
  
  
  print("Resampling in progress...")
  ## getting slopes and 95% CI by re sampling
  preds <- sapply(X=1:5000, FUN=function(i){
    coefs <- mvrnorm(n=1, mu=Fixed_effects_Trait, Sigma=V_Cov_sub)
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
  
  preds_Trait <- data.frame(Median=apply(X=preds, MARGIN=1, FUN=median),
                             Upper=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.975),
                             Lower=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.025))
  
  preds_Trait <- cbind(preds_Trait, Newdata)
  Col <- which(grepl(Trait, colnames(preds_Trait)))
  preds_Trait <- preds_Trait[, c(1:5)]
  preds_Trait <- unique(preds_Trait)
  preds_Trait <- subset(preds_Trait, Median!=0 & Lower!=0 & Upper!=0)
  preds_Trait$Class <- Class
  preds_Trait$Trait <- Trait
  preds_Trait$RefLine <- preds_Trait$Median[preds_Trait$LandUseGrouped=="Primary vegetation" & preds_Trait$Use_intensity=="Minimal use"]
  
  return(preds_Trait)
  
}

## loading models
Birds <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Birds_Diet.rds")
Mammals <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Mammals_Diet.rds")
Amphibians <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Amphibians_Diet.rds") ## no primary diet
Reptiles <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Reptiles_Diet.rds")

###############################################################################################################################
###############################################################################################################################

## Natural habitat specialisation
BirdsSpe <- Get_effects_resampling(Birds, PD=TRUE, "Natural habitat specialist", "Birds")
AmphibiansSpe <- Get_effects_resampling(Amphibians, PD=FALSE, "Natural habitat specialist", "Amphibians")
MammalsSpe <- Get_effects_resampling(Mammals, PD=TRUE, "Natural habitat specialist", "Mammals")
ReptilesSpe <- Get_effects_resampling(Reptiles, PD=FALSE, "Natural habitat specialist", "Reptiles")
SpecialisationEffect <- rbind(BirdsSpe, AmphibiansSpe, MammalsSpe, ReptilesSpe)

# # setting no error bars for minimally-used primary vegetation
# SpecialisationEffect$Lower[SpecialisationEffect$Use_intensity=="Minimal use" & SpecialisationEffect$LandUseGrouped=="Primary vegetation"] <- NA
# SpecialisationEffect$Upper[SpecialisationEffect$Use_intensity=="Minimal use" & SpecialisationEffect$LandUseGrouped=="Primary vegetation"] <- NA

# Amphibians remove effects for urban - errors bars too large. null effects
SpecialisationEffect$Median[SpecialisationEffect$Class=="Amphibians" & SpecialisationEffect$LandUseGrouped=="Urban"] <- NA
SpecialisationEffect$Lower[SpecialisationEffect$Class=="Amphibians" & SpecialisationEffect$LandUseGrouped=="Urban"] <- NA
SpecialisationEffect$Upper[SpecialisationEffect$Class=="Amphibians" & SpecialisationEffect$LandUseGrouped=="Urban"] <- NA

SpecialisationEffect$Class <- factor(SpecialisationEffect$Class, labels = c("(a) Amphibians", "(b) Birds", "(c) Mammals", "(d) Reptiles"))

Plot_Spe <- 
  ggplot(SpecialisationEffect,
         aes(LandUseGrouped, Median, shape=Use_intensity, col=Use_intensity, ymin = Lower, ymax = Upper)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  xlab("") +
  #geom_hline(aes(yintercept=RefLine, group=Use_intensity),linetype="dashed") +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2.5, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary:\nnat. hab. spec.", 
                            "Secondary:\nnat. hab. spec.",
                            "Plantation:\nnat. hab. spec.", 
                            "Agricultural:\nnat. hab. spec.",
                            "Urban:\nnat. hab. spec.",
                            "Light use:\nnat. hab. spec.")) +
  GGPoptions +
  ylab("Effect (+/- 95% confidence interval)")  +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust=1.15)) +
  ggtitle("Natural-habitat specialisation: effects on occurrence probability (log-odds); \nReference: artificial habitat user")+
  viridis::scale_color_viridis(discrete = TRUE, option = "B",  end = 0.7, name="Use intensity") +
  facet_wrap(~Class, scales="free_y")  +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
                                              strip.text.y = element_text(size = 12, face = "bold")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 12, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
  theme(panel.spacing.y = unit(0, "lines"))  +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") 

Plot_Spe

###############################################################################################################################
###############################################################################################################################

BirdsDA <- Get_effects_resampling(Birds, PD=TRUE, "Other", "Birds")
AmphibiansDA <- Get_effects_resampling(Amphibians, PD=FALSE, "Other", "Amphibians")
MammalsDA <- Get_effects_resampling(Mammals, PD=TRUE, "Other", "Mammals")
ReptilesDA <- Get_effects_resampling(Reptiles, PD=FALSE, "Other", "Reptiles")
DAEffect <- rbind(BirdsDA, AmphibiansDA, MammalsDA, ReptilesDA)

# # setting no error bars for minimally-used primary vegetation
# DAEffect$Lower[DAEffect$Use_intensity=="Minimal use" & DAEffect$LandUseGrouped=="Primary vegetation"] <- NA
# DAEffect$Upper[DAEffect$Use_intensity=="Minimal use" & DAEffect$LandUseGrouped=="Primary vegetation"] <- NA

# # Amphibians and mammals: remove effects for urban - errors bars too large. null effects
DAEffect$Median[DAEffect$Class=="Amphibians" & DAEffect$LandUseGrouped=="Urban"] <- NA
DAEffect$Lower[DAEffect$Class=="Amphibians" & DAEffect$LandUseGrouped=="Urban"] <- NA
DAEffect$Upper[DAEffect$Class=="Amphibians" & DAEffect$LandUseGrouped=="Urban"] <- NA
DAEffect$Median[DAEffect$Class=="Mammals" & DAEffect$LandUseGrouped=="Urban"] <- NA
DAEffect$Lower[DAEffect$Class=="Mammals" & DAEffect$LandUseGrouped=="Urban"] <- NA
DAEffect$Upper[DAEffect$Class=="Mammals" & DAEffect$LandUseGrouped=="Urban"] <- NA

DAEffect$Class <- factor(DAEffect$Class, labels = c("(a) Amphibians", "(b) Birds", "(c) Mammals", "(d) Reptiles"))

Plot_DA <- 
  ggplot(DAEffect,
         aes(LandUseGrouped, Median, shape=Use_intensity, col=Use_intensity, ymin = Lower, ymax = Upper)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  xlab("") +
  #geom_hline(aes(yintercept=RefLine, group=Use_intensity),linetype="dashed") +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2.5, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary:\nnon-nocturnal", 
                            "Secondary:\nnon-nocturnal",
                            "Plantation:\nnon-nocturnal", 
                            "Agricultural:\nnon-nocturnal",
                            "Urban:\nnon-nocturnal",
                            "Light use:\nnon-nocturnal")) +
  GGPoptions +
  ylab("Effect (+/- 95% confidence interval)")  +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust=1.15)) +
  ggtitle("Non-nocturnality: effects on occurrence probability (log-odds); \nReference: nocturnal")+
  viridis::scale_color_viridis(discrete = TRUE, option = "B",  end = 0.7, name="Use intensity") +
  facet_wrap(~Class, scales="free_y")  +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 12, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
  theme(panel.spacing.y = unit(0, "lines"))  +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") 

Plot_DA


###############################################################################################################################
###############################################################################################################################
## Primary diet (mammals and birds only)

## birds and mammals -- levels ref = FR|NE

BirdsIN <- Get_effects_resampling(Birds, TRUE, "IN", "Birds")
BirdsIN$Diet <- "IN"
BirdsVE <- Get_effects_resampling(Birds, TRUE, "VE", "Birds")
BirdsVE$Diet <- "VE"
BirdsOM <- Get_effects_resampling(Birds, TRUE, "OM", "Birds")
BirdsOM$Diet <- "OM"
BirdsPLSE <- Get_effects_resampling(Birds, TRUE, "PL|SE", "Birds")
BirdsPLSE$Diet <- "PL|SE"

MammalsIN <- Get_effects_resampling(Mammals, TRUE, "IN", "Mammals")
MammalsIN$Diet <- "IN"
MammalsVE <- Get_effects_resampling(Mammals, TRUE, "VE", "Mammals")
MammalsVE$Diet <- "VE"
MammalsOM <- Get_effects_resampling(Mammals, TRUE, "OM", "Mammals")
MammalsOM$Diet <- "OM"
MammalsPLSE <- Get_effects_resampling(Mammals, TRUE, "PL|SE", "Mammals")
MammalsPLSE$Diet <- "PL|SE"

BirdsPD <- rbind(BirdsIN, BirdsVE, BirdsOM, BirdsPLSE)
BirdsPD$Effect <- rep(c("PV:Diet", "SV:Diet", "PF:Diet", "AGR:Diet", "UR:Diet"), 4)
BirdsPD$Effect <- factor(BirdsPD$Effect, levels = c("PV:Diet", "SV:Diet", "PF:Diet", "AGR:Diet", "UR:Diet"))

MammalsPD <- rbind(MammalsIN, MammalsVE, MammalsOM, MammalsPLSE)
MammalsPD$Effect <- rep(c("PV:Diet", "SV:Diet", "PF:Diet", "AGR:Diet", "UR:Diet"),12)

## remove mammal omnivores in urban -- error bars too large
MammalsPD$Median[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"] <- NA
MammalsPD$Lower[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"] <- NA
MammalsPD$Upper[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"] <- NA

PDeffects <- rbind(BirdsPD, MammalsPD)
PDeffects$Diet <- factor(PDeffects$Diet,
                         levels=c("IN", "VE", "PL|SE", "OM"),
                         labels=c("Invertebrates", "Vertebrates", "Plants/seeds", "Omnivore"))
Plot_PD <-
ggplot(PDeffects,
       aes(LandUseGrouped, Median, shape=Use_intensity, col=Use_intensity, ymin = Lower, ymax = Upper)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  xlab("") +
  geom_hline(yintercept = 0,linetype="dashed") +
  #geom_hline(aes(yintercept=RefLine, group=Use_intensity),linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2.5, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary:diet", 
                            "Secondary:diet",
                            "Plantation:diet", 
                            "Agricultural:diet",
                            "Urban:diet")) +
  GGPoptions +
  ylab("Effect (+/- 95% confidence interval)")  +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust=1.15)) +
  ggtitle("Primary diet: effects on occurrence probability (log-odds); \nReference: fruit/nectar eater")+
  viridis::scale_color_viridis(discrete = TRUE, option = "B",  end = 0.7, name="Use intensity") +
  facet_grid(Class~Diet, scales="free_y")  +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 1, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 12, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
  theme(panel.spacing = unit(0, "lines"))  +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") 

Plot_PD



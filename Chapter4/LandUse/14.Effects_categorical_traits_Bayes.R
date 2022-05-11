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


Get_effects_cat <- function(Model, Trait, Class) {
  
  Summary <- as.data.frame(summary(Model)$solutions)
  
  # subset the summary for the effect of the trait on ocurrence probability (log-odds)
  Sub <- Summary[grepl(Trait, rownames(Summary)),]
  # Sub$CI_up <- Sub$Estimate + Sub$`Std. Error`*1.96
  # Sub$CI_low <- Sub$Estimate - Sub$`Std. Error`*1.96
  
  # function for backtransforming
  # logit2prob <- function(logit){
  #   return(exp(logit)/(1+exp(logit)))
  # }
  # 
  # Sub$Estimate <- logit2prob(Sub$Estimate)
  # Sub$CI_up <- logit2prob(Sub$CI_up)
  # Sub$CI_low <- logit2prob(Sub$CI_low)
  
  Sub$Class <- Class
  Sub$Effect <- outer(c("PV", "SV", "PF", "AGR", "UR", "Light use", "Intense use"), paste0(":", Trait), FUN="paste0")
  Sub$Effect<- factor(Sub$Effect, levels=outer(c("PV", "SV", "PF", "AGR", "UR", "Light use", "Intense use"), paste0(":", Trait), FUN="paste0"))
  
  return(Sub)
}

## loading models
gc()
memory.limit(size=900000000000)
Birds <- readRDS("../Results/MCMC_glmm_models/Class_specific/Birds_Diet.rds")
Mammals <- readRDS("../Results/MCMC_glmm_models/Class_specific/Mammals_Diet.rds")
Amphibians <- readRDS("../Results/MCMC_glmm_models/Class_specific/Amphibians_Diet.rds") ## no primary diet
Reptiles <- readRDS("../Results/MCMC_glmm_models/Class_specific/Reptiles_Diet.rds")


## looking at effects from summary
summary(Amphibians)
summary(Reptiles)

###############################################################################################################################

## Natural habitat specialisation
BirdsSpe <- Get_effects_cat(Birds, "Natural habitat specialist", "Birds")
AmphibiansSpe <- Get_effects_cat(Amphibians, "Natural habitat specialist", "Amphibians")
# AmphibiansSpe$CI_up[5] <- NA
# AmphibiansSpe$CI_low[5] <- NA
# AmphibiansSpe$Estimate[5] <- NA
MammalsSpe <- Get_effects_cat(Mammals, "Natural habitat specialist", "Mammals")
ReptilesSpe <- Get_effects_cat(Reptiles, "Natural habitat specialist", "Reptiles")
SpecialisationEffect <- rbind(BirdsSpe, AmphibiansSpe, MammalsSpe, ReptilesSpe)
colnames(SpecialisationEffect)[1:3] <- c("post.mean", "CI_low", "CI_up")

Plot_Spe <- 
  ggplot(SpecialisationEffect,
         aes(Effect, post.mean, shape=Class, group=Class, col=Class, ymin = CI_low, ymax = CI_up)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  #geom_rect(xmin=5.5, xmax=6.5,ymin=-Inf,ymax=Inf, fill="beige", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=5.5, xmax=5.52,ymin=-Inf,ymax=Inf, col="black", linetype="solid", size=1) +
  xlab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2.5, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary:\nnat. hab. spec.", 
                            "Secondary:\nnat. hab. spec.",
                            "Plantation:\nnat. hab. spec.", 
                            "Agricultural:\nnat. hab. spec.",
                            "Urban:\nnat. hab. spec.",
                            "Light use:\nnat. hab. spec.",
                            "Intense use:\nnat. hab. spec.")) +
  GGPoptions +
  ylab("Posterior mean (+/- 95% credible interval)")  +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust=1.15)) +
  ggtitle("(a) Natural-habitat specialisation: \neffects on occurrence probability (log-odds); \nReference: artificial habitat user")+
  viridis::scale_color_viridis(discrete = TRUE,option = "B",  end = 0.7)
Plot_Spe

## Diel activity
BirdsDa <- Get_effects_cat(Birds, "Other", "Birds")
AmphibiansDa <- Get_effects_cat(Amphibians, "Other", "Amphibians")
# AmphibiansDa$CI_up[5] <- NA
# AmphibiansDa$CI_low[5] <- NA
# AmphibiansDa$Estimate[5] <- NA
MammalsDa <- Get_effects_cat(Mammals, "Other", "Mammals")
# MammalsDa$CI_up[5] <- NA
# MammalsDa$CI_low[5] <- NA
# MammalsDa$Estimate[5] <- NA
ReptilesDa <- Get_effects_cat(Reptiles, "Other", "Reptiles")
DaEffect <- rbind(BirdsDa, AmphibiansDa, MammalsDa, ReptilesDa)
colnames(DaEffect)[1:3] <- c("post.mean", "CI_low", "CI_up")

Plot_DA <- 
  ggplot(DaEffect,
         aes(Effect, post.mean, shape=Class, group=Class, ymin = CI_low, ymax = CI_up, col=Class)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=5.5, xmax=5.52,ymin=-Inf,ymax=Inf, col="black", linetype="solid", size=1) +
  xlab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2.5, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary:\nnon-nocturnal", 
                            "Secondary:\nnon-nocturnal",
                            "Plantation:\nnon-nocturnal", 
                            "Agricultural:\nnon-nocturnal",
                            "Urban:\nnon-nocturnal",
                            "Light use:\nnon-nocturnal",
                            "Intense use:\nnon-nocturnal")) +
  GGPoptions +
  ylab("Posterior mean (+/- 95% credible interval)")  +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.15, hjust=1.15)) +
  ggtitle("(b) Non-nocturnality: \neffects on occurrence probability (log-odds); \nReference: nocturnal") +
  viridis::scale_color_viridis(discrete = TRUE,option = "B",  end = 0.7)
Plot_DA

## Primary diet (mammals and birds only)

## birds and mammals -- levels ref = FR|NE
# Birds@frame$Primary_diet %>%  unique
# Mammals@frame$Primary_diet %>%  unique

BirdsIN<- Get_effects_cat(Birds, "IN", "Birds")
BirdsIN$Diet <- "IN"
BirdsVE<- Get_effects_cat(Birds, "VE", "Birds")
BirdsVE$Diet <- "VE"
BirdsOM<- Get_effects_cat(Birds, "OM", "Birds")
BirdsOM$Diet <- "OM"
BirdsPLSE<- Get_effects_cat(Birds, "PL|SE", "Birds")
BirdsPLSE$Diet <- "PL|SE"

MammalsIN<- Get_effects_cat(Mammals, "IN", "Mammals")
MammalsIN$Diet <- "IN"
MammalsVE<- Get_effects_cat(Mammals, "VE", "Mammals")
MammalsVE$Diet <- "VE"
MammalsOM<- Get_effects_cat(Mammals, "OM", "Mammals")
MammalsOM$Diet <- "OM"
MammalsPLSE<- Get_effects_cat(Mammals, "PL|SE", "Mammals")
MammalsPLSE$Diet <- "PL|SE"

BirdsPD <- rbind(BirdsIN, BirdsVE, BirdsOM, BirdsPLSE)
BirdsPD$Effect <- rep(c("PV:Diet", "SV:Diet", "PF:Diet", "AGR:Diet", "UR:Diet", "Light use:Diet", "Intense use:Diet"), 4)
BirdsPD$Effect <- factor(BirdsPD$Effect, levels = c("PV:Diet", "SV:Diet", "PF:Diet", "AGR:Diet", "UR:Diet", "Light use:Diet", "Intense use:Diet"))

MammalsPD <- rbind(MammalsIN, MammalsVE, MammalsOM, MammalsPLSE)
MammalsPD$Effect <- rep(c("PV:Diet", "SV:Diet", "PF:Diet", "AGR:Diet", "UR:Diet", "Light use:Diet", "Intense use:Diet"), 4)
MammalsPD$Effect <- factor(MammalsPD$Effect, levels = c("PV:Diet", "SV:Diet", "PF:Diet", "AGR:Diet", "UR:Diet", "Light use:Diet", "Intense use:Diet"))

## remove mammal omnivores in urban -- error bars too large
# MammalsPD$Estimate[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"]
# MammalsPD$CI_up[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"]
# MammalsPD$CI_low[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"]
# 
# MammalsPD$Estimate[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"] <- NA
# MammalsPD$CI_up[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"] <- NA
# MammalsPD$CI_low[MammalsPD$Effect=="UR:Diet" & MammalsPD$Diet=="OM"] <- NA


PDeffects <- rbind(BirdsPD, MammalsPD)
PDeffects$Diet <- factor(PDeffects$Diet,
                         levels=c("IN", "VE", "PL|SE", "OM"),
                         labels=c("Invertebrates", "Vertebrates", "Plants/seeds", "Omnivore"))
colnames(PDeffects)[1:3] <- c("post.mean", "CI_low", "CI_up")


Plot_PD <- 
  ggplot(PDeffects,
         aes(Effect, post.mean, shape=Diet, group=Diet, col=Diet, ymin = CI_low, ymax = CI_up)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=8,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=5.5, xmax=5.52,ymin=-Inf,ymax=Inf, col="black", linetype="solid", size=1) +
  xlab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, stat="identity",position=position_dodge(width = 0.6)) +
  geom_point(size=2.5, position=position_dodge(width = 0.6)) +
  scale_x_discrete(labels=c("Primary:\ndiet", 
                            "Secondary:\ndiet",
                            "Plantation:\ndiet", 
                            "Agricultural:\ndiet",
                            "Urban:\ndiet",
                            "Light use:\ndiet",
                            "Intense use:\ndiet")) +
  GGPoptions +
  ylab("Posterior mean (+/- 95% credible interval)")  +
  theme(axis.text.x = element_text(angle = 60, vjust = 1.05, hjust=1)) +
  ggtitle("(c) Primary diet for birds and mammals: \neffects on occurrence probability (log-odds); \nReference: Fruit/nectar eater") +
  facet_wrap(~Class, nrow=2, scales = "free_y") + 
  theme(panel.spacing =  unit(0, "lines")) +
  scale_colour_colorblind() +
  theme(legend.position = "top")

Plot_PD

############################################################################################


## Arranging plots for manuscript

Plot <-
  ggarrange(
    ggarrange(Plot_Spe, Plot_DA, common.legend = TRUE, nrow=2),
    Plot_PD
  )


ggsave(Plot, filename="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/MCMCglmm_Land_use_categorical_traits.pdf",
       width=12, height = 11)


# probs <- seq(-5,1,0.01)
# plot(probs, logit2prob(probs), type="line")






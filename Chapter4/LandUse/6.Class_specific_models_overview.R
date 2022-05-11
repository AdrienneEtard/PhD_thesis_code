## Analysing GLMER models

library(lme4)
library(performance)
library(DHARMa)
library(dplyr)
library(ggpubr)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=11, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=11))

AnovaVar <- function(model) {
  #browser()
  ## anova on model terms
  AnovaVar <- anova(model)
  Anova_ss <- AnovaVar$"Sum Sq"
  Anova_ss <- cbind(AnovaVar, PctExp=Anova_ss/sum(Anova_ss)*100)
  Anova_ss <- Anova_ss[order(Anova_ss$PctExp),]
  return(Anova_ss)
}

#setwd("E:/3.Explanatory_traits/Code/")

## loading model for amphibians
Amphibians <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Amphibians_Diet.rds")
Birds <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Birds_Diet.rds")
Mammals <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Mammals_Diet.rds")
Reptiles <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Reptiles_Diet.rds")

Amphibians@frame$Best_guess_binomial %>%  unique() %>%  length
Birds@frame$Best_guess_binomial %>%  unique() %>%  length
Mammals@frame$Best_guess_binomial %>%  unique() %>%  length
Reptiles@frame$Best_guess_binomial %>%  unique() %>%  length

################################################################################################################################
# model diagnostics
r2(Amphibians)
r2(Birds)
r2(Mammals)
r2(Reptiles)

# diagnostic plots with DHARMa
simulationOutputA <- simulateResiduals(fittedModel = Amphibians)
simulationOutputB <- simulateResiduals(fittedModel = Birds)
simulationOutputM <- simulateResiduals(fittedModel = Mammals)
simulationOutputR <- simulateResiduals(fittedModel = Reptiles)
plot(simulationOutputA)
plot(simulationOutputB)
plot(simulationOutputM)
plot(simulationOutputR)


################################################################################################################################
## plotting sample sizes

ss_amphibians <- Amphibians@frame %>% 
  group_by(LandUseGrouped, Use_intensity, SSBS) %>% 
  summarise(C=n()) %>% 
  group_by(LandUseGrouped, Use_intensity) %>% 
  summarise(C=n()) %>%
  mutate(Class="(a) Amphibians")

ss_birds <- Birds@frame %>% 
  group_by(LandUseGrouped, Use_intensity, SSBS) %>% 
  summarise(C=n()) %>% 
  group_by(LandUseGrouped, Use_intensity) %>% 
  summarise(C=n()) %>%
  mutate(Class="(b) Birds")

ss_reptiles <- Reptiles@frame %>% 
  group_by(LandUseGrouped, Use_intensity, SSBS) %>% 
  summarise(C=n()) %>% 
  group_by(LandUseGrouped, Use_intensity) %>% 
  summarise(C=n()) %>%
  mutate(Class="(c) Reptiles")

ss_mammals <- Mammals@frame %>% 
  group_by(LandUseGrouped, Use_intensity, SSBS) %>% 
  summarise(C=n()) %>% 
  group_by(LandUseGrouped, Use_intensity) %>% 
  summarise(C=n()) %>%
  mutate(Class="(d) Mammals")

ss <- rbind(ss_amphibians, ss_birds, ss_mammals, ss_reptiles)

Plot_SS <- 
    ggplot(ss, aes(LandUseGrouped, C, group=Use_intensity, fill=Use_intensity)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values=c("#000000", "coral", "coral4"), name="Use intensity") +
    GGPoptions+
    scale_x_discrete(limits=c("Primary vegetation", 
                              "Secondary vegetation", 
                              "Plantation forest", 
                              "AGR", 
                              "Urban"),
                     labels=c("Primary",
                              "Secondary",
                              "Plantation",
                              "Agricultural",
                              "Urban"))+ 
    xlab("") + ylab("Number of sites") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
    facet_wrap(~Class, scales="free") +
    geom_text(aes(label=C), position=position_dodge(width=0.9), vjust=-0.25, size=3) +
  theme(legend.position = "top")
Plot_SS
ggsave(Plot_SS, 
       filename="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/Sample_size_figure.pdf",
       width=9, height = 10)

################################################################################################################################

## phylogenetic signal in models' residuals

Phylo_Mammals <- read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Mammals_complete_TreeAnnotator.nex")  
Phylo_Birds <- read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Birds_TreeAnnotator.nex") 
Phylo_Amphibians <- read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Amphibians_TreeAnnotator.nex")  
Phylo_Reptiles <- read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Reptiles_TreeAnnotator.nex")  


PhySignal <- function(residuals, Names, Tree, simulations, N) {
  
  .Format_tiplabels <- function (Phylogeny) {
    Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
    return(Phylogeny)
  }
  
  Residuals <- residuals
  names(Residuals) <- Names
  
  Res_lamba <- c()
  Tree <- .Format_tiplabels(Tree)
  Signal <- phytools::phylosig(Tree, Residuals, method="lambda", test = TRUE) %>% 
    unlist()
  if(!simulations){return(Signal)} 
  
  if(simulations){
    Res_lambda <- Signal$lambda
    for(i in 1:N) {
      # randomise residuals
      Residuals_sim <- sample(residuals, size = length(residuals), replace = FALSE)
      names(Residuals_sim) <- Names
      Signalsim <- phytools::phylosig(Tree, Residuals_sim, method="lambda", test = FALSE) %>% 
        unlist()
      Res_lamba <- c(Res_lamba, Signalsim$lambda)
      print(i)
    }
    return(Res_lamba)
  }
}

Run <- function(Model, Phylo){
  Residuals <- residuals(Model)
  Names <- Model@frame$Best_guess_binomial
  S <- Sys.time()
  lambda_residuals <- PhySignal(Residuals, Names,Tree= Phylo, simulations = FALSE, N=NULL)
  E <- Sys.time()
  print(S-E)
  return(lambda_residuals)

}

LambdaMammals <- Run(Mammals, Phylo_Mammals)
LambdaAmphibians <- Run(Amphibians, Phylo_Amphibians)
LambdaReptiles <- Run(Reptiles, Phylo_Reptiles)
LambdaBirds <- Run(Birds, Phylo_Birds)

################################################################################################################################
# variation explained by the different predictors

Anova_var_amphibians <- AnovaVar(Amphibians)
Anova_var_amphibians$Class <- "Amphibians"

Anova_var_birds <- AnovaVar(Birds)
Anova_var_birds$Class <- "Birds"

Anova_var_mammals <- AnovaVar(Mammals)
Anova_var_mammals$Class <- "Mammals"

Anova_var_reptiles <- AnovaVar(Reptiles)
Anova_var_reptiles$Class <- "Reptiles"

Anova_classes <- rbind(Anova_var_amphibians, Anova_var_birds, Anova_var_mammals, Anova_var_reptiles)
Anova_classes$Effects <- rownames(Anova_classes)

More_than_10 <- Anova_classes %>% 
  group_by(Class) %>%
  filter(PctExp>=5) %>% 
  dplyr::select(Class, PctExp, Effects)

Less_than_10 <- Anova_classes %>% 
  group_by(Class) %>%
  filter(PctExp<5) %>% 
  summarise(PctExp=sum(PctExp)) %>% 
  mutate(Effects="Other effects")

PctExpClasses <- rbind(More_than_10, Less_than_10)
# PctExpClasses$Effects <- c("Land use:Litter/clutch size",
#                            "Land use:Range area",
#                            "Land use:Habitat breadth",
#                            "Land use:Primary diet",
#                            "Land use:Range area",
#                            "Land use:Habitat breadth",
#                            "Land use:Primary diet",
#                            "Use intensity:Body mass",
#                            "Land use:Range area",
#                            "Land use:Diet breadth",
#                            "Land use:Lifespan proxy",
#                            rep("Other effects",4))

PctExpClasses$Effects <- c("Land-use intensity:Artificial habitat use",
                           "Land-use intensity:Habitat breadth",
                           "Land use",
                           "Land use:Land-use intensity",
                           "Land use:Litter/clutch size",
                           "Land use:Range area",
                           "Land use:Habitat breadth",
                           "Land use:Land-use intensity",
                           "Land use:Primary diet",
                           "Land use:Range area",
                           "Land use:Habitat breadth",
                           "Land-use intensity:Diel activity",
                           "Land use:Land-use intensity",
                           "Land use:Range area",
                           "Primary diet",
                           "Land-use intensity:Primary diet",
                           "Land use:Diet breadth",
                           "Land use:Primary diet",
                           "Land-use intensity:Body mass",   
                           "Land use:Body mass",
                           "Body mass",
                           "Land use:Diel activity",
                           "Land-use intensity:Diet breadth",
                           "Land-use intensity:Habitat breadth",
                           "Land use:Habitat breadth",
                           "Land use",
                           "Land use:Range area",
                           "Land use:Diet breadth",
                           "Land use:Lifespan proxy",
                           rep("Other effects",4))


# PctExpClasses$Effects <- factor(PctExpClasses$Effects, levels=c("Land use:Range area",
#                                                                 "Land use:Litter/clutch size",
#                                                                 "Land use:Diet breadth",
#                                                                 "Land use:Primary diet",     
#                                                                 "Use intensity:Body mass",
#                                                                 "Land use:Lifespan proxy",
#                                                                 "Land use:Habitat breadth",
#                                                                 "Other effects"))

PctExpClasses$Effects <- factor(PctExpClasses$Effects, levels=unique(PctExpClasses$Effects))

PctExpClasses <- PctExpClasses[order(PctExpClasses$PctExp, decreasing = FALSE),]

PctExpClasses2 <- PctExpClasses[order(PctExpClasses$PctExp, decreasing = FALSE),]
PctExpClasses2 <- PctExpClasses2 %>%  group_by(Class) %>% mutate(position=rank(PctExp))

## plotting % variation explained by the factors (barplot)
VarLU_plot <- 
ggplot(PctExpClasses2, aes(Effects, PctExp, fill=Class, group=position)) + 
  geom_bar(stat="identity", position="dodge", width=0.7) +
  GGPoptions +  
  coord_flip() + ylab("Variance explained (%)") + xlab("") +
  scale_fill_viridis_d() + 
  geom_hline(yintercept=10, lty="dashed") + facet_grid(~Class) + 
  theme(legend.position = "none")+
  theme(panel.spacing =  unit(0, "lines"))+
  scale_x_discrete(limits = rev(levels(PctExpClasses2$Effects))) +
  ggtitle("(a) % variance explained by main effects (land-use models, factoring out residual variation)") +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
                    strip.text.y = element_text(size = 12, face = "bold"))

CheckLevels <- unique(PctExpClasses2$Effects)
CheckLevels[order(CheckLevels)]

## combining with variance explained from climate-change models

VarSensi <- readRDS("../Results/Results_climate_sensitivity/Variance_explained_withoutres.rds")
#VarSensi <- readRDS("E:/3.Explanatory_traits/Results/Results_climate_sensitivity/Variance_explained_withres.rds")

VarSensi <- VarSensi[order(VarSensi$PctExp, decreasing = FALSE),]
VarSensi <- VarSensi %>%  group_by(Class) %>% mutate(position=rank(PctExp))
VarSensi$Effects[VarSensi$Effects=="Other effects"] <-
  "                              Other effects"

## plotting % variation explained by the factors (barplot)
VarCl_plot <- 
ggplot(VarSensi, aes(Effects, PctExp, fill=Class, group=position)) + 
  geom_bar(stat="identity", position="dodge", width=0.7) +
  GGPoptions +  
  coord_flip() + ylab("Variance explained (%)") + xlab("") +
  scale_fill_viridis_d() + 
  geom_hline(yintercept=10, lty="dashed") + facet_grid(~Class) + 
  theme(legend.position = "none")+
  theme(panel.spacing =  unit(0, "lines"))+
  scale_x_discrete(limits = rev(levels(VarSensi$Effects))) +
  ggtitle("(b) % variance explained by main effects (PGLS models, factoring out residual variation)") +
  facet_grid(~Class) + scale_y_continuous(limits = c(0, 99), breaks = seq(0, 99, by = 20))+
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold"))


## wihout residual var and without range area
VarSensi2 <- readRDS("../Results/Results_climate_sensitivity/Variance_explained_withoutres_withoutRA.rds")
VarSensi2 <- VarSensi2[order(VarSensi2$PctExp, decreasing = FALSE),]
VarSensi2 <- VarSensi2 %>%  group_by(Class) %>% mutate(position=rank(PctExp))
VarSensi2$Effects[VarSensi2$Effects=="Other effects"] <-
  "                             Other effects"

VarCl_plot2 <- 
  ggplot(VarSensi2, aes(Effects, PctExp, fill=Class, group=position)) + 
  geom_bar(stat="identity", position="dodge", width=0.7) +
  GGPoptions +  
  coord_flip() + ylab("Variance explained (%)") + xlab("") +
  scale_fill_viridis_d() + 
  geom_hline(yintercept=10, lty="dashed") + facet_grid(~Class) + 
  theme(legend.position = "none")+
  theme(panel.spacing =  unit(0, "lines"))+
  scale_x_discrete(limits = rev(levels(VarSensi$Effects))) +
  ggtitle("(c) % variance explained by main effects (PGLS models, factoring out residual variation and variation attributable to range area)") +
  facet_grid(~Class) + scale_y_continuous(limits = c(0, 99), breaks = seq(0, 99, by = 20))+
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold"))


library(ggpubr)
P_Variance <- 
ggarrange(
VarLU_plot, 
VarCl_plot,
VarCl_plot2,
nrow=3, heights = c(0.55, 0.2, 0.25))
P_Variance

ggsave(P_Variance, 
       filename = "../Results/Results_plots/Variance_breakdown.pdf",
       width=12.2, height=11)

################################################################################################################################
# plotting effects: GLMER vs mcmcglmm


################################################################################################################################
# Plotting maps of PREDICTS sites
library("ggmap")
library("rnaturalearth")
library("rnaturalearthdata")
library("patchwork")

Predicts <- readRDS("../Results/Predicts_merged_sites.rds")

Draw_map <- function(Model, Class){
  
  Sites <- Predicts[Predicts$SSBS %in% Model@frame$SSBS,]
  print(length(unique(Sites$SSBS)))
  
  world <- ne_countries(scale = "medium", returnclass = "sf")

  # ggplot2::ggplot(data = world) + xlab("") + ylab("") +
  # geom_sf(fill="lightgrey", colour="lightgrey") +
  # geom_point(data=Sites[!is.na(Sites$Longitude),],
  #            aes(y=Latitude, x=Longitude), size=1.5) +
  # theme_classic() + GGPoptions + theme(panel.grid.major = element_line(colour = 'transparent')) +
  # theme(legend.title = element_blank())# + ggtitle(Class)
  # 
  }

Draw_map(Amphibians, "Amphibians") +
Draw_map(Birds, "Birds") +
Draw_map(Mammals, "Mammals") +
Draw_map(Reptiles, "Reptiles")






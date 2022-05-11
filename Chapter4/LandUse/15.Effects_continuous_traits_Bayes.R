## continous traits: effects of land use & use intensity on slopes of relationships between traits and occurrence probability

##############################################################################################################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(StatisticalModels)
library(scales)
library(viridis)
library(dplyr)
library(lme4)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

## function to get slopes by resampling in the models posteriors

GetSlopes <- function(Sols, Trait){
  
  SolTrait <- Sols[, grepl(Trait, colnames(Sols))] 
  
  ## add columns together -- estimates for minimal use
  Minimal_uses <- data.frame(SolTrait[,1], SolTrait[,2:(ncol(SolTrait)-2)] + SolTrait[,1])
  
  ## for light uses
  Light_uses <- Minimal_uses + SolTrait[,ncol(SolTrait)-1]
  
  ## for intense uses
  Intense_uses <- Minimal_uses + SolTrait[,ncol(SolTrait)]
  
  ## get median and 95%CI for each column
  Qt <- apply(Minimal_uses, 2, quantile, probs = c(0.025, 0.975)) %>%
    t %>% 
    as.data.frame()
  Qt$Median <- apply(Minimal_uses, 2, median)
  Qt$Use_intensity <- "Minimal use"
  
  Qt2 <- apply(Light_uses, 2, quantile, probs = c(0.025, 0.975)) %>%
    t %>% 
    as.data.frame()
  Qt2$Median <- apply(Light_uses, 2, median)
  Qt2$Use_intensity <- "Light use"
  
  Qt3 <- apply(Intense_uses, 2, quantile, probs = c(0.025, 0.975)) %>%
    t %>% 
    as.data.frame()
  Qt3$Median <- apply(Intense_uses, 2, median)
  Qt3$Use_intensity <- "Intense use"
  
  
  Q <- rbind(Qt, Qt2, Qt3)
  
  return(Q)
}

Run_on_traits <- function(Model, Traits, Class) {
  
  Sols <- Model$Sol
  
  ResList <- list()
  
  for (t in Traits){
    cat("Running for", t, "\n", "\n")
    res <- GetSlopes(Sols, Trait = t)
    res$Trait <- t
    res$LandUseGrouped <- rep(c("Primary", "Secondary", "Plantation", "AGR", "Urban"),3)
    res$Class <-Class
    ResList[[t]] <- res
  }
  
  return(data.table::rbindlist(ResList))
  
}

## CONTINUOUS TRAITS:

Cont_traits <- c("log10_Body_mass_g",
                 "log10_Litter_size",
                 "log10_Lifespan_proxy",
                 "log10_Range_area",
                 "sqrt_Habitat_breadth_IUCN",
                 "sqrt_Diet_breadth")


## loading models
gc()
memory.limit(size=90000000000)
Birds <- readRDS("../Results/MCMC_glmm_models/Class_specific/Birds_Diet.rds")
Mammals <- readRDS("../Results/MCMC_glmm_models/Class_specific/Mammals_Diet.rds")
Amphibians <- readRDS("../Results/MCMC_glmm_models/Class_specific/Amphibians_Diet.rds") ## no primary diet

## SLOPES for mammals
SlopesMammals <- Run_on_traits(Model = Mammals, Traits = Cont_traits, Class = "Mammals")
## SLOPES for amphibians
SlopesAmphibians <- Run_on_traits(Model = Amphibians, Traits = Cont_traits, Class="Amphibians")
## SLOPES for birds
SlopesBirds <- Run_on_traits(Model = Birds, Traits = Cont_traits, Class="Birds")

rm(Mammals, Amphibians)

Reptiles <- readRDS("../Results/MCMC_glmm_models/Class_specific/Reptiles_Diet.rds")
## SLOPES for reptiles
SlopesReptiles <- Run_on_traits(Model = Reptiles, Traits = Cont_traits, Class="Reptiles")
rm(Reptiles)

##############################################################################################################################

## all slopes
Slopes <- rbind(SlopesAmphibians, SlopesBirds, 
                SlopesMammals, SlopesReptiles)

colnames(Slopes)[1:2] <- c("Lower", "Upper")

## setting CI for reference to NA ## minimally used PV
Slopes$Lower[Slopes$LandUseGrouped=="Primary" & Slopes$Use_intensity=="Minimal use"] <- NA
Slopes$Upper[Slopes$LandUseGrouped=="Primary" & Slopes$Use_intensity=="Minimal use"] <- NA


# ## removing certain effects before plotting because of massive error bars
 
## removing all traits, amphibians, urban -- all null effects
Slopes$Median[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

# ## removing diet breadth urban mammals
Slopes$Median[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Mammals"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Mammals"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Mammals"] <- NA

## removing all traits, reptiles, urban -- all null effects
Slopes$Median[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Reptiles
Slopes$Median[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA

Slopes$Median[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA

Slopes$Median[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA

Slopes$Median[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA

Slopes$Median[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Reptiles"] <- NA


############################################# Plotting slopes


Plot_slopes_fun_Allfacet <- function(SlopesDF) {
  
  #browser()
  
  SlopesDF$H_line <- NA
  
  for(comb in levels(SlopesDF$Comb)){
    
    ref <- SlopesDF$Median[SlopesDF$Comb==comb & 
                             SlopesDF$LandUseGrouped=="Primary" &
                             SlopesDF$Use_intensity=="Minimal use"]
    SlopesDF$H_line[SlopesDF$Comb==comb] <- ref
  }
  
  # ## assigning lines of reference
  # Vals <- SlopesDF %>%
  #   dplyr::group_by(Comb) %>%
  #   dplyr::filter(LandUseGrouped=="Primary") %>%
  #   dplyr::filter(Use_intensity=="Minimal use") %>%
  #   dplyr::distinct(Median) %>%
  #   as.data.frame()
  # 
  # #Vals <- Vals[order(match(Vals$Comb, levels(SlopesDF$Comb))),]
  # 
  # Vals <- Vals$Median
  # 
  # Mult <- SlopesDF %>%
  #   group_by(Comb) %>%
  #   summarise(C=n())
  # 
  # Mult <- Mult$C
  # 
  # SlopesDF$H_line <- NA
  # 
  # Seq <- 1
  # for(i in 1:length(Mult)){
  #   a <- Seq
  #   b <- Seq + Mult[i]-1
  #   Seq <- Seq + Mult[i]
  #   SlopesDF$H_line[a:b] <- Vals[i]
  # 
  # }
  
  browser()
  
  ggplot(SlopesDF,
         aes(LandUseGrouped, Median, ymin = Lower, ymax = Upper, group=Use_intensity, shape=Use_intensity, col=Use_intensity)) +
    geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=8.5, xmax=9.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    ylab("") + xlab("") +
    geom_hline(aes(yintercept=H_line, group=Use_intensity),linetype="dashed") +
    geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
    geom_point(size=2, position=position_dodge(width = 0.7)) +
    GGPoptions +
    ggtitle("") +
    theme(panel.spacing = unit(0, "lines")) +
    theme( strip.text.x = element_text(size = 12, face = "bold"),
           strip.text.y = element_text(size = 12, face = "bold")) +
    scale_colour_manual(values=c("black", "#FF7F0E", "#D55E00"), name="Use intensity")+
    theme(legend.position = "right") +
    scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
    ylab("Slope estimate (+/- 95% CI)") +
    theme(legend.position = "top") + facet_wrap(~Comb, scale="free_y", ncol=4)+
    scale_x_discrete(labels=c("Primary", "Secondary","Plantation", "Agricultural", "Urban")) +
    theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
          axis.text.x = element_text(size = 10, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")))
  
}

Slopes$Comb <- paste(Slopes$Class, Slopes$Trait)
Slopes$Comb <- as.factor(Slopes$Comb)
Slopes$Comb %>%  levels()

Names <- expand.grid(c("BM", "LP","LCS","RA", "DB", "HB"), c("Amphibians", "Birds","Mammals", "Reptiles"))
Names <- paste0(Names$Var2, ": ",Names$Var1)
Levels  <- expand.grid(c("Amphibians", "Birds", "Mammals", "Reptiles"), c("BM", "LP","LCS","RA", "DB", "HB"))
Levels <- paste0(Levels$Var1, ": ",Levels$Var2)
Slopes$Comb <- factor(Slopes$Comb, labels=Names)
Slopes$Comb <- factor(Slopes$Comb, levels=Levels)

Slopes$Use_intensity <- factor(Slopes$Use_intensity, levels=c("Minimal use", "Light use", "Intense use"))
Slopes$LandUseGrouped <- factor(Slopes$LandUseGrouped,
                                levels=c("Primary", 
                                         "Secondary", 
                                         "Plantation",
                                         "AGR", 
                                         "Urban"))

p <- Plot_slopes_fun_Allfacet(Slopes)
p


ggsave(p, 
       filename="E:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Figures/Slopes_LU_mcmcglmm.pdf",
       width = 8, height = 10)


















## Looking at the models -- filtered for species whose RS > threshold (100km2)

library(dplyr)
library(caper)

#setwd("F:/PhD/PhD_R_projects/5.Climatic_niche_space/Code/")

## mammals, with diet
MammalsModel <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Mammals_diet_v2.rds")
MammalsModel$call
MammalsModel$fitted %>%  nrow() ## sample size: 4712
MammalsModel$data$data %>%  nrow() ## sample size: 4712
as.data.frame(summary(MammalsModel)$coefficients)

## reptiles, with diet
ReptilesModel <- readRDS("../Results/5.PGLS_models_results/Species_over_100km/Reptiles_diet_v2.rds")
ReptilesModel$call
ReptilesModel$data$data$Primary_diet %>% table
ReptilesModel$data$data$sqrt_Diet_breadth %>% table
ReptilesModel$fitted %>%  nrow() ## sample size: 7330
ReptilesModel$data$data %>%  nrow() ## sample size: 7330

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

##############################################################################
## summaries
summary(AmphibiansModel)
summary(MammalsModel)
summary(BirdsModel)
summary(ReptilesModel)

##############################################################################
## plotting distribution of log10 sensitivity for Figure 1C of manuscript

S_Mammals <- data.frame(Sensitivitylog10=log10(MammalsModel$data$data$sensitivity)) 
S_Mammals$Class <- "Mammals"
S_Birds <- data.frame(Sensitivitylog10=log10(BirdsModel$data$data$sensitivity))
S_Birds$Class <- "Birds"
S_Amphibians <- data.frame(Sensitivitylog10=log10(AmphibiansModel$data$data$sensitivity))
S_Amphibians$Class <- "Amphibians"
S_Reptiles <- data.frame(Sensitivitylog10=log10(ReptilesModel$data$data$sensitivity))
S_Reptiles$Class <- "Reptiles"

Sdf <- rbind(S_Mammals, S_Birds, S_Amphibians, S_Reptiles)

library(ggplot2)
library(viridis)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=15, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=15), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=15),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=15))

ggplot(Sdf, aes(Sensitivitylog10)) + geom_histogram() + GGPoptions +
  ylab("Number of species") + xlab("Climate change sensitivity (log10)") +
  facet_wrap(~Class)

  ggplot(Sdf, aes(Sensitivitylog10)) + geom_histogram() + GGPoptions +
    ylab("Number of species") + xlab("Climate change sensitivity (log10)") +  
    scale_fill_viridis_d()+ scale_color_viridis_d() + facet_wrap(~Class, ncol=4)+
    theme(panel.margin = unit(0, "lines")) 

ggplot(Sdf, aes(Sensitivitylog10, col=Class, fill=Class)) +
    geom_density(alpha=0.3, adjust=4) + GGPoptions +
    xlab("Climate change sensitivity (log10)") + ylab("Density of species") +
  scale_fill_viridis_d() + scale_color_viridis_d()
  
  
##############################################################################
## Variance explained by main effects  -- with residual variation

AnovaVar <- function(anova_model) {
  #anova_model <- anova_model[c(1:(nrow(anova_model)-1)),]
  Anova_ss <- anova_model$"Sum Sq"
  Anova_ss <- cbind(anova_model, PctExp=Anova_ss/sum(Anova_ss)*100)
  Anova_ss <- Anova_ss[order(Anova_ss$PctExp),]
  return(Anova_ss)
}

Anova_Amphibians <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Amphibians_ANOVA.rds")
Anova_Birds <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Birds_ANOVA.rds")
Anova_Mammals <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Mammals_ANOVA.rds")
Anova_Reptiles <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Reptiles_ANOVA.rds")

Anova_Amphibians <- AnovaVar(Anova_Amphibians)
Anova_Birds <- AnovaVar(Anova_Birds)
Anova_Mammals <- AnovaVar(Anova_Mammals)
Anova_Reptiles <- AnovaVar(Anova_Reptiles)

Anova_Amphibians$Class <- "Amphibians"
Anova_Birds$Class <- "Birds"
Anova_Mammals$Class <- "Mammals"
Anova_Reptiles$Class <- "Reptiles"

Anova_classes <- rbind(Anova_Amphibians, Anova_Birds, Anova_Mammals, Anova_Reptiles)
Anova_classes$Effects <- rownames(Anova_classes)

More_than_5 <- Anova_classes %>%
  group_by(Class) %>%
  filter(PctExp>=5) %>%
  dplyr::select(Class, PctExp, Effects)

Less_than_5 <- Anova_classes %>%
  group_by(Class) %>%
  filter(PctExp<5) %>%
  summarise(PctExp=sum(PctExp)) %>%
  mutate(Effects="Other effects")

PctExpClasses <- rbind(More_than_5, Less_than_5)
PctExpClasses$Effects <- c(rep(c("Residuals", "Range area"), 4), rep("Other effects", 4))

saveRDS(PctExpClasses, "E:/3.Explanatory_traits/Results/Results_climate_sensitivity/Variance_explained_withres.rds")

## factoring out range area and residual variation

Without_RA <- Anova_classes[!grepl("Range_area", rownames(Anova_classes)),]
Without_RA <- Without_RA[!grepl("Residuals", rownames(Without_RA)),]

colnames(Without_RA)[2] <- "SumSq"
Without_RA <- Without_RA %>% 
  group_by(Class) %>% 
  mutate(PctExp=SumSq/sum(SumSq)*100)

More_than_5 <- Without_RA %>%
  group_by(Class) %>%
  filter(PctExp>=5) %>%
  dplyr::select(Class, PctExp, Effects)

Less_than_5 <- Without_RA %>%
  group_by(Class) %>%
  filter(PctExp<5) %>%
  summarise(PctExp=sum(PctExp)) %>%
  mutate(Effects="Other effects")

PctExpClasses <- rbind(More_than_5, Less_than_5)
PctExpClasses$Effects <- c("Lifespan proxy", "Habitat breadth", "Body mass", "Litter/clutch size", 
                           "Primary diet", "Body mass", "Artificial habitat use", "Litter/clutch size",
                           "Primary diet", "Habitat breadth", "Litter/clutch size", "Body mass",
                           "Lifespan proxy", "Artificial habitat use", "Body mass", rep("Other effects",4))

saveRDS(PctExpClasses, "D:/3.Explanatory_traits/Results/Results_climate_sensitivity/Variance_explained_withoutres_withoutRA.rds")


## Variance explained by main effects  -- without residual variation

AnovaVar <- function(anova_model) {
  anova_model <- anova_model[c(1:(nrow(anova_model)-1)),]
  Anova_ss <- anova_model$"Sum Sq"
  Anova_ss <- cbind(anova_model, PctExp=Anova_ss/sum(Anova_ss)*100)
  Anova_ss <- Anova_ss[order(Anova_ss$PctExp),]
  return(Anova_ss)
}

Anova_Amphibians <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Amphibians_ANOVA.rds")
Anova_Birds <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Birds_ANOVA.rds")
Anova_Mammals <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Mammals_ANOVA.rds")
Anova_Reptiles <- readRDS("../Results/5.PGLS_models_results/ANOVA_diet_models/Reptiles_ANOVA.rds")

write.csv(as.data.frame(Anova_Amphibians), "../Results/5.PGLS_models_results/Species_over_100km/ANOVA_table_amphibians.csv", row.names = TRUE)
write.csv(as.data.frame(Anova_Birds), "../Results/5.PGLS_models_results/Species_over_100km/ANOVA_table_birds.csv", row.names = TRUE)
write.csv(as.data.frame(Anova_Mammals), "../Results/5.PGLS_models_results/Species_over_100km/ANOVA_table_mammals.csv", row.names = TRUE)
write.csv(as.data.frame(Anova_Reptiles), "../Results/5.PGLS_models_results/Species_over_100km/ANOVA_table_reptiles.csv", row.names = TRUE)

Anova_Amphibians <- AnovaVar(Anova_Amphibians)
Anova_Birds <- AnovaVar(Anova_Birds)
Anova_Mammals <- AnovaVar(Anova_Mammals)
Anova_Reptiles <- AnovaVar(Anova_Reptiles)

Anova_Amphibians$Class <- "Amphibians"
Anova_Birds$Class <- "Birds"
Anova_Mammals$Class <- "Mammals"
Anova_Reptiles$Class <- "Reptiles"

Anova_classes <- rbind(Anova_Amphibians, Anova_Birds, Anova_Mammals, Anova_Reptiles)
Anova_classes$Effects <- rownames(Anova_classes)

More_than_5 <- Anova_classes %>% 
  group_by(Class) %>%
  filter(PctExp>=5) %>% 
  dplyr::select(Class, PctExp, Effects)

Less_than_5 <- Anova_classes %>% 
  group_by(Class) %>%
  filter(PctExp<5) %>% 
  summarise(PctExp=sum(PctExp)) %>% 
  mutate(Effects="Other effects")

PctExpClasses <- rbind(More_than_5, Less_than_5)
PctExpClasses$Effects <- c(rep("Range area", 3), "Body mass", "Range area", rep("Other effects", 4))

saveRDS(PctExpClasses, "E:/3.Explanatory_traits/Results/Results_climate_sensitivity/Variance_explained_withoutres.rds")




















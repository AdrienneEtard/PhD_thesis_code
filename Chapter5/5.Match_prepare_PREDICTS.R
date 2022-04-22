## Match PREDICTS and metabolic rate data and prepare PREDICTS data

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

## for plotting

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

######################################################################################################################

## loading trait and resting metabolic rates data

Traits <- read.csv("../Results/4.BMR_data_residuals.csv")
Predicts <- readRDS("../Data/PredictsVertebrates.rds")

# add trophic levels to data
TL <- readRDS("../Data/Imputed_traits.rds")[[8]]
TL <- rbind(TL$M$Imputed.Dataset[c("Best_guess_binomial", "Trophic_level")],
            TL$B$Imputed.Dataset[c("Best_guess_binomial", "Trophic_level")],
            TL$R$Imputed.Dataset[c("Best_guess_binomial", "Trophic_level")],
            TL$A$Imputed.Dataset[c("Best_guess_binomial", "Trophic_level")])
colnames(TL)[1] <- "Species"
Traits <- left_join(Traits, TL)

# checking species number
length(unique(Traits$Species))
length(unique(Predicts$Best_guess_binomial))
colnames(Traits)[12] <- "Best_guess_binomial"

# subsetting data for relevant columns
TraitsSub <- Traits %>% 
  dplyr::select(Best_guess_binomial,
                log_BMR,
                Thermoregulation, 
                Trophic_level, 
                Body_mass_g,
                Residual_BMR_log_log)

# checking herbivore ectotherms (small sample sizes)
#Ecto_herb <- TraitsSub %>%  filter(Thermoregulation=="Ectotherms", Trophic_level=="Herbivore")

######################################################################################################################

## PREDICTS data

# All measurements in PREDICTS as presence/absence data
Predicts$Occurrence <- ifelse(Predicts$Measurement > 0, 1, 0)

# joining PREDICTS and traits
PredictsTraits <- left_join(Predicts, TraitsSub, by="Best_guess_binomial")
any(is.na(PredictsTraits$Trophic_level))
any(is.na(PredictsTraits$Thermoregulation))
any(is.na(PredictsTraits$log_BMR_standardised))

# checking RMR distribution
PredictsTraits$log_BMR %>%  hist()
exp(PredictsTraits$log_BMR) %>%  hist()
PredictsTraits$log_BMR %>%  median()

# checking residual RMR distribution
PredictsTraits$Residual_BMR_log_log %>%  hist()
PredictsTraits$Residual_BMR_log_log %>%  mean()
PredictsTraits$Residual_BMR_log_log %>%  median()
PredictsTraits$Residual_BMR_log_log %>%  range()

# PredictsTraits$Residual_BMR_log_log_rescaled <- scale(PredictsTraits$Residual_BMR_log_log, center = TRUE, scale = TRUE)
# PredictsTraits$Residual_BMR_log_log_rescaled %>%  hist()
# PredictsTraits$Residual_BMR_log_log %>%  hist()

######################################################################################################################

## adding annual mean temperature to PREDICTS

# load climate variables, at 2.5 arc-minute resolution (~4.6 km sq at the equator), and select those that are not too collinear
WCfiles <- paste("../Data/wc2.1_2.5m_bio", 
                 list.files(path = "../Data/wc2.1_2.5m_bio"), 
                 sep="/")
MeanAnnualTemp <- raster(WCfiles[[1]])
names(MeanAnnualTemp)

## mean annual temperature at PREDICTS sites
xy <- PredictsTraits %>% 
  dplyr::select(Longitude, Latitude)

result <- extract(MeanAnnualTemp, xy, cellnumbers = TRUE) %>% 
  as.data.frame()
min(result$wc2.1_2.5m_bio_1, na.rm=TRUE)
max(result$wc2.1_2.5m_bio_1,  na.rm=TRUE)
PredictsTraits$Annual_mean_temperature <- result$wc2.1_2.5m_bio_1

######################################################################################################################

## plotting RMR against body mass and residual RMRs against body mass (for conceptual figure)

Predicts_unique_species <- unique(PredictsTraits[,c("Residual_BMR_log_log",
                                                    "log_BMR",
                                                    "Body_mass_g",
                                                    "Trophic_level",
                                                    "Class", 
                                                    "Family", 
                                                    "Best_guess_binomial",
                                                    "Thermoregulation")])

Predicts_unique_species$Class <- factor(Predicts_unique_species$Class,
                                        levels = c("Amphibia",
                                                   "Aves",
                                                   "Mammalia",
                                                   "Reptilia"),
                                        labels=c("Amphibians",
                                                 "Birds",
                                                 "Mammals",
                                                 "Reptiles"))

## plot for residual RMRs against log body mass
p_residuals <- ggplot(Predicts_unique_species, 
                      aes(log(Body_mass_g), Residual_BMR_log_log)) +
  geom_point() + 
  GGPoptions + 
  geom_hline(yintercept = 0, lty="dashed", col="darkgrey", size=1.3) +
  xlim(c(-1,30)) +
  ylab("Residual RMR (mL O2/hour, log)") + xlab("Body mass (log, g)") +
  #scale_color_viridis_d(end=.7) + 
  theme(legend.position = "bottom") +
  geom_segment(x = 16, y = 0, xend = 16, yend = 2.2, col="red", size=1.3,
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_segment(x = 16, y = 0, xend = 16, yend = -2.2, col="blue", size=1.3,
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_text(x=17, y=1, label="Positive deviations:
            \n RMR higher than expected \n from body mass", col="red", hjust=0)+
  geom_text(x=17, y=-1, label="Negative deviations:
            \n RMR lower than expected \n from body mass", col="blue", hjust=0)

# ggplot(Predicts_unique_species,
#        aes(log(Body_mass_g), Residual_BMR_log_log)) +
#   geom_point() + 
#   GGPoptions + facet_wrap(~Class)


## plot for residual RMRs against log body mass
p_rmr <-ggplot(Predicts_unique_species, aes(log(Body_mass_g), log_BMR, col=Class)) +
  geom_point() + GGPoptions + 
  ylab("RMR (mL O2/hour, log)") + xlab("Body mass (log, g)") +
  scale_color_viridis_d(end=.8) + 
  theme(legend.position=c(0.835,0.20),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = "grey"),
        legend.title = element_blank())

library(patchwork)
#p_framework <- ggarrange(p_rmr, p_residuals)
#p_rmr + plot_spacer() + p_residuals

p_framework <- (p_rmr + theme(plot.margin = unit(c(0,60,0,0), "pt"))) +
  (p_residuals + theme(plot.margin = unit(c(0,0,0,60), "pt")))

ggsave(p_framework, filename="../Results/p_framework.png", width=12, height=4)
ggsave(p_framework, filename="../Results/p_framework.pdf", width=12, height=4)


# plot(NULL, xlim=c(-1,1), 
#      ylim=c(-1,1), yaxt='n', xaxt='n', 
#      xlab="Residual RMR", ylab="Occurrence probability")
# abline(a=0,b=-1, col="red")
# abline(a=0,b=-0.2, col="blue")


## looking at the species with the most extreme values in residual RMR
Predicts_unique_species %>%
  group_by(Class) %>% 
  summarise(Min=min(Residual_BMR_log_log), 
            Max=max(Residual_BMR_log_log))


GetMin <- function(Predicts_species, VClass){
  PredSub <- subset(Predicts_species, Predicts_species$Class==VClass)
  Min <- min(PredSub$Residual_BMR_log_log)
  cat("Residual RMR is ", Min, "\n")
  BM <- PredSub$Body_mass_g[PredSub$Residual_BMR_log_log==Min]
  cat("Body mass is ", BM, "g\n")
  return(PredSub$Best_guess_binomial[PredSub$Residual_BMR_log_log==Min])
}

GetMax <- function(Predicts_species, VClass){
  PredSub <- subset(Predicts_species, Predicts_species$Class==VClass)
  Max <- max(PredSub$Residual_BMR_log_log)
  cat("Residual RMR is ", Max, "\n")
  BM <- PredSub$Body_mass_g[PredSub$Residual_BMR_log_log==Max]
  cat("Body mass is ", BM, "g\n")
  return(PredSub$Best_guess_binomial[PredSub$Residual_BMR_log_log==Max])
}

GetMin(Predicts_unique_species, "Amphibians")
GetMax(Predicts_unique_species, "Amphibians")

GetMin(Predicts_unique_species, VClass="Birds")
GetMax(Predicts_unique_species, "Birds")

GetMin(Predicts_unique_species, "Mammals") # "Gracilinanus microtarsus"
GetMax(Predicts_unique_species, "Mammals")

GetMin(Predicts_unique_species, "Reptiles") 
GetMax(Predicts_unique_species, "Reptiles")


Predicts_unique_species$Best_guess_binomial[Predicts_unique_species$Residual_BMR_log_log==min(Predicts_unique_species$Residual_BMR_log_log)]
Predicts_unique_species$Best_guess_binomial[Predicts_unique_species$Residual_BMR_log_log==max(Predicts_unique_species$Residual_BMR_log_log)]

## min = "Oophaga pumilio"
## max = "Phyllobates lugubris"


#################################################################################################################################

## setting levels in PREDICTS for land use and land-use intensity

# use intensity
PredictsTraits$Use_intensity <- as.character(PredictsTraits$Use_intensity)
PredictsTraits$Use_intensity %>%  unique()
PredictsTraits$Use_intensity <- ifelse(PredictsTraits$Use_intensity=="Cannot decide", NA, PredictsTraits$Use_intensity)
PredictsTraits$Use_intensity <- factor(PredictsTraits$Use_intensity, levels=c("Minimal use", "Light use", "Intense use"))
PredictsTraits$Use_intensity %>%  levels()

# predominant land use with secondary vegetation grouped
PredictsTraits$Predominant_land_use <- as.character(PredictsTraits$Predominant_land_use)
PredictsTraits$LandUse <- PredictsTraits$Predominant_land_use 
PredictsTraits$LandUse[PredictsTraits$LandUse=="Cannot decide"] <- NA 

PredictsTraits$LandUse %>%  unique()
PredictsTraits$LandUse[grepl("econdary vegetation", PredictsTraits$LandUse)] <- "Secondary vegetation"
PredictsTraits$LandUse <- factor(PredictsTraits$LandUse,
                                              levels=c("Primary vegetation",
                                                       "Secondary vegetation",
                                                       "Plantation forest",
                                                       "Pasture",
                                                       "Cropland",
                                                       "Urban"))
PredictsTraits$LandUse %>%  levels()

# predominant land use without secondary vegetation grouped
PredictsTraits$Predominant_land_use %>%  unique
PredictsTraits$Predominant_land_use[PredictsTraits$Predominant_land_use=="Cannot decide"] <- NA 

PredictsTraits$Predominant_land_use <- factor(PredictsTraits$Predominant_land_use,
                                 levels=c("Primary vegetation",
                                          "Mature secondary vegetation",
                                          "Intermediate secondary vegetation",
                                          "Young secondary vegetation",
                                          "Plantation forest",
                                          "Pasture",
                                          "Cropland",
                                          "Urban"))

PredictsTraits$Predominant_land_use %>%  levels()

saveRDS(PredictsTraits, "../Results/PredictsTraits.rds")

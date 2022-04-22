## Model validations (for the models with class as an effect)

library(StatisticalModels)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(viridis)
library(scales)
source("Functions.R")

##############################################################################################################
## Within vertebrate classes: robustness to variation in imputed values

Reptiles_8 <- readRDS("../../Results/dbFD_indices/8_setsReptiles.rds")
Mammals_8 <- readRDS("../../Results/dbFD_indices/8_sets_Mammals.rds")
Amphibians_8 <- readRDS("../../Results/dbFD_indices/8_setsAmphibians.rds")
Birds_8 <- readRDS("../../Results/dbFD_indices/8_setsBirds.rds")

Reptiles_8 <- lapply(Reptiles_8, Clean)
Mammals_8 <- lapply(Mammals_8, Clean)
Amphibians_8 <- lapply(Amphibians_8, Clean)
Birds_8 <- lapply(Birds_8, Clean)

Reptiles_8 <- lapply(Reptiles_8, function(x){
  x$FRic[x$SR==1] <- NA
  x$FDis[x$SR==1] <- NA
  # tranform values
  x$asin_sqrt_FRic <- asin(sqrt(x$FRic))
  x$asin_sqrt_FDis <- asin(sqrt(x$FDis))
  # # filter out NA use intensity
  # x <- x %>%
  #   filter(!is.na(Use_intensity))
  # group land uses
  x$Predominant_land_use <- as.character(x$Predominant_land_use)
  x$Predominant_land_use[x$Predominant_land_use %in% c("Mature secondary vegetation",
                                                       "Intermediate secondary vegetation",
                                                       "Young secondary vegetation")]  <- "Secondary vegetation"
  x$Predominant_land_use[x$Predominant_land_use %in% c("Cropland", "Pasture")]  <- "Agricultural"
  x$Predominant_land_use %>%  unique()
  x$Predominant_land_use<- factor(x$Predominant_land_use, levels=c("Primary vegetation",
                                                                   "Secondary vegetation",
                                                                   "Plantation forest",
                                                                   "Agricultural",
                                                                   "Urban"))


  return(x)
})
Mammals_8 <-  lapply(Mammals_8, function(x){
  x$FRic[x$SR==1] <- NA
  x$FDis[x$SR==1] <- NA
  # tranform values
  x$asin_sqrt_FRic <- asin(sqrt(x$FRic))
  x$asin_sqrt_FDis <- asin(sqrt(x$FDis))
  # # filter out NA use intensity
  # x <- x %>%
  #   filter(!is.na(Use_intensity))
  # group land uses
  x$Predominant_land_use <- as.character(x$Predominant_land_use)
  x$Predominant_land_use[x$Predominant_land_use %in% c("Mature secondary vegetation",
                                                       "Intermediate secondary vegetation",
                                                       "Young secondary vegetation")]  <- "Secondary vegetation"
  x$Predominant_land_use[x$Predominant_land_use %in% c("Cropland", "Pasture")]  <- "Agricultural"
  x$Predominant_land_use %>%  unique()
  x$Predominant_land_use<- factor(x$Predominant_land_use, levels=c("Primary vegetation",
                                                                   "Secondary vegetation",
                                                                   "Plantation forest",
                                                                   "Agricultural",
                                                                   "Urban"))


  return(x)
})
Amphibians_8 <- lapply(Amphibians_8, function(x){
  x$FRic[x$SR==1] <- NA
  x$FDis[x$SR==1] <- NA
  # tranform values
  x$asin_sqrt_FRic <- asin(sqrt(x$FRic))
  x$asin_sqrt_FDis <- asin(sqrt(x$FDis))
  # # filter out NA use intensity
  # x <- x %>%
  #   filter(!is.na(Use_intensity))
  # group land uses
  x$Predominant_land_use <- as.character(x$Predominant_land_use)
  x$Predominant_land_use[x$Predominant_land_use %in% c("Mature secondary vegetation",
                                                       "Intermediate secondary vegetation",
                                                       "Young secondary vegetation")]  <- "Secondary vegetation"
  x$Predominant_land_use[x$Predominant_land_use %in% c("Cropland", "Pasture")]  <- "Agricultural"
  x$Predominant_land_use %>%  unique()
  x$Predominant_land_use<- factor(x$Predominant_land_use, levels=c("Primary vegetation",
                                                                   "Secondary vegetation",
                                                                   "Plantation forest",
                                                                   "Agricultural",
                                                                   "Urban"))


  return(x)
})
Birds_8 <- lapply(Birds_8, function(x){
  x$FRic[x$SR==1] <- NA
  x$FDis[x$SR==1] <- NA
  # tranform values
  x$asin_sqrt_FRic <- asin(sqrt(x$FRic))
  x$asin_sqrt_FDis <- asin(sqrt(x$FDis))
  # # filter out NA use intensity
  # x <- x %>%
  #   filter(!is.na(Use_intensity))
  # group land uses
  x$Predominant_land_use <- as.character(x$Predominant_land_use)
  x$Predominant_land_use[x$Predominant_land_use %in% c("Mature secondary vegetation",
                                                       "Intermediate secondary vegetation",
                                                       "Young secondary vegetation")]  <- "Secondary vegetation"
  x$Predominant_land_use[x$Predominant_land_use %in% c("Cropland", "Pasture")]  <- "Agricultural"
  x$Predominant_land_use %>%  unique()
  x$Predominant_land_use<- factor(x$Predominant_land_use, levels=c("Primary vegetation",
                                                                   "Secondary vegetation",
                                                                   "Plantation forest",
                                                                   "Agricultural",
                                                                   "Urban"))


  return(x)
})

# check levels
Reptiles_8[[1]]$Predominant_land_use
Reptiles_8[[1]]$Use_intensity

# add tropical versus temperate divide
# define realms as tropical versus temperate (depending on difference in latitude)
PredictsRealmBiome <- readRDS("../../Results/Predicts_merged_sites.rds") %>%
  dplyr::select(SSBS, Realm, Biome, Longitude, Latitude)
PredictsRealmBiome <- unique(PredictsRealmBiome)
PredictsRealmBiome$Biome <- ifelse(abs(PredictsRealmBiome$Latitude)<=23.5, "Tropical", "Temperate")
# a few sites don't have longitude or latitude - resolve manually
PredictsRealmBiome$Biome[PredictsRealmBiome$SSBS %in%
                           c("DB1_2010__Garden 1  1",
                             "DB1_2010__Garden 1  10",
                             "DB1_2010__Garden 1  29")] <- "Tropical"

PredictsRealmBiome$Biome[PredictsRealmBiome$SSBS %in%
                           c("SE1_2011__Rosselli 1  1",
                             "SE1_2011__Rosselli 1  4",
                             "SE1_2011__Rosselli 1  2")] <- "Tropical"

PredictsRealmBiome$Biome[PredictsRealmBiome$SSBS %in%
                           c("SC1_2008__Eigenbrod 1  1",
                             "SC1_2008__Eigenbrod 1  2",
                             "SC1_2008__Eigenbrod 1  29")] <- "Temperate"


Reptiles_8 <- lapply(Reptiles_8, function(x){
  Y <- match(x$SSBS, PredictsRealmBiome$SSBS)
  x$Realm <- PredictsRealmBiome$Realm[Y]
  x$Biome <- PredictsRealmBiome$Biome[Y]
  return(x)
})
Amphibians_8 <- lapply(Amphibians_8, function(x){
  Y <- match(x$SSBS, PredictsRealmBiome$SSBS)
  x$Realm <- PredictsRealmBiome$Realm[Y]
  x$Biome <- PredictsRealmBiome$Biome[Y]
  return(x)
})
Mammals_8 <- lapply(Mammals_8, function(x){
  Y <- match(x$SSBS, PredictsRealmBiome$SSBS)
  x$Realm <- PredictsRealmBiome$Realm[Y]
  x$Biome <- PredictsRealmBiome$Biome[Y]
  return(x)
})
Birds_8 <- lapply(Birds_8, function(x){
  Y <- match(x$SSBS, PredictsRealmBiome$SSBS)
  x$Realm <- PredictsRealmBiome$Realm[Y]
  x$Biome <- PredictsRealmBiome$Biome[Y]
  return(x)
})

## binding together
Classes_set <- list()
for(i in 1:8){
  Amphibians_8[[i]]$Class <- "Amphibians"
  Birds_8[[i]]$Class <- "Birds"
  Mammals_8[[i]]$Class <- "Mammals"
  Reptiles_8[[i]]$Class <- "Reptiles"
  Classes_set[[i]] <- rbind(Amphibians_8[[i]], Birds_8[[i]], Mammals_8[[i]], Reptiles_8[[i]])
}


## running models 2a and 2b and predictions for each set of metrics

Models_2 <- function(set, N) {

  model2a <- lme4::lmer(asin_sqrt_FRic~Predominant_land_use+Use_intensity+Biome+Class+
                          Predominant_land_use:Use_intensity+
                          Predominant_land_use:Class+
                          Use_intensity:Biome+
                          Class:Biome+
                          (1|SS/SSB), data=set)

  model2b <- lme4::lmer(asin_sqrt_FDis~Predominant_land_use+Use_intensity+Class+
                          Predominant_land_use:Use_intensity+
                          Predominant_land_use:Class+
                          Class:Use_intensity+
                          (1|SS/SSB), data=set)

  print("predictions")

  Newdata <- expand.grid(Predominant_land_use=levels(set$Predominant_land_use),
                         Use_intensity=levels(set$Use_intensity),
                         Biome=c("Temperate", "Tropical"),
                         Class=c("Amphibians", "Birds", "Mammals", "Reptiles"),
                         asin_sqrt_FDis=0,
                         asin_sqrt_FRic=0)

  Predict_effects <- function(Newdata, Model, rescale) {

    preds <- sapply(X=1:1000, FUN=function(i){

      coefs <- mvrnorm(n = 1, mu = fixef(object=Model), Sigma = vcov(object=Model))
      mm <- model.matrix(terms(Model), Newdata)

      # drop coefs that couldn't be estimated
      print(setdiff(colnames(mm), names(coefs)))
      to_drop <- print(setdiff(colnames(mm), names(coefs)))
      mm <- as.data.frame(mm)
      mm <- mm[, -which(colnames(mm) %in% to_drop)]
      mm <- as.matrix(mm)

      y <- mm%*%coefs

      # backtranforming
      y <- sin(y)^2

      # rescaling
      if(rescale){
        # initialisation
        seq <- 1:5
        y[seq] <- y[seq]/y[seq[1]]*100
        # for loop to rescale all values
        for(i in 1:(nrow(Newdata)/5-1)){
          seq <- seq + 5
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

  # predictions for FRic
  preds_fric <- Predict_effects(Newdata, model2a, TRUE)
  preds_fric$Median[preds_fric$Predominant_land_use_2=="Urban" & preds_fric$Class=="Reptiles"] <- NA
  preds_fric$Lower[preds_fric$Predominant_land_use_2=="Urban" & preds_fric$Class=="Reptiles"] <- NA
  preds_fric$Upper[preds_fric$Predominant_land_use_2=="Urban" & preds_fric$Class=="Reptiles"] <- NA

  preds_fric$Median <- preds_fric$Median-100
  preds_fric$Lower <- preds_fric$Lower-100
  preds_fric$Upper <- preds_fric$Upper-100

  # predictions for FDis
  preds_fdis <- Predict_effects(Newdata, model2b, TRUE)
  preds_fdis$Median[preds_fdis$Predominant_land_use_2=="Urban" & preds_fdis$Class=="Reptiles"] <- NA
  preds_fdis$Lower[preds_fdis$Predominant_land_use_2=="Urban" & preds_fdis$Class=="Reptiles"] <- NA
  preds_fdis$Upper[preds_fdis$Predominant_land_use_2=="Urban" & preds_fdis$Class=="Reptiles"] <- NA

  preds_fdis$Median <- preds_fdis$Median-100
  preds_fdis$Lower <- preds_fdis$Lower-100
  preds_fdis$Upper <- preds_fdis$Upper-100

  preds_fric$Metric <- "FRic"
  preds_fdis$Metric <- "FDis"

  preds_fric$set <- N
  preds_fdis$set <- N

  return(rbind(preds_fric, preds_fdis))
}

Results_8_models <- list()

for(i in 1:8) {
  Results_8_models[[i]] <- Models_2(Classes_set[[i]],i)
}

Results_8_models <- data.table::rbindlist(Results_8_models)


## plotting (Supplementary Information, Figure S23)
Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Agricultural",
            "Urban")

Labels=c("PV", "SV", "PF", "AGR", "UR")

pFRic <- ggplot(Results_8_models[Results_8_models$Metric=="FRic",],
                aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, col=Use_intensity, group=interaction(set,Use_intensity))) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.8), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.8)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=c("black", "red", "orange")) +
  GGPoptions +
  facet_grid(Biome~Class)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

pFDis <- ggplot(Results_8_models[Results_8_models$Metric=="FDis" & Results_8_models$Biome=="Temperate",],
                aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, col=Use_intensity, group=interaction(set,Use_intensity))) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.8), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.8)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=c("black", "red", "orange")) +
  GGPoptions +
  facet_grid(~Class)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(b) FDis (Temperate; similar effects for tropical)") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

p_fig_SI_23 <- ggarrange(pFRic, pFDis, nrow = 2, heights = c(0.60, 0.40))
ggsave(p_fig_SI_23, filename="../Revisions/PDF_figures/Figure_SI_23.pdf",
       height=8, width=10)

rm(pFRic, pFDis)

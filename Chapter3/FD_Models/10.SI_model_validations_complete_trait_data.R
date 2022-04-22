#####################################################################################################
## trait data and gower distance matrices for the validations

source("Functions.R")
library(dplyr)
library(FD)
library(usedist)
library(beepr)

## across all vertebrates: complete trait data  + complete trait data and geographical range size
## also by class and for herptiles

Amphibians <- read.csv("../../Data/GlobalGaps_traitdata/Amphibians.csv")
Birds <- read.csv("../../Data/GlobalGaps_traitdata/Birds.csv")
Mammals <- read.csv("../../Data/GlobalGaps_traitdata/Mammals.csv")
Reptiles <- read.csv("../../Data/GlobalGaps_traitdata/Reptiles.csv")

## select the traits in each class
Amphibians <- Amphibians %>%
  select(Best_guess_binomial,
         Body_length_mm,
         Litter_size,
         Trophic_level,
         Habitat_breadth_IUCN,
         Specialisation,
         Diel_activity,
         Maturity_d)

colnames(Amphibians)[c(2, 8)] <- c("Body_mass_g", "Generation_length_d")

Mammals <- Mammals %>%
  select(Best_guess_binomial,
         Body_mass_g,
         Litter_size,
         Trophic_level,
         Habitat_breadth_IUCN,
         Specialisation,
         Diel_activity,
         Generation_length_d)

Birds <- Birds %>%
  select(Best_guess_binomial,
         Body_mass_g,
         Litter_size,
         Trophic_level,
         Habitat_breadth_IUCN,
         Specialisation,
         Diel_activity,
         Generation_length_d)

Reptiles <- Reptiles %>%
  select(Best_guess_binomial,
         Body_mass_g,
         Litter_size,
         Trophic_level,
         Habitat_breadth_IUCN,
         Specialisation,
         Diel_activity,
         Longevity_d)

colnames(Reptiles)[c(8)] <- "Generation_length_d"


## add geographical range sizes
Ranges_sizes <- read.csv("../../Data/RangeSizes.csv")
Ranges_sizes$RangeSize_sqkm_AfterCuts[Ranges_sizes$RangeSize_sqkm_AfterCuts==0] <- NA
colnames(Ranges_sizes)[2] <- "Best_guess_binomial"
Ranges_sizes <- Ranges_sizes %>%
  dplyr::select(Best_guess_binomial, RangeSize_sqkm_AfterCuts)

Amphibians <- left_join(Amphibians, Ranges_sizes)
Birds <- left_join(Birds, Ranges_sizes)
Mammals <- left_join(Mammals, Ranges_sizes)
Reptiles <- left_join(Reptiles, Ranges_sizes)

## get species with complete trait data

Completeness <- function(data){

  data1 <- data %>%
    select(-Best_guess_binomial, -RangeSize_sqkm_AfterCuts)

  data2 <- data %>%
    select(-Best_guess_binomial)

  data$completeness_no_rs <- apply(data1, 1, FUN = function(x){
    return(length(which(!is.na(x)))/length(x)*100)
  })

  data$completeness_rs <- apply(data2, 1, FUN = function(x){
    return(length(which(!is.na(x)))/length(x)*100)
  })
  return(data)
}

Amphibians <- Completeness(Amphibians)
Birds <- Completeness(Birds)
Mammals <- Completeness(Mammals)
Reptiles <- Completeness(Reptiles)

## across vertebrates

complete_trait_data_no_rs <- rbind(Amphibians[Amphibians$completeness_no_rs==100,],
                                   Birds[Birds$completeness_no_rs==100,],
                                   Mammals[Mammals$completeness_no_rs==100,],
                                   Reptiles[Reptiles$completeness_no_rs==100,])

complete_trait_data_rs <- rbind(Amphibians[Amphibians$completeness_rs==100,],
                                   Birds[Birds$completeness_rs==100,],
                                   Mammals[Mammals$completeness_rs==100,],
                                   Reptiles[Reptiles$completeness_rs==100,])

## for each class
Complete_amphibians <- Amphibians %>%
  filter(completeness_no_rs==100)

Complete_birds <- Birds %>%
  filter(completeness_no_rs==100)

Complete_mammals <- Mammals %>%
  filter(completeness_no_rs==100)

Complete_reptiles <- Reptiles %>%
  filter(completeness_no_rs==100)

Complete_herptiles <- rbind(Complete_amphibians, Complete_reptiles)

## log10 transform and rescale

## all verts, no range size
complete_trait_data_no_rs$log10_Body_mass_g <- scale(log10(complete_trait_data_no_rs$Body_mass_g), center = TRUE, scale = TRUE) %>% as.numeric()
complete_trait_data_no_rs$log10_Litter_size <- scale(log10(complete_trait_data_no_rs$Litter_size), center = TRUE, scale = TRUE)%>% as.numeric()
complete_trait_data_no_rs$log10_Lifespan <- scale(log10(complete_trait_data_no_rs$Generation_length_d), center = TRUE, scale = TRUE)%>% as.numeric()
complete_trait_data_no_rs$sqrt_Habitat_breadth_IUCN <- scale(sqrt(complete_trait_data_no_rs$Habitat_breadth_IUCN), center = TRUE, scale = TRUE)%>% as.numeric()

## all verts, range sizes
complete_trait_data_rs$log10_Body_mass_g <- scale(log10(complete_trait_data_rs$Body_mass_g), center = TRUE, scale = TRUE)%>% as.numeric()
complete_trait_data_rs$log10_Litter_size <- scale(log10(complete_trait_data_rs$Litter_size), center = TRUE, scale = TRUE)%>% as.numeric()
complete_trait_data_rs$log10_Lifespan <- scale(log10(complete_trait_data_rs$Generation_length_d), center = TRUE, scale = TRUE)%>% as.numeric()
complete_trait_data_rs$sqrt_Habitat_breadth_IUCN <- scale(sqrt(complete_trait_data_rs$Habitat_breadth_IUCN), center = TRUE, scale = TRUE)%>% as.numeric()
complete_trait_data_rs$log10_rangesize <- scale(log10(complete_trait_data_rs$RangeSize_sqkm_AfterCuts), center = TRUE, scale = TRUE)%>% as.numeric()

## amphibians
Complete_amphibians$log10_Body_mass_g <- scale(log10(Complete_amphibians$Body_mass_g), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_amphibians$log10_Litter_size <- scale(log10(Complete_amphibians$Litter_size), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_amphibians$log10_Lifespan <- scale(log10(Complete_amphibians$Generation_length_d), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_amphibians$sqrt_Habitat_breadth_IUCN <- scale(sqrt(Complete_amphibians$Habitat_breadth_IUCN), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_amphibians$log10_rangesize <- scale(log10(Complete_amphibians$RangeSize_sqkm_AfterCuts), center = TRUE, scale = TRUE)%>% as.numeric()

## birds
Complete_birds$log10_Body_mass_g <- scale(log10(Complete_birds$Body_mass_g), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_birds$log10_Litter_size <- scale(log10(Complete_birds$Litter_size), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_birds$log10_Lifespan <- scale(log10(Complete_birds$Generation_length_d), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_birds$sqrt_Habitat_breadth_IUCN <- scale(sqrt(Complete_birds$Habitat_breadth_IUCN), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_birds$log10_rangesize <- scale(log10(Complete_birds$RangeSize_sqkm_AfterCuts), center = TRUE, scale = TRUE)%>% as.numeric()

## mammals
Complete_mammals$log10_Body_mass_g <- scale(log10(Complete_mammals$Body_mass_g), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_mammals$log10_Litter_size <- scale(log10(Complete_mammals$Litter_size), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_mammals$log10_Lifespan <- scale(log10(Complete_mammals$Generation_length_d), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_mammals$sqrt_Habitat_breadth_IUCN <- scale(sqrt(Complete_mammals$Habitat_breadth_IUCN), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_mammals$log10_rangesize <- scale(log10(Complete_mammals$RangeSize_sqkm_AfterCuts), center = TRUE, scale = TRUE)%>% as.numeric()

## reptiles
Complete_reptiles$log10_Body_mass_g <- scale(log10(Complete_reptiles$Body_mass_g), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_reptiles$log10_Litter_size <- scale(log10(Complete_reptiles$Litter_size), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_reptiles$log10_Lifespan <- scale(log10(Complete_reptiles$Generation_length_d), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_reptiles$sqrt_Habitat_breadth_IUCN <- scale(sqrt(Complete_reptiles$Habitat_breadth_IUCN), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_reptiles$log10_rangesize <- scale(log10(Complete_reptiles$RangeSize_sqkm_AfterCuts), center = TRUE, scale = TRUE)%>% as.numeric()

## herptiles
Complete_herptiles$log10_Body_mass_g <- scale(log10(Complete_herptiles$Body_mass_g), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_herptiles$log10_Litter_size <- scale(log10(Complete_herptiles$Litter_size), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_herptiles$log10_Lifespan <- scale(log10(Complete_herptiles$Generation_length_d), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_herptiles$sqrt_Habitat_breadth_IUCN <- scale(sqrt(Complete_herptiles$Habitat_breadth_IUCN), center = TRUE, scale = TRUE)%>% as.numeric()
Complete_herptiles$log10_rangesize <- scale(log10(Complete_herptiles$RangeSize_sqkm_AfterCuts), center = TRUE, scale = TRUE)%>% as.numeric()


## Getting Gower distances (subsetted for PREDICTS species)

Predicts <- readRDS("../../Results/Predicts_merged_sites.rds")
verts <- unique(Predicts$Best_guess_binomial)
amph <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Amphibia"])
rep <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Reptilia"])
bir <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Aves"])
mam <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Mammalia"])

## filter out PREDICTS species that are not represented in the complete trait data
verts_no__rs <- verts[verts %in% complete_trait_data_no_rs$Best_guess_binomial]
verts_rs <- verts[verts %in% complete_trait_data_rs$Best_guess_binomial]

amph <- amph[amph %in% Complete_amphibians$Best_guess_binomial]
rep <- rep[rep %in% Complete_reptiles$Best_guess_binomial]
bir <- bir[bir %in% Complete_birds$Best_guess_binomial]
mam <- mam[mam %in% Complete_mammals$Best_guess_binomial]

gc()
memory.limit(size=50000)

Selected_traits <- c("log10_Body_mass_g",
                     "log10_Lifespan",
                     "log10_Litter_size",
                     "sqrt_Habitat_breadth_IUCN",
                     "Specialisation",
                     "Diel_activity",
                     "Trophic_level")

Selected_traits_rs <- c("log10_Body_mass_g",
                     "log10_Lifespan",
                     "log10_Litter_size",
                     "sqrt_Habitat_breadth_IUCN",
                     "Specialisation",
                     "Diel_activity",
                     "Trophic_level",
                     "log10_rangesize")

## Gower across complete verts without rs
rownames(complete_trait_data_no_rs) <- as.character(complete_trait_data_no_rs$Best_guess_binomial)
Gower_complete_verts <- gowdis(complete_trait_data_no_rs[,Selected_traits])
Gower_complete_verts <- dist_subset(Gower_complete_verts, verts_no__rs)
saveRDS(Gower_complete_verts, "../../Results/Gower_distances/Complete_vertebrates.rds")

## Gower across complete verts with rs
rownames(complete_trait_data_rs) <- as.character(complete_trait_data_rs$Best_guess_binomial)
Gower_complete_verts_rs <- gowdis(complete_trait_data_rs[,Selected_traits_rs])
Gower_complete_verts_rs <- dist_subset(Gower_complete_verts_rs, verts_rs)
saveRDS(Gower_complete_verts_rs, "../../Results/Gower_distances/Complete_vertebrates_with_rangesize.rds")

## Gower across complete amphibians
rownames(Complete_amphibians) <- as.character(Complete_amphibians$Best_guess_binomial)
Gower_complete_amph <- gowdis(Complete_amphibians[,Selected_traits])
Gower_complete_amph <- dist_subset(Gower_complete_amph, amph)
saveRDS(Gower_complete_amph, "../../Results/Gower_distances/Complete_amphibians.rds")

## Gower across complete birds
rownames(Complete_birds) <- as.character(Complete_birds$Best_guess_binomial)
Gower_complete_bir <- gowdis(Complete_birds[,Selected_traits])
Gower_complete_bir <- dist_subset(Gower_complete_bir, bir)
saveRDS(Gower_complete_bir, "../../Results/Gower_distances/Complete_birds.rds")

## Gower across complete mammals
rownames(Complete_mammals) <- as.character(Complete_mammals$Best_guess_binomial)
Gower_complete_mam <- gowdis(Complete_mammals[,Selected_traits])
Gower_complete_mam <- dist_subset(Gower_complete_mam, mam)
saveRDS(Gower_complete_mam, "../../Results/Gower_distances/Complete_mammals.rds")

## Gower across complete reptiles
rownames(Complete_reptiles) <- as.character(Complete_reptiles$Best_guess_binomial)
Gower_complete_rep <- gowdis(Complete_reptiles[,Selected_traits])
Gower_complete_rep <- dist_subset(Gower_complete_rep, rep)
saveRDS(Gower_complete_rep, "../../Results/Gower_distances/Complete_reptiles.rds")

## Gower across complete herptiles
rownames(Complete_herptiles) <- as.character(Complete_herptiles$Best_guess_binomial)
Gower_complete_herp <- gowdis(Complete_herptiles[,Selected_traits])
Gower_complete_herp <- dist_subset(Gower_complete_herp, c(amph,rep))
saveRDS(Gower_complete_herp, "../../Results/Gower_distances/Complete_herptiles.rds")


#####################################################################################################
## compute FD for thesetaxonomic groups with complete trait data
source("Functions.R")
source("Function_to_prepare_Predicts_communities.R")
library(dplyr)
library(FD)
library(usedist)
library(beepr)

Predicts <- readRDS("../../Results/Predicts_merged_sites.rds")
Predicts_site_info <- read.csv("../../Results/Predicts_site_info_ER.csv")
Reg <- read.csv("../../Results/dbFD_indices_for_analysis/dbFD_region_std.csv")

#####################################################################################################
## validation 1: across vertebrates with complete trait data

Gower_complete_verts <- readRDS("../../Results/Gower_distances/Complete_vertebrates.rds")
verts <- intersect(labels(Gower_complete_verts), Predicts$Best_guess_binomial)
Predicts_complete_verts <- subset(Predicts, Best_guess_binomial %in% verts)

# get community matrices
Complete_verts_comm <-  Generate_communities(Predicts_complete_verts, FALSE)

# compute FD
Complete_verts <- To_run(List_of_communities = na.omit.list(Complete_verts_comm),
                  Gowerdist = Gower_complete_verts,
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = TRUE,
                  Predicts_site_info=Predicts_site_info)

Complete_verts <- Extract(Complete_verts, Predicts_site_info)
write.csv(Complete_verts, "../../Results/dbFD_indices/Complete_verts_validation.csv", row.names = FALSE)

Complete_verts <- read.csv( "../../Results/dbFD_indices/Complete_verts_validation.csv")


## fitting models 1a and 1b

Complete_verts <- Clean(Complete_verts)
Complete_verts$FRic[Complete_verts$SR==1] <- NA
Complete_verts$FDis[Complete_verts$SR==1] <- NA
# check levels
Complete_verts$Predominant_land_use
Complete_verts$Use_intensity
# add tropical versus temperate divide
Temperate <- Reg$SS[Reg$Biome=="Temperate"] %>%  unique()
Tropical <- Reg$SS[Reg$Biome=="Tropical"] %>%  unique()
Complete_verts$Biome <- NA
Complete_verts$Biome[Complete_verts$SS %in% Temperate] <- "Temperate"
Complete_verts$Biome[Complete_verts$SS %in% Tropical] <- "Tropical"

# tranform values
Complete_verts$asin_sqrt_FRic <- asin(sqrt(Complete_verts$FRic))
Complete_verts$asin_sqrt_FDis <- asin(sqrt(Complete_verts$FDis))

model1a <- lme4::lmer(asin_sqrt_FRic~Predominant_land_use+Use_intensity+Biome+
                        Predominant_land_use:Use_intensity+
                        Predominant_land_use:Biome+
                        (1|SS/SSB), data=Complete_verts)

model1b <- lme4::lmer(asin_sqrt_FDis~Predominant_land_use+Use_intensity+Biome+
                        Predominant_land_use:Use_intensity+
                        (1|SS/SSB), data=Complete_verts)

## predictions

Newdata <- expand.grid(Predominant_land_use=levels(Complete_verts$Predominant_land_use),
                       Use_intensity=levels(Complete_verts$Use_intensity),
                       Biome=c("Temperate", "Tropical"),
                       asin_sqrt_FDis=0,
                       asin_sqrt_FRic=0)

levels(Newdata$Predominant_land_use)

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
      y[1:8] <- y[1:8]/y[1]*100
      y[9:16] <- y[9:16]/y[9]*100
      y[17:24] <- y[17:24]/y[17]*100
      y[25:32] <- y[25:32]/y[25]*100
      y[33:40] <- y[33:40]/y[33]*100
      y[41:48] <- y[41:48]/y[41]*100
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
preds_fric <- Predict_effects(Newdata, model1a, TRUE)
preds_fric$Median[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA
preds_fric$Lower[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA
preds_fric$Upper[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA

preds_fric$Median <- preds_fric$Median-100
preds_fric$Lower <- preds_fric$Lower-100
preds_fric$Upper <- preds_fric$Upper-100

# predictions for FDis
preds_fdis <- Predict_effects(Newdata, model1b, TRUE)
preds_fdis$Median[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA
preds_fdis$Lower[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA
preds_fdis$Upper[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA

preds_fdis$Median <- preds_fdis$Median-100
preds_fdis$Lower <- preds_fdis$Lower-100
preds_fdis$Upper <- preds_fdis$Upper-100

## plotting
Limits <- c("Primary vegetation",
            "Mature secondary vegetation",
            "Intermediate secondary vegetation",
            "Young secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")

Labels=c("PV", "MSV", "ISV", "YSV", "PF", "PA", "CR", "UR")

cols <- scales::show_col(viridis(option = "plasma", n=8))
cbPalette <- c("#000000", viridis(option = "plasma", n=8)[1:7])
show_col(cbPalette)

p1 <- ggplot(preds_fric, aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  GGPoptions +
  facet_wrap(~Biome)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

p2 <- ggplot(preds_fdis[preds_fdis$Biome=="Temperate",], aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  scale_colour_manual(values=cbPalette) +
  GGPoptions + ggtitle("(b) FDis (Temperate; similar effects for tropical)")+
  guides(colour=FALSE, shape=FALSE)+
  ylab("Relative effect (%)")
# facet_wrap(~Biome, labeller = labeller(Biome =  c("Temperate" = "Temperate; similar effects for tropical"))) +
# theme( strip.text.x = element_text(size = 12, face = "bold"),
#        strip.text.y = element_text(size = 12, face = "bold"))

p2 <- ggarrange(p2, ggplot() + theme_void(), common.legend = TRUE)

## plotting figure 2
fig_validation_verts1 <- ggarrange(p1, p2, common.legend=TRUE, nrow=2)

## saving
ggsave(fig_validation_verts1, filename="../Revisions/PDF_figures/SI_Figure18.pdf",
       height=6, width=10)


#####################################################################################################
## validation 2: across vertebrates with complete trait data, with range size as an additional trait

rm(p1, p2, Gower_complete_verts, verts, Predicts_complete_verts, Complete_verts_comm, Complete_verts,
   preds_fdis, preds_fric)

Gower_complete_verts <- readRDS("../../Results/Gower_distances/Complete_vertebrates_with_rangesize.rds")
verts <- intersect(labels(Gower_complete_verts), Predicts$Best_guess_binomial)
Predicts_complete_verts <- subset(Predicts, Best_guess_binomial %in% verts)

# get community matrices
Complete_verts_comm <-  Generate_communities(Predicts_complete_verts, FALSE)

# compute FD
Complete_verts <- To_run(List_of_communities = na.omit.list(Complete_verts_comm),
                         Gowerdist = Gower_complete_verts,
                         Abundance = FALSE,
                         SimTRUE = FALSE,
                         nsim=NULL,
                         sim_pool = NULL,
                         global_pool_character = NULL,
                         Std.FRic = TRUE,
                         Predicts_site_info=Predicts_site_info)

Complete_verts <- Extract(Complete_verts, Predicts_site_info)
write.csv(Complete_verts, "../../Results/dbFD_indices/Complete_verts_validation_Range_size.csv", row.names = FALSE)

## fitting models 1a and 1b

Complete_verts <- Clean(Complete_verts)
Complete_verts$FRic[Complete_verts$SR==1] <- NA
Complete_verts$FDis[Complete_verts$SR==1] <- NA
# check levels
Complete_verts$Predominant_land_use
Complete_verts$Use_intensity
# add tropical versus temperate divide
Temperate <- Reg$SS[Reg$Biome=="Temperate"] %>%  unique()
Tropical <- Reg$SS[Reg$Biome=="Tropical"] %>%  unique()
Complete_verts$Biome <- NA
Complete_verts$Biome[Complete_verts$SS %in% Temperate] <- "Temperate"
Complete_verts$Biome[Complete_verts$SS %in% Tropical] <- "Tropical"

# tranform values
Complete_verts$asin_sqrt_FRic <- asin(sqrt(Complete_verts$FRic))
Complete_verts$asin_sqrt_FDis <- asin(sqrt(Complete_verts$FDis))

model1a <- lme4::lmer(asin_sqrt_FRic~Predominant_land_use+Use_intensity+Biome+
                        Predominant_land_use:Use_intensity+
                        Predominant_land_use:Biome+
                        (1|SS/SSB), data=Complete_verts)

model1b <- lme4::lmer(asin_sqrt_FDis~Predominant_land_use+Use_intensity+Biome+
                        Predominant_land_use:Use_intensity+
                        (1|SS/SSB), data=Complete_verts)

## predictions

Newdata <- expand.grid(Predominant_land_use=levels(Complete_verts$Predominant_land_use),
                       Use_intensity=levels(Complete_verts$Use_intensity),
                       Biome=c("Temperate", "Tropical"),
                       asin_sqrt_FDis=0,
                       asin_sqrt_FRic=0)

levels(Newdata$Predominant_land_use)

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
      y[1:8] <- y[1:8]/y[1]*100
      y[9:16] <- y[9:16]/y[9]*100
      y[17:24] <- y[17:24]/y[17]*100
      y[25:32] <- y[25:32]/y[25]*100
      y[33:40] <- y[33:40]/y[33]*100
      y[41:48] <- y[41:48]/y[41]*100
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
preds_fric <- Predict_effects(Newdata, model1a, TRUE)
preds_fric$Median[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA
preds_fric$Lower[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA
preds_fric$Upper[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA

preds_fric$Median <- preds_fric$Median-100
preds_fric$Lower <- preds_fric$Lower-100
preds_fric$Upper <- preds_fric$Upper-100

# predictions for FDis
preds_fdis <- Predict_effects(Newdata, model1b, TRUE)
preds_fdis$Median[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA
preds_fdis$Lower[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA
preds_fdis$Upper[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA

preds_fdis$Median <- preds_fdis$Median-100
preds_fdis$Lower <- preds_fdis$Lower-100
preds_fdis$Upper <- preds_fdis$Upper-100

## plotting
Limits <- c("Primary vegetation",
            "Mature secondary vegetation",
            "Intermediate secondary vegetation",
            "Young secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")

Labels=c("PV", "MSV", "ISV", "YSV", "PF", "PA", "CR", "UR")

cols <- scales::show_col(viridis(option = "plasma", n=8))
cbPalette <- c("#000000", viridis(option = "plasma", n=8)[1:7])
show_col(cbPalette)

p1 <- ggplot(preds_fric, aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  GGPoptions +
  facet_wrap(~Biome)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

p2 <- ggplot(preds_fdis[preds_fdis$Biome=="Temperate",], aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  scale_colour_manual(values=cbPalette) +
  GGPoptions + ggtitle("(b) FDis (Temperate; similar effects for tropical)")+
  guides(colour=FALSE, shape=FALSE)+
  ylab("Relative effect (%)")
# facet_wrap(~Biome, labeller = labeller(Biome =  c("Temperate" = "Temperate; similar effects for tropical"))) +
# theme( strip.text.x = element_text(size = 12, face = "bold"),
#        strip.text.y = element_text(size = 12, face = "bold"))

p2 <- ggarrange(p2, ggplot() + theme_void(), common.legend = TRUE)

## plotting figure 2
fig_validation_verts2 <- ggarrange(p1, p2, common.legend=TRUE, nrow=2)

## saving
ggsave(fig_validation_verts1, filename="../Revisions/PDF_figures/SI_Figure19.pdf",
       height=6, width=10)





#####################################################################################################
## validation 3: across each classes with complete trait data

rm(p1, p2, Gower_complete_verts, verts, Predicts_complete_verts, Complete_verts_comm, Complete_verts,
   preds_fdis, preds_fric)

Gower_complete_mammals <- readRDS("../../Results/Gower_distances/Complete_mammals.rds")
Gower_complete_birds <- readRDS("../../Results/Gower_distances/Complete_birds.rds")
Gower_complete_reptiles <- readRDS("../../Results/Gower_distances/Complete_reptiles.rds")
Gower_complete_herptiles <- readRDS("../../Results/Gower_distances/Complete_herptiles.rds")
Gower_complete_amphibians <- readRDS("../../Results/Gower_distances/Complete_amphibians.rds")

mammals <- intersect(labels(Gower_complete_mammals), Predicts$Best_guess_binomial)
Predicts_complete_mammals <- subset(Predicts, Best_guess_binomial %in% mammals)
birds <- intersect(labels(Gower_complete_birds), Predicts$Best_guess_binomial)
Predicts_complete_birds <- subset(Predicts, Best_guess_binomial %in% birds)
reptiles <- intersect(labels(Gower_complete_reptiles), Predicts$Best_guess_binomial)
Predicts_complete_reptiles <- subset(Predicts, Best_guess_binomial %in% reptiles)
herptiles <- intersect(labels(Gower_complete_herptiles), Predicts$Best_guess_binomial)
Predicts_complete_herptiles <- subset(Predicts, Best_guess_binomial %in% herptiles)
amphibians <- intersect(labels(Gower_complete_amphibians), Predicts$Best_guess_binomial)
Predicts_complete_amphibians <- subset(Predicts, Best_guess_binomial %in% amphibians)

# get community matrices
Complete_mammals_comm <-  Generate_communities(Predicts_complete_mammals, FALSE)
Complete_birds_comm <-  Generate_communities(Predicts_complete_birds, FALSE)
Complete_reptiles_comm <-  Generate_communities(Predicts_complete_reptiles, FALSE)
Complete_herptiles_comm <-  Generate_communities(Predicts_complete_herptiles, FALSE)
Complete_amphibians_comm <-  Generate_communities(Predicts_complete_amphibians, FALSE)

# compute FD
Complete_mammals <- To_run(List_of_communities = na.omit.list(Complete_mammals_comm),
                         Gowerdist = Gower_complete_mammals,
                         Abundance = FALSE,
                         SimTRUE = FALSE,
                         nsim=NULL,
                         sim_pool = NULL,
                         global_pool_character = NULL,
                         Std.FRic = TRUE,
                         Predicts_site_info=Predicts_site_info)

Complete_birds <- To_run(List_of_communities = na.omit.list(Complete_birds_comm),
                           Gowerdist = Gower_complete_birds,
                           Abundance = FALSE,
                           SimTRUE = FALSE,
                           nsim=NULL,
                           sim_pool = NULL,
                           global_pool_character = NULL,
                           Std.FRic = TRUE,
                           Predicts_site_info=Predicts_site_info)

Complete_reptiles <- To_run(List_of_communities = na.omit.list(Complete_reptiles_comm),
                         Gowerdist = Gower_complete_reptiles,
                         Abundance = FALSE,
                         SimTRUE = FALSE,
                         nsim=NULL,
                         sim_pool = NULL,
                         global_pool_character = NULL,
                         Std.FRic = TRUE,
                         Predicts_site_info=Predicts_site_info)

Complete_herptiles <- To_run(List_of_communities = na.omit.list(Complete_herptiles_comm),
                            Gowerdist = Gower_complete_herptiles,
                            Abundance = FALSE,
                            SimTRUE = FALSE,
                            nsim=NULL,
                            sim_pool = NULL,
                            global_pool_character = NULL,
                            Std.FRic = TRUE,
                            Predicts_site_info=Predicts_site_info)

Complete_amphibians <- To_run(List_of_communities = na.omit.list(Complete_amphibians_comm),
                             Gowerdist = Gower_complete_amphibians,
                             Abundance = FALSE,
                             SimTRUE = FALSE,
                             nsim=NULL,
                             sim_pool = NULL,
                             global_pool_character = NULL,
                             Std.FRic = TRUE,
                             Predicts_site_info=Predicts_site_info)

Complete_mammals <- Extract(Complete_mammals, Predicts_site_info)
Complete_mammals$Class <- "Mammals"
Complete_birds <- Extract(Complete_birds, Predicts_site_info)
Complete_birds$Class <- "Birds"
Complete_reptiles <- Extract(Complete_reptiles, Predicts_site_info)
Complete_reptiles$Class <- "Reptiles"
Complete_herptiles <- Extract(Complete_herptiles, Predicts_site_info)
Complete_herptiles$Class <- "Herptiles"
Complete_amphibians <- Extract(Complete_amphibians, Predicts_site_info)
Complete_amphibians$Class <- "Amphibians"

Complete_classes <- rbind(Complete_mammals, Complete_birds, Complete_reptiles, Complete_herptiles, Complete_amphibians)
write.csv(Complete_classes, "../../Results/dbFD_indices/Complete_classes_validation.csv", row.names = FALSE)

Complete_classes <- read.csv("../../Results/dbFD_indices/Complete_classes_validation.csv")
unique(Complete_classes$Class)

Complete_classes$Class %>%  levels()

Complete_classes <- Complete_classes %>%
  filter(Class!="Amphibians")
Complete_classes$Class <- droplevels(Complete_classes$Class)

## fitting models 2a and ab

Complete_classes <- Clean(Complete_classes)
Complete_classes$FRic[Complete_classes$SR==1] <- NA
Complete_classes$FDis[Complete_classes$SR==1] <- NA
# check levels
Complete_classes$Predominant_land_use
Complete_classes$Use_intensity
# add tropical versus temperate divide
Temperate <- Reg$SS[Reg$Biome=="Temperate"] %>%  unique()
Tropical <- Reg$SS[Reg$Biome=="Tropical"] %>%  unique()
Complete_classes$Biome <- NA
Complete_classes$Biome[Complete_classes$SS %in% Temperate] <- "Temperate"
Complete_classes$Biome[Complete_classes$SS %in% Tropical] <- "Tropical"

# tranform values
Complete_classes$asin_sqrt_FRic <- asin(sqrt(Complete_classes$FRic))
Complete_classes$asin_sqrt_FDis <- asin(sqrt(Complete_classes$FDis))

# change land uses
Complete_classes$Predominant_land_use <- as.character(Complete_classes$Predominant_land_use)
Complete_classes$Predominant_land_use[Complete_classes$Predominant_land_use %in% c("Mature secondary vegetation",
                                                                     "Intermediate secondary vegetation",
                                                                     "Secondary vegetation (indeterminate age)",
                                                                     "Young secondary vegetation")]  <- "Secondary vegetation"
Complete_classes$Predominant_land_use[Complete_classes$Predominant_land_use %in% c("Cropland",
                                                                     "Pasture")]  <- "Agricultural"
Complete_classes$Predominant_land_use %>%  unique()
Complete_classes$Predominant_land_use <- factor(Complete_classes$Predominant_land_use, levels=c("Primary vegetation",
                                                                                 "Secondary vegetation",
                                                                                 "Plantation forest",
                                                                                 "Agricultural",
                                                                                 "Urban"))


## fitting models
model2a <- lme4::lmer(asin_sqrt_FRic~Predominant_land_use+Use_intensity+Class+Biome+
                        Predominant_land_use:Use_intensity+
                        Predominant_land_use:Class+
                        Use_intensity:Biome+
                        Biome:Class+
                        (1|SS/SSB), data=Complete_classes)

model2b <- lme4::lmer(asin_sqrt_FDis~Predominant_land_use+Use_intensity+Class+
                        Predominant_land_use:Use_intensity+
                        Predominant_land_use:Class+
                        Use_intensity:Class+
                        (1|SS/SSB), data=Complete_classes)

## predictions
Newdata <- expand.grid(Predominant_land_use=levels(Complete_classes$Predominant_land_use),
                       Use_intensity=levels(Complete_classes$Use_intensity),
                       Class=c("Birds", "Mammals", "Herptiles", "Reptiles"),
                       Biome=c("Temperate", "Tropical"),
                       asin_sqrt_FDis=0,
                       asin_sqrt_FRic=0)

levels(Newdata$Predominant_land_use)

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
preds_fric$Median[preds_fric$Predominant_land_use=="Urban" & preds_fric$Class=="Reptiles"] <- NA
preds_fric$Lower[preds_fric$Predominant_land_use=="Urban" & preds_fric$Class=="Reptiles"] <- NA
preds_fric$Upper[preds_fric$Predominant_land_use=="Urban" & preds_fric$Class=="Reptiles"] <- NA

preds_fric$Median <- preds_fric$Median-100
preds_fric$Lower <- preds_fric$Lower-100
preds_fric$Upper <- preds_fric$Upper-100

# predictions for FDis
preds_fdis <- Predict_effects(Newdata, model2b, TRUE)
preds_fdis$Median[preds_fdis$Predominant_land_use=="Urban" & preds_fdis$Class=="Reptiles"] <- NA
preds_fdis$Lower[preds_fdis$Predominant_land_use=="Urban" & preds_fdis$Class=="Reptiles"] <- NA
preds_fdis$Upper[preds_fdis$Predominant_land_use=="Urban" & preds_fdis$Class=="Reptiles"] <- NA

preds_fdis$Median <- preds_fdis$Median-100
preds_fdis$Lower <- preds_fdis$Lower-100
preds_fdis$Upper <- preds_fdis$Upper-100

## plotting
Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Agricultural",
            "Urban")

Labels=c("PV", "SV", "PF", "AGR", "UR")

cols <- scales::show_col(viridis(option = "plasma", n=8))
cbPalette <- c("#000000", viridis(option = "plasma", n=8)[c(2,4,5,6)])
show_col(cbPalette)

preds_fric$Class <- factor(preds_fric$Class, levels = c("Amphibians","Birds", "Mammals", "Reptiles", "Herptiles"))
preds_fdis$Class <- factor(preds_fdis$Class, levels = c("Amphibians","Birds", "Mammals", "Reptiles", "Herptiles"))

p1 <- ggplot(preds_fric, aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  GGPoptions +
  facet_grid(Biome~Class)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

p2 <- ggplot(preds_fdis[preds_fdis$Biome=="Temperate",], aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  scale_colour_manual(values=cbPalette) +
  GGPoptions + ggtitle("(b) FDis (Temperate; similar effects for tropical)")+
  guides(colour=FALSE, shape=FALSE)+
  ylab("Relative effect (%)") +
  facet_grid(~Class)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold"))
# facet_wrap(~Biome, labeller = labeller(Biome =  c("Temperate" = "Temperate; similar effects for tropical"))) +
# theme( strip.text.x = element_text(size = 12, face = "bold"),
#        strip.text.y = element_text(size = 12, face = "bold"))

## plotting figure 2
fig_validation_verts_classes <- ggarrange(p1, p2, common.legend=TRUE, nrow=2, heights = c(0.6, 0.4))

## saving
ggsave(fig_validation_verts_classes,
       filename="../Revisions/PDF_figures/SI_Figure22.pdf",
       height=7, width=10)

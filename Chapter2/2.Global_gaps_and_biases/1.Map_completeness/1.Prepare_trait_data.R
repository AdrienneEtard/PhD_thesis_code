## Prepares trait data for mapping
# calculates completeness

library(dplyr)
library(stringr)
library(ggplot2)
#library(ggthemes)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
`%nin%` <- Negate(`%in%`)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13)) 

# function to id species with index
Link_index <- function(traits, indices) {
  x <- match(traits$Best_guess_binomial, indices$Binomial)
  print(paste(length(x[!is.na(x)]), "species matches:", round(length(x[!is.na(x)])/nrow(traits), digits=3)*100, "% species"))
  traits$Index <- indices$Index[x]
  return(traits)
}

# function to compute trait completeness
Completeness <- function(traits) {
  
  if(any(grepl("Family", colnames(traits)))) {
      Col <- which(colnames(traits)=="Family")
      Col2 <- which(colnames(traits)=="Best_guess_binomial")
      Col3 <- which(colnames(traits)=="Index")
      Col4 <- which(colnames(traits)=="Genus")
      Col5  <- which(colnames(traits)=="Order")
      traits$completeness <- apply(traits[, -c(Col, Col2, Col3, Col4, Col5)], 1, function(x)
    {return (length(x[!is.na(x)])/length(x)*100)})
  }
  
  else {
    traits$completeness <- apply(traits, 1, function(x)
    {return (length(x[!is.na(x)])/length(x)*100)})
  }
  
  return(traits)
  }


# # # # # # # # #

## load trait data
Amphibians <- read.csv("../../Data/Trait_data/Amphibians.csv")
Birds <- read.csv("../../Data/Trait_data/Birds.csv")
Mammals <- read.csv("../../Data/Trait_data/Mammals.csv")
Reptiles <- read.csv("../../Data/Trait_data/Reptiles.csv")

is.unsorted(Amphibians$Best_guess_binomial)
is.unsorted(Birds$Best_guess_binomial)
is.unsorted(Mammals$Best_guess_binomial)
is.unsorted(Reptiles$Best_guess_binomial)


## load species list for which spatial analyses can be run, based on ranges
ListBeforeCuts <- read.csv("../../Data/Distribution_maps/SpeciesBeforeCuts_List.csv")
ListAfterCuts <- read.csv("../../Data/Distribution_maps/SpeciesAfterCuts_List.csv")


##  species indices
Ind <- read.csv("../../Data/Distribution_maps/SpeciesIndices.csv")

## Correlations among traits

# amphibians: body mass - body length
cor(log(Amphibians$Body_mass_g), log(Amphibians$Body_length_mm), use = "complete.obs", method="pearson")

p1 <-  ggplot(Amphibians, aes(log(Body_mass_g), log(Body_length_mm), col=Order)) +
  geom_point() + GGPoptions + ylab("log(Body length) (mm)") + xlab("log(Body mass) (g)") +
  scale_color_brewer(palette="Dark2") + scale_color_manual(limits=c("ANURA", "CAUDATA", "GYMNOPHIONA"), 
                                                           labels=c("Anura", "Caudata", "Gymnophiona"),
                                                           values=c("#1B9E77", "#D95F02", "#7570B3"))
  
# longevity measures
plot(log(Reptiles$Maturity_d), log(Reptiles$Max_longevity_d))
cor(log(Reptiles$Maturity_d), log(Reptiles$Max_longevity_d), use = "complete.obs", method="pearson")

p2 <- ggplot(Amphibians, aes(log(Maturity_d), log(Longevity_d), col=Order)) +
  geom_point() + GGPoptions + ylab("log(Age at sexual maturity) (d)") + xlab("log(Longevity) (d)")+
  scale_color_manual(limits=c("ANURA", "CAUDATA", "GYMNOPHIONA"), 
                       labels=c("Anura", "Caudata", "Gymnophiona"),
                       values=c("#1B9E77", "#D95F02", "#7570B3"))

p <- ggpubr::ggarrange(p1+ggtitle("(a)"), p2+ggtitle("(b)"), common.legend=TRUE)

cor(log(Amphibians$Maturity_d), log(Amphibians$Max_longevity_d), use = "complete.obs", method="pearson")
ggsave(p, filename="../../Results/Correlations_traits/MatLong_BMBL_amphibians.pdf", width = 8, height = 3.5)
rm(p)

plot(log(Mammals$Generation_length_d), log(Mammals$Longevity_d))
cor(log(Mammals$Generation_length_d), log(Mammals$Longevity_d), use = "complete.obs", method="pearson")
cor(log(Birds$Generation_length_d), log(Birds$Longevity_d), use = "complete.obs", method="pearson")

p1 <- ggplot(Mammals, aes(log(Generation_length_d), log(Longevity_d))) + geom_point() + GGPoptions + ylab("log(Generation length) (d)") + xlab("log(Longevity) (d)")
p2 <- ggplot(Birds, aes(log(Generation_length_d), log(Longevity_d))) + geom_point() + GGPoptions + ylab("log(Generation length) (d)") + xlab("log(Longevity) (d)")
p <- ggpubr::ggarrange(p1+ggtitle("(a) Mammals"), p2+ggtitle("(b) Birds"), common.legend=TRUE)

ggsave(p, filename="../../Results/Correlations_traits/GL_Longe_birds_mammals.pdf", width = 8, height = 3.5)
rm(p)



## Identify species with index

## species match with Index

Amphibians <- Link_index(Amphibians, Ind)
Birds <- Link_index(Birds, Ind)
Mammals <- Link_index(Mammals, Ind)
Reptiles <- Link_index(Reptiles, Ind)

## Select traits and put trait data together in a list
# mammals and birds: BM, GL; TL; DA; HB; Spe; LCS; (RS)
# reptiles: BM, longevity; TL; DA; HB; Spe; LCS; (RS)
# amphibians: BL; maturity;...

## Amphibians

Amphibians.tomap <- Amphibians %>%
  dplyr::select("Body_length_mm", 
         "Maturity_d",
         "Litter_size",
         "Trophic_level",
         "Diel_activity",
         "Habitat_breadth_IUCN",
         "Specialisation", 
         "Index",
         "Family",
         "Best_guess_binomial") %>%
  dplyr::filter(!is.na(Index)) %>%
  mutate(Index=paste("sp", Index, sep=""))
rownames(Amphibians.tomap) <- Amphibians.tomap$Index
colnames(Amphibians.tomap)[c(1,2)] <- c("Body_size", "Life_span_proxy")

Amphibians.completeness <- Amphibians %>%
  dplyr::select("Body_length_mm", 
                "Maturity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation", 
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Amphibians.completeness)[c(1,2)] <- c("Body_size", "Life_span_proxy")

## Birds

Birds.tomap <- Birds %>%
  dplyr::select("Body_mass_g", 
         "Generation_length_d",
         "Litter_size",
         "Trophic_level",
         "Diel_activity",
         "Habitat_breadth_IUCN",
         "Specialisation", 
         "Index",
         "Family",
         "Best_guess_binomial") %>%
  dplyr::filter(!is.na(Index)) %>%
  mutate(Index=paste("sp", Index, sep=""))
rownames(Birds.tomap) <- Birds.tomap$Index
colnames(Birds.tomap)[c(1,2)] <- c("Body_size", "Life_span_proxy")

Birds.completeness <- Birds %>%
  dplyr::select("Body_mass_g", 
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation", 
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Birds.completeness)[c(1,2)] <- c("Body_size", "Life_span_proxy")


## Mammals

Mammals.tomap <- Mammals %>%
  dplyr::select("Body_mass_g", 
         "Generation_length_d",
         "Litter_size",
         "Trophic_level",
         "Diel_activity",
         "Habitat_breadth_IUCN",
         "Specialisation",
         "Index",
         "Family",
         "Best_guess_binomial") %>%
  dplyr::filter(!is.na(Index)) %>%
  mutate(Index=paste("sp", Index, sep=""))
rownames(Mammals.tomap) <- Mammals.tomap$Index
colnames(Mammals.tomap)[c(1,2)] <- c("Body_size", "Life_span_proxy")

Mammals.completeness <- Mammals %>%
  dplyr::select("Body_mass_g", 
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Mammals.completeness)[c(1,2)] <- c("Body_size", "Life_span_proxy")


## Reptiles

Reptiles.tomap <- Reptiles %>%
  dplyr::select("Body_mass_g", 
         "Max_longevity_d",
         "Litter_size",
         "Trophic_level",
         "Diel_activity",
         "Habitat_breadth_IUCN",
         "Specialisation",
         "Index",
         "Family",
         "Best_guess_binomial") %>%
  dplyr::filter(!is.na(Index)) %>%
  mutate(Index=paste("sp", Index, sep=""))
rownames(Reptiles.tomap) <- Reptiles.tomap$Index
colnames(Reptiles.tomap)[c(1,2)] <- c("Body_size", "Life_span_proxy")


Reptiles.completeness <- Reptiles %>%
  dplyr::select("Body_mass_g", 
                "Max_longevity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Reptiles.completeness)[c(1,2)] <- c("Body_size", "Life_span_proxy")


## Put files for plotting taxonomic & phylogenetic biases together and calculate completeness for these.
traits.completeness <- list(Amphibians.completeness, Birds.completeness, Mammals.completeness, Reptiles.completeness)
names(traits.completeness) <- c("Amphibians", "Birds", "Mammals", "Reptiles")
traits.completeness <- lapply(traits.completeness, Completeness)
saveRDS(traits.completeness, "../../Results/Traits_to_map/traits_completeness_V2.rds")
# NB: traits.completeness: these files include traits + family + order + genus + best guess binomial for all species 


## Put files together for mapping spatial biases and also calculate completeness for these

## For working with distribution maps before cutting by altitudinal limits
Amphibians.tomap_B <- Amphibians.tomap %>% 
  filter(Index %in% ListBeforeCuts$IndexSp)

Amphibians.tomap_A <- Amphibians.tomap %>% 
  filter(Index %in% ListAfterCuts$IndexSp)

Birds.tomap_B <- Birds.tomap %>% 
  filter(Index %in% ListBeforeCuts$IndexSp) 

Birds.tomap_A <- Birds.tomap %>% 
  filter(Index %in% ListAfterCuts$IndexSp)

Mammals.tomap_B <- Mammals.tomap %>% 
  filter(Index %in% ListBeforeCuts$IndexSp)

Mammals.tomap_A <- Mammals.tomap %>% 
  filter(Index %in% ListAfterCuts$IndexSp)

Reptiles.tomap_B <- Reptiles.tomap %>% 
  filter(Index %in% ListBeforeCuts$IndexSp)

Reptiles.tomap_A <- Reptiles.tomap %>% 
  filter(Index %in% ListAfterCuts$IndexSp) 


traits.tomap_B <- list(Amphibians.tomap_B, Birds.tomap_B, Mammals.tomap_B, Reptiles.tomap_B)
names(traits.tomap_B) <- c("Amphibians", "Birds", "Mammals", "Reptiles")
traits.tomap_B <- lapply(traits.tomap_B, Completeness)
saveRDS(traits.tomap_B, "../../Results/Traits_to_map/traits_tomap_SpeciesBeforeCuts.rds")


traits.tomap_A <- list(Amphibians.tomap_A, 
                       Birds.tomap_A,
                       Mammals.tomap_A,
                       Reptiles.tomap_A)
names(traits.tomap_A) <- c("Amphibians", "Birds", "Mammals", "Reptiles")
traits.tomap_A <- lapply(traits.tomap_A, Completeness)
saveRDS(traits.tomap_A, "../../Results/Traits_to_map/traits_tomap_SpeciesAfterCuts.rds")

lapply(traits.tomap_A, nrow)


# NB: traits.tomap: these files include traits + family  + best guess binomial for species which have a range map file that I can work with (species in ListAfter and ListBefore).


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


## using the trait datasets uncorrected for taxonomy for taxonomic biases.

Amphibians_unc <- read.csv("../../Data/Trait_data_UNCORRECTED_TAXONOMY/Amphibians.csv")
Birds_unc <- read.csv("../../Data/Trait_data_UNCORRECTED_TAXONOMY/Birds.csv")
Mammals_unc <- read.csv("../../Data/Trait_data_UNCORRECTED_TAXONOMY/Mammals.csv")
Reptiles_unc <- read.csv("../../Data/Trait_data_UNCORRECTED_TAXONOMY/Reptiles.csv")


Amphibians.completeness_unc <- Amphibians_unc %>%
  dplyr::select("Body_length_mm", 
                "Maturity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation", 
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Amphibians.completeness_unc)[c(1,2)] <- c("Body_size", "Life_span_proxy")

Birds.completeness_unc <- Birds_unc %>%
  dplyr::select("Body_mass_g", 
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation", 
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Birds.completeness_unc)[c(1,2)] <- c("Body_size", "Life_span_proxy")

Mammals.completeness_unc <- Mammals_unc %>%
  dplyr::select("Body_mass_g", 
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Mammals.completeness_unc)[c(1,2)] <- c("Body_size", "Life_span_proxy")


Reptiles.completeness_unc <- Reptiles_unc %>%
  dplyr::select("Body_mass_g", 
                "Max_longevity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                "Family",
                "Order",
                "Genus",
                "Best_guess_binomial")
colnames(Reptiles.completeness_unc)[c(1,2)] <- c("Body_size", "Life_span_proxy")


## Put files for plotting taxonomic & phylogenetic biases together and calculate completeness for these.
traits.completeness_unc <- list(Amphibians.completeness_unc, Birds.completeness_unc, Mammals.completeness_unc, Reptiles.completeness_unc)
names(traits.completeness_unc) <- c("Amphibians", "Birds", "Mammals", "Reptiles")
traits.completeness_unc <- lapply(traits.completeness_unc, Completeness)
saveRDS(traits.completeness_unc, "../../Results/Traits_to_map/traits_completeness_V2_UNCORRECTED_TAXONOMY.rds")
# NB: traits.completeness: these files include traits + family + order + genus + best guess binomial for all species 



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Tim asked for a map with % of habitat specialist species
# ## For habitat specialisation - uncorrected for taxonomy
# 
# ## Load trait data
# Amphibians_un <- read.csv("../../Data/Trait_data/Uncorrected_taxonomy/Amphibians.csv")
# Birds_un <- read.csv("../../Data/Trait_data/Uncorrected_taxonomy/Birds.csv")
# Mammals_un <- read.csv("../../Data/Trait_data/Uncorrected_taxonomy/Mammals.csv")
# Reptiles_un <- read.csv("../../Data/Trait_data/Uncorrected_taxonomy/Reptiles.csv")
# 
# is.unsorted(Amphibians_un$Best_guess_binomial)
# is.unsorted(Birds_un$Best_guess_binomial)
# is.unsorted(Mammals_un$Best_guess_binomial)
# is.unsorted(Reptiles_un$Best_guess_binomial)
# 
# Amphibians_un <- Amphibians_un[order(Amphibians_un$Best_guess_binomial),]
# Birds_un <- Birds_un[order(Birds_un$Best_guess_binomial),]
# Mammals_un <- Mammals_un[order(Mammals_un$Best_guess_binomial),]
# Reptiles_un <- Reptiles_un[order(Reptiles_un$Best_guess_binomial),]
# 
# ## Identify species with index
# Amphibians_un <- Link_index(Amphibians_un, Ind)
# Birds_un <- Link_index(Birds_un, Ind)
# Mammals_un <- Link_index(Mammals_un, Ind)
# Reptiles_un <- Link_index(Reptiles_un, Ind)
# 
# ## select for mapping proportion of habitat specialist (without taxonomic corrections)
# Amphibians_un.tomap <- Amphibians_un %>%
#   dplyr::select("Specialisation", 
#                 "Index",
#                 "Best_guess_binomial") %>%
#   dplyr::filter(!is.na(Index)) %>%
#   mutate(Index=paste("sp", Index, sep=""))
# rownames(Amphibians_un.tomap) <- Amphibians_un.tomap$Index
# Amphibians_un.tomap <- Amphibians_un.tomap %>% dplyr::select(-Index)
# 
# Mammals_un.tomap <- Mammals_un %>%
#   dplyr::select("Specialisation", 
#                 "Index",
#                 "Best_guess_binomial") %>%
#   dplyr::filter(!is.na(Index)) %>%
#   mutate(Index=paste("sp", Index, sep=""))
# rownames(Mammals_un.tomap) <- Mammals_un.tomap$Index
# Mammals_un.tomap <- Mammals_un.tomap %>% dplyr::select(-Index)
# 
# Birds_un.tomap <- Birds_un %>%
#   dplyr::select("Specialisation", 
#                 "Index",
#                 "Best_guess_binomial") %>%
#   dplyr::filter(!is.na(Index)) %>%
#   mutate(Index=paste("sp", Index, sep=""))
# rownames(Birds_un.tomap) <- Birds_un.tomap$Index
# Birds_un.tomap <- Birds_un.tomap %>% dplyr::select(-Index)
# 
# Reptiles_un.tomap <- Reptiles_un %>%
#   dplyr::select("Specialisation", 
#                 "Index",
#                 "Best_guess_binomial") %>%
#   dplyr::filter(!is.na(Index)) %>%
#   mutate(Index=paste("sp", Index, sep=""))
# rownames(Reptiles_un.tomap) <- Reptiles_un.tomap$Index
# Reptiles_un.tomap <- Reptiles_un.tomap %>% dplyr::select(-Index)
# 
# Un_hab_spe <- list(Amphibians_un.tomap, Birds_un.tomap, Mammals_un.tomap, Reptiles_un.tomap)
# names(Un_hab_spe) <- c("Amphibians", "Birds", "Mammals", "Reptiles")
# saveRDS(Un_hab_spe, "../../Results/Traits_to_map/Habitat_uncorrected.rds")



## Trait data formatting for release

library(dplyr)

Amphibians <- read.csv("../Data/Trait_data/Amphibians.csv")
Habitat <- colnames(Amphibians)[22:34]
Uncor.Amphibians <- read.csv("../Data/Trait_data_UNCORRECTED_TAXONOMY/Amphibians.csv")

Amphibians <- Amphibians %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Body_length_mm", 
                "Body_mass_g",
                "Maturity_d",
                "Max_longevity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note") %>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))

Uncor.Amphibians <- Uncor.Amphibians %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Body_length_mm", 
                "Body_mass_g",
                "Maturity_d",
                "Max_longevity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note") %>% 
  dplyr::filter(!Best_guess_binomial %in% c("Rana cafe", "Rana hoja", "Rana hojarasca", "Rana verde", "Arthroleptis sp")) %>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))

colnames(Amphibians)[9] <- "Litter_clutch_size"
colnames(Uncor.Amphibians)[9] <- "Litter_clutch_size"
colnames(Amphibians)[13] <- "Artificial_habitat_use"
colnames(Uncor.Amphibians)[13] <- "Artificial_habitat_use"

write.csv(Amphibians, "../Data/Trait_data/DataForSharing/Final/Amphibians.csv", row.names = FALSE)
write.csv(Uncor.Amphibians, "../Data/Trait_data/DataForSharing/Final/Taxonomy_uncorrected_Amphibians.csv", row.names = FALSE)


Birds <- read.csv("../Data/Trait_data/Birds.csv")
Uncor.Birds <- read.csv("../Data/Trait_data_UNCORRECTED_TAXONOMY/Birds.csv")

Birds <- Birds %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Body_mass_g",
                "Adult_svl_cm",
                "Maturity_d",
                "Max_longevity_d",
                "Longevity_d",
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note")%>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))

Uncor.Birds <- Uncor.Birds %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Body_mass_g",
                "Adult_svl_cm", 
                "Maturity_d",
                "Max_longevity_d",
                "Longevity_d",
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note") %>% 
  dplyr::filter(!is.na(Family))%>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))

colnames(Birds)[11] <- "Litter_clutch_size"
colnames(Uncor.Birds)[11] <- "Litter_clutch_size"
colnames(Birds)[15] <- "Artificial_habitat_use"
colnames(Uncor.Birds)[15] <- "Artificial_habitat_use"

write.csv(Birds, "../Data/Trait_data/DataForSharing/Final/Birds.csv", row.names = FALSE)
write.csv(Uncor.Birds, "../Data/Trait_data/DataForSharing/Final/Taxonomy_uncorrected_Birds.csv", row.names = FALSE)


Mammals <- read.csv("../Data/Trait_data/Mammals.csv")
Uncor.Mammals <- read.csv("../Data/Trait_data_UNCORRECTED_TAXONOMY/Mammals.csv")


Mammals <- Mammals %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Body_mass_g",
                "Adult_svl_cm",
                "Maturity_d",
                "Max_longevity_d",
                "Longevity_d",
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note")%>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))


Uncor.Mammals <- Uncor.Mammals %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Body_mass_g",
                "Adult_svl_cm", 
                "Maturity_d",
                "Max_longevity_d",
                "Longevity_d",
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note") %>% 
  dplyr::filter(!is.na(Family))%>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))


colnames(Mammals)[11] <- "Litter_clutch_size"
colnames(Uncor.Mammals)[11] <- "Litter_clutch_size"
colnames(Mammals)[15] <- "Artificial_habitat_use"
colnames(Uncor.Mammals)[15] <- "Artificial_habitat_use"

write.csv(Mammals, "../Data/Trait_data/DataForSharing/Final/Mammals.csv", row.names = FALSE)
write.csv(Uncor.Mammals, "../Data/Trait_data/DataForSharing/Final/Taxonomy_uncorrected_Mammals.csv", row.names = FALSE)

Reptiles <- read.csv("../Data/Trait_data/Reptiles.csv")
Uncor.Reptiles <- read.csv("../Data/Trait_data_UNCORRECTED_TAXONOMY/Reptiles.csv")

Reptiles <- Reptiles %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Adult_svl_cm", 
                "Body_mass_g",
                "Maturity_d",
                "Max_longevity_d",
                "Longevity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note") %>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))


# ## checks
# nrow(Reptiles[!is.na(Reptiles$Longevity_d),])
# nrow(Reptiles[!is.na(Reptiles$Max_longevity_d),])
# Reptiles %>%  filter(!is.na(Longevity_d)) %>%  filter(!is.na(Max_longevity_d)) %>%  nrow() +
# Reptiles %>%  filter(is.na(Longevity_d)) %>%  filter(!is.na(Max_longevity_d)) %>%  nrow() +
# Reptiles %>%  filter(is.na(Max_longevity_d)) %>%  filter(!is.na(Longevity_d)) %>%  nrow()


Uncor.Reptiles <- Uncor.Reptiles %>%
  dplyr::select("Order", 
                "Family", 
                "Genus", 
                "Best_guess_binomial",
                "Adult_svl_cm", 
                "Body_mass_g",
                "Maturity_d",
                "Max_longevity_d",
                "Longevity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                all_of(Habitat),
                "Note") %>% 
  dplyr::filter(!is.na(Family)) %>% 
  mutate(Specialisation=ifelse(Specialisation=="Natural habitat specialist", FALSE, TRUE))

colnames(Reptiles)[10] <- "Litter_clutch_size"
colnames(Uncor.Reptiles)[10] <- "Litter_clutch_size"
colnames(Reptiles)[14] <- "Artificial_habitat_use"
colnames(Uncor.Reptiles)[14] <- "Artificial_habitat_use"


write.csv(Reptiles, "../Data/Trait_data/DataForSharing/Final/Reptiles.csv", row.names = FALSE)
write.csv(Uncor.Reptiles, "../Data/Trait_data/DataForSharing/Final/Taxonomy_uncorrected_Reptiles.csv", row.names = FALSE)



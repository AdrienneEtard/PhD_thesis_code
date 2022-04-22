## prepare data for FigShare

setwd("../../Gaps_biases_trait_data/Results/Traits_to_map")


library(dplyr)
`%nin%` <- Negate(`%in%`)

x <- readRDS("traits_completeness.rds")

#######################################

## Amphibians
Amphibians <- x[["Amphibians"]] %>%
  mutate(Class="Amphibians") %>% 
  select(-completeness)

Amphibians <- Amphibians[, c(10,8,9, 1:7)]
colnames(Amphibians)[4] <- "Body_size_mm"
colnames(Amphibians)[5] <- "Age_sexual_maturity_d"
colnames(Amphibians)[6] <- "Litter_Clutch_size"
colnames(Amphibians)[9] <- "Habitat_breadth"

is.unsorted(Amphibians$Best_guess_binomial)

Amp <- read.csv("../../Data/Trait_data/Amphibians.csv")
Amp <- Amp[order(Amp$Best_guess_binomial), ]

Amphibians$Order <- Amp$Order
Amphibians <- Amphibians[, c(1,11, 2:10 )]

########################################

## Birds

Birds <- x[["Birds"]] %>%
  mutate(Class="Birds") %>% 
  select(-completeness)

Birds <- Birds[, c(10,8,9, 1:7)]
colnames(Birds)[4] <- "Body_mass_g"
colnames(Birds)[5] <- "Generation_length_d"
colnames(Birds)[6] <- "Litter_Clutch_size"
colnames(Birds)[9] <- "Habitat_breadth"

is.unsorted(Birds$Best_guess_binomial)

Bi <- read.csv("../../Data/Trait_data/Birds.csv")
Bi <- Bi[order(Bi$Best_guess_binomial), ]

Birds$Order <- Bi$Order
Birds <- Birds[, c(1,11, 2:10 )]

#########################################

## Mammals

Mammals <- x[["Mammals"]] %>%
  mutate(Class="Mammals") %>% 
  select(-completeness)

Mammals <- Mammals[, c(10,8,9, 1:7)]
colnames(Mammals)[4] <- "Body_mass_g"
colnames(Mammals)[5] <- "Generation_length_d"
colnames(Mammals)[6] <- "Litter_Clutch_size"
colnames(Mammals)[9] <- "Habitat_breadth"

is.unsorted(Mammals$Best_guess_binomial)

Mam <- read.csv("../../Data/Trait_data/Mammals.csv")
Mam <- Mam[order(Mam$Best_guess_binomial), ]

Marine_mammals1 <- Mam$Best_guess_binomial[Mam$Order=="SIRENIA"]
Marine_mammals2 <- Mam$Best_guess_binomial[Mam$Family %in% c("OTARIIDAE", "PHOCIDAE", "ODOBENIDAE", 
                                                                     "BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", 
                                                                     "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE",
                                                                     "ESCHRICHTIIDAE", "INIIDAE", "PHYSETERIDAE","LIPOTIDAE",
                                                                     "PHOCOENIDAE", "PLATANISTIDAE", "PONTOPORIIDAE")]
Marine_mammals <- c(as.character(Marine_mammals1), as.character(Marine_mammals2)) %>% unique()
rm(Marine_mammals1, Marine_mammals2)

Mam <- Mam %>%
  filter(Order %nin% c("SIRENIA")) %>%
  filter(Family %nin% c("OTARIIDAE", "PHOCIDAE", "ODOBENIDAE", 
                        "BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", 
                        "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE",
                        "ESCHRICHTIIDAE", "INIIDAE", "PHYSETERIDAE","LIPOTIDAE",
                        "PHOCOENIDAE", "PLATANISTIDAE", "PONTOPORIIDAE"))

Mammals$Order <- Mam$Order
Mammals <- Mammals[, c(1,11, 2:10 )]

#########################################

## Reptiles

Reptiles <- x[["Reptiles"]] %>%
  mutate(Class="Reptiles") %>% 
  select(-completeness)

Reptiles <- Reptiles[, c(10,8,9, 1:7)]
colnames(Reptiles)[4] <- "Body_mass_g"
colnames(Reptiles)[5] <- "Longevity_d"
colnames(Reptiles)[6] <- "Litter_Clutch_size"
colnames(Reptiles)[9] <- "Habitat_breadth"

is.unsorted(Reptiles$Best_guess_binomial)

Re <- read.csv("../../Data/Trait_data/Reptiles.csv")
Re <- Re[order(Re$Best_guess_binomial), ]

Reptiles$Order <- Re$Order
Reptiles <- Reptiles[, c(1,11, 2:10 )]

###########################################


Trait_data_to_share <- list(Amphibians, Mammals, Birds, Reptiles)
saveRDS(Trait_data_to_share, file = "../../Results/Trait_data_for_sharing/TerVerTraits.rds")

write.csv(Amphibians, "../../Results/Trait_data_for_sharing/Amphibians.csv", row.names = FALSE)
write.csv(Birds, "../../Results/Trait_data_for_sharing/Birds.csv", row.names = FALSE)
write.csv(Mammals, "../../Results/Trait_data_for_sharing/Mammals.csv", row.names = FALSE)
write.csv(Reptiles, "../../Results/Trait_data_for_sharing/Reptiles.csv", row.names = FALSE)





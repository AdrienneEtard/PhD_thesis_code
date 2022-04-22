# #  Processing synonym datasets  # # 


# #  P R E A M B L E  # #  

X <- c("dplyr", "taxize", "phytools", "stringr", "rredlist", "stringdist", "plyr", "pbapply", "GlobalOptions", "data.table", "ngram", "curl", "RCurl")
lapply(X, library, character.only=TRUE); rm(X)
`%nin%` = Negate(`%in%`)

opt <- options(iucn_redlist_key="ba30954f38deda075bd9b52495d65092ccf1b220b0c7c67a41465646e50ef72c")

source("Resolve_taxonomy_functions.R")

# Load files
Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Mammals.csv") %>% ToChar()
Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Reptiles.csv") %>% ToChar()
Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Amphibians.csv") %>% ToChar()
Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Birds.csv") %>% ToChar()


#  #  #  #  #


## Check if species names have 2 words (original, corrected and accepted)
Check_binomial(Mammals)
Check_binomial(Reptiles)
Check_binomial(Amphibians)
Check_binomial(Birds)

## For accepted names that have only one word, take corrected name instead; for names that have more than two, take first two words
Mammals$Accepted <- word(Mammals$Accepted, 1,2)
Reptiles$Accepted <- word(Reptiles$Accepted, 1,2)
Amphibians$Accepted <- word(Amphibians$Accepted, 1,2)
Birds$Accepted <- word(Birds$Accepted, 1,2)

One.M <- which(sapply(strsplit(Mammals$Accepted, " "), length)==1)
Mammals$Accepted[One.M] <- Mammals$CorrectedTypos[One.M]

One.R <- which(sapply(strsplit(Reptiles$Accepted, " "), length)==1)
Reptiles$Accepted[One.R] <- Reptiles$CorrectedTypos[One.R]

One.A <- which(sapply(strsplit(Amphibians$Accepted, " "), length)==1)
Amphibians$Accepted[One.A] <- Amphibians$CorrectedTypos[One.A]

One.B <- which(sapply(strsplit(Birds$Accepted, " "), length)==1)
Birds$Accepted[One.B] <- Birds$CorrectedTypos[One.B]
One.B <- which(sapply(strsplit(Birds$Accepted, " "), length)==1)
Birds$Accepted[One.B] <- Birds$Original[One.B]


## For species that are entered as Genus sp. or Genus cf.: genus level
Mammals <- ToGenus(Mammals)
Reptiles <- ToGenus(Reptiles)
Amphibians <- ToGenus(Amphibians)
Birds <- ToGenus(Birds)


## Delta species
length(unique(Reptiles$Original[Reptiles$Accepted!=""])) - length(unique(Reptiles$Accepted[Reptiles$Accepted!=""]))
length(unique(Mammals$Original[Mammals$Accepted!=""])) - length(unique(Mammals$Accepted[Mammals$Accepted!=""]))
length(unique(Amphibians$Original[Amphibians$Accepted!=""])) - length(unique(Amphibians$Accepted[Amphibians$Accepted!=""]))
length(unique(Birds$Original[Birds$Accepted!=""])) - length(unique(Birds$Accepted[Birds$Accepted!=""]))


## Genus column
Reptiles$Genus <- word(Reptiles$Accepted, 1,1)
Reptiles$Genus[Reptiles$Genus==""] <- word(Reptiles$CorrectedTypos[Reptiles$Genus==""],1,1)
Mammals$Genus <- word(Mammals$Accepted, 1, 1)
Mammals$Genus[Mammals$Genus==""] <- word(Mammals$CorrectedTypos[Mammals$Genus==""],1,1)
Amphibians$Genus <- word(Amphibians$Accepted, 1, 1)
Amphibians$Genus[Amphibians$Genus==""] <- word(Amphibians$CorrectedTypos[Amphibians$Genus==""],1,1)
Birds$Genus <- word(Birds$Accepted, 1, 1)
Birds$Genus[Birds$Genus==""] <- word(Birds$CorrectedTypos[Birds$Genus==""],1,1)


## Family and Order information that is missing: fill with informaiton extracted from GBIF

## Reptiles
# ForGBIF_reptiles <- Order_family_info(Reptiles)$ForGBIF
ToComplement_Reptiles <- Order_family_info(Reptiles)$TaxInfo %>% filter(!is.na(Order))

## Mammals
# ForGBIF_mammals <- Order_family_info(Mammals)$ForGBIF
ToComplement_Mammals <- Order_family_info(Mammals)$TaxInfo %>% filter(!is.na(Order))

## Birds
# ForGBIF_birds <- Order_family_info(Birds)$ForGBIF
ToComplement_Birds <- Order_family_info(Birds)$TaxInfo %>% filter(!is.na(Order))

## Amphibians
# ForGBIF_amphibians <- Order_family_info(Amphibians)$ForGBIF
ToComplement_Amphibians <- Order_family_info(Amphibians)$TaxInfo %>% filter(!is.na(Order))


# ## Write GBIF info to extract
# write.csv(ForGBIF_mammals, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Mammals.csv", row.names = FALSE)
# write.csv(ForGBIF_reptiles, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Reptiles.csv", row.names = FALSE)
# write.csv(ForGBIF_amphibians, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Amphibians.csv", row.names = FALSE)
# write.csv(ForGBIF_birds, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Birds.csv", row.names = FALSE)


## Open GBIF outputs, and filter (those are complements)

Complements_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Mammals.csv") %>%
  select(c("genus", "family", "order", "verbatimScientificName")) %>% 
  mutate(family=toupper(family), order=toupper(order))
# Manual input: one for mammal
Complements_Mammals$genus <- as.character(Complements_Mammals$genus)
Complements_Mammals$family[Complements_Mammals$verbatimScientificName=="Galerella nigrata"] <- "HERPESTIDAE"
Complements_Mammals$genus[Complements_Mammals$verbatimScientificName=="Galerella nigrata"] <- as.character("Galerella")
Complements_Mammals$order[Complements_Mammals$verbatimScientificName=="Galerella nigrata"] <- "CARNIVORA"

Complements_Mammals <- Complements_Mammals %>%
  select(-verbatimScientificName) %>%
  setNames(., c("Genus", "Family", "Order"))

Complements_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Reptiles.csv") %>%
  select(c("genus", "family", "order"))%>% 
  mutate(family=toupper(family), order=toupper(order)) %>%
  setNames(., c("Genus", "Family", "Order"))

Complements_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Amphibians.csv") %>%
  select(c("genus", "family", "order"))%>% 
  mutate(family=toupper(family), order=toupper(order)) %>%
  setNames(., c("Genus", "Family", "Order"))

Complements_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Birds.csv") %>%
  select(c("genus", "family", "order"))%>%
  mutate(family=toupper(family), order=toupper(order))%>%
  setNames(., c("Genus", "Family", "Order")) %>%
  filter(Order!="")


## Merge
Mammals_additional_taxInfo <- rbind(ToComplement_Mammals, Complements_Mammals)
Mammals_additional_taxInfo <- Mammals_additional_taxInfo[order(Mammals_additional_taxInfo$Genus),]

Reptiles_additional_taxInfo <- rbind(ToComplement_Reptiles, Complements_Reptiles)
Reptiles_additional_taxInfo <- Reptiles_additional_taxInfo[order(Reptiles_additional_taxInfo$Genus),]

Amphibians_additional_taxInfo <- rbind(ToComplement_Amphibians, Complements_Amphibians)
Amphibians_additional_taxInfo <- Amphibians_additional_taxInfo[order(Amphibians_additional_taxInfo$Genus),]

Birds_additional_taxInfo <- rbind(ToComplement_Birds, Complements_Birds)
Birds_additional_taxInfo <- Birds_additional_taxInfo[order(Birds_additional_taxInfo$Genus),]


## Add additional taxonomic information to main synonym datasets
Mammals <- AddFamilyOrder(Mammals, Mammals_additional_taxInfo)$Syn
# GBIF_2_Mammals <- AddFamilyOrder(Mammals, Mammals_additional_taxInfo)$Missing %>% as.data.frame()%>% setNames(., "scientificName")

Reptiles <- AddFamilyOrder(Reptiles, Reptiles_additional_taxInfo)$Syn
# GBIF_2_Reptiles <- AddFamilyOrder(Reptiles, Reptiles_additional_taxInfo)$Missing %>% as.data.frame()%>% setNames(., "scientificName")

Amphibians <- AddFamilyOrder(Amphibians, Amphibians_additional_taxInfo)$Syn
# GBIF_2_Amphibians <- AddFamilyOrder(Amphibians, Amphibians_additional_taxInfo)$Missing %>% as.data.frame()%>% setNames(., "scientificName")

# Correct one typo 
Birds$CorrectedTypos[Birds$Original=="Brachypteryx cruralis"] <- "Brachypteryx cruralis"
Birds$Accepted[Birds$Original=="Brachypteryx cruralis"] <- "Brachypteryx cruralis"
Birds <- AddFamilyOrder(Birds, Birds_additional_taxInfo)$Syn
GBIF_2_Birds <- AddFamilyOrder(Birds, Birds_additional_taxInfo)$Missing %>% as.data.frame()%>% setNames(., "scientificName")

# write.csv(GBIF_2_Mammals, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Mammals2.csv", row.names = FALSE)
# write.csv(GBIF_2_Reptiles, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Reptiles2.csv", row.names = FALSE)
# write.csv(GBIF_2_Amphibians, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Amphibians2.csv", row.names = FALSE)
# write.csv(GBIF_2_Birds, "../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/ToGet/Birds2.csv", row.names = FALSE)


## Second round of GBIF addition: load files
MammalsGBIF2 <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Mammals2.csv")
ReptilesGBIF2 <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Reptiles2.csv")
AmphibiansGBIF2 <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Amphibians2.csv")
BirdsGBIF2 <- read.csv("../../Results/0.Data_resolved_taxonomy/GBIF_Extra_family_order_information/GBIFOutputs/Mammals2.csv")

Mammals_final <- GBIF2(Mammals, MammalsGBIF2)
Reptiles_final <- GBIF2(Reptiles, ReptilesGBIF2)
Amphibians_final <- GBIF2(Amphibians, AmphibiansGBIF2)
Birds_final <- GBIF2(Birds, BirdsGBIF2)


## Some species still without order or family information for reptiles: MANUAL ENTRIES - CHECKS AT GENUS LEVEL
Reptiles_left_to_check <- subset(Reptiles_final, is.na(Order))
Manual_reptile_additions <- unique(Reptiles_left_to_check$CorrectedTypos) %>%
  as.data.frame() %>%
  setNames(., "Corrected_name")
Manual_reptile_additions$Genus <- word(Manual_reptile_additions$Corrected_name,1,1)  #38 genuses to check manually
Manual_reptile_additions$Family <- toupper(c("Amphisbaenidae", rep("Typhlopidae",2), "Amphisbaenidae", "Typhlopidae", "Scincidae", "Viperidae", "Gekkonidae",
                                     rep("Amphisbaenidae",3), "Scincidae", "Colubridae", rep("Amphisbaenidae",2), "Scincidae", rep("Dactyloidae", 4),
                                     rep("Colubridae",2), "Kinosternidae", rep("Dactyloidae", 4), rep("Colubridae",2), "Gekkonidae", "Scincidae",
                                     "Elapidae", rep("Amphisbaenidae",2), "Elapidae", rep("Typhlopidae",2), rep("Dactyloidae", 22), "Colubridae",
                                     "Gymnophthalmidae", rep("Colubridae",2), rep("Gymnophthalmidae",2), "Colubridae", rep("Scincidae",2),
                                     rep("Colubridae",2), "Xenodermidae", "Typhlopidae", "Elapidae", rep("Teiidae",5), "Colubridae"))
Manual_reptile_additions$Order <-toupper(c(rep("Squamata",22), "Testudines", rep("Squamata",56)))

for(i in 1:nrow(Manual_reptile_additions)) {
  Reptiles_final$Order[Reptiles_final$CorrectedTypos==Manual_reptile_additions$Corrected_name[i]] <- Manual_reptile_additions$Order[i]
  Reptiles_final$Family[Reptiles_final$CorrectedTypos==Manual_reptile_additions$Corrected_name[i]] <- Manual_reptile_additions$Family[i]
}


## Some species still without order or family information for birds: MANUAL ENTRIES - CHECKS AT GENUS LEVEL
Birds_left_to_check <- subset(Birds_final, is.na(Order))
Manual_bird_additions <- unique(Birds_left_to_check$CorrectedTypos) %>%
  as.data.frame() %>%
  setNames(., "Corrected_name")
Manual_bird_additions$Genus <- word(Manual_bird_additions$Corrected_name,1,1)  #38 genuses to check manually
Manual_bird_additions$Family <- toupper(c("Accipitridae", rep("Muscicapidae",2), rep("Thraupidae",2), "Megalaimidae", "Acanthizidae", rep("Muscicapidae",5), "Dicaeidae",
                                        "Platysteiridae", rep("Campephagidae", 4), "Fringillidae", rep("Phasianidae",2), "Muscicapidae", rep("Psittacidae",3), "Apodidae",
                                        "Muscicapidae", rep("Laridae",4), rep("Megalapterygidae",2), "Tityridae", "Timaliidae", "Apodidae", "Columbidae", "Thraupidae", "Psittacidae"))

Manual_bird_additions$Order <-toupper(c("Accipitriformes", rep("Passeriformes",4), "Piciformes", rep("Passeriformes",13), rep("Galliformes",2), "Passeriformes", rep("Psittaciformes",3),
                                      "Caprimulgiformes", "Passeriformes", rep("Charadriiformes",4), rep("Dinornithiformes",2), rep("Passeriformes",2), "Apodiformes", "Columbiformes",
                                      "Passeriformes","Psittaciformes"))

for(i in 1:nrow(Manual_bird_additions)) {
  Birds_final$Order[Birds_final$CorrectedTypos==Manual_bird_additions$Corrected_name[i]] <- Manual_bird_additions$Order[i]
  Birds_final$Family[Birds_final$CorrectedTypos==Manual_bird_additions$Corrected_name[i]] <- Manual_bird_additions$Family[i]
}


## Reorganise columns
Mammals_final <- Mammals_final[, c("Original", "CorrectedTypos", "IsCorrected","IsAccepted", 
                                   "IsSynonym", "Accepted", "Synonyms", "InRedList", "InITIS", "FuzzyMatch",
                                   "Genus_level", "Genus", "Family", "Order")]

Amphibians_final <- Amphibians_final[, c("Original", "CorrectedTypos", "IsCorrected","IsAccepted", 
                                         "IsSynonym", "Accepted", "Synonyms", "InRedList", "InITIS", "FuzzyMatch",
                                         "Genus_level", "Genus", "Family", "Order")]


Reptiles_final <- Reptiles_final[, c("Original", "CorrectedTypos", "IsCorrected","IsAccepted", 
                                         "IsSynonym", "Accepted", "Synonyms", "InRedList", "InITIS", "FuzzyMatch",
                                         "Genus_level", "Genus", "Family", "Order")]


Birds_final <- Birds_final[, c("Original", "CorrectedTypos", "IsCorrected","IsAccepted", 
                                     "IsSynonym", "Accepted", "Synonyms", "InRedList", "InITIS", "FuzzyMatch",
                                     "Genus_level", "Genus", "Family", "Order")]


## Write final datasets: amphibians and mammals are final
write.csv(Mammals_final, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals.csv", row.names = FALSE)
write.csv(Amphibians_final, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians.csv", row.names = FALSE)
write.csv(Reptiles_final, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles.csv", row.names = FALSE)
write.csv(Birds_final, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds.csv", row.names = FALSE)



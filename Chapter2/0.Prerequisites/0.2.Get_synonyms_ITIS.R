## Here: For species that did not appear in the Red List, extract synonyms from ITIS

# #  P R E A M B L E  # #  

X <- c("dplyr", "taxize", "phytools", "stringr", "rredlist", "stringdist", "plyr", "pbapply", "GlobalOptions", "data.table", "ngram", "curl", "RCurl")
lapply(X, library, character.only=TRUE); rm(X)
`%nin%` = Negate(`%in%`)

opt <- options(iucn_redlist_key="ba30954f38deda075bd9b52495d65092ccf1b220b0c7c67a41465646e50ef72c")

source("Resolve_taxonomy_functions.R")

# #   Extracting synonyms - from ITIS  # # 

# Load data (list of synonyms)
Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_mammals_V2.csv")
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_amphibians_V2.csv")
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_birds_V2.csv")
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_reptiles_V2.csv")

# # Check that corrected typos all have 2 words, otherwise take original name (Amphibians, mammals and reptiles are OK, birds need to be corrected)
Syn_Birds$CorrectedTypos <- as.character(Syn_Birds$CorrectedTypos)
Syn_Birds$N_words <- sapply(strsplit(Syn_Birds$CorrectedTypos, " "), length)
Syn_Birds$CorrectedTypos[Syn_Birds$N_words==1] <- as.character(Syn_Birds$Original[Syn_Birds$N_words==1])

# Run the function
Syn_Mammals_ITIS <- Complement_ITIS(Syn_Mammals, Split = FALSE,1)
Syn_Birds_ITIS <- Complement_ITIS(Syn_Birds, Split = TRUE, 12)
Syn_Reptiles_ITIS <- Complement_ITIS(Syn_Reptiles, Split = FALSE,1)

Syn_Amphibians$Original <- as.character(Syn_Amphibians$Original)
Syn_Amphibians$CorrectedTypos <- as.character(Syn_Amphibians$CorrectedTypos)
Syn_Amphibians$CorrectedTypos[Syn_Amphibians$Original=="Xenopus [Silurana]"] <- "Xenopus tropicalis"
Syn_Amphibians <- Complement_ITIS(Syn_Amphibians)

# Save file after addition of ITIS synonyms
write.csv(Syn_Mammals_ITIS, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Mammals.csv", row.names = FALSE)
write.csv(Syn_Amphibians, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Amphibians.csv", row.names = FALSE)
write.csv(Syn_Reptiles_ITIS, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Reptiles.csv", row.names = FALSE)
write.csv(Syn_Birds_ITIS, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Birds.csv", row.names = FALSE)

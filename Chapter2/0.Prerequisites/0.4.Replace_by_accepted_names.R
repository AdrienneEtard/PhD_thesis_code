## Correct trait datasets, phylogenies and PREDICTS for taxonomy.

# #  P R E A M B L E  # #  

X <- c("dplyr", "taxize", "phytools", "stringr", "rredlist", "stringdist", "plyr", "pbapply", "GlobalOptions", "data.table", "ngram")
lapply(X, library, character.only=TRUE); rm(X)
`%nin%` = Negate(`%in%`)

opt <- options(iucn_redlist_key="ba30954f38deda075bd9b52495d65092ccf1b220b0c7c67a41465646e50ef72c")

source("Resolve_taxonomy_functions.R")

## Load synonym datasets
Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals.csv") 
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds.csv")
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles.csv")
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians.csv")


## Load datasets
Predicts <-  readRDS("../../Data/PREDICTS_database.rds")
Predicts <- subset(Predicts, Class %in% c("Aves", "Amphibia", "Mammalia", "Reptilia"))
Predicts$Best_guess_binomial <- as.character(Predicts$Best_guess_binomial)

# Non identifiable species
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Duellmanohyla eutisanota"] <- ""

# Phylogenies
PhyloMammal <- read.newick("../../Data/Mammals/Phylogenies/TTOL_mammals_smoothed_interpolated.nwk") %>% .Format_tiplabels()
PhyloAmphibian <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()
PhyloBird <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk") %>% .Format_tiplabels()
PhyloReptile <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()

# Traits
Myhrvold <- read.csv("../../Data/Amniotes_Myhrvold_2015/Amniote_Database_Aug_2015.csv", sep=",")
Myhrvold$Binomial_name <- paste(Myhrvold$genus, Myhrvold$species, sep=" ")
MyhrvoldMammal <- subset(Myhrvold, class=="Mammalia")
MyhrvoldBird <- subset(Myhrvold, class=="Aves")
MyhrvoldReptile <- subset(Myhrvold, class=="Reptilia")

Pantheria <- read.csv("../../Data/Mammals/PanTHERIA/Pantheria_1_0_WR05_Aug2008.csv", sep=",")
Pacifici <- read.csv("../../Data/Mammals/PacificiMammals.csv", sep=",")

## Put mammalian diet information together
Kissling <- read.csv("../../Data/Mammals/Kissling_Mammal_diet_2014.csv", sep=",")
Kissling$Binomial_name <- paste(Kissling$Genus, Kissling$Species, sep=" ")
Kissling <- subset(Kissling, Binomial_name != "Mico sp. nov.")
Kissling <- Kissling %>% select(-TaxonID, -TaxonomicNote, -FillCode)
Kissling <- Kissling[, c(28, 1:27)]

MammalDIET2 <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET2.csv")[, c(1:28)]
colnames(MammalDIET2)[28] <- "DataSource"

MammalDIET2_supp <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET_2_supp.csv")
MammalDIET2_supp$TrophicLevel <- as.character(MammalDIET2_supp$TrophicLevel)
MammalDIET2_supp <- MammalDIET2_supp %>%
  mutate(TrophicLevel=ifelse(TrophicLevel %in% c("Carnivore", "Omnivore", "Herbivore"), TrophicLevel, NA)) %>%
  filter(!is.na(TrophicLevel)) %>%
  mutate(Mammal=ifelse(Mammal %in% c(0:3), Mammal, NA)) %>%
  mutate(MammalEater=ifelse(MammalEater %in% c(0:3), MammalEater, NA)) %>%
  mutate(Insectivore=ifelse(Insectivore %in% c(0:3), Insectivore, NA))%>%
  mutate(Frugivore=ifelse(Frugivore %in% c(0:3), Frugivore, NA))%>%
  mutate(Granivore=ifelse(Granivore %in% c(0:3), Granivore, NA)) %>%
  mutate(Folivore=ifelse(Folivore %in% c(0:3), Folivore, NA))
MammalDIET2_supp[, c(6:27)] <- droplevels(MammalDIET2_supp[, c(6:27)])
  
MammalDIET2 <- rbind(MammalDIET2, MammalDIET2_supp)
MammalDIET2$Binomial <- as.character(MammalDIET2$Binomial)
MammalDIET2$Binomial[MammalDIET2$Binomial=="Rhinolophus hildebrandtii"] <- "Rhinolophus hildebrandti"
MammalDIET2$Binomial[MammalDIET2$Binomial=="Tadarida bivittatus"] <- "Tadarida bivittata"
colnames(MammalDIET2)[1] <- "Binomial_name"

Y <- intersect(MammalDIET2$Binomial_name, Kissling$Binomial_name)
Kissling <- Kissling %>% filter(Binomial_name %nin% Y)
MammalDIET <- rbind(Kissling, MammalDIET2) # dataset to use

Butchart_BM <- read.csv("../../Data/Birds/Butchart_BM.csv", sep=",")
Butchart_GL <- read.csv("../../Data/Birds/ButchartGenerationLength.csv", sep=",")
Butchart_GL <- subset(Butchart_GL, Binomial !="")
Sekercioglu_Diet <- read.csv("../../Data/Birds/SekerciogluDiet.csv", sep=",")

Scharf <- read.csv("../../Data/Reptiles/Scharf.csv", sep=",")
Meiri_0 <- read.csv("../../Data/Reptiles/MeiriReptileMasses.csv", sep=",")
Meiri_0 <- subset(Meiri_0, Rank=="Species")

# Reptiles traits to add, a posteriori extract synonyms for species for which it has not been done before
Vidan <- read.csv("../../Data/Reptiles/Vidan2017_Dielactivity.csv")

Stark<- read.csv("../../Data/Reptiles/Stark2018_GEB_longevity.csv")
colnames(Stark)[1] <- "species"
Stark_longevity <- Stark %>% 
  filter(species !="")

Schwarz <- read.csv("../../Data/Reptiles/Schwarz_Meiri_GEB_2017.csv")
Novosolov <- read.csv("../../Data/Reptiles/Novosolov_2017_GEB.csv") %>%
  filter(Taxonomic.group!="Mammals") %>%
  filter(Taxonomic.group!="Birds")
Novosolov$Taxonomic.group <- droplevels(Novosolov$Taxonomic.group)
Novosolov_2 <- read.csv("../../Data/Reptiles/Novosolov_GEB_2013.csv")
Slavenko <- read.csv("../../Data/Reptiles/Body_sizes_of_all_extant_reptiles_Slavenko_2016_GEB.csv") 
Meiri <- read.csv("../../Data/Reptiles/Meiri_2015_Evolutionary_Biology.csv")

## reviewers suggested adding Meiri (2018, GEB) and Feldman et al (GEB, 2016) to the trait compilation
Meiri2018 <- read.csv("../../Data/Reptiles/MeiriGEB2018.csv")
Feldman2016 <- read.csv("../../Data/Reptiles/FeldmanGEB2016.csv")



Amphibio <- read.csv("../../Data/Amphibians/AmphiBIO_v1.csv", sep=",")
Amphibio$Species <- as.character(Amphibio$Species)
Amphibio$Species[Amphibio$Species=="Duttaphrynus pariet alis"] <- "Duttaphrynus parietalis"
Cooper <- read.csv("../../Data/Amphibians/Cooper2008.csv", sep=",")
Cooper <- subset(Cooper, Binomial!="")
Senior <- read.csv("../../Data/Amphibians/Senior_svl_data.csv", sep=",")
Senior <- subset(Senior, Rank=="Species")
Bickford <- read.csv("../../Data/Amphibians/Bickford.csv", sep=",")
Bickford$Binomial_name <-  paste(Bickford$Genus, Bickford$Species, sep=" ")

Mammal_range <- read.csv("../../Data/Range_sizes/mammal_range_areas.csv", sep=",")
Amphibian_range <- read.csv("../../Data/Range_sizes/amphibian_range_areas.csv", sep=",")
Bird_range <- read.csv("../../Data/Range_sizes/bird_range_areas.csv", sep=",")
Reptile_range <- read.csv("../../Data/Range_sizes/reptile_range_areas.csv", sep=",")
Reptile_range <- subset(Reptile_range, Binomial_name!="Chelonia mydas Hawaiian subpopulation")

IUCN_mammal <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Mammals_20180726.csv", sep=",")
IUCN_amphibian <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Amphibians_20180726.csv", sep=",")
IUCN_bird <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Birds_20180726.csv", sep=",")
IUCN_reptile <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Reptiles_20180726.csv", sep=",")

# Elton traits
Elton_birds <- read.csv("../../Data/Birds/EltonTraits_Birds.csv")
Elton_mammals <- read.csv("../../Data/Mammals/EltonTraits_Mammals.csv")

## Replace vernacular names by scientific names

# Amphibians
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana cafe"] <- "Eleutherodactylus maurus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana hoja"] <- "Noblella lochites"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana hojarasca"] <- "Craugastor angelicus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana verde"] <- "Lithobates palmipes"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Arthroleptis sp"] <- ""

# Reptiles
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Serp cafe"] <- "Ninia sebae"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Serpiente cafe"] <- "Ninia sebae"

# Mammals
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Samango monkey"] <- "Cercopithecus albogularis"

# Birds
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Black backed starling"] <- "Acridotheres melanopterus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Black backed weaver"] <- "Ploceus bicolor"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Blue rock pigeon"] <- "Columba livia"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Blue whistling thrush"] <- "Myophonus caeruleus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Brown crested tit"] <- "Lophophanes cristatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Chestnut beilled nuthatch"] <- "Sitta cinnamoventris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Large yellow nape"] <- "Amazona auropalliata"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Common stone chat"] <- "Saxicola torquatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Rufous backed shrike"] <- "Lanius schach"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="White breasted forktail"] <- "Enicurus leschenaulti"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="White breasted kingfisher"] <- "Halcyon smyrnensis"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Jungle babbler"] <- "Turdoides striata"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Grey bulbul"] <- "Pycnonotus cyaniventris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Pied wagtail"] <- "Motacilla alba"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian rufous babbler"] <-"Argya subrufa"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Large wagtail"] <- "Motacilla madaraspatensis"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Common swallow"] <- "Hirundo rustica"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Eurasian tree creeper"] <- "Certhia familiaris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Red headed tit"] <- "Aegithalos iredalei"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Red jungle fowl"] <- "Gallus gallus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian myna"] <-"Acridotheres tristis"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian robin"] <- "Saxicoloides fulicatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian tree pie"] <-"Dendrocitta vagabunda"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Sparrow hawk"] <-"Accipiter nisus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Creasted bunting"] <-"Emberiza lathami"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Yellow eyed warbler"] <- "Chrysomma sinense"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Canary blue flycatcher"] <- "Eumyias thalassinus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Hawk booted eagle"] <- "Hieraaetus pennatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Yellow backed greenbul"] <- "Chlorocichla flaviventris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="White piegon"] <- "Columba livia"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Simla crested tit"] <-"Lophophanes dichrous"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Cacomantis esculena"] <- "Collocalia esculenta"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Bush lark"] <- ""
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Imperial pigeon"] <- ""
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Tailor bird"] <- ""
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Paradise flycatcher"] <- ""


# ## Correct a few noticeable taxonomic errors
Slavenko$Binomial <- Slavenko$Binomial %>% as.character()
Slavenko$Binomial[Slavenko$Binomial=="Chelonoidis nigra (complex)"] <- "Chelonoidis nigra"

# Amphibians
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Kaloula pleurostigma"] <- "Kalophrynus pleurostigma"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Pristimantis bransfordii"] <- "Craugastor bransfordii"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Pristimantis fitzingerii"] <-  "Craugastor fitzingeri"
# Reptiles
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Mabuya elegans"] <- "Lygosoma punctata"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Cryptoblepharus carnabyi"] <- "Cryptoblepharus australis"
# Birds
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Dicaeum maculatus"] <- "Prionochilus maculatus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Dicaeum percussus"] <- "Prionochilus percussus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Picoides leucotos"] <- "Dendrocopos leucotos"

# Mammals
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="sus scrofa"] <- "Sus scrofa"

# Select data where species known by binomial names
Predicts <- subset(Predicts, Best_guess_binomial != "")



 #  #  #  #  #  #  #   REPLACE NAMES   #  #  #  #  #  #  # 

## NB: minor problem: some names do not appear in the synonym dataset (weird!!??). Those are
# for mammals: Desmodus draculae (from MammalDiet)
# plus one taxonomic error 


# Trait datasets
print("Myhrvold")
MyhrvoldMammal <- Replace_by_accepted_name(Syn_Mammals, MyhrvoldMammal, "Binomial_name")
MyhrvoldBird <- Replace_by_accepted_name(Syn_Birds, MyhrvoldBird, "Binomial_name")
MyhrvoldReptile <- Replace_by_accepted_name(Syn_Reptiles, MyhrvoldReptile, "Binomial_name")
Myhrvold <- rbind(MyhrvoldBird, MyhrvoldMammal, MyhrvoldReptile)

print("Pantheria")
Pantheria <- Replace_by_accepted_name(Syn_Mammals, Pantheria, "MSW05_Binomial")
print("Pacifici")
Pacifici <- Replace_by_accepted_name(Syn_Mammals, Pacifici, "Scientific_name")
#Kissling <- Replace_by_accepted_name(Syn_Mammals, Kissling, "Binomial_name")
print("Mammal diet: Kissling")
MammalDIET <- Replace_by_accepted_name(Syn_Mammals, MammalDIET, "Binomial_name")
print("Elton mammals")
Elton_mammals <- Replace_by_accepted_name(Syn_Mammals, Elton_mammals, "Scientific")

print("Butchart BM")
Butchart_BM <- Replace_by_accepted_name(Syn_Birds, Butchart_BM, "Sci.name")
print("Butchart GL")
Butchart_GL <- Replace_by_accepted_name(Syn_Birds, Butchart_GL, "Binomial")
print("Sekercioglu")
Sekercioglu_Diet <- Replace_by_accepted_name(Syn_Birds, Sekercioglu_Diet, "Latin")
print("Elton birds")
Elton_birds <- Elton_birds %>% filter(Scientific!="")
Elton_birds <- Replace_by_accepted_name(Syn_Birds, Elton_birds, "Scientific")

print("Scharf")
Scharf <- Replace_by_accepted_name(Syn_Reptiles, Scharf, "species")
print("Meiri_0")
Meiri_0 <- Replace_by_accepted_name(Syn_Reptiles, Meiri_0, "Taxon")
print("Meiri")
Meiri <- Replace_by_accepted_name(Syn_Reptiles, Meiri, "species")
print("Vidan")
Vidan <- Replace_by_accepted_name(Syn_Reptiles, Vidan, "Species")
print("Stark")
Stark <- Stark %>% filter(species!="")
Stark <- Replace_by_accepted_name(Syn_Reptiles, Stark, "species")
print("Novosolov")
Novosolov <- Replace_by_accepted_name(Syn_Reptiles, Novosolov, "Binomial")
print("Novosolov 2")
Novosolov_2 <- Replace_by_accepted_name(Syn_Reptiles, Novosolov_2, "species")
print("Schwarz")
Schwarz <- Replace_by_accepted_name(Syn_Reptiles, Schwarz, "Species")
print("Slavenko")
Slavenko <- Replace_by_accepted_name(Syn_Reptiles, Slavenko, "Binomial")


## NEW ADDITIONS, suggested by peer-reviewers - tax. information added a posterior, work is done here
print("Meiri2018GEB")
Meiri2018 <- Replace_by_accepted_name(Syn_Reptiles, Meiri2018, "Binomial")
Meiri2018c <- Meiri2018$Data
Meiri_NOMATCH <- Meiri2018$No_match

print("Feldman")
Feldman2016 <- Feldman2016 %>% filter(binomial!="")
Feldman2016 <- Feldman2016 %>% filter(valid!="")
Feldman2016 <- Replace_by_accepted_name(Syn_Reptiles, Feldman2016, "binomial")
Feldman2016c <- Feldman2016$Data 
Feldman_NOMATCH <- Feldman2016$No_match

No_match <- unique(c(as.character(Meiri_NOMATCH), as.character(Feldman_NOMATCH)))
No_match <- No_match[order(No_match)] %>% 
  as.data.frame() %>% 
  setNames(., "CorrectedTypos")
No_match$CorrectedTypos <- as.character(No_match$CorrectedTypos)

RL_no_match <- RunSyn(No_match)
RL_no_match$Original <- as.character(RL_no_match$CorrectedTypos)
RL_no_match_ITIS <- Complement_ITIS(RL_no_match, FALSE, 1)
write.csv(RL_no_match_ITIS, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Complements_reptiles.csv", row.names=FALSE)

RL_no_match_ITIS <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Complements_reptiles.csv")
RL_no_match_ITIS$Accepted <- as.character(RL_no_match_ITIS$Accepted)
RL_no_match_ITIS$Accepted[is.na(RL_no_match_ITIS$Accepted)] <- RL_no_match_ITIS$CorrectedTypos[is.na(RL_no_match_ITIS$Accepted)] %>% as.character()
RL_no_match_ITIS$Genus <- word(RL_no_match_ITIS$Accepted, 1)

# add family and order information
unique(RL_no_match_ITIS$Genus) %>%  length()
Order_Family_info <- unique(Syn_Reptiles[, c("Order", "Family", "Genus")]) %>% 
  filter(Genus %in% unique(RL_no_match_ITIS$Genus) )

To_add <- setdiff(unique(RL_no_match_ITIS$Genus), Order_Family_info$Genus) %>% 
  as.data.frame() %>% 
  setNames(., "Genus")# tax. info to be added manually for these generq (order and family information)
# [1] "Andinosaura"  "Bolyeria" "Brachyseps"    "Cachryx"       "Carinascincus" "Flexiseps"     "Gelanesaurus"  "Glaucomastix"  "Gowidon"      
# [9] "Grandidierina" "Diploderma"    "Lubuya"        "Oreosaurus"    "Phaeoscincus"  "Psilops"       "Rondonops"

To_add$Order <- c("Squamata", "Squamata", "Squamata", "Squamata", "Squamata", "Squamata", "Squamata", "Squamata", 
                  "Squamata", "Squamata", "Squamata", "Squamata", "Squamata", "Squamata", "Squamata", "Squamata") %>%
  toupper()

To_add$Family <- c("Gymnophthalmidae", "Bolyeriidae","Scincidae", "Iguanidae", "Scincidae", "Scincidae", "Gymnophthalmidae", "Teiidae", "Agamidae",
                   "Scincidae", "Iguania", "Scincidae", "Gymnophthalmidae", "Scincidae", "Gymnophthalmidae", "Gymnophthalmidae") %>% 
  toupper

Order_Family_info <- rbind(Order_Family_info, To_add)
RL_no_match_ITIS <- RL_no_match_ITIS %>% 
  dplyr::select(-Family, -Order)
RL_ITIS_additions <- join(RL_no_match_ITIS, Order_Family_info, by="Genus")
RL_ITIS_additions$Manual_edits <- NA
RL_ITIS_additions$Notes <- "Added a posteriori - additions of two  new sources in the compilation"
RL_ITIS_additions$IsCorrected <- NA # here names were not corrected for typos
RL_ITIS_additions$Genus_level <- FALSE 

# lapply(RL_ITIS_additions[, "Accepted"] %>% as.list, wordcount) %>%  unlist() %>% unique
# (this line checks that no species are entered at the genus level)

Syn_Reptiles2 <- rbind(Syn_Reptiles, RL_ITIS_additions)
write.csv(Syn_Reptiles2, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles2.csv")
# read again and harmonize taxonomy
Meiri2018 <- read.csv("../../Data/Reptiles/MeiriGEB2018.csv")
Feldman2016 <- read.csv("../../Data/Reptiles/FeldmanGEB2016.csv")
Feldman2016 <- Feldman2016 %>% filter(binomial!="")
Feldman2016 <- Feldman2016 %>% filter(valid!="")

Meiri2018 <- Replace_by_accepted_name(Syn_Reptiles2, Meiri2018, "Binomial")
Feldman2016 <- Replace_by_accepted_name(Syn_Reptiles2, Feldman2016, "binomial")



##########
print("Amphibio")
Amphibio <- Replace_by_accepted_name(Syn_Amphibians, Amphibio, "Species")
print("Cooper")
Cooper <- Replace_by_accepted_name(Syn_Amphibians, Cooper, "Binomial")
print("Senior")
Senior <- Replace_by_accepted_name(Syn_Amphibians, Senior, "Taxon")
print("Bickford")
Bickford <- Replace_by_accepted_name(Syn_Amphibians, Bickford, "Binomial_name")

print("Ranges")
Mammal_range <- Replace_by_accepted_name(Syn_Mammals, Mammal_range, "Species")
Amphibian_range <- Replace_by_accepted_name(Syn_Amphibians, Amphibian_range, "Species")
Bird_range <- Replace_by_accepted_name(Syn_Birds, Bird_range, "Species")
Reptile_range <- Replace_by_accepted_name(Syn_Reptiles, Reptile_range, "Binomial_name")

print("IUCN data")
IUCN_amphibian <-  Replace_by_accepted_name(Syn_Amphibians, IUCN_amphibian, "binomial")
IUCN_bird <-  Replace_by_accepted_name(Syn_Birds, IUCN_bird, "binomial")
IUCN_reptile <-  Replace_by_accepted_name(Syn_Reptiles, IUCN_reptile, "binomial")
IUCN_mammal <-  Replace_by_accepted_name(Syn_Mammals, IUCN_mammal, "binomial")

# Predicts dataset ## problem wiht species replacement here
PredictsGenusLevel <- subset(Predicts, Best_guess_binomial == "")

PredictsAmphibians <- subset(Predicts, Class=="Amphibia")
PredictsReptiles <- subset(Predicts, Class=="Reptilia")
PredictsMammals <- subset(Predicts, Class=="Mammalia")
PredictsBirds <- subset(Predicts, Class=="Aves")

print("PREDICTS")
PredictsAmphibians <- ReplacePredicts(PredictsAmphibians, Syn_Amphibians)
PredictsReptiles <- ReplacePredicts(PredictsReptiles, Syn_Reptiles)
PredictsBirds <- ReplacePredicts(PredictsBirds, Syn_Birds)
PredictsMammals <- ReplacePredicts(PredictsMammals,Syn_Mammals)

Predicts_Vertebrates <- rbind(PredictsAmphibians, PredictsBirds, PredictsMammals, PredictsReptiles, PredictsGenusLevel)
Predicts_Vertebrates <- Predicts_Vertebrates[order(Predicts_Vertebrates$SSBS),]

# Phylogenies
print("Phylogenies")
PhyloMammal <- Replace_by_accepted_name(Syn_Mammals, PhyloMammal)
PhyloBird <- Replace_by_accepted_name(Syn_Birds, PhyloBird)
PhyloReptile <- Replace_by_accepted_name(Syn_Reptiles, PhyloReptile)
PhyloAmphibian <- Replace_by_accepted_name(Syn_Amphibians, PhyloAmphibian)


## Saving processed datasets
write.csv(Meiri2018, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Meiri2018GEB.csv", row.names=FALSE)
write.csv(Feldman2016, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Feldman2016.csv", row.names=FALSE)
write.csv(Myhrvold, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Myhrvold.csv", row.names=FALSE)
write.csv(Pantheria, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Pantheria.csv", row.names=FALSE)
write.csv(Pacifici, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Pacifici.csv", row.names=FALSE)
write.csv(MammalDIET, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/MammalDIET.csv", row.names=FALSE)
write.csv(Butchart_BM, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Butchart_BM.csv", row.names=FALSE)
write.csv(Butchart_GL, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Butchart_GL.csv", row.names=FALSE)
write.csv(Sekercioglu_Diet, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Sekercioglu_Diet.csv", row.names=FALSE)
write.csv(Scharf, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Scharf.csv", row.names=FALSE)
write.csv(Meiri_0, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Meiri_0.csv", row.names=FALSE)
write.csv(Meiri, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Meiri.csv", row.names=FALSE)
write.csv(Vidan, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Vidan.csv", row.names=FALSE)
write.csv(Stark, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Stark.csv", row.names=FALSE)
write.csv(Novosolov, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Novosolov.csv", row.names=FALSE)
write.csv(Novosolov_2, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Novosolov_2.csv", row.names=FALSE)
write.csv(Schwarz, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Schwarz.csv", row.names=FALSE)
write.csv(Slavenko, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Slavenko.csv", row.names=FALSE)
write.csv(Amphibio, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Amphibio.csv", row.names=FALSE)
write.csv(Cooper, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Cooper.csv", row.names=FALSE)
write.csv(Senior, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Senior.csv", row.names=FALSE)
write.csv(Bickford, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Bickford.csv", row.names=FALSE)
write.csv(Mammal_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Mammal_range.csv", row.names=FALSE)
write.csv(Amphibian_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Amphibian_range.csv", row.names=FALSE)
write.csv(Bird_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Bird_range.csv", row.names=FALSE)
write.csv(Reptile_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Reptile_range.csv", row.names=FALSE)
write.csv(IUCN_amphibian, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_amphibians.csv", row.names=FALSE)
write.csv(IUCN_bird, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_birds.csv", row.names=FALSE)
write.csv(IUCN_mammal, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_mammals.csv", row.names=FALSE)
write.csv(IUCN_reptile, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_reptiles.csv", row.names=FALSE)
write.csv(Elton_birds, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Elton_birds.csv", row.names=FALSE)
write.csv(Elton_mammals, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Elton_mammals.csv", row.names=FALSE)

saveRDS(Predicts_Vertebrates, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

write.tree(PhyloAmphibian, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloAmphibians.nwk")
write.tree(PhyloBird, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloBirds.nwk")
write.tree(PhyloMammal, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloMammals.nwk")
write.tree(PhyloReptile, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloReptiles.nwk")

## phylogeny for all vertebrates
PredictsCor <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")
AllSyn <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/SingleFinalDataset.csv")
AllVert <-  read.tree("../../Data/Phylogenies/TTOL_animals_unsmoothed.nwk") %>% .Format_tiplabels()
AllVert$tip.label[AllVert$tip.label=="Xenopus "] <- "Xenopus sp."
PhyloAll <- Replace_by_accepted_name(AllSyn, AllVert)

write.tree(PhyloAll, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloVertebrates.nwk")


# to_drop <- dplyr::setdiff(AllVert$tip.label, PredictsCor$Best_guess_binomial) %>% unlist()
# Phylo <- ape::drop.tip(AllVert, to_drop)
# length(Phylo$tip.label)
# length(setdiff(unique(PredictsCor$Best_guess_binomial), Phylo$tip.label)) # 1301 species in PREDICTS are not represented in the verterbate phylogeny

## after correcting
# to_drop <- dplyr::setdiff(PhyloAll$tip.label, PredictsCor$Best_guess_binomial) %>% unlist()
# Phylo <- ape::drop.tip(PhyloAll, to_drop)
# length(Phylo$tip.label)
# length(setdiff(unique(PredictsCor$Best_guess_binomial), Phylo$tip.label)) # 863 species in PREDICTS are not represented in the verterbate phylogeny


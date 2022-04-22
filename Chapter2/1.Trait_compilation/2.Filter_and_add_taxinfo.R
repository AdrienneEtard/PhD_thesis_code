## Add family / order information to compiled trait datasets.
## Filter out marine mammals
## Check which species have completeness = 0%


## PREAMBLE
X <- c("dplyr", "stringr")
lapply(X, library, character.only=TRUE); rm(X)
`%nin%` = Negate(`%in%`)

Completeness0 <- function(TraitDF, Traits) {
  # completeness
  TraitDF$Percent <- apply(TraitDF[,Traits], 1, function(y) sum(!is.na(y)))
  TraitDF$Percent <- TraitDF$Percent / length(Traits) * 100 
  
  Completeness_0 <- TraitDF$Best_guess_binomial[TraitDF$Percent==0]
  print(length(Completeness_0))
  return(Completeness_0)
}

## Traits
Traits <- c("Body_mass_g",
            "Longevity_d",
            "Litter_size", 
            "Habitat_breadth_IUCN",
            "Specialisation",
            "Trophic_level",
            "Diel_activity")
            #"Primary_diet",
            #"Diet_breadth")

TraitsReptiles <- c(Traits[1:7], "Adult_svl_cm", "Maturity_d", "Max_longevity_d")
TraitsAmphibians <- c(Traits[-which(Traits=="Longevity_d")], "Body_length_mm", "Max_longevity_d", "Maturity_d")
TraitsBirds <- c(Traits, "Generation_length_d", "Adult_svl_cm", "Maturity_d", "Max_longevity_d")
TraitsMammals <- TraitsBirds

## Load data 

# Corrected
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Reptiles.csv")

# Uncorrected
U_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Mammals.csv")
U_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Birds.csv")
U_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Amphibians.csv")
U_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Reptiles.csv")

# for manuscript
ms_Mammals <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Mammals.csv")
ms_Birds <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Birds.csv")
ms_Amphibians <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Amphibians.csv")
ms_Reptiles <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Reptiles.csv")
Ums_Mammals <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Mammals.csv")
Ums_Birds <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Birds.csv")
Ums_Amphibians <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Amphibians.csv")
Ums_Reptiles <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Reptiles.csv")


# Synonym datasets
Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals_final.csv")
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians_final.csV")
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds_final.csv")
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles_final.csv")


# Add family and order to the data

Add_taxinfo <- function(TraitDF, SynDF, Cor) {
 
  if(Cor){
    colnames(SynDF)[colnames(SynDF)=="Accepted"] <- "Best_guess_binomial"
  }
  
  else {
    colnames(SynDF)[colnames(SynDF)=="Original"] <- "Best_guess_binomial"
  }
  
  SynDF <- SynDF[, c("Order", "Family", "Genus", "Best_guess_binomial")]
  SynDF <- unique(SynDF)
  Out <- dplyr::left_join(TraitDF, SynDF, by="Best_guess_binomial")
  Out <- Out[order(Out$Best_guess_binomial),]
  
  Out1 <- Out[,c("Order", "Family", "Genus", "Best_guess_binomial")]
  Out2 <- Out %>% dplyr::select(-Order, -Family, -Genus, -Best_guess_binomial)

  return(cbind(Out1, Out2))
  }

# Add family and order information to the datasets
Mammals <- Add_taxinfo(Mammals, Syn_Mammals, TRUE)
Mammals$Genus[Mammals$Best_guess_binomial=="Desmodus draculae"] <- "Desmodus"
Mammals$Family[Mammals$Best_guess_binomial=="Desmodus draculae"] <- "PHYLLOSTOMIDAE"
Mammals$Order[Mammals$Best_guess_binomial=="Desmodus draculae"] <- "CHIROPTERA"

Amphibians <- Add_taxinfo(Amphibians, Syn_Amphibians, TRUE)
Reptiles <- Add_taxinfo(Reptiles, Syn_Reptiles, TRUE)
Birds <- Add_taxinfo(Birds, Syn_Birds, TRUE)

U_Mammals <- Add_taxinfo(U_Mammals, Syn_Mammals, FALSE)
U_Birds <- Add_taxinfo(U_Birds, Syn_Birds, FALSE) %>% filter(Best_guess_binomial != "")
U_Amphibians <- Add_taxinfo(U_Amphibians, Syn_Amphibians, FALSE)
U_Reptiles <- Add_taxinfo(U_Reptiles, Syn_Reptiles, FALSE)

## for ms
ms_Mammals <- Add_taxinfo(ms_Mammals, Syn_Mammals, TRUE)
ms_Mammals$Genus[ms_Mammals$Best_guess_binomial=="Desmodus draculae"] <- "Desmodus"
ms_Mammals$Family[ms_Mammals$Best_guess_binomial=="Desmodus draculae"] <- "PHYLLOSTOMIDAE"
ms_Mammals$Order[ms_Mammals$Best_guess_binomial=="Desmodus draculae"] <- "CHIROPTERA"

ms_Amphibians <- Add_taxinfo(ms_Amphibians, Syn_Amphibians, TRUE)
ms_Reptiles <- Add_taxinfo(ms_Reptiles, Syn_Reptiles, TRUE)
ms_Birds <- Add_taxinfo(ms_Birds, Syn_Birds, TRUE)

Ums_Mammals <- Add_taxinfo(Ums_Mammals, Syn_Mammals, FALSE)
Ums_Birds <- Add_taxinfo(Ums_Birds, Syn_Birds, FALSE) %>% filter(Best_guess_binomial != "")
Ums_Amphibians <- Add_taxinfo(Ums_Amphibians, Syn_Amphibians, FALSE)
Ums_Reptiles <- Add_taxinfo(Ums_Reptiles, Syn_Reptiles, FALSE)


## mammals: filter non-terrestrial species out
Mammals <- Mammals %>%
  filter(Order %nin% c("SIRENIA")) %>%
  filter(Family %nin% c("OTARIIDAE", "PHOCIDAE", "ODOBENIDAE", 
                        "BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", 
                        "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE",
                        "ESCHRICHTIIDAE", "INIIDAE", "PHYSETERIDAE","LIPOTIDAE",
                        "PHOCOENIDAE", "PLATANISTIDAE", "PONTOPORIIDAE"))

U_Mammals <- U_Mammals %>%
  filter(Order %nin% c("SIRENIA")) %>%
  filter(Family %nin% c("OTARIIDAE", "PHOCIDAE", "ODOBENIDAE", 
                        "BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", 
                        "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE",
                        "ESCHRICHTIIDAE", "INIIDAE", "PHYSETERIDAE","LIPOTIDAE",
                        "PHOCOENIDAE", "PLATANISTIDAE", "PONTOPORIIDAE"))

ms_Mammals <- ms_Mammals %>%
  filter(Order %nin% c("SIRENIA")) %>%
  filter(Family %nin% c("OTARIIDAE", "PHOCIDAE", "ODOBENIDAE", 
                        "BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", 
                        "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE",
                        "ESCHRICHTIIDAE", "INIIDAE", "PHYSETERIDAE","LIPOTIDAE",
                        "PHOCOENIDAE", "PLATANISTIDAE", "PONTOPORIIDAE"))

Ums_Mammals <- Ums_Mammals %>%
  filter(Order %nin% c("SIRENIA")) %>%
  filter(Family %nin% c("OTARIIDAE", "PHOCIDAE", "ODOBENIDAE", 
                        "BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", 
                        "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE",
                        "ESCHRICHTIIDAE", "INIIDAE", "PHYSETERIDAE","LIPOTIDAE",
                        "PHOCOENIDAE", "PLATANISTIDAE", "PONTOPORIIDAE"))



# Check which species have 0% completeness

Completeness0(ms_Mammals, TraitsMammals)
# [1] "Amblyrhiza inundata"      "Clidomys osborni"         "Nesophontes longirostris"
# [4] "Nesophontes submicrus"    "Nesophontes superstes"    "Quemisia gravis"

Completeness0(ms_Birds, TraitsBirds)
# [1] Ara atwoodi                Ara erythrocephala         Ara gossei                
# [4] Ara guadeloupensis         Argusianus bipunctatus     Dysmoropelia dekarchiskos 
# [7] Falco buboisi              Fulica newtoni             Gallicolumba norfolciensis

Completeness0(ms_Amphibians, TraitsAmphibians)
# [1] "Brachytarsophrys platyparietus" "Bufo kabischi"                 
# [3] "Bufo minshanicus"               "Bufo wolongensis"              
# [5] "Dryaderces inframaculata"       "Eleutherodactylus dixoni"      
# [7] "Fejervarya frithii"             "Fejervarya multistriata"       
# [9] "Hyla altipotens"                "Hyla angustilineata"           
# [11] "Hyla auraria"                   "Hyla debilis"                  
# [13] "Hyla fimbrimembra"              "Hyla hazelae"                  
# [15] "Hyla miliaria"                  "Hyla pellita"                  
# [17] "Hyla picadoi"                   "Hyla vigilans"                 
# [19] "Leptobrachium echinatum"        "Leptodactylus hallowelli"      
# [21] "Lithobates subaquavocalis"      "Phaeognathus ainsworthi"       
# [23] "Philautus pleurotaenia"         "Philautus semiruber"           
# [25] "Philautus wynaadensis"          "Platymantis rhipiphalcus"      
# [27] "Pristimantis kelephus"          "Proceratophrys fryi"           
# [29] "Pseudepidalea pewzowi"          "Pseudopaludicola mirandae"     
# [31] "Pyxicephalus cordofanus"        "Quasipaa delacouri"            
# [33] "Rana attigua"                   "Rana aurata"                   
# [35] "Rana charlesdarwini"            "Rana chitwanensis"             
# [37] "Rana emeljanovi"                "Rana faber"                    
# [39] "Rana garoensis"                 "Rana humeralis"                
# [41] "Rana lungshengensis"            "Rana margariana"               
# [43] "Rana minima"                    "Rana oatesii"                  
# [45] "Rana sangzhiensis"              "Rana spinulosa"                
# [47] "Rana swinhoana"                 "Rana tientaiensis"             
# [49] "Rana tytleri"                   "Rana volkerjane"               
# [51] "Rana weiningensis"              "Rana wuchuanensis"             
# [53] "Rhacophorus gongshanensis"      "Rhacophorus taronensis"        
# [55] "Rhacophorus zhaojuensis"        "Rhinella beebei"               
# [57] "Rhinella sima"                  "Scutiger ruginosus"            
# [59] "Sphaerotheca maskeyi"           "Sphaerotheca swani"            
# [61] "Theloderma albopunctatum"       "Zachaenus roseus"   

Completeness0(ms_Reptiles, TraitsReptiles)
# [1] Anolis eewi                      Anolis ricordi                  
# [3] Chabanaudia boulengeri           Cnemaspis nigridius             
# [5] Colobosaura kraepelini           Cyrtodactylus aravallensis      
# [7] Geckoella yakhuna                Latastia carinata               
# [9] Liolaemus filiorum               Rhoptropus braconnieri          
# [11] Sphenomorphus amblyplacodes      Trachylepis hildae              
# [13] Trioceros tremperi               Alopoglossus embera             
# [15] Amphisbaena acrobeles            Amphisbaena bilabialata         
# [17] Amphisbaena filiformis           Amphisbaena metallurga          
# [19] Andinosaura crypta               Andinosaura hyposticta          
# [21] Andinosaura laevis               Andinosaura petrorum            
# [23] Andinosaura vieta                Anolis demissus                 
# [25] Anolis maia                      Anolis peruensis                
# [27] Aprasia wicherina                Brachymeles dalawangdaliri      
# [29] Brachymeles ilocandia            Brachymeles ligtas              
# [31] Brachyseps anosyensis            Brachyseps splendidus           
# [33] Cachryx alfredschmidti           Calotes nigriplicatus           
# [35] Calumma linotum                  Carlia insularis                
# [37] Carlia isostriacantha            Chondrodactylus pulitzerae      
# [39] Cnemaspis aceh                   Cnemaspis minang                
# [41] Cnemaspis pagai                  Cnemaspis rajakarunai           
# [43] Ctenotus superciliaris           Cyrtodactylus equestris         
# [45] Cyrtodactylus hitchi             Cyrtodactylus lenya             
# [47] Cyrtodactylus payarhtanensis     Cyrtodactylus pharbaungensis    
# [49] Cyrtodactylus rex                Dibamus floweri                 
# [51] Diploderma ngoclinensis          Eremiascincus rubiginosus       
# [53] Flexiseps alluaudi               Flexiseps andranovahensis       
# [55] Flexiseps crenni                 Flexiseps decaryi               
# [57] Flexiseps elongatus              Flexiseps johannae              
# [59] Flexiseps mandokava              Flexiseps meva                  
# [61] Flexiseps tsaratananensis        Gehyra granulum                 
# [63] Gehyra paranana                  Gehyra pluraporosa              
# [65] Gehyra pseudopunctata            Gelanesaurus flavogularis       
# [67] Gerrhonotus lazcanoi             Glaucomastix cyanura            
# [69] Goggia matzikamaensis            Grandidierina lineata           
# [71] Grandidierina petiti             Grandidierina rubrocaudata      
# [73] Hemidactylus fragilis            Hemidactylus kangerensis        
# [75] Hemidactylus pauciporosus        Hemiphyllodactylus hongkongensis
# [77] Hemiphyllodactylus linnwayensis  Hemiphyllodactylus tonywhitteni 
# [79] Ichnotropis tanganicana          Kinyongia msuyae                
# [81] Kinyongia mulyai                 Larutia penangensis             
# [83] Lepidoblepharis emberawoundule   Lepidoblepharis rufigularis     
# [85] Lepidoblepharis victormartinezi  Leposternon bagual              
# [87] Leposternon cerradensis          Leposternon kisteumacheri       
# [89] Leposternon maximus              Leposternon mineiro             
# [91] Leposternon octostegum           Leposternon scutigerum          
# [93] Lerista hobsoni                  Lerista vanderduysi             
# [95] Lipinia inconspicua              Loxopholis caparensis           
# [97] Loxopholis ioanna                Lygosoma samajaya               
# [99] Lygosoma tabonorum               Mesaspis cuchumatanus           
# [101] Nannoscincus fuscus              Nannoscincus koniambo           
# [103] Nessia gansi                     Oedura bella                    
# [105] Oedura cincta                    Oedura fimbria                  
# [107] Ophiodes enso                    Ophiodes luciae                 
# [109] Ophiomorus kardesi               Oreosaurus luctuosus            
# [111] Oreosaurus serranus              Pachydactylus macrolepis        
# [113] Panaspis duruarum                Panaspis seydeli                
# [115] Papuascincus buergersi           Paracontias ampijoroensis       
# [117] Paracontias mahamavo             Petracola angustisoma           
# [119] Phaeoscincus ouinensis           Phaeoscincus taomensis          
# [121] Pholidobolus dicrus              Pholidoscelis turukaeraensis    
# [123] Pholidoscelis umbratilis         Plestiodon lotus                
# [125] Proctoporus rahmi                Proctoporus spinalis            
# [127] Pseudogekko isapa                Ptyodactylus togensis           
# [129] Scolecoseps broadleyi            Sigaloseps conditus             
# [131] Sigaloseps ferrugicauda          Sigaloseps pisinnus             
# [133] Sphenomorphus malaisei           Sphenomorphus rufus             
# [135] Trachylepis gonwouoi             Tribolonotus choiseulensis      
# [137] Tribolonotus parkeri             Tropidophorus sebi              
# [139] Tupinambis zuliensis             Tytthoscincus langkawiensis     
# [141] Tytthoscincus martae             Tytthoscincus textus            
# [143] Voeltzkowia yamagishii           Xenosaurus fractus              
# [145] Zygaspis maraisi                 Amphiglossus stylus             
# [147] Atractus ayeush                  Coniophanes sarae               
# [149] Dipsas petersi                   Epictia rubrolineata            
# [151] Hebius nicobariense              Leiopython montanus             
# [153] Liotyphlops caissara             Lycodon bibonius                
# [155] Oligodon wagneri                 Phisalixella iarakaensis        
# [157] Synophis calamitus               Tantilla sertula                
# [159] Typhlops cariei                  Sphenomorphus arborens          
# [161] Sphenomorphus cherriei   

## examples of species with 0% completeness
r0 <- Completeness0(ms_Reptiles, TraitsReptiles)
R0 <- ms_Reptiles %>% filter(Best_guess_binomial %in% r0)


## write files
write.csv(Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Mammals.csv", row.names = FALSE)
write.csv(Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Amphibians.csv", row.names = FALSE)
write.csv(Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Reptiles.csv", row.names = FALSE)
write.csv(Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Birds.csv", row.names=FALSE)

write.csv(U_Mammals, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Mammals.csv", row.names = FALSE)
write.csv(U_Amphibians, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Amphibians.csv", row.names = FALSE)
write.csv(U_Reptiles, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Reptiles.csv", row.names = FALSE)
write.csv(U_Birds, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Birds.csv", row.names=FALSE)

write.csv(Ums_Mammals, "../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Mammals_final.csv", row.names = FALSE)
write.csv(Ums_Amphibians, "../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Amphibians_final.csv", row.names = FALSE)
write.csv(Ums_Reptiles, "../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Reptiles_final.csv", row.names = FALSE)
write.csv(Ums_Birds, "../../Results/Trait_data_for_ms_GLOBALGAPS/UN_Birds_final.csv", row.names=FALSE)


write.csv(ms_Mammals, "../../Results/Trait_data_for_ms_GLOBALGAPS/Mammals_final.csv", row.names = FALSE)
write.csv(ms_Amphibians, "../../Results/Trait_data_for_ms_GLOBALGAPS/Amphibians_final.csv", row.names = FALSE)
write.csv(ms_Reptiles, "../../Results/Trait_data_for_ms_GLOBALGAPS/Reptiles_final.csv", row.names = FALSE)
write.csv(ms_Birds, "../../Results/Trait_data_for_ms_GLOBALGAPS/Birds_final.csv", row.names=FALSE)


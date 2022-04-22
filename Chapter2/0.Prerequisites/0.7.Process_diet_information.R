## Process Diet information from Kissling (Mammals) ("MammalDiet"), Amphibio (Amphibians) and Elton traits (Birds and mammals)
## Sekercioglu not used 

## Adopting a binary classification for food items. The variables are: VE, IN, FR, SE, PL, NE, SCV
## Define 5 categories:
# 1: plant/seed
# 2: fruit/nectar
# 3: vertebrates (includes carrion)
# 4: invertebrates
# 5: omnivores


## Preamble

X <- c("data.table", "plyr", "dplyr", "tidyr", "magrittr", "reshape", "reshape2", "stringr", "stringi") 
invisible(lapply(X, library, character.only=TRUE)); rm(X)
`%nin%` = Negate(`%in%`)

# Corrected data
MammalDIET_C <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/MammalDIET.csv")
Amphibio_C <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Amphibio.csv")
Elton_birds_C <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Elton_birds.csv")
Elton_mammals_C <-  read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Elton_mammals.csv")

# Uncorrected data
Amphibio_UN <- read.csv("../../Data/Amphibians/AmphiBIO_v1.csv")
colnames(Amphibio_UN)[5] <- "Best_guess_binomial"
Elton_birds_UN <- read.csv("../../Data/Birds/EltonTraits_Birds.csv")
Elton_mammals_UN <-  read.csv("../../Data/Mammals/EltonTraits_Mammals.csv")

## Change Diet.5Cat in Elton birds
Elton_birds_C <- Elton_birds_C %>%
  mutate(Diet.5Cat=ifelse(Diet.5Cat=="PlantSeed", "PL|SE", 
                          ifelse(Diet.5Cat=="Omnivore", "OM", 
                                 ifelse(Diet.5Cat=="FruiNect", "FR|NE", 
                                        ifelse(Diet.5Cat=="Invertebrate", "IN",
                                               ifelse(Diet.5Cat=="VertFishScav", "VE", NA))))))

Elton_birds_UN <- Elton_birds_UN %>%
  mutate(Diet.5Cat=ifelse(Diet.5Cat=="PlantSeed", "PL|SE", 
                          ifelse(Diet.5Cat=="Omnivore", "OM", 
                                 ifelse(Diet.5Cat=="FruiNect", "FR|NE", 
                                        ifelse(Diet.5Cat=="Invertebrate", "IN",
                                               ifelse(Diet.5Cat=="VertFishScav", "VE", NA))))))


# Put mammalian diet information together
Kissling <- read.csv("../../Data/Mammals/Kissling_Mammal_diet_2014.csv", sep=",")
Kissling$Binomial_name <- paste(Kissling$Genus, Kissling$Species, sep=" ")
Kissling <- subset(Kissling, Binomial_name != "Mico sp. nov.")
Kissling <- Kissling %>% select(-TaxonID, -TaxonomicNote, -FillCode)
Kissling <- Kissling[, c(28, 1:27)]

MammalDIET2 <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET2.csv")[, c(1:28)]
colnames(MammalDIET2)[28] <- "DataSource"

MammalDIET2_supp <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET_2_supp.csv")
MammalDIET2_supp$TrophicLevel %<>% as.character()
MammalDIET2_supp <- MammalDIET2_supp %>%
  mutate(TrophicLevel=ifelse(TrophicLevel %in% c("Carnivore", "Omnivore", "Herbivore"), TrophicLevel, NA)) %>%
  filter(!is.na(TrophicLevel)) %>%
  mutate(Mammal=ifelse(Mammal %in% c(0:3), Mammal, NA)) %>%
  mutate(MammalEater=ifelse(MammalEater %in% c(0:3), MammalEater, NA)) %>%
  mutate(Insectivore=ifelse(Insectivore %in% c(0:3), Insectivore, NA))%>%
  mutate(Frugivore=ifelse(Frugivore %in% c(0:3), Frugivore, NA))%>%
  mutate(Granivore=ifelse(Granivore %in% c(0:3), Granivore, NA)) %>%
  mutate(Folivore=ifelse(Folivore %in% c(0:3), Folivore, NA))
MammalDIET2_supp[, c(6:27)] %<>% droplevels()

MammalDIET2 <- rbind(MammalDIET2, MammalDIET2_supp)
MammalDIET2$Binomial <- as.character(MammalDIET2$Binomial)
MammalDIET2$Binomial[MammalDIET2$Binomial=="Rhinolophus hildebrandtii"] <- "Rhinolophus hildebrandti"
MammalDIET2$Binomial[MammalDIET2$Binomial=="Tadarida bivittatus"] <- "Tadarida bivittata"
colnames(MammalDIET2)[1] <- "Binomial_name"

Y <- intersect(MammalDIET2$Binomial_name, Kissling$Binomial_name)
Kissling <- Kissling %>% filter(Binomial_name %nin% Y)
MammalDIET_UN <- rbind(Kissling, MammalDIET2) # dataset to use

rm(Kissling, MammalDIET2, MammalDIET2_supp, Y)

write.csv(Elton_mammals_C, "../../../3.Explanatory_traits/Data/EltonMammalsDiet_corrected.csv", row.names = FALSE)
write.csv(Elton_birds_C, "../../../3.Explanatory_traits/Data/EltonBirdsDiet_corrected.csv", row.names = FALSE)

# # Functions

## 1. AmphiBIO: changing the column names to PL, NE, SE, FR, IN, VE

For.Amphibio <- function(AmphiData) {
  
  Func <- function(X) {ifelse(is.na(X), 0, 1)}
  
  # Food items as binary variables (0 or 1)
  colnames(AmphiData)[c(10:15)] <- c("PL", "NE", "SE", "FR", "IN", "VE")
  AmphiData <- AmphiData %>%
    mutate_at(.vars = vars(PL, NE, SE, FR, IN, VE),
              .funs = Func)
  
  # Diet breadth
  AmphiData$Diet_breadth <- apply(AmphiData[, c("PL", "NE", "SE", "FR", "IN", "VE")], 1, sum, na.rm=T)
  AmphiData$Diet_breadth[AmphiData$Diet_breadth==0] <- NA
  
  # Trophic levels
  AmphiData$Trophic_level <- NA
  AmphiData <- AmphiData %>%
    mutate(Trophic_level=ifelse( (PL==1|NE==1|SE==1|FR==1) & (IN==0 & VE==0), "Herbivore",
                                 ifelse( (PL==0 & NE==0 & SE==0 & FR==0) & (IN==1|VE==1), "Carnivore",
                                         ifelse( (PL==1|NE==1|SE==1|FR==1) & (IN==1| VE==1), "Omnivore", 
                                                 ifelse( (PL==0 & NE==0 & SE==0 & FR==0 & IN==0 & VE==0), NA, Trophic_level)))))
  
  # Primary diet
  AmphiData$Primary_diet <- NA
  AmphiData <- AmphiData %>%
    mutate(Primary_diet=ifelse((PL==1|SE==1) & (NE==0 & FR==0 & IN==0 & VE==0), "PL|SE", 
                               ifelse((FR==1|NE==1) & (PL==0 & SE==0 & IN==0 & VE==0), "FR|NE", 
                                      ifelse(VE==1 & (PL==0 & SE==0 & NE==0 & FR==0 & IN==0), "VE", 
                                              ifelse(IN==1 & (PL==0 & SE==0 & NE==0 & FR==0 & VE==0), "IN", 
                                                     ifelse(IN==0 & PL==0 & SE==0 & NE==0 & FR==0 & VE==0, NA, "OM"))))))
  
  return(AmphiData)
  
  
}


## 2. Kissling (MammalDIET): 1 are primary food items, 0 are items that score less than 50%

For.MammalDiet <- function(MammalDiet) {
  
  # Correct
  MammalDiet[MammalDiet=="Could not find data"] <- NA
  MammalDiet$TrophicLevel <- as.character(MammalDiet$TrophicLevel)
  MammalDiet <- MammalDiet %>%
    mutate(TrophicLevel=ifelse(TrophicLevel=="NotAssigned", NA, TrophicLevel))
  
  # Create diet binary variables
  MammalDiet <- MammalDiet %>% 
    mutate(VE=ifelse(Vertebrate==1|Mammal==1|Bird==1|Herptile==1|Fish==1, 1, 0)) %>% 
    mutate(IN=ifelse(Invertebrate==1, 1, 0)) %>% 
    mutate(SE=ifelse(Seed==1, 1, 0)) %>% 
    mutate(FR=ifelse(Fruit==1, 1, 0)) %>% 
    mutate(NE=ifelse(Nectar==1, 1, 0)) %>% 
    mutate(PL=ifelse((Plant==1 & Seed==0 & Fruit==0 & Nectar==0)|Root==1|Leaf==1|Woody==1|Herbaceous==1, 1, 0))
  
  # Diet breadth
  MammalDiet$Diet_breadth <- apply(MammalDiet[, c("PL", "NE", "SE", "FR", "IN", "VE")], 1, sum, na.rm=T)
  MammalDiet$Diet_breadth[MammalDiet$Diet_breadth==0] <- NA
  
  # Primary diet
  MammalDiet$Primary_diet <- NA
  MammalDiet <- MammalDiet %>%
    mutate(Primary_diet=ifelse((PL==1|SE==1) & (NE==0 & FR==0 & IN==0 & VE==0), "PL|SE", 
                               ifelse((FR==1|NE==1) & (PL==0 & SE==0 & IN==0 & VE==0), "FR|NE", 
                                      ifelse(VE==1 & (PL==0 & SE==0 & NE==0 & FR==0 & IN==0), "VE", 
                                             ifelse(IN==1 & (PL==0 & SE==0 & NE==0 & FR==0 & VE==0), "IN", 
                                                    ifelse(IN==0 & PL==0 & SE==0 & NE==0 & FR==0 & VE==0, NA, "OM"))))))
  
  
  return(MammalDiet)
  
}

## 3. Elton birds and mammals

# For Elton mammals, need to define the primary diet based on food items that score more than 50 percent
# Primary diet is already present in Elton birds

Elton.Mammals <- function(Elton) {
  
  # Vertebrate items as one
  Elton$Diet.Vert <- Elton$Diet.Vend + Elton$Diet.Vect + Elton$Diet.Vfish + Elton$Diet.Vunk + Elton$Diet.Scav
  
  # Primary diet based on predominant food items (here, cases where all < 50 or one >50)
  Elton <- Elton %>% 
    mutate(Diet.5Cat=ifelse((Diet.Fruit<50 & Diet.Inv<50 & Diet.Nect<50 & Diet.Seed<50 & Diet.PlantO<50 & Diet.Vert <50), "OM",  # cases where all below 50
                            # cases where one food item is strictly above 50 percent    
                            ifelse(Diet.Fruit>50 | Diet.Nect>50, "FR|NE",
                                             ifelse(Diet.PlantO>50 | Diet.Seed>50, "PL|SE",
                                                    ifelse(Diet.Vert>50, "VE",
                                                           ifelse(Diet.Inv>50, "IN", 
                                                                  # cases where one above 50 
                                                                  ifelse(Diet.Fruit==50 & Diet.Nect<50 & Diet.Inv <50 & Diet.PlantO <50 & Diet.Seed<50 & Diet.Vert<50, "FR|NE", 
                                                                         ifelse(Diet.Fruit<50 & Diet.Nect==50 & Diet.Inv <50 & Diet.PlantO <50 & Diet.Seed<50 & Diet.Vert<50, "FR|NE",
                                                                                ifelse(Diet.Fruit<50 & Diet.Nect<50 & Diet.Inv==50 & Diet.PlantO <50 & Diet.Seed<50 & Diet.Vert<50, "IN",
                                                                                       ifelse(Diet.Fruit<50 & Diet.Nect<50 & Diet.Inv <50 & Diet.PlantO==50 & Diet.Seed<50 & Diet.Vert<50, "PL|SE",
                                                                                              ifelse(Diet.Fruit<50 & Diet.Nect<50 & Diet.Inv <50 & Diet.PlantO <50 & Diet.Seed==50 & Diet.Vert<50, "PL|SE",
                                                                                                     ifelse(Diet.Fruit<50 & Diet.Nect<50 & Diet.Inv <50 & Diet.PlantO <50 & Diet.Seed<50 & Diet.Vert==50, "VE",
                                                                                                            #cases where two food items have equal percent use
                                                                                                            ifelse(Diet.Fruit==50 & Diet.Nect==50, "FR|NE", 
                                                                                                                   ifelse(Diet.PlantO==50 & Diet.Seed==50, "PL|SE", "OM"))))))))))))))
  return(Elton)
}

For.Elton <- function(Elton) {
  
  # Change column name of primary diet
  colnames(Elton)[colnames(Elton)=="Diet.5Cat"] <- "Primary_diet"
  
  # Food items as binary variables
  Elton <- Elton %>%
    mutate(VE=ifelse(Diet.Vunk>=50|Diet.Vect>=50|Diet.Vend>=50|Diet.Vfish>=50|Diet.Scav>=50, 1, 0)) %>%
    mutate(IN=ifelse(Diet.Inv>=50, 1, 0)) %>%
    mutate(FR=ifelse(Diet.Fruit>=50, 1, 0)) %>%
    mutate(NE=ifelse(Diet.Nect>=50, 1, 0)) %>%
    mutate(SE=ifelse(Diet.Seed>=50,1, 0)) %>%
    mutate(PL=ifelse(Diet.PlantO>=50, 1, 0))
  
  # Trophic levels
  Elton$Trophic_level <- NA
  Elton <- Elton %>%
    mutate(Trophic_level=ifelse( (PL==1|NE==1|SE==1|FR==1) & (IN==0 & VE==0), "Herbivore",
                                 ifelse( (PL==0 & NE==0 & SE==0 & FR==0) & (IN==1|VE==1), "Carnivore",
                                         ifelse( (PL==1|NE==1|SE==1|FR==1) & (IN==1| VE==1), "Omnivore", 
                                                 ifelse( (PL==0 & NE==0 & SE==0 & FR==0 & IN==0 & VE==0), NA, Trophic_level)))))
  
  # Diet breadth
  Elton$Diet_breadth <- apply(Elton[, c("IN", "VE", "SE", "PL", "FR", "NE")], 1, sum, na.rm=T)
  Elton$Diet_breadth[Elton$Diet_breadth==0] <- NA
  
  # Remove some columns
  Elton <- Elton %>%
    select(-Diet.Vect, -Diet.Vend, -Diet.Vunk, -Diet.Vfish, -Diet.Scav, -Diet.Inv, -Diet.Fruit, -Diet.Nect, -Diet.Seed, -Diet.PlantO) 
  
  return(Elton)
}


# # Runs

Amphibio_C <- For.Amphibio(Amphibio_C)
Amphibio_UN <- For.Amphibio(Amphibio_UN)

MammalDIET_C <- For.MammalDiet(MammalDIET_C)
MammalDIET_UN <- For.MammalDiet(MammalDIET_UN)

Elton_mammals_C <- Elton.Mammals(Elton_mammals_C)
Elton_mammals_UN <- Elton.Mammals(Elton_mammals_UN)

Elton_mammals_C <- For.Elton(Elton_mammals_C)
Elton_mammals_UN <- For.Elton(Elton_mammals_UN)

Elton_birds_C <- For.Elton(Elton_birds_C)
Elton_birds_UN <- For.Elton(Elton_birds_UN)



# Save files --------------------------------------------------------------
write.csv(Amphibio_C, "../../Results/0.Processed_diet_datasets/Amphibio_processed_diet.csv", row.names = F)
write.csv(MammalDIET_C, "../../Results/0.Processed_diet_datasets/MammalDIET_processed_diet.csv", row.names = F)
write.csv(Elton_birds_C, "../../Results/0.Processed_diet_datasets/Elton_birds_processed.csv", row.names = F)
write.csv(Elton_mammals_C, "../../Results/0.Processed_diet_datasets/Elton_mammals_processed.csv", row.names = F)

write.csv(Amphibio_UN, "../../Data/Amphibians/Amphibio_processed_diet.csv", row.names = F)
write.csv(MammalDIET_UN, "../../Data/Mammals/MammalDIET_processed_diet.csv", row.names = F)
write.csv(Elton_birds_UN, "../../Data/Birds/Elton_birds_processed.csv", row.names = F)
write.csv(Elton_mammals_UN, "../../Data/Mammals/Elton_mammals_processed.csv", row.names = F)




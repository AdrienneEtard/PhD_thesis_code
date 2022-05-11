## SCRIPT AIM:
## Add diet data compiled by Rhiannon Osborne to the trait datasets for amphibians and reptiles
## And for mammals and birds, run some additional checks as data has already been compiled
## then, define Primary diet for use in further models

## HERE ALTERNATIVE SCRIPT WITH DIFFERENT GROUPINGS FOR PRIMARY DIET

## some clarifications about mammal and bird diet:
# for mammals, the trophic levels are obtained from the Kissling dataset (MammalDiet)
# From Elton traits for mammals we get: PD, DB, Diet as binary (FR, NE, SE, PL, IN, VE)
# for birds, we get Diet, PD, DB from Elton Traits; and we define the trophic levels on the basis of the Diet.

# for Elton, a 0 indicates that the food item makes up less than 50% of the diet because the definitions are based on primarily consumed items

library(dplyr)
`%nin%` <- Negate(`%in%`)

## AMPHIBIANS #####

AmphDiet <- read.csv("../../Data/Amphibian_Diet.csv") ## Rhiannon's compilation
AmphDiet[!is.na(AmphDiet$Primary_diet),] %>%  nrow
AmphDiet$Primary_diet[!is.na(AmphDiet$Primary_diet)] %>%  table
AmphDiet$Paper[AmphDiet$Paper!=""] %>% unique %>% length() ## 26 papers
levels(AmphDiet$Paper)

# cleaning up amphibian diet data
unique(AmphDiet$Primary_diet)
unique(AmphDiet$Trophic_level)
AmphDiet$Primary_diet[AmphDiet$Primary_diet==""] <- NA
AmphDiet$Trophic_level[AmphDiet$Trophic_level==""] <- NA
colnames(AmphDiet)[1] <- "Best_guess_binomial"

# Amphibian trait data
Amphibians <- read.csv("../../Data/GlobalGaps_traitdata/Amphibians.csv")

# how much information from Amphibio?
Amphibians$DB <- apply(Amphibians[, c("IN", "FR", "NE", "VE", "PL", "SE")], 1, sum)
nrow(Amphibians[Amphibians$DB!=0,])

# all names intersect between collected data and Amphibio
intersect(AmphDiet$Best_guess_binomial, Amphibians$Best_guess_binomial) %>% length()

# for these species in PREDICTS
To.complete <- subset(Amphibians, Best_guess_binomial %in% intersect(AmphDiet$Best_guess_binomial, Amphibians$Best_guess_binomial))
To.complete <- subset(To.complete, is.na(Trophic_level))
Species.to.complete <- To.complete$Best_guess_binomial # 244

# retrieve from AmphDiet the species for which information is needed and clean up data
Data.to.complete <- AmphDiet[AmphDiet$Best_guess_binomial %in% Species.to.complete,]
unique(Data.to.complete$Primary_diet)
unique(Data.to.complete$Trophic_level)
unique(Data.to.complete$Invertebrate)
unique(Data.to.complete$Vertebrate)
unique(Data.to.complete$Invertebrate)
unique(Data.to.complete$Invertebrate)
unique(Data.to.complete$Invertebrate)

Data.to.complete$Trophic_level[Data.to.complete$Primary_diet=="IN"] <- "Carnivore"

Data.to.complete <- subset(Data.to.complete, !is.na(Primary_diet))
Data.to.complete <- Data.to.complete %>% select(Best_guess_binomial, Primary_diet, Invertebrate, Vertebrate, Fruit, Nectar, Seed, Plant, Trophic_level, Paper)
colnames(Data.to.complete) <- c("Best_guess_binomial", "Primary_diet", "IN", "VE", "FR", "NE", "SE", "PL", "Trophic_level", "Source")

## how many species do we add from the literature?
Data.to.complete$Best_guess_binomial %>% unique  %>% length()
Data.to.complete$Source  %>% unique()
Data.to.complete$Primary_diet  %>% unique()

## adding diet data for 108 PREDICTS amphibians in the global trait dataset
Amphibians$Best_guess_binomial <- as.character(Amphibians$Best_guess_binomial)
Data.to.complete$Best_guess_binomial <- as.character(Data.to.complete$Best_guess_binomial)

Amphibians$Trophic_level[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$Trophic_level
Amphibians$Primary_diet[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$Primary_diet
Amphibians$IN[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$IN
Amphibians$VE[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$VE
Amphibians$FR[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$FR
Amphibians$NE[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$NE
Amphibians$SE[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$SE
Amphibians$PL[Amphibians$Best_guess_binomial %in% Data.to.complete$Best_guess_binomial] <- Data.to.complete$PL

## cleaning up diet data
Amphibians[is.na(Amphibians$Primary_diet), c( "IN", "VE", "FR", "NE", "SE", "PL")] <- NA
Amphibians$Diet_breadth <- apply(Amphibians[,c( "IN", "VE", "FR", "NE", "SE", "PL")], 1, sum)
Amphibians$IN[Amphibians$Diet_breadth==0] <- NA
Amphibians$VE[Amphibians$Diet_breadth==0] <- NA
Amphibians$FR[Amphibians$Diet_breadth==0] <- NA
Amphibians$NE[Amphibians$Diet_breadth==0] <- NA
Amphibians$SE[Amphibians$Diet_breadth==0] <- NA
Amphibians$PL[Amphibians$Diet_breadth==0] <- NA
Amphibians$Diet_breadth <- apply(Amphibians[,c( "IN", "VE", "FR", "NE", "SE", "PL")], 1, sum)
unique(Amphibians$Diet_breadth)

## redefining Primary diet
Amphibians <- Amphibians %>% 
  mutate(Primary_diet=ifelse(VE==1 & SE==0 & PL==0 & IN==0 & FR==0 & NE==0, "VE", 
                             ifelse(IN==1 & SE==0 & PL==0 & VE==0 & FR==0 & NE==0, "IN", 
                                    ifelse(IN==0 & SE==0 & PL==0 & VE==0 & FR==1 & NE==0, "FR", 
                                           ifelse(IN==0 & SE==0 & PL==0 & VE==0 & FR==0 & NE==1, "NE",
                                                  ifelse(IN==0 & SE==0 & PL==0 & VE==0 & FR==1 & NE==1, "FRNE", 
                                                         ifelse(IN==0 & SE==0 & PL==1 & VE==0 & FR==0 & NE==0, "PL",
                                                                ifelse(IN==0 & SE==1 & PL==0 & VE==0 & FR==0 & NE==0, "SE",
                                                                       ifelse(IN==0 & SE==1 & PL==1 & VE==0 & FR==0 & NE==0, "PLSE",
                                                                              ifelse(is.na(IN) & is.na(VE) & is.na(SE) & is.na(PL) & is.na(FR) & is.na(NE), NA, "OM"))))))))))




unique(Amphibians$Primary_diet)
write.csv(Amphibians, "../../Data/TraitsWithDiet/Diet_redefined/Amphibians.csv", row.names = FALSE)

EltAm <- read.csv("../../Data/TraitsWithDiet/Amphibians.csv")
table(EltAm$Primary_diet)
table(Amphibians$Primary_diet) 


## REPTILES #####

ReptDiet <- read.csv("../../Data/Reptile_Diet.csv")
colnames(ReptDiet)[1] <- "Best_guess_binomial"

Papers <- unique(ReptDiet$Paper)
Papers[order(Papers)]

# reptile trait data
Reptiles <- read.csv("../../Data/GlobalGaps_traitdata/Reptiles.csv")

# all names intersect
intersect(ReptDiet$Best_guess_binomial, Reptiles$Best_guess_binomial) %>% length()

# adding diet for those 329 PREDICTS species -- we delete the Funghi columns - this only affects one species and does not change the primary diet for that species
ReptDiet <- ReptDiet %>% 
  dplyr::select(Best_guess_binomial, Seed, Fruit, Plant, Nectar, Invertebrate, Vertebrate_gen) %>% 
  setNames(., c("Best_guess_binomial", "SE", "FR", "PL", "NE", "IN", "VE"))

ReptDiet$Best_guess_binomial <- as.character(ReptDiet$Best_guess_binomial)
ReptDiet$Best_guess_binomial %>%  length()

Reptiles$Best_guess_binomial <- as.character(Reptiles$Best_guess_binomial)
Reptiles <- left_join(Reptiles, ReptDiet, by="Best_guess_binomial")

## verifications
apply(ReptDiet[,-1], 1, sum) %>%  unique()
Reptiles$Diet_breadth <- apply(Reptiles[,c( "IN", "VE", "FR", "NE", "SE", "PL")], 1, sum)


## redefining Primary diet
Reptiles <- Reptiles %>% 
  mutate(Primary_diet=ifelse(VE==1 & SE==0 & PL==0 & IN==0 & FR==0 & NE==0, "VE", 
                             ifelse(IN==1 & SE==0 & PL==0 & VE==0 & FR==0 & NE==0, "IN", 
                                    ifelse(IN==0 & SE==0 & PL==0 & VE==0 & FR==1 & NE==0, "FR", 
                                           ifelse(IN==0 & SE==0 & PL==0 & VE==0 & FR==0 & NE==1, "NE",
                                                  ifelse(IN==0 & SE==0 & PL==0 & VE==0 & FR==1 & NE==1, "FRNE", 
                                                         ifelse(IN==0 & SE==0 & PL==1 & VE==0 & FR==0 & NE==0, "PL",
                                                                ifelse(IN==0 & SE==1 & PL==0 & VE==0 & FR==0 & NE==0, "SE",
                                                                       ifelse(IN==0 & SE==1 & PL==1 & VE==0 & FR==0 & NE==0, "PLSE",
                                                                              ifelse(is.na(IN) & is.na(VE) & is.na(SE) & is.na(PL) & is.na(FR) & is.na(NE), NA, "OM"))))))))))


write.csv(Reptiles, "../../Data/TraitsWithDiet/Diet_redefined/Reptiles.csv", row.names = FALSE)

EltRep <- read.csv("../../Data/TraitsWithDiet/Reptiles.csv")
table(EltRep$Primary_diet)
table(Reptiles$Primary_diet)



## COMPILING DIET FOR MAMMALS AND BIRDS #####

## sources for diet (NB: for trophic levels, mammals: using Kissling's MammalDiet)
EltonMammals <- read.csv("../../Data/EltonMammalsDiet_corrected.csv")
EltonBirds <- read.csv("../../Data/EltonBirdsDiet_corrected.csv") 

EltonBirds <- EltonBirds %>% 
  dplyr::select(Best_guess_binomial, Scientific, Diet.Inv, Diet.Vend, Diet.Vect, Diet.Vfish, Diet.Vunk, Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO) %>% 
  mutate(Diet.Vert= Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk + Diet.Scav) %>% 
  mutate(DietPL=Diet.PlantO, DietSE=Diet.Seed, DietFR=Diet.Fruit, DietNE=Diet.Nect) %>% 
  mutate(DietPLSE=Diet.PlantO+Diet.Seed, DietFRNE=Diet.Fruit+Diet.Nect) %>% 
  
  mutate(Diet.5Cat=ifelse(Diet.Inv<=50 & Diet.Vert<=50 & DietPLSE <=50 & DietFRNE <=50, "OM", 
                          ifelse(Diet.Inv>50, "IN", 
                                 ifelse(Diet.Vert>50, "VE",
                                        ifelse(DietFR>50, "FR", 
                                               ifelse(DietNE>50, "NE", 
                                                      ifelse(DietPL>50, "PL", 
                                                             ifelse(DietSE>50, "SE",
                                                                    ifelse(DietFRNE>50 & DietFR<=50 & DietNE<=50, "FR|NE", 
                                                                           ifelse(DietPLSE>50  & DietPL<=50 & DietSE<=50, "PL|SE","OM")))))))))) %>% 
  mutate(VE=ifelse(Diet.Vert!=0, 1, 0)) %>%
  mutate(IN=ifelse(Diet.Inv!=0, 1, 0)) %>%
  mutate(FR=ifelse(Diet.Fruit!=0, 1, 0)) %>%
  mutate(NE=ifelse(Diet.Nect!=0, 1, 0)) %>%
  mutate(SE=ifelse(Diet.Seed!=0,1, 0)) %>%
  mutate(PL=ifelse(Diet.PlantO!=0, 1, 0)) %>% 
  mutate(DB=IN+VE+FR+NE+SE+PL) %>% 
  mutate(Trophic_level.Elton=ifelse(Diet.5Cat %in% c("FR", "NE", "FR|NE", "PL", "SE", "PL|SE"), "Herbivore",
                                    ifelse(Diet.5Cat %in% c("IN", "VE"), "Carnivore",
                                           ifelse(DB==3 & ((NE==1&PL==1&FR==1)|(NE==1&PL==1&SE==1)|(NE==1&FR==1&SE==1)|(FR==1&PL==1&SE==1)), "Herbivore",
                                                  ifelse(DB==2 & ((NE==1&PL==1)|(NE==1&FR==1)|(NE==1&SE==1)|(PL==1&FR==1)|(PL==1&SE==1)|(FR==1&SE==1)), "Herbivore",
                                                         "Omnivore"))))) %>% 
  dplyr::select(VE, IN, SE, NE, PL, FR, DB, Diet.5Cat, Best_guess_binomial, Scientific, Trophic_level.Elton)

colnames(EltonBirds)[c(7,8)] <- c("Diet_breadth", "Primary_diet")

## processing diet for mammals
EltonMammals <- EltonMammals %>% 
  dplyr::select(Best_guess_binomial, Scientific, Diet.Inv, Diet.Vend, Diet.Vect, Diet.Vfish, Diet.Vunk, Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO) %>% 
  mutate(Diet.Vert= Diet.Vend + Diet.Vect + Diet.Vfish + Diet.Vunk + Diet.Scav) %>% 
  mutate(DietPL=Diet.PlantO, DietSE=Diet.Seed, DietFR=Diet.Fruit, DietNE=Diet.Nect) %>% 
  mutate(DietPLSE=Diet.PlantO+Diet.Seed, DietFRNE=Diet.Fruit+Diet.Nect) %>% 
  
  mutate(Diet.5Cat=ifelse(Diet.Inv<=50 & Diet.Vert<=50 & DietPLSE <=50 & DietFRNE <=50, "OM", 
                          ifelse(Diet.Inv>50, "IN", 
                                 ifelse(Diet.Vert>50, "VE",
                                        ifelse(DietFR>50, "FR", 
                                               ifelse(DietNE>50, "NE", 
                                                      ifelse(DietPL>50, "PL", 
                                                             ifelse(DietSE>50, "SE",
                                                                    ifelse(DietFRNE>50 & DietFR<=50 & DietNE<=50, "FR|NE", 
                                                                           ifelse(DietPLSE>50  & DietPL<=50 & DietSE<=50, "PL|SE","OM")))))))))) %>% 
  mutate(VE=ifelse(Diet.Vert!=0, 1, 0)) %>%
  mutate(IN=ifelse(Diet.Inv!=0, 1, 0)) %>%
  mutate(FR=ifelse(Diet.Fruit!=0, 1, 0)) %>%
  mutate(NE=ifelse(Diet.Nect!=0, 1, 0)) %>%
  mutate(SE=ifelse(Diet.Seed!=0,1, 0)) %>%
  mutate(PL=ifelse(Diet.PlantO!=0, 1, 0)) %>% 
  mutate(DB=IN+VE+FR+NE+SE+PL) %>% 
  mutate(Trophic_level.Elton=ifelse(Diet.5Cat %in% c("FR", "NE", "FR|NE", "PL", "SE", "PL|SE"), "Herbivore",
                                    ifelse(Diet.5Cat %in% c("IN", "VE"), "Carnivore",
                                           ifelse(DB==3 & ((NE==1&PL==1&FR==1)|(NE==1&PL==1&SE==1)|(NE==1&FR==1&SE==1)|(FR==1&PL==1&SE==1)), "Herbivore",
                                                  ifelse(DB==2 & ((NE==1&PL==1)|(NE==1&FR==1)|(NE==1&SE==1)|(PL==1&FR==1)|(PL==1&SE==1)|(FR==1&SE==1)), "Herbivore",
                                                         "Omnivore"))))) %>% 
  dplyr::select(VE, IN, SE, NE, PL, FR, DB, Diet.5Cat, Best_guess_binomial, Scientific, Trophic_level.Elton)

colnames(EltonMammals)[c(7,8)] <- c("Diet_breadth", "Primary_diet")

apply(EltonMammals[,c( "IN", "VE", "FR", "NE", "SE", "PL")], 1, sum) %>%  unique
EltonMammals[EltonMammals$Diet_breadth==0,c( "IN", "VE", "FR", "NE", "SE", "PL")] <- NA
EltonMammals$Diet_breadth <- apply(EltonMammals[,c( "IN", "VE", "FR", "NE", "SE", "PL")], 1, sum) 
unique(EltonMammals$Diet_breadth)
EltonMammals$Trophic_level.Elton[is.na(EltonMammals$Diet_breadth)] <- NA

apply(EltonBirds[,c( "IN", "VE", "FR", "NE", "SE", "PL")], 1, sum) %>%  unique

nrow(distinct(EltonMammals[,-10]))
nrow(distinct(EltonBirds[,-10]))

length(unique(EltonMammals$Best_guess_binomial))
length(unique(EltonBirds$Best_guess_binomial))

EltonMammals <- distinct(EltonMammals[, -10])
EltonBirds <- distinct(EltonBirds[, -10])

# species for which there is still a duplicate for mammals
DM <- EltonMammals %>%  
  group_by(Best_guess_binomial) %>% 
  summarise(C=n()) %>% 
  filter(C>1) %>% 
  dplyr::select(Best_guess_binomial)
DM <- DM$Best_guess_binomial %>%  as.character()

DM <- EltonMammals %>% filter(Best_guess_binomial %in% DM)
DM <- DM %>%  group_by(Best_guess_binomial) %>% 
  filter(Diet_breadth==max(Diet_breadth)) %>% 
  slice(1) # this is an approximation for a small number of species

# reassemble for mammals
EltonMammals <- rbind(EltonMammals[EltonMammals$Best_guess_binomial %nin% DM$Best_guess_binomial,], as.data.frame(DM))

# species for which there is still a duplicate for birds

DB <- EltonBirds %>%  
  group_by(Best_guess_binomial) %>% 
  summarise(C=n()) %>% 
  filter(C>1) %>% 
  dplyr::select(Best_guess_binomial)
DB <- DB$Best_guess_binomial %>%  as.character()

DB <- EltonBirds %>% filter(Best_guess_binomial %in% DB)
DB <- DB %>%  group_by(Best_guess_binomial) %>% 
  filter(Diet_breadth==max(Diet_breadth)) %>% 
  slice(1) # this is an approximation for a small number of species

# reassemble for birds
EltonBirds <- rbind(EltonBirds[EltonBirds$Best_guess_binomial %nin% DB$Best_guess_binomial,], as.data.frame(DB))



## Assembling datasets
Mammals <- read.csv("../../Data/GlobalGaps_traitdata/Mammals.csv") %>% 
  dplyr::select(-VE, -IN, -SE, -NE, -PL, -FR, -Diet_breadth, -Primary_diet)
Mammals$Best_guess_binomial <- as.character(Mammals$Best_guess_binomial)

M <- left_join(Mammals, EltonMammals, by="Best_guess_binomial")
table(M$Primary_diet)

# comparing alternative diet with Elton's classification for mammals
EltMa <- read.csv("../../Data/TraitsWithDiet/Mammals.csv")
table(EltMa$Primary_diet)

## the numbers are congruent

# comparing TL from Kissling and TL from Elton for mammals
MTL <- M %>%  dplyr::select(Trophic_level, Trophic_level.Elton)
nrow(MTL[MTL$Trophic_level==MTL$Trophic_level.Elton,])/nrow(Mammals)*100 # 70% congruent

M <- M %>% 
  dplyr::select(-Trophic_level.Elton)

write.csv(M, "../../Data/TraitsWithDiet/Diet_redefined/Mammals.csv", row.names = FALSE)

## Assembling for birds
Birds <- read.csv("../../Data/GlobalGaps_traitdata/Birds.csv") %>% 
  dplyr::select(-VE, -IN, -SE, -NE, -PL, -FR, -Diet_breadth, -Primary_diet, -Trophic_level)
Birds$Best_guess_binomial <- as.character(Birds$Best_guess_binomial)

Birds2 <- read.csv("../../Data/GlobalGaps_traitdata/Birds.csv") %>% 
  dplyr::select(-VE, -IN, -SE, -NE, -PL, -FR, -Diet_breadth, -Primary_diet)
Birds2$Best_guess_binomial <- as.character(Birds2$Best_guess_binomial)

B <- left_join(Birds, EltonBirds, by="Best_guess_binomial")
B2 <- left_join(Birds2, EltonBirds, by="Best_guess_binomial")

# Comp <- B2$Trophic_level==B2$Trophic_level.Elton 
# table(Comp)
# 369/length(Comp)*100
# comparing alternative diet with Elton's classification for birds

EltBi <- read.csv("../../Data/TraitsWithDiet/Birds.csv")

table(EltBi$Primary_diet)
table(B$Primary_diet)

## moslty congruent in terms of numbers with very few differences 

# om1 <- EltBi[EltBi$Primary_diet=="IN",]
# om2 <- B[B$Primary_diet=="IN",]
# s <- setdiff(om2$Best_guess_binomial, om1$Best_guess_binomial)
# om2 <- B %>%  filter(Best_guess_binomial %in% s)
# om1 <- om1 %>%  filter(Best_guess_binomial %in% s)
# sub <- EltBi %>% filter(Best_guess_binomial %in% s)
# EltonBirds<- read.csv("../../Data/EltonBirdsDiet_corrected.csv") 
# sube <- EltonBirds %>% filter(Best_guess_binomial %in% s)


write.csv(B, "../../Data/TraitsWithDiet/Diet_redefined/Birds.csv", row.names = FALSE)


#####################################################
## Analysing coverage for diet variables

rm(list=ls())
library(dplyr)
library(ggplot2)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


Mammals <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Mammals.csv") %>% 
  mutate(Class="Mammals")
Birds <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Birds.csv") %>% 
  mutate(Class="Birds")
colnames(Birds)[37] <- "Trophic_level"
Amphibians <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Amphibians.csv") %>% 
  mutate(Class="Amphibians")
Reptiles <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Reptiles.csv") %>% 
  mutate(Class="Reptiles")

Dietdata <- rbind(Mammals[, c("Class", "Primary_diet", "Trophic_level", "Best_guess_binomial")], 
                  Birds[, c("Class", "Primary_diet", "Trophic_level", "Best_guess_binomial")], 
                  Amphibians[, c("Class", "Primary_diet", "Trophic_level", "Best_guess_binomial")], 
                  Reptiles[, c("Class", "Primary_diet", "Trophic_level", "Best_guess_binomial")])


## Primary diet and trophic levels - all species 

Dietdata$Primary_diet <- factor(Dietdata$Primary_diet, levels=c("VE", "IN", "FR|NE", "FR", "NE", "PL|SE", "PL", "SE","OM"))

## percentage of species
pPD <- ggplot(Dietdata, aes(x=Primary_diet, y= ..prop.., group = 1), stat = "count") + geom_bar() + 
  scale_y_continuous(labels = scales::percent_format()) + 
  facet_wrap(~Class) + GGPoptions + xlab("") + ylab("Percentage of species") +
  #scale_x_discrete(labels=c("Vertivore", "Invertivore", "Frugivore/Nectarivore", "Plant/seed eater", "Omnivore" ,"NA"))+
  theme(axis.text.x = element_text(angle = 60,  hjust =1))
ggsave(pPD, filename = "../../Results/DietDataCoverage/diet_redefined/PD_all.pdf", width = 7, height= 6)


## number of species
Dietdata$Primary_diet
Dietdata$Primary_diet <- as.character(Dietdata$Primary_diet)
Dietdata$Primary_diet[Dietdata$Primary_diet=="SE"] <- "PL|SE"
Dietdata$Primary_diet[Dietdata$Primary_diet=="PL"] <- "PL|SE"
Dietdata$Primary_diet[Dietdata$Primary_diet=="NE"] <- "FR|NE"
Dietdata$Primary_diet[Dietdata$Primary_diet=="FR"] <- "FR|NE"


## plot for thesis / manuscript
Dietdata$Class <- factor(Dietdata$Class, levels = c("Birds", "Mammals", "Amphibians", "Reptiles"))
pPD2 <- ggplot(Dietdata, aes(x=Primary_diet), stat = "count") + geom_bar(fill="blue") + 
  #scale_y_continuous(labels = scales::percent_format()) + 
  facet_wrap(~Class) + GGPoptions + xlab("") + ylab("Number of species") +
  #scale_x_discrete(labels=c("Vertivore", "Invertivore", "Frugivore/Nectarivore", "Plant/seed eater", "Omnivore" ,"NA"))+
  theme(axis.text.x = element_text(angle = 60,  hjust =1)) + 
  geom_text( stat='count', aes(label=..count..), vjust=-1) + ylim(0, 12000) +
  ggtitle("Distribution of diets among vertebrate classes (collected data)") +
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold"))

 ggplot(Dietdata, aes(x=Primary_diet), stat = "count") + geom_bar(fill="darkblue") + 
  #scale_y_continuous(labels = scales::percent_format()) + 
  GGPoptions + xlab("") + ylab("Number of species") +
  #scale_x_discrete(labels=c("Vertivore", "Invertivore", "Frugivore/Nectarivore", "Plant/seed eater", "Omnivore" ,"NA"))+
  theme(axis.text.x = element_text(angle = 60,  hjust =1)) + 
  geom_text( stat='count', aes(label=..count..), vjust=-1) + ylim(0, 20000) +
  ggtitle("Distribution of diet data among vertebrates") +
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold"))
  
ggsave(pPD2, filename = "../../Results/DietDataCoverage/diet_redefined/Diet_distribution.pdf", width = 7, height= 6)


pTL <- ggplot(Dietdata, aes(x=Trophic_level, y= ..prop.., group = 1), stat = "count") + geom_bar(fill="darkgrey") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  facet_wrap(~Class) + GGPoptions + xlab("") + ylab("Percentage of species")
ggsave(pTL, filename = "../../Results/DietDataCoverage/diet_redefined/TL_all.pdf", width = 7, height= 4)


## Primary diet and TL - for PREDICTS species
Predicts <- readRDS("../../Data/PredictsVertebrates.rds")
Psp <- unique(Predicts$Best_guess_binomial)
Dietdata_Predicts <- Dietdata %>% 
  filter(Best_guess_binomial %in% Psp)

pPD_Predicts <- ggplot(Dietdata_Predicts, aes(x=Primary_diet, y= ..prop.., group = 1), stat = "count") + geom_bar() + 
  scale_y_continuous(labels = scales::percent_format()) + 
  facet_wrap(~Class) + GGPoptions + xlab("") + ylab("Percentage of species") +
  #scale_x_discrete(labels=c("Vertivore", "Invertivore", "Frugivore/Nectarivore", "Plant/seed eater", "Omnivore" ,"NA"))+
  theme(axis.text.x = element_text(angle = 60,  hjust =1))
ggsave(pPD_Predicts, filename = "../../Results/DietDataCoverage/diet_redefined/PD_Predicts.pdf", width = 7, height= 6)

pPD_Predicts$data

pTL_Predicts <- ggplot(Dietdata_Predicts, aes(x=Trophic_level, y= ..prop.., group = 1), stat = "count") + geom_bar(fill="darkgrey") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  facet_wrap(~Class) + GGPoptions + xlab("") + ylab("Percentage of species")
ggsave(pTL_Predicts, filename = "../../Results/DietDataCoverage/diet_redefined/TL_Predicts.pdf", width = 7, height= 4)


####################################################################################################################################

## trait coverage for all traits

Mammals <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Mammals.csv") %>% 
  mutate(Class="Mammals")
Birds <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Birds.csv") %>% 
  mutate(Class="Birds")
colnames(Birds)[37] <- "Trophic_level"
Amphibians <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Amphibians.csv") %>% 
  mutate(Class="Amphibians")
Reptiles <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Reptiles.csv") %>% 
  mutate(Class="Reptiles")

colnames(Mammals)[10] <- "Lifespan_proxy"
colnames(Birds)[11] <- "Lifespan_proxy"
colnames(Amphibians)[10] <- "Lifespan_proxy"
colnames(Reptiles)[9] <- "Lifespan_proxy"

colnames(Amphibians)[8] <- "Body_mass_g_true"
colnames(Amphibians)[5] <- "Body_mass_g"


Traits_cat <- c("Diel_activity", "Habitat_breadth_IUCN","Specialisation", "Primary_diet", "Diet_breadth")
Traits_cont <- c("Body_mass_g", "Lifespan_proxy", "Litter_size")
Traits <- c(Traits_cont, Traits_cat)

Dietdata <- rbind(Mammals[, c("Class", "Best_guess_binomial", Traits)], 
                  Birds[, c("Class", "Best_guess_binomial", Traits)], 
                  Amphibians[, c("Class", "Best_guess_binomial", Traits)], 
                  Reptiles[, c("Class", "Best_guess_binomial", Traits)])


Plot_coverage_V2 <- function(traits, FontSize){
  
  GGPoptions <- theme_classic()+ theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=FontSize, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=FontSize), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=FontSize),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=FontSize)) 
  
  Get_cov <- function(df, vclass) {
  
    df <- subset(df, Class==vclass) %>% 
      dplyr::select(-Best_guess_binomial, -Class)
    
    Coverage <- colnames(df) %>%
      as.data.frame() %>%
      setNames(., "Original")
    
    Coverage$FP <- c("Body size",
                     "Lifespan proxy",
                     "Litter/clutch size",
                     "Diel activity", 
                     "Habitat breadth", 
                     "Use of artificial habitats",
                     "Primary diet",
                     "Diet breadth")
    
    Coverage$Cov <- apply(df, 2,  function(y) sum(!is.na(y))) 
    Coverage$Cov <- Coverage$Cov/nrow(df)*100
    Coverage$Class <- vclass
    
    return(Coverage)
  }
  
  MammalCov <- Get_cov(traits, vclass="Mammals")
  BirdCov <- Get_cov(traits, vclass="Birds")
  AmphibianCov <- Get_cov(traits, vclass="Amphibians")
  ReptileCov <- Get_cov(traits, vclass="Reptiles")
  
  Coverage <- rbind(MammalCov, BirdCov, AmphibianCov, ReptileCov)
  
  # order by highest to lowest mean coverage
  Order <- Coverage %>%  group_by(FP) %>% 
    summarise(Mean=mean(Cov)) %>% 
    as.data.frame()
  Order <- Order[order(Order$Mean, decreasing = FALSE),]
  Order <- Order$FP
  
  Plotdata <- Coverage %>%  group_by(FP) %>% mutate(position=rank(Cov))
  
  # plot for all classes together
  p <- ggplot(Plotdata, aes(FP, Cov, fill=Class, group=position)) + 
    geom_bar(stat="identity", position="dodge", width=0.7) +
    GGPoptions +  
    scale_x_discrete(limits=Order)+
    coord_flip() + ylab("Coverage (%)") + xlab("") +
    scale_fill_viridis_d() + 
    geom_hline(yintercept=50, lty="dashed")
  
  return(p)
}

Coverageplot <- Plot_coverage_V2(Dietdata, 13)
ggsave(Coverageplot, filename="D:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/Coverage_all_species.pdf",
       width=8, height = 4)

####################################################################################################################################

## Phylogenetic signal in primary diet

rm(list=ls())

## Code is parallelised
library(parallel)
library(dplyr)
library(phytools)
source("Functions_Borges_et_al_2018_delta_statistic.R")

## START CLUSTER
Cluster <- makeCluster(6)

## EXCECUTE ANY PRE PROCESSING CODE NECESSARY
clusterEvalQ(Cluster, {
  library(dplyr)
  library(phytools)
  library(picante)
  library(geiger)
  library(ape)
  library(pbmcapply)
  library(pbapply)
})

# Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  return(Phylogeny)
}


# Function to apply for categorical traits: PhySignal_Cat. Based on Borges et al 2018, Bioinformatics
PhySignal_Cat <- function(Traitdata, Names, Phylo, n) {
  
  # browser()
  
  ## The Borges function does not work when branches have 0 length
  ## in that case, add a very small number to these branches
  #Phylo$edge.length[Phylo$edge.length==0] <- 10e-10
  
  # n = number of simulations
  ## For the current trait, match and prune phylogeny
  Traitdata <- as.character(Traitdata)
  names(Traitdata) <- Names
  Match <- picante::match.phylo.data(Phylo, Traitdata)
  Phylo <- Match$phy
  Trait <- Match$data
  rm(Match)
  
  ## Priors and parameters
  lambda0 <- 0.1   #rate parameter of the proposal 
  se      <- 0.5   #standard deviation of the proposal
  sim     <- 10000 #number of iterations
  thin    <- 10    #we kept only each 10th iterate 
  burn    <- 100   #100 iterates are burned-in
  
  ## Run the function (Borges et al, 2018, Bioinformatics: Measuring phylogenetic signal between categorical traits and phylogenies)
  ## with trycatch to avoid the process crashing when encountering errors.
  Delta <- tryCatch(expr={delta(Trait, Phylo, lambda0, se, sim, thin, burn)}, error = function(e) {NA})
  
  ## If the signal is not NA, then calculate a null distribution of delta values for the trait
  if(!is.na(Delta)) {
    
    # Generate randomised trait values - n times to generate a null distribution of delta -- stored in a list
    Func <- function(Trait){
      N <- length(Trait)
      L <- levels(as.factor(Trait))
      return(sample(L, size=N, replace=TRUE))
    }
    ListRandom <- lapply(rep(list(Trait), n), Func)
    
    Func_delta_toapply <- function(trait, tree, lambda0, se, sim, thin, burn) {
      Result <- tryCatch(expr = {delta(trait, tree, lambda0, se, sim, thin, burn)}, 
                         error=function(e){NA})
      return(Result)
    }
    
    Random_Delta <- pbmapply(FUN=Func_delta_toapply,
                             trait=ListRandom,
                             tree=rep(list(Phylo), n),
                             lambda0=rep(list(lambda0), n),
                             se=rep(list(se), n),
                             sim=rep(list(sim), n),
                             thin=rep(list(thin), n),
                             burn=rep(list(burn), n)) %>%
      as.data.frame()
    return(list(Delta=Delta, Delta0=Random_Delta))
    
  }
  
  else{return(Delta)}
}


## Load phylogenies, the ones used in Global Gaps article (consensus trees).
Phylo_Mammals <- ape::read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Mammals_complete_TreeAnnotator.nex") %>% .Format_tiplabels()
Phylo_Amphibians <- ape::read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Amphibians_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Birds <- ape::read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Birds_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Reptiles <- ape::read.nexus("../../Data/Phylogenies_GlobalGaps_Consensus_Trees_TreeAnnotator/Reptiles_TreeAnnotator.nex")  %>% .Format_tiplabels()

## Load data
Mammals <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Mammals.csv") %>% 
  mutate(Class="Mammals")
Birds <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Birds.csv") %>% 
  mutate(Class="Birds")
Amphibians <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Amphibians.csv") %>% 
  mutate(Class="Amphibians")
Reptiles <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Reptiles.csv") %>% 
  mutate(Class="Reptiles")


## measuring delta

# Names
Names.Mammals <- Mammals$Best_guess_binomial
Names.Birds <- Birds$Best_guess_binomial
Names.Reptiles <- Reptiles$Best_guess_binomial
Names.Amphibians <- Amphibians$Best_guess_binomial

clusterExport(cl=Cluster, list(".Format_tiplabels", "PhySignal_Cat",
                               "Mammals", "Birds", "Reptiles", "Amphibians",
                               "Phylo_Birds", "Phylo_Reptiles", "Phylo_Amphibians","Phylo_Mammals", 
                               "Names.Mammals", "Names.Birds", "Names.Reptiles", "Names.Amphibians",
                               "nentropy", "lpalpha", "lpbeta", "mhalpha", "mhbeta", "emcmc", "ratematrix", "rtrait", "delta"), envir=environment())


# 23 hours for birds
print(Sys.time())
start_time <- Sys.time()
delta_Birds_PD <- parApply(Cluster, Birds[, c("Primary_diet", "Trophic_level.Elton")], 2, PhySignal_Cat, Names=Names.Birds, Phylo=Phylo_Birds, n=50)
end_time <- Sys.time()
saveRDS(delta_Birds_PD, "../../Results/Phylogenetic_signal_diet/Diet_redefined/Birds.rds")

# 36 mins for mammals
print(Sys.time())
start_time <- Sys.time()
delta_Mammals_PD <- parApply(Cluster, Mammals[, c("Primary_diet", "Trophic_level")], 2, PhySignal_Cat, Names=Names.Mammals, Phylo=Phylo_Mammals, n=50)
end_time <- Sys.time()
saveRDS(delta_Mammals_PD, "../../Results/Phylogenetic_signal_diet/Diet_redefined/Mammals.rds")
print(start_time-end_time)

# 1.2 hours for reptiles
print(Sys.time())
start_time <- Sys.time()
delta_Reptiles_PD <- parApply(Cluster, Reptiles[, c("Primary_diet", "Trophic_level")], 2, PhySignal_Cat, Names=Names.Reptiles, Phylo=Phylo_Reptiles, n=50)
end_time <- Sys.time()
saveRDS(delta_Reptiles_PD, "../../Results/Phylogenetic_signal_diet/Diet_redefined/Reptiles.rds")

# 3.6 hours for amphibians
print(Sys.time())
start_time <- Sys.time()
delta_Amphibians_PD <- parApply(Cluster, Amphibians[, c("Primary_diet", "Trophic_level")], 2, PhySignal_Cat, Names=Names.Amphibians, Phylo=Phylo_Amphibians, n=50)
end_time <- Sys.time()
saveRDS(delta_Amphibians_PD, "../../Results/Phylogenetic_signal_diet/Diet_redefined/Amphibians.rds")

## DESTROY CLUSTER
stopCluster(Cluster)


### Loading phylogenetic signals

DeltaA <- readRDS("../../Results/Phylogenetic_signal_diet/Diet_redefined/Amphibians.rds")
wilcox.test(mu=DeltaA$Trophic_level$Delta, DeltaA$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=DeltaA$Primary_diet$Delta, DeltaA$Primary_diet$Delta0$., alternative = "less")

DeltaR <- readRDS("../../Results/Phylogenetic_signal_diet/Diet_redefined/Reptiles.rds")
wilcox.test(mu=DeltaR$Trophic_level$Delta, DeltaR$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=DeltaR$Primary_diet$Delta, DeltaR$Primary_diet$Delta0$., alternative = "less")

DeltaM <- readRDS("../../Results/Phylogenetic_signal_diet/Diet_redefined/Mammals.rds")
wilcox.test(mu=DeltaM$Trophic_level$Delta, DeltaM$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=DeltaM$Primary_diet$Delta, DeltaM$Primary_diet$Delta0$., alternative = "less")

DeltaM <- readRDS("../../Results/Phylogenetic_signal_diet/Mammals.rds")
wilcox.test(mu=DeltaM$Primary_diet$Delta, DeltaM$Primary_diet$Delta0$., alternative = "less")

DeltaB <- readRDS("../../Results/Phylogenetic_signal_diet/Diet_redefined/Birds.rds")
wilcox.test(mu=DeltaB$Trophic_level$Delta, DeltaB$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=DeltaB$Primary_diet$Delta, DeltaB$Primary_diet$Delta0$., alternative = "less")



# DeltaA <- readRDS("../../Results/Phylogenetic_signal_diet/Diet_redefined/Amphibians.rds")
# wilcox.test(mu=DeltaA$Trophic_level$Delta, DeltaA$Trophic_level$Delta0$., alternative = "less")
# wilcox.test(mu=DeltaA$Primary_diet$Delta, DeltaA$Primary_diet$Delta0$., alternative = "less")
# 
# DeltaR <- readRDS("../../Results/Phylogenetic_signal_diet/Reptiles.rds")
# wilcox.test(mu=DeltaR$Trophic_level$Delta, DeltaR$Trophic_level$Delta0$., alternative = "less")
# wilcox.test(mu=DeltaR$Primary_diet$Delta, DeltaR$Primary_diet$Delta0$., alternative = "less")
# 
# DeltaM <- readRDS("../../Results/Phylogenetic_signal_diet/Diet_redefined/Mammals.rds")
# wilcox.test(mu=DeltaM$Trophic_level$Delta, DeltaM$Trophic_level$Delta0$., alternative = "less")
# wilcox.test(mu=DeltaM$Primary_diet$Delta, DeltaM$Primary_diet$Delta0$., alternative = "less")
# 
# DeltaM <- readRDS("../../Results/Phylogenetic_signal_diet/Mammals.rds")
# wilcox.test(mu=DeltaM$Primary_diet$Delta, DeltaM$Primary_diet$Delta0$., alternative = "less")
# 
# DeltaB <- readRDS("../../Results/Phylogenetic_signal_diet/Birds.rds")
# wilcox.test(mu=DeltaB$Trophic_level$Delta, DeltaB$Trophic_level$Delta0$., alternative = "less")
# wilcox.test(mu=DeltaB$Primary_diet$Delta, DeltaB$Primary_diet$Delta0$., alternative = "less")


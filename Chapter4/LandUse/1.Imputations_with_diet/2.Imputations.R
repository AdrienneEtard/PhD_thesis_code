# Imputations using missForest, with  phylogenetic eigenvectors as predictors
# TAKES ~14 HOURS
# Code is parallelised
# https://www.jottr.org/2018/06/23/future.apply_1.0.0/

library(parallel)
library(dplyr)

setwd("D:/3.Explanatory_traits/Code/1.Imputations_with_diet/")

## START CLUSTER
Cluster <- makeCluster(detectCores()-4)

## EXCECUTE ANY PRE PROCESSING CODE NECESSARY
clusterEvalQ(Cluster, {
  library(dplyr)
  library(phytools)
  library(missForest)
  #library(pbmcapply)
  library(pbapply)
})

## Preamble
`%nin%` <- Negate(`%in%`)
source("Functions.R")

## Load trait data with phylogenetic imformation as eigenvectors
Mammals <- read.csv("../../Results/Traits_with_phy_eigenvectors/Mammals.csv") %>% 
   dplyr::select(-Trophic_level.Elton)
Birds <- read.csv("../../Results/Traits_with_phy_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/Traits_with_phy_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/Traits_with_phy_eigenvectors/Reptiles.csv")

colnames(Birds)[colnames(Birds)=="Trophic_level.Elton"] <- "Trophic_level"

## Define variables for imputations

Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")

Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
             "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")

Traits_cat <- c("Diel_activity","Trophic_level", "Specialisation", Habitat, Diet, "Primary_diet")

# Family and Genus have too many factor levels, so only include order information.
Taxinfo <- c("Order") #, "Family", "Genus")

MammalsCont <- colnames(Mammals)[c(5,6,7,10,11)]
BirdsCont <- colnames(Birds)[c(5,6,7,8, 11)]
ReptilesCont <- colnames(Reptiles)[c(5,6,7,8,9,10)]
AmphibiansCont <- colnames(Amphibians)[c(5,6,7,8,9,10)]

## add diet breadth as a continuous traits to impute.
Mammals$Diet_breadth <- apply(Mammals[, Diet], 1, sum)
Amphibians$Diet_breadth <- apply(Amphibians[, Diet], 1, sum)
Reptiles$Diet_breadth <- apply(Reptiles[, Diet], 1, sum)
Birds$Diet_breadth <- apply(Birds[, Diet], 1, sum)

Mammals$Diet_breadth <- as.factor(Mammals$Diet_breadth)
Amphibians$Diet_breadth <- as.factor(Amphibians$Diet_breadth)
Birds$Diet_breadth <- as.factor(Birds$Diet_breadth)
Reptiles$Diet_breadth <- as.factor(Reptiles$Diet_breadth)

## add habitat breadth as a continuous traits to impute.
Mammals$Habitat_breadth_IUCN <- apply(Mammals[, Habitat], 1, sum)
Amphibians$Habitat_breadth_IUCN <- apply(Amphibians[, Habitat], 1, sum)
Reptiles$Habitat_breadth_IUCN <- apply(Reptiles[, Habitat], 1, sum)
Birds$Habitat_breadth_IUCN <- apply(Birds[, Habitat], 1, sum)

# set as NA when 0
Mammals$Habitat_breadth_IUCN[Mammals$Habitat_breadth_IUCN==0] <- NA
Amphibians$Habitat_breadth_IUCN[Amphibians$Habitat_breadth_IUCN==0] <- NA
Reptiles$Habitat_breadth_IUCN[Reptiles$Habitat_breadth_IUCN==0] <- NA
Birds$Habitat_breadth_IUCN[Birds$Habitat_breadth_IUCN==0] <- NA

# set as factor
Mammals$Habitat_breadth_IUCN <- as.factor(Mammals$Habitat_breadth_IUCN)
Amphibians$Habitat_breadth_IUCN <- as.factor(Amphibians$Habitat_breadth_IUCN)
Birds$Habitat_breadth_IUCN <- as.factor(Birds$Habitat_breadth_IUCN)
Reptiles$Habitat_breadth_IUCN <- as.factor(Reptiles$Habitat_breadth_IUCN)

## Traits categorical
Traits_cat <-c(Traits_cat, "Diet_breadth", "Habitat_breadth_IUCN")

## phylogenetic eigenvectors
EV <- c("EV_1","EV_2", "EV_3", "EV_4", "EV_5", "EV_6", "EV_7", "EV_8", "EV_9", "EV_10")

## Function arguments as lists, nested into one bigger list - each of these list elements are agurments for the function Imputations_missForest
DF.TraitsList <- list(M=Mammals, B=Birds, R=Reptiles, A=Amphibians)
Taxinfo.List <- list(M="Order", B="Order", R="Order", A="Order")
Cont.TraitsList <- list(M=MammalsCont, B=BirdsCont, R=ReptilesCont, A=AmphibiansCont)
Cat.TraitsList <- list(M=Traits_cat, B=Traits_cat, R=Traits_cat, A=Traits_cat)
EV.List <- list(M=EV, B=EV, R=EV, A=EV)
ErrorTrue.List <- list(M=TRUE, B=TRUE, R=TRUE, A=TRUE)
DietTRUE.List <- list(M=TRUE, B=TRUE, R=TRUE, A=TRUE)

# List of function arguments. This list  will be replicated 8 times (number of cores) for parallel imputations.
# On each cluster, imputation of 4 datasets (one for each class). 
# ArgumentsList_st <- list(TraitDF=DF.TraitsList_st,
# 					            Taxinfo=Taxinfo.List,
# 					            Traits_cont=Cont.TraitsList_st,
# 					            Traits_cat=Cat.TraitsList,
# 					            EV=EV.List,
# 					            ErrorTrue=ErrorTrue.List,
# 					            DietTRUE=DietTRUE.List, 
# 					            std=std.List_st)

ArgumentsList <- list(TraitDF=DF.TraitsList,
                      Taxinfo=Taxinfo.List,
                      Traits_cont=Cont.TraitsList,
                      Traits_cat=Cat.TraitsList,
                      EV=EV.List,
                      ErrorTrue=ErrorTrue.List,
                      DietTRUE=DietTRUE.List)

# Replicate this list N times so that: list with N elements, each of these are ArgumentsLists
N <- 8
To_impute_parallel <- rep(list(ArgumentsList), N)

rm(Mammals, Birds, Reptiles, Amphibians,
   Habitat, Taxinfo, Traits_cat, 
   MammalsCont, BirdsCont, ReptilesCont, AmphibiansCont,
   DF.TraitsList, Taxinfo.List, Cont.TraitsList, Cat.TraitsList,
   EV.List, ErrorTrue.List, DietTRUE.List, N,
   ArgumentsList)

## Export variables in all clusters
clusterExport(cl=Cluster, list("Imputations_missForest",
                               "To_apply_parallel_imputations",
                               "To_impute_parallel",
                               "%nin%"),
              envir=environment())

## Parallel imputations on 8 cores (takes about 1.2 day)
# print("Imputations on tranformed and standardised data")
# system.time(Imputed_sets_st <- parLapply(cl=Cluster,
#                           X=To_impute_parallel_st,
#                           fun=To_apply_parallel_imputations))

print("Imputations with diet (For species-level analyses")

Start <- Sys.time()
print(Start)
Imputed_sets <- parLapply(cl=Cluster,
                          X=To_impute_parallel,
                          fun=To_apply_parallel_imputations)

End <- Sys.time()
print(Start-End)
## Save results
saveRDS(Imputed_sets, "../../Results/Imputed_traits/List_of_8_sets.rds")

## STOP CLUSTER
stopCluster(Cluster)

table(Imputed_sets[[1]]$M$Imputed.Dataset$Primary_diet)
table(Imputed_sets[[1]]$M$Imputed.Dataset$Habitat_breadth_IUCN)
table(Imputed_sets[[1]]$M$Imputed.Dataset$Diet_breadth)


# ## test 
# TestFunction <- Imputations_missForest(Amphibians[1:100,],
#                        Taxinfo = "Order", 
#                        Traits_cont = AmphibiansCont, 
#                        Traits_cat=Traits_cat,
#                        EV = EV, 
#                        ErrorTrue = TRUE, 
#                        DietTRUE = TRUE)
# 
# TestFunction[[1]]$Imputed.Dataset$Habitat_breadth_IUCN %>%  levels()
# TestFunction[[1]]$Imputed.Dataset$Diet_breadth %>%  levels()

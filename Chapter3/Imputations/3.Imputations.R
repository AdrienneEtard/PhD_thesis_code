# Imputations of missing trait values using missForest, with  phylogenetic eigenvectors as predictors
# TAKES ~14 HOURS
# Code is parallelised
# https://www.jottr.org/2018/06/23/future.apply_1.0.0/
# traits need not be transformed as missForest is a non parametric approach

library(parallel)

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
Mammals <- read.csv("../../Results/Datasets_FD_imputations/Mammals.csv")
Birds <- read.csv("../../Results/Datasets_FD_imputations/Birds.csv")
Amphibians <- read.csv("../../Results/Datasets_FD_imputations/Amphibians.csv")
Reptiles <- read.csv("../../Results/Datasets_FD_imputations/Reptiles.csv")

## Define variables for imputations

# Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")

Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
             "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")

Traits_cat <- c("Diel_activity","Trophic_level", "Specialisation", Habitat)

# Family and Genus have too many factor levels to be used with missForest, so only include order information.
Taxinfo <- c("Order") #, "Family", "Genus")

# names of categorical traits
MammalsCont <- colnames(Mammals)[c(5,6,7,10)]
BirdsCont <- colnames(Birds)[c(5,6,7,10)]
ReptilesCont <- colnames(Reptiles)[c(5,6,7,8,9,12)]
AmphibiansCont <- colnames(Amphibians)[c(5,6,7,8,9,12)]

# eigenvectors
EV <- c("EV_1","EV_2", "EV_3", "EV_4", "EV_5", "EV_6", "EV_7", "EV_8", "EV_9", "EV_10")

## Function arguments as lists, nested into one bigger list - each of these list elements are agurments for the function Imputations_missForest
DF.TraitsList <- list(M=Mammals, B=Birds, R=Reptiles, A=Amphibians)
Taxinfo.List <- list(M="Order", B="Order", R="Order", A="Order")
Cont.TraitsList <- list(M=MammalsCont, B=BirdsCont, R=ReptilesCont, A=AmphibiansCont)
Cat.TraitsList <- list(M=Traits_cat, B=Traits_cat, R=Traits_cat, A=Traits_cat)
EV.List <- list(M=EV, B=EV, R=EV, A=EV)
ErrorTrue.List <- list(M=TRUE, B=TRUE, R=TRUE, A=TRUE)
DietTRUE.List <- list(M=FALSE, B=FALSE, R=FALSE, A=FALSE)
std.List <- list(M=FALSE, B=FALSE, R=FALSE, A=FALSE)

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
                      DietTRUE=DietTRUE.List,
                      std=std.List)

# Replicate this list N times so that: list with N elements, each of these are ArgumentsLists
N <- 8
To_impute_parallel <- rep(list(ArgumentsList), N)

rm(Mammals, Birds, Reptiles, Amphibians,
   Habitat, Taxinfo, Traits_cat,
   MammalsCont, BirdsCont, ReptilesCont, AmphibiansCont,
   DF.TraitsList, Taxinfo.List, Cont.TraitsList, Cat.TraitsList,
   EV.List, ErrorTrue.List, DietTRUE.List, N, std.List,
   ArgumentsList)

## Export variables in all clusters
clusterExport(cl=Cluster, list("Imputations_missForest",
                               "To_apply_parallel_imputations",
                               "To_impute_parallel",
                               "%nin%"),
              envir=environment())

## Parallel imputations on 8 cores (15 hours)
# print("Imputations on tranformed and standardised data")
# system.time(Imputed_sets_st <- parLapply(cl=Cluster,
#                           X=To_impute_parallel_st,
#                           fun=To_apply_parallel_imputations))

print("Imputations without diet (For functional diversity analyses")

Imputed_sets <- parLapply(cl=Cluster,
                          X=To_impute_parallel,
                          fun=To_apply_parallel_imputations)

## Save results
saveRDS(Imputed_sets, "../../Results/Imputed_sets/List_of_8_sets.rds")

## STOP CLUSTER
stopCluster(Cluster)

## compute dbFD indices -- FRic, FDis, FDiv, RaoQ and compute null expectations

# # # # # # # DO NOT PARALLELISE THESE BECAUSE OF THE "VERT.TXT" FILE GENERATED DURING THE PROCESS # # # # # # # # #

library(dplyr)
library(FD)
library(vegan)
library(swfscMisc)
library(funrar)
library(lme4)
library(StatisticalModels)
source("Functions.R")
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

To_run <- function(List_of_communities,
                   Gowerdist,
                   Abundance,
                   SimTRUE,
                   nsim,
                   sim_pool,
                   global_pool_character,
                   Std.FRic,
                   Predicts_site_info)
  {

  List <- list()

  for (i in 1:length(List_of_communities)) {

    # print(Sys.time())
    # print(paste("Starting study", i, "of", length(List_of_communities)))

    List[[i]] <- dbFD_to_apply(Distmatrix=Gowerdist,
                                 community_matrix=List_of_communities[[i]],
                                 Abundance_weighted=Abundance,
                                 Simulations=SimTRUE,
                                 nsim=nsim,
                                 sim_pool=sim_pool,
                                 global_pool_character=global_pool_character,
                                 Std.FRic=Std.FRic,
                                 Predicts_site_info=Predicts_site_info)

    print(paste("Finished study", i, "of", length(List_of_communities)))
    #print(Sys.time())
  }
  return(List)
}


## Load data: Gower distance matrices for PREDICTS vertebrates and community matrices (presence-abscence)
Gower_8 <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets.rds")
Communities_Occu <- readRDS("../../Data/Presence_absence.rds") %>% na.omit.list

## Read Predicts site-level information with ecoregion information
Predicts_site_info <- read.csv("../../Results/Predicts_site_info_ER.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Metrics without simulations, using all of 8 gower distances from the set of 8 imputed trait datasets

## With standardisation: 43 minutes
Results_8_st <- list()

t1 <- Sys.time()

for (i in 1:8) {

  print(paste("Processing Gower data #", i))

  Res_i <- To_run(List_of_communities = Communities_Occu,
                  Gowerdist = Gower_8[[i]],
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = TRUE,
                  Predicts_site_info=Predicts_site_info)


  Results_8_st[[i]] <- Extract(Res_i, Predicts_site_info)
}

t2 <- Sys.time()

saveRDS(Results_8_st, "../../Results/dbFD_indices/8_sets.rds")

## across vertebrates, main results - with 8th set of imputed traits
Verts_res <- To_run(List_of_communities = Communities_Occu,
                Gowerdist = Gower_8[[8]],
                Abundance = FALSE,
                SimTRUE = FALSE,
                nsim=NULL,
                sim_pool = NULL,
                global_pool_character = NULL,
                Std.FRic = TRUE,
                Predicts_site_info=Predicts_site_info)

Verts_res <- Extract(Verts_res, Predicts_site_info)
write.csv(Verts_res, "../../Results/dbFD_indices/Verts_main_results.csv", row.names = FALSE)


## without standardisation

Results_8 <- list()

for (i in 1:8) {

  print(paste("Processing Gower data #", i))

  t1 <- Sys.time()
  Res_i <- To_run(List_of_communities = Communities_Occu,
                  Gowerdist = Gower_8[[i]],
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = FALSE,
                  Predicts_site_info=Predicts_site_info)
  t2 <- Sys.time()

  Results_8[[i]] <- Extract(Res_i, Predicts_site_info)
}

saveRDS(Results_8, "../../Results/dbFD_indices/8_sets_nst.rds")

## metrics by class, without simulations

Communities_Mammals <- readRDS("../../Data/Presence_absence_Mammals.rds") %>% na.omit.list
Communities_Birds <- readRDS("../../Data/Presence_absence_Birds.rds") %>% na.omit.list
Communities_Reptiles <- readRDS("../../Data/Presence_absence_Reptiles.rds") %>% na.omit.list
Communities_Amphibians <- readRDS("../../Data/Presence_absence_Amphibians.rds") %>% na.omit.list

GowerDist_Mammals <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Mammals.rds")
GowerDist_Birds <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Birds.rds")
GowerDist_Reptiles <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Reptiles.rds")
GowerDist_Amphibians <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Amphibians.rds")


## for mammals

Results_8_st_mammals <- list()

t1 <- Sys.time()

for (i in 1:8) {

  print(paste("Processing Gower data #", i))

  Res_i <- To_run(List_of_communities = Communities_Mammals,
                  Gowerdist = GowerDist_Mammals[[i]],
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = TRUE,
                  Predicts_site_info=Predicts_site_info)

  Results_8_st_mammals[[i]] <- Extract(Res_i, Predicts_site_info)
}

t2 <- Sys.time()

saveRDS(Results_8_st_mammals, "../../Results/dbFD_indices/8_sets_Mammals.rds")

## for birds

Results_8_st_birds <- list()

t1 <- Sys.time()

for (i in 1:8) {

  print(paste("Processing Gower data #", i))

  Res_i <- To_run(List_of_communities = Communities_Birds,
                  Gowerdist = GowerDist_Birds[[i]],
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = TRUE,
                  Predicts_site_info=Predicts_site_info)

  Results_8_st_birds[[i]] <- Extract(Res_i, Predicts_site_info)
}

t2 <- Sys.time()

saveRDS(Results_8_st_birds, "../../Results/dbFD_indices/8_setsBirds.rds")


## for amphibians

Results_8_st_amphibians <- list()

t1 <- Sys.time()

for (i in 1:8) {

  print(paste("Processing Gower data #", i))

  Res_i <- To_run(List_of_communities = Communities_Amphibians,
                  Gowerdist = GowerDist_Amphibians[[i]],
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = TRUE,
                  Predicts_site_info=Predicts_site_info)

  Results_8_st_amphibians[[i]] <- Extract(Res_i, Predicts_site_info)
}

t2 <- Sys.time()

saveRDS(Results_8_st_amphibians, "../../Results/dbFD_indices/8_setsAmphibians.rds")


## for reptiles

Results_8_st_reptiles <- list()

t1 <- Sys.time()

for (i in 1:8) {

  print(paste("Processing Gower data #", i))

  Res_i <- To_run(List_of_communities = Communities_Reptiles,
                  Gowerdist = GowerDist_Reptiles[[i]],
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = TRUE,
                  Predicts_site_info=Predicts_site_info)

  Results_8_st_reptiles[[i]] <- Extract(Res_i, Predicts_site_info)
}

t2 <- Sys.time()

saveRDS(Results_8_st_reptiles, "../../Results/dbFD_indices/8_setsReptiles.rds")

### main results datasets for each class

Mammals <- To_run(List_of_communities = Communities_Mammals,
       Gowerdist = GowerDist_Mammals[[8]],
       Abundance = FALSE,
       SimTRUE = FALSE,
       nsim=NULL,
       sim_pool = NULL,
       global_pool_character = NULL,
       Std.FRic = TRUE,
       Predicts_site_info=Predicts_site_info)
Mammals <- Extract(Mammals, Predicts_site_info)
Mammals$Class <- "Mammals"
write.csv(Mammals, "../../Results/dbFD_indices/Mammals_main_results.csv", row.names = FALSE)

Birds <- To_run(List_of_communities = Communities_Birds,
                  Gowerdist = GowerDist_Birds[[8]],
                  Abundance = FALSE,
                  SimTRUE = FALSE,
                  nsim=NULL,
                  sim_pool = NULL,
                  global_pool_character = NULL,
                  Std.FRic = TRUE,
                  Predicts_site_info=Predicts_site_info)
Birds <- Extract(Birds, Predicts_site_info)
Birds$Class <- "Birds"
write.csv(Birds, "../../Results/dbFD_indices/Birds_main_results.csv", row.names = FALSE)

Amphibians <- To_run(List_of_communities = Communities_Amphibians,
                Gowerdist = GowerDist_Amphibians[[8]],
                Abundance = FALSE,
                SimTRUE = FALSE,
                nsim=NULL,
                sim_pool = NULL,
                global_pool_character = NULL,
                Std.FRic = TRUE,
                Predicts_site_info=Predicts_site_info)

Amphibians <- Extract(Amphibians, Predicts_site_info)
Amphibians$Class <- "Amphibians"
write.csv(Amphibians, "../../Results/dbFD_indices/Amphibians_main_results.csv", row.names = FALSE)

Reptiles <- To_run(List_of_communities = Communities_Reptiles,
                Gowerdist = GowerDist_Reptiles[[8]],
                Abundance = FALSE,
                SimTRUE = FALSE,
                nsim=NULL,
                sim_pool = NULL,
                global_pool_character = NULL,
                Std.FRic = TRUE,
                Predicts_site_info=Predicts_site_info)
Reptiles <- Extract(Reptiles, Predicts_site_info)
Reptiles$Class <- "Reptiles"
write.csv(Reptiles, "../../Results/dbFD_indices/Reptiles_main_results.csv", row.names = FALSE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # WITH SIMULATIONS: USING THE 8TH SET OF IMPUTED TRAITS (8TH GOWER DISTANCE MATRIX) # # # # # #

Gower_distances <- Gower_8[[8]]

To_run_Regional <- function(Gowerdist,
                            List_of_communities,
                            SimTRUE,
                            nsim,
                            sim_pool,
                            Std.FRic,
                            Predicts_site_info)
{

  List <- list()

  for (i in 1:length(List_of_communities)) {

    List[[i]] <- dbFD_to_apply_Regional(Distmatrix = Gower_distances,
                                        community_matrix=List_of_communities[[i]],
                                        Simulations=SimTRUE,
                                        nsim=nsim,
                                        sim_pool=sim_pool,
                                        Std.FRic=Std.FRic,
                                        Predicts_site_info=Predicts_site_info)

    print(paste("Finished study", i, "of", length(List_of_communities)))
    #print(Sys.time())
  }
  return(List)
}

## with FRic standardised ~ 6 hours for 100 simulations, or 1.17 days for 500

t1 = Sys.time()
dbFD_Regional.std <- To_run_Regional(List_of_communities = Communities_Occu,
                               Gowerdist = Gower_distances,
                               SimTRUE = TRUE,
                               nsim=500,
                               sim_pool = "regional",
                               Std.FRic = TRUE,
                               Predicts_site_info=Predicts_site_info)

t2 = Sys.time()
print(t1-t2)

saveRDS(dbFD_Regional.std, "../../Results/dbFD_indices/Sim_Region/All_results_std.rds")
ResultsReg.std <- Extract(dbFD_Regional.std, Predicts_site_info)
write.csv(ResultsReg.std, "../../Results/dbFD_indices/Sim_Region/dbFD_Predicts_sites_std.csv", row.names=FALSE)

## with FRic not standardised

t1 = Sys.time()
dbFD_Regional <- To_run_Regional(List_of_communities = Communities_Occu,
                                     Gowerdist = Gower_distances,
                                     SimTRUE = TRUE,
                                     nsim=100,
                                     sim_pool = "regional",
                                     Std.FRic = FALSE,
                                     Predicts_site_info=Predicts_site_info)

t2 = Sys.time()
print(t1-t2)

saveRDS(dbFD_Regional, "../../Results/dbFD_indices/Sim_Region/All_results.rds")
ResultsReg <- Extract(dbFD_Regional, Predicts_site_info)
write.csv(ResultsReg, "../../Results/dbFD_indices/Sim_Region/dbFD_Predicts_sites.csv", row.names=FALSE)

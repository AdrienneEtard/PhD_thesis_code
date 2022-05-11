## Implementing CENFA (Climate & Ecological Niche Factor Analysis) to estimate species' sensitivity to CC

## allows an estimation of sensitivity to climate change that combines 
## a marginality metric
## a specialisation metric
## those metrics are based on climatic niche properties for a species - compared to the overall, global climate niche (or that of an area of reference) 

## here test at 5 km resolution

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(dplyr)
library(CENFA)
library(car)
library(pedometrics)
library(virtualspecies)
library(stringr)
library(pbapply)
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

#setwd("C:/Users/Adrienne/OneDrive - University College London/PhD/PhD_R_work/5.Climatic_niche_space/Code/")
behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
setwd("F:/PhD/PhD_R_projects/5.Climatic_niche_space/Code/")
# Index to species name dataframe
Index <- read.csv("../../Range_maps_work/Results/7_GenerateSpeciesIndices/SpeciesIndices.csv")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## load climate variables, at 2.5 arc-minute resolution (~4.6 km sq at the equator), and select those that are not too collinear
WCfiles <- paste("../Data/wc2.1_2.5m_bio", list.files(path = "../Data/wc2.1_2.5m_bio"), sep="/")
Climatevars <- stack(WCfiles)
names(Climatevars) <- c("Annual_mean_temperature",
                        "Mean_temperature_warmest_quarter",
                        "Mean_temperature_coldest_quarter",
                        "Annual_precipitation",
                        "Precipitation_wettest_month",
                        "Precipitation_driest_month",
                        "Precipitation_seasonality",
                        "Precipitation_wettest_quarter",
                        "Precipitation_driest_quarter",
                        "Precipitation_warmest_quarter",
                        "Precipitation_coldest_quarter",
                        "Mean_diurnal_range",
                        "Isothermality",
                        "Temperature_seasonality",
                        "Max_temperature_warmest_month",
                        "Min_temperature_coldest_month",
                        "Temperature_annual_range",
                        "Mean_temperature_wettest_quarter",
                        "Mean_temperature_driest_quarter")

# load distribution for a species -- at 5km resolution, in Behrman projection 
Mammals_ex <- raster("../../Range_maps_work/Results/11_AggregateRasters_5km/Mammals0/sp17890")
crs(Mammals_ex) <- behrCRS

# realign projections of climate rasters to that of species distribution files (5 km square, Behrman)
Climatevars_Behrman_5 <- projectRaster(from=Climatevars, to=Mammals_ex, method = "bilinear")
crs(Climatevars_Behrman_5) <- behrCRS

# # select climatic variables that are not too correlated
par(family="serif")
Selection <- removeCollinearity(Climatevars_Behrman_5,
                                plot = TRUE,
                                sample.points = FALSE,
                                method = "spearman",
                                multicollinearity.cutoff = 0.65)


Selected_Climvars_5 <- Climatevars_Behrman_5[[c("Max_temperature_warmest_month", 
                                                  "Mean_diurnal_range",
                                                  "Annual_mean_temperature",
                                                  "Precipitation_seasonality",
                                                  "Annual_precipitation",
                                                  "Precipitation_coldest_quarter")]]

plot(Selected_Climvars_5)

writeRaster(Selected_Climvars_5, filename="../Results/Selected_Climvars_5.grd")
saveRDS(Selected_Climvars_5, "../Results/Selected_Climvars_5.rds")

# Calculate global covariance matrix at 5 km
Selected_Climvars_5 <- readRDS("../Results/Selected_Climvars_5.rds")

Global_cov_5 <- GLcenfa(Selected_Climvars_5,
                         progress=TRUE, 
                         center=TRUE,
                         scale=TRUE, 
                         filename="../Results/CENFA/GlobalClimCov5",
                         overwrite=TRUE)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## FUNCTIONS

## functions for running the CENFA analysis

CENFA_analysis <- function(GlobalCovMatrix, X, IndexDF, resolution, Class) {
  
  ## load species distribution raster
  Distrib <- raster::raster(X)
  behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
  crs(Distrib) <- behrCRS
  
  if(resolution==50) {
    Distrib <- aggregate(Distrib, fact=5)
  }
  
  ## run CENFA analysis
  Cenfa <- tryCatch(expr={cnfa(x = GlobalCovMatrix,
                               s.dat = Distrib,
                               progress = TRUE,
                               scale = FALSE)},
                    error = function(e){print("error")})
  
  if(class(Cenfa)=="character") {
    print("error in CENFA analysis")
  } else {
    
    ## save CENFA results
    if(resolution==50){
      pathsave <- paste0("../Results/CENFA/50_km/", Class, "/", names(Distrib), ".rds")
      saveRDS(Cenfa, pathsave)
    }
    
    if(resolution==10){
      pathsave <- paste0("../Results/CENFA/10_km/", Class, "/" ,names(Distrib), ".rds")
      saveRDS(Cenfa, pathsave)
    }
    if(resolution==5){
      pathsave <- paste0("../Results/CENFA/5_km/", Class, "/" ,names(Distrib), ".rds")
      saveRDS(Cenfa, pathsave)
    }
    
  }
}

Run_CENFA <- function(GlobalCovMatrix, ListPaths, IndexDF, resolution, Class) {
  for(i in 1:length(ListPaths)) {
    CENFA_analysis(GlobalCovMatrix, X=ListPaths[i], IndexDF, resolution, Class)
    print(paste(i, "over", length(ListPaths)))
  }
}

## extrating results and storing in a dataframe 
Extract_results <- function(path) {
  species <- str_split(path, pattern="/") %>% 
    unlist()
  species <- species[length(species)]
  species <- gsub(species, pattern=".rds", replacement = "")
  res <- readRDS(path)
  Res_df <- data.frame(sensitivity=res@sensitivity, marginality=res@marginality, Species=species)
  return(Res_df)
} 

## FUNCTION TO get species for which CENFA was run
Get_remaining_species <- function(Distrib_files, Class) {
  
  pathDone <- paste0("../Results/CENFA/5_km/", Class, "/")
  Done_files <- dir(pathDone, recursive = TRUE, full.names = TRUE)
  
  print(paste(length(Distrib_files) - length(Done_files), "species remaining"))
  
  ## process done files
  Done_files <- lapply(Done_files, FUN = function(x) {
    x <- strsplit(x, "/") %>%
      unlist()
    return(x[length(x)]) 
  }) %>% 
    unlist()
  Done_files <- gsub(Done_files, pattern=".rds", replacement = "")
  
  ## process all files
  Dist_files <- lapply(Distrib_files, FUN = function(x) {
    x <- strsplit(x, "/") %>%
      unlist()
    return(x[length(x)]) 
  }) %>% 
    unlist()
  Dist_files <- gsub(Dist_files, pattern=".rds", replacement = "")
  
  ## intersect
  Left_to_run <- setdiff(Dist_files, Done_files)
  
  ## get the paths for these species
  FT <- Dist_files %in% Left_to_run
  Paths_Left <- Distrib_files[FT]
  return(Paths_Left)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## List paths to all species distribution files
files <- dir("../../Range_maps_work/Results/11_AggregateRasters_5km/", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("info", files)]
files <- lapply(files, FUN = function(x) {
  x <- strsplit(x, "/") %>%
    unlist()
  return(paste0(x[1:(length(x)-1)], collapse = "/"))
})
files <- unique(unlist(files))
files <- files[grepl("sp", files)]

filesAmphibians <- files[grepl("Amphibians", files)]
filesBirds <- files[grepl("Birds", files)]
filesMammals <- files[grepl("Mammals", files)]
filesReptiles <- files[grepl("Reptiles", files)]

# ## test
# T1 <- Sys.time()
# Distrib <- raster(filesMammals[29])
# crs(Distrib) <- behrCRS
# cnfa(x = Global_cov_5,
#      s.dat = Distrib,
#      progress = TRUE,
#      scale = FALSE)
# T2 <- Sys.time()
# print(T2-T1)

## running CENFA at 5 km square resolution

## mammals
# filesMammalsLeft <- Get_remaining_species(filesMammals, "Mammals")
# length(filesMammalsLeft)
# gc(verbose = TRUE)
# Run_CENFA(Global_cov_5, ListPaths = filesMammalsLeft, IndexDF = Index, resolution = 5, Class = "Mammals")

# ## amphibians
# filesAmphibiansLeft <- Get_remaining_species(filesAmphibians, "Amphibians")
# length(filesAmphibiansLeft)
# gc(verbose = TRUE)
# Run_CENFA(Global_cov_5, ListPaths = filesAmphibiansLeft, IndexDF = Index, resolution = 5, Class = "Amphibians")

## birds filesBirdsLeft <- Get_remaining_species(filesBirds, "Birds")
#length(filesBirdsLeft) gc() Run_CENFA(Global_cov_5, ListPaths = filesBirdsLeft,
#IndexDF = Index, resolution = 5, Class = "Birds") #66 species at least

## reptiles
filesReptilesLeft <- Get_remaining_species(filesReptiles, "Reptiles")
length(filesReptilesLeft)
gc(verbose = TRUE)
Run_CENFA(Global_cov_5, ListPaths = rev(filesReptilesLeft), IndexDF = Index, resolution = 5, Class = "Reptiles")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# results at 5 km square resolution

Mammals_5 <- dir("../Results/CENFA/5_km/Mammals/", recursive = TRUE, full.names = TRUE)
length(Mammals_5)
Results_Mammals_5 <- pblapply(as.list(Mammals_5), Extract_results)
Results_Mammals_5 <- data.table::rbindlist(Results_Mammals_5)
Results_Mammals_5$Class <- "Mammals"
write.csv(Results_Mammals_5, "../Results/CENFA_summary_dataframes/Mammals_5km.csv", row.names = FALSE)

Amphibians_5 <- dir("../Results/CENFA/5_km/Amphibians/", recursive = TRUE, full.names = TRUE)
length(Amphibians_5)
Results_Amphibians_5 <- pblapply(as.list(Amphibians_5), Extract_results)
Results_Amphibians_5 <- data.table::rbindlist(Results_Amphibians_5)
Results_Amphibians_5$Class <- "Amphibians"
write.csv(Results_Amphibians_5, "../Results/CENFA_summary_dataframes/Amphibians_5km.csv", row.names = FALSE)

Birds_5 <- dir("../Results/CENFA/5_km/Birds/", recursive = TRUE, full.names = TRUE)
length(Birds_5)
Results_Birds_5 <- pblapply(as.list(Birds_5), Extract_results)
Results_Birds_5 <- data.table::rbindlist(Results_Birds_5)
Results_Birds_5$Class <- "Birds"
write.csv(Results_Birds_5, "../Results/CENFA_summary_dataframes/Birds_5km.csv", row.names = FALSE)

Reptiles_5 <- dir("../Results/CENFA/5_km/Reptiles/", recursive = TRUE, full.names = TRUE)
length(Reptiles_5)
Results_Reptiles_5<- pblapply(as.list(Reptiles_5), Extract_results)
Results_Reptiles_5 <- data.table::rbindlist(Results_Reptiles_5)
Results_Reptiles_5$Class <- "Reptiles"
write.csv(Results_Reptiles_5, "../Results/CENFA_summary_dataframes/Reptiles_5km.csv", row.names = FALSE)


Results_Amphibians_5 <- "../Results/CENFA_summary_dataframes/Amphibians_5km.csv" %>%  read.csv()
Results_Birds_5 <- "../Results/CENFA_summary_dataframes/Birds_5km.csv" %>%  read.csv()
Results_Mammals_5 <- "../Results/CENFA_summary_dataframes/Mammals_5km.csv" %>%  read.csv()
Results_Reptiles_5 <- "../Results/CENFA_summary_dataframes/Reptiles_5km.csv" %>%  read.csv()

Vertebrates_5 <- rbind(Results_Amphibians_5, Results_Birds_5, Results_Mammals_5, Results_Reptiles_5)
write.csv(Vertebrates_5, "../Results/CENFA_summary_dataframes/Vertebrates_5km.csv", row.names = FALSE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# 
# ## getting sensitivity maps (spatially explicit) from CENFA analysis
# library(CENFA)
# 
# Function_sensitivity_maps <- function(path) {
#   
#   # get species and paths for saving
#   SplitPath <-  unlist(strsplit(path, "/")) 
#   Species <- SplitPath[c((length(SplitPath)-1):length(SplitPath))]
#   Species <- paste(Species, collapse = "/")
#   
#   # read cenfa saved output, computed at 5km
#   cenfafile <- readRDS(path)
#   
#   # compute and save sensitivity map (spatially explicit) for that species
#   s.map <- sensitivity_map(cenfafile)
#   # pathsave_sensitivitymap <- paste0("../Results/1.CENFA/sensitivity_maps_species/", Species)
#   # saveRDS(s.map, pathsave_sensitivitymap)
#   return(s.map)
# }
# 
# Mammals_paths <- dir("../Results/1.CENFA/5_km/Mammals/", recursive = TRUE, full.names = TRUE)
# length(Mammals_paths)
# 
# S <- Sys.time()
# test <- Function_sensitivity_maps(Mammals_paths[1])
# E <- Sys.time()
# print(S-E)
# 
# # parallelise
# nCores <- parallel::detectCores()-4
# gc()
# cl <- parallel::makeCluster(nCores)
# parallel::clusterEvalQ(cl, {library(CENFA)})
# parallel::clusterExport(cl = cl,
#                         varlist = c("Mammals_paths", "Function_sensitivity_maps"),
#                         envir = environment())
# ParallelLogger::clusterApply(cluster = cl,
#                              fun = Function_sensitivity_maps,
#                              x = as.list(Mammals_paths[2:5]),
#                              progressBar = TRUE)
# stopCluster(cl)
# 
# 

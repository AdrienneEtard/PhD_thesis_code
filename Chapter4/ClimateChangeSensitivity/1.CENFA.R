## Implementing CENFA (Climate & Ecological Niche Factor Analysis) to estimate species' sensitivity to CC

   ## allows an estimation of sensitivity to climate change that combines 
        ## a marginality metric
        ## a specialisation metric
  ## those metrics are based on climatic niche properties for a species - compared to the overall, global climate niche (or that of an area of reference) 


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

setwd("E:/PhD/PhD_R_work/Climatic_niche_space/Code/")
behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

# Index to species name dataframe
Index <- read.csv("../../Range_maps_work/Results/7_GenerateSpeciesIndices/SpeciesIndices.csv")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## load climate variables, at 5 arc-minute resolution (~8.3kmsq at the equator), and select those that are not too collinear
WCfiles <- paste("../Data/wc2.1_5m_bio", list.files(path = "../Data/wc2.1_5m_bio"), sep="/")
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

# load distribution for a species -- at 10km resolution, in Behrman projection 
Mammals_ex <- raster("../../Range_maps_work/Results/11_AggregateRasters/Mammals0/sp17890")
crs(Mammals_ex) <- behrCRS

# realign projections of climate rasters to that of species distribution files (10 km square, Behrman)
Climatevars_Behrman_10 <- projectRaster(from=Climatevars, to=Mammals_ex, method = "bilinear")
crs(Climatevars_Behrman_10) <- behrCRS

# select climatic variables that are not too correlated
Selection <- removeCollinearity(Climatevars_Behrman_10, 
                                plot = TRUE, 
                                sample.points = FALSE,
                                method = "spearman", 
                                multicollinearity.cutoff = 0.65)

Selected_Climvars_10 <- Climatevars_Behrman_10[[c("Max_temperature_warmest_month", 
                                                "Mean_diurnal_range",
                                                "Annual_mean_temperature",
                                                "Precipitation_seasonality",
                                                "Annual_precipitation",
                                                "Precipitation_coldest_quarter")]]

plot(Selected_Climvars_10)

# # also chack correlation coefficients
# cors <- layerStats(Climatevars_Behrman_10[[c("Max_temperature_warmest_month",
#                                              "Temperature_seasonality",
#                                              "Mean_diurnal_range",
#                                              "Annual_mean_temperature",
#                                              "Precipitation_seasonality",
#                                              "Annual_precipitation",
#                                              "Precipitation_coldest_quarter")]],
#                    'pearson', na.rm=TRUE)
# corr_matrix <- cors$'pearson correlation coefficient'

# aggregate at 50 km resolution for a calculation of CENFA also at 50 km
Selected_Climvars_50 <- aggregate(Selected_Climvars_10, fact=5)

# Climatevars_Behrman_50 <- aggregate(Climatevars_Behrman_10, fact=5)
# Selection50 <- removeCollinearity(Climatevars_Behrman_50, plot = TRUE, sample.points = FALSE)

# Calculate global covariance matrix at 10 and 50 km

Global_cov_10 <- GLcenfa(Selected_Climvars_10,
                      progress=TRUE, 
                      center=TRUE,
                      scale=TRUE, 
                      filename="../Results/CENFA/GlobalClimCov6",
                      overwrite=TRUE)

Global_cov_50 <- GLcenfa(Selected_Climvars_50,
                        progress=TRUE, 
                        center=TRUE,
                        scale=TRUE, 
                        filename="../Results/CENFA/GlobalClimCov10",
                        overwrite=TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
  }
}

Run_CENFA <- function(GlobalCovMatrix, ListPaths, IndexDF, resolution, Class) {
  for(i in 1:length(ListPaths)) {
    CENFA_analysis(GlobalCovMatrix, X=ListPaths[i], IndexDF, resolution, Class)
    print(paste(i, "over", length(ListPaths)))
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## List paths to all species distribution files
files <- dir("../../Range_maps_work/Results/11_AggregateRasters", recursive = TRUE, full.names = TRUE)
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
# Distrib <- raster(filesReptiles[391])
# crs(Distrib) <- behrCRS
# cnfa(x = Global_cov_10,
#      s.dat = Distrib,
#      progress = TRUE,
#      scale = FALSE)

## running CENFA at 50 km square resolution
Run_CENFA(Global_cov_50, ListPaths = filesMammals, IndexDF = Index, resolution = 50, Class = "Mammals")
Run_CENFA(Global_cov_50, ListPaths = filesAmphibians, IndexDF = Index, resolution = 50, Class = "Amphibians")
Run_CENFA(Global_cov_50, ListPaths = filesBirds, IndexDF = Index, resolution = 50, Class = "Birds")
Run_CENFA(Global_cov_50, ListPaths = filesReptiles, IndexDF = Index, resolution = 50, Class = "Reptiles")

## running CENFA at 10 km square resolution
Run_CENFA(Global_cov_10, ListPaths = filesMammals, IndexDF = Index, resolution = 10, Class = "Mammals")
Run_CENFA(Global_cov_10, ListPaths = filesAmphibians, IndexDF = Index, resolution = 10, Class = "Amphibians")
Run_CENFA(Global_cov_10, ListPaths = filesBirds, IndexDF = Index, resolution = 10, Class = "Birds")
Run_CENFA(Global_cov_10, ListPaths = filesReptiles[1087:length(filesReptiles)], IndexDF = Index, resolution = 10, Class = "Reptiles")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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

# results at 50 km square resolution

Mammals_50 <- dir("../Results/CENFA/50_km/Mammals/", recursive = TRUE, full.names = TRUE)
length(Mammals_50)
Results_Mammals_50 <- pblapply(as.list(Mammals_50), Extract_results)
Results_Mammals_50 <- data.table::rbindlist(Results_Mammals_50)
Results_Mammals_50$Class <- "Mammals"

Amphibians_50 <- dir("../Results/CENFA/50_km/Amphibians/", recursive = TRUE, full.names = TRUE)
length(Amphibians_50)
Results_Amphibians_50 <- pblapply(as.list(Amphibians_50), Extract_results)
Results_Amphibians_50 <- data.table::rbindlist(Results_Amphibians_50)
Results_Amphibians_50$Class <- "Amphibians"

Birds_50 <- dir("../Results/CENFA/50_km/Birds/", recursive = TRUE, full.names = TRUE)
length(Birds_50)
Results_Birds_50 <- pblapply(as.list(Birds_50), Extract_results)
Results_Birds_50 <- data.table::rbindlist(Results_Birds_50)
Results_Birds_50$Class <- "Birds"

Reptiles_50 <- dir("../Results/CENFA/50_km/Reptiles/", recursive = TRUE, full.names = TRUE)
length(Reptiles_50)
Results_Reptiles_50 <- pblapply(as.list(Reptiles_50), Extract_results)
Results_Reptiles_50 <- data.table::rbindlist(Results_Reptiles_50)
Results_Reptiles_50$Class <- "Reptiles"

Vertebrates_50 <- rbind(Results_Amphibians_50, Results_Birds_50, Results_Mammals_50, Results_Reptiles_50)
write.csv(Vertebrates_50, "../Results/CENFA_summary_dataframes/Vertebrates_50km.csv", row.names = FALSE)

# results at 10 km square resolution

Mammals_10 <- dir("../Results/CENFA/10_km/Mammals/", recursive = TRUE, full.names = TRUE)
length(Mammals_10)
Results_Mammals_10 <- pblapply(as.list(Mammals_10), Extract_results)
Results_Mammals_10 <- data.table::rbindlist(Results_Mammals_10)
Results_Mammals_10$Class <- "Mammals"

Amphibians_10 <- dir("../Results/CENFA/10_km/Amphibians/", recursive = TRUE, full.names = TRUE)
length(Amphibians_10)
Results_Amphibians_10 <- pblapply(as.list(Amphibians_10), Extract_results)
Results_Amphibians_10 <- data.table::rbindlist(Results_Amphibians_10)
Results_Amphibians_10$Class <- "Amphibians"

Birds_10 <- dir("../Results/CENFA/10_km/Birds/", recursive = TRUE, full.names = TRUE)
length(Birds_10)
Results_Birds_10 <- pblapply(as.list(Birds_10), Extract_results)
Results_Birds_10 <- data.table::rbindlist(Results_Birds_10)
Results_Birds_10$Class <- "Birds"

Reptiles_10 <- dir("../Results/CENFA/10_km/Reptiles/", recursive = TRUE, full.names = TRUE)
length(Reptiles_10)
Results_Reptiles_10 <- pblapply(as.list(Reptiles_10), Extract_results)
Results_Reptiles_10 <- data.table::rbindlist(Results_Reptiles_10)
Results_Reptiles_10$Class <- "Reptiles"

Vertebrates_10 <- rbind(Results_Amphibians_10, Results_Birds_10, Results_Mammals_10, Results_Reptiles_10)
write.csv(Vertebrates_10, "../Results/CENFA_summary_dataframes/Vertebrates_10km.csv", row.names = FALSE)


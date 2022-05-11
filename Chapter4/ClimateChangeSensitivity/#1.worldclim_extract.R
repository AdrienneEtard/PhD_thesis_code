# Species distribution rasters are projected with Behrman projection, at 1km square resolution.
# WorldClim data rasters have a Long Lat projection with a resolution of 30 arc-second res (~0.86 km sq at the Equator)

# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month

# Before extracting range-wide WorldClim variables for a species, need to reproject WorldClim into Behrman equal area projection

library(raster)
library(dplyr)

# The species distributions rasters have Behrman equal area projection  
# behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
# crs(spdist) <- behrCRS


## Function to get range-wide Worldclim variables and some summary stats for a species
Extract_rangewide_WorldClim <- function(sp_raster_path, WorldClim_raster) {
  
  # load raster file for species distribution   
  sp <- raster(sp_raster_path)
  behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
  crs(sp) <- behrCRS
  
  # align projections
  WClim_reproj <- projectRaster(from=WorldClim_raster, to=sp, method = "bilinear")
  
  # get range-wide WorldClim values for the species
  Inter <- mask(x=WClim_reproj, mask = sp)
  
  # store results 
  # species index and name
  x <- strsplit(sp_raster_path, "/") %>%
    unlist()
  spIndex <- x[length(x)]
  spName <- Index$Binomial[Index$Index==sub("sp", spIndex, replacement = "")]

  ResultsDF <- data.frame(SpeciesIndex=spIndex,
                          SpeciesName=spName,
                          Min=minValue(Inter),
                          Max=maxValue(Inter),
                          Mean=cellStats(Inter, stat='mean', na.rm=TRUE, asSample=NULL),
                          Sd=cellStats(Inter, stat='sd', na.rm=TRUE, asSample=FALSE),
                          Median=median(values(Inter)[!is.na(values(Inter))]))


  return(ResultsDF)
  
}


####################################################################################################

## Load WorldClim rasters
WClim_5 <- raster("../Data/wc2.1_30s_bio/wc2.1_30s_bio_5.tif")

## List paths to all species distribution files
files <- dir("../../Range_maps_work/Results/10_CutByElevationalRanges", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("info", files)]
files <- lapply(files, FUN = function(x) {
  x <- strsplit(x, "/") %>%
    unlist()
  return(paste0(x[1:(length(x)-1)], collapse = "/"))
})
files <- unique(unlist(files))
files <- files[grepl("sp", files)]

## species index
Index <- read.csv("../../Range_maps_work/Results/7_GenerateSpeciesIndices/SpeciesIndices.csv")


####################################################################################################

## Applying the function for BIO5 = Max Temperature of Warmest Month

filesAmphibians <- files[grepl("Amphibians", files)]

Iterate_over_species <- function(fileslist, WCraster) {
  
  Results <- list()
  
  for(i in 1:length(fileslist)) {
  
  Results[[i]] <- Extract_rangewide_WorldClim(fileslist[i], WCraster)
  print(paste("species",i, "on", length(fileslist)))
}

  return(Results)
  
} 

memory.limit(size=5000000)
start <- Sys.time()
Bio5_Amphibians <- Iterate_over_species(filesAmphibians, WClim_5)
end <- Sys.time()
print(end-start)
Bio5_Amphibians <- data.table::rbindlist(Bio5_Amphibians) %>% 
  as.data.frame()


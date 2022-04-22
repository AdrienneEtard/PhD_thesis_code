## Mapping trait completeness

# # # # # # # # # # # # # 
#       preamble        #
# # # # # # # # # # # # #

## things to consider: - masks and rasters in different projections ; 

## NB: All raster files are ALREADY projected using Behrman but crs appears NA.
 # behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

library(raster)
library(RColorBrewer)
library(stringr)
library(parallel)
library(ParallelLogger)
source("Functions_map_completeness.R")

# # # # # Set paths to range map files # # # # # # 

# Relevant species:
SpeciesAfter <- read.csv("../../Data/Distribution_maps/SpeciesAfterCuts_List.csv")

rangesDir <- "D:/Range_maps_work/Results/11_AggregateRasters/"

rangeMapFiles <- paste(rangesDir, dir(path = rangesDir, recursive = TRUE), sep="")
rangeMapFiles <- rangeMapFiles[!grepl("aux",rangeMapFiles)]
rangeMapFiles <- rangeMapFiles[!grepl("info",rangeMapFiles)]
rangeMapFiles <- rangeMapFiles[!grepl("log",rangeMapFiles)]
rangeMapFiles <- rangeMapFiles[!grepl("SppDidntWork",rangeMapFiles)]
rangeMapFiles <- lapply(rangeMapFiles, function(x){
  y <- strsplit(x, "/") %>%  unlist
  y <- y[1:length(y)-1]
  return(paste0(y, collapse = "/"))
  })
rangeMapFiles <- unlist(rangeMapFiles) %>%
  unique %>% 
  as.data.frame %>% 
  setNames(., "Path")

rangeMapFiles$Index <- apply(rangeMapFiles, 1, function(x) {
  y <- strsplit(x, "/") %>%  unlist
  return(y[length(y)])})
  
## Filter for relevant species 

## after cuts by elevational limits
rangeMaps_pathsAfter <- rangeMapFiles$Path[rangeMapFiles$Index %in% SpeciesAfter$IndexSp] %>% 
  as.character()

amphibianFilesAfter <- rangeMaps_pathsAfter[grep("Amphibians",rangeMaps_pathsAfter)]
birdFilesAfter <- rangeMaps_pathsAfter[grep("Birds",rangeMaps_pathsAfter)]
reptileFilesAfter <- rangeMaps_pathsAfter[grep("Reptiles",rangeMaps_pathsAfter)]
mammalFilesAfter <- rangeMaps_pathsAfter[grep("Mammals",rangeMaps_pathsAfter)]


# all paths to range files (after cutting by elevation limits)
allFilesAfter <- list(amphibianFilesAfter,birdFilesAfter,mammalFilesAfter,reptileFilesAfter)
saveRDS(allFilesAfter, "../../Results/SpatialAnalyses/Paths_to_distribution_maps/PathsToMapsAfter.rds")

rm(SpeciesAfter, 
   allFilesAfter,
   amphibianFilesAfter,
   birdFilesAfter,
   mammalFilesAfter, 
   reptileFilesAfter,
   rangeMaps_pathsAfter, 
   rangesDir, 
   rangeMapFiles)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

allFiles <- readRDS("../../Results/SpatialAnalyses/Paths_to_distribution_maps/PathsToMapsAfter.rds")
names(allFiles) <- c("Amphibians", "Birds", "Mammals", "Reptiles")
lapply(allFiles, length)

mask <- raster("../../Data/Distribution_maps/mask10k")
values(mask)[!is.na(values(mask))] <- 1

nCores <- parallel::detectCores()-2

# need to increase memory size otherwise the script will crash
gc()
memory.limit(900000)

# load range sizes files
RangeSizes <- read.csv("../../Results/Range_sizes/AfterCuts.csv")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                     
#
#       Generate species richness rasters at 10 km resolution for each class              #      
#                                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


## to run to calculate species richness in each grid cell (stacking species range files), to generate both SR by class and TotalRichness

# Birds species richness
Birds_SR <- ApplyToChunks(rangeMapFiles=allFiles[["Birds"]], mask=mask, nCores=nCores, ChunkSize = 500)
saveRDS(Birds_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Birds_10k.rds")

# Amphibians species richness
Amphibians_SR <- ApplyToChunks(rangeMapFiles=allFiles[["Amphibians"]], mask=mask, nCores=nCores, ChunkSize = 500)
values(Amphibians_SR)[values(Amphibians_SR)>10e7] <- NA 
saveRDS(Amphibians_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Amphibians_10k.rds")

# Reptiles species richness
Reptiles_SR <- ApplyToChunks(rangeMapFiles=allFiles[["Reptiles"]], mask=mask, nCores=nCores, ChunkSize = 500)
saveRDS(Reptiles_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Reptiles_10k.rds")

# Mammals species richness
Mammals_SR <- ApplyToChunks(rangeMapFiles=allFiles[["Mammals"]], mask=mask, nCores=nCores, ChunkSize = 500)
saveRDS(Mammals_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Mammals_10k.rds")
plot(Mammals_SR)

## to run to calculate RANGE-WEIGHTED species richness in each grid cell (1/range size of the species)
## to generate both range-weigthed SR by class and TotalRichness -- 10 km resolution

# Birds species richness
Birds_SR <- ApplyToChunks_rangeweighted(rangeMapFiles=allFiles[["Birds"]], mask=mask, nCores=nCores, ChunkSize = 500, rangeSizes = RangeSizes)
saveRDS(Birds_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Birds_10k_range_weighted.rds")
plot(log(Birds_SR))

gc()
# Reptiles species richness
Reptiles_SR <- ApplyToChunks_rangeweighted(rangeMapFiles=allFiles[["Reptiles"]], mask=mask, nCores=nCores, ChunkSize = 500, rangeSizes = RangeSizes)
saveRDS(Reptiles_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Reptiles_10k_range_weighted.rds")
plot(log(Reptiles_SR))
gc()
# Mammals species richness
Mammals_SR <- ApplyToChunks_rangeweighted(rangeMapFiles=allFiles[["Mammals"]], mask=mask, nCores=nCores, ChunkSize = 500, rangeSizes = RangeSizes)
saveRDS(Mammals_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Mammals_10k_range_weighted.rds")
plot(log(Mammals_SR))

gc()
# Amphibians species richness
Amphibians_SR <- ApplyToChunks_rangeweighted(rangeMapFiles=allFiles[["Amphibians"]], mask=mask, nCores=nCores, ChunkSize = 500, rangeSizes = RangeSizes)
plot(log(Amphibians_SR))
#values(Amphibians_SR)[values(Amphibians_SR)>10e7] <- NA 
saveRDS(Amphibians_SR, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Amphibians_10k_range_weighted.rds")


## checking
Birds_SR <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Birds_10k_range_weighted.rds")
Reptiles_SR <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Reptiles_10k_range_weighted.rds")
Mammals_SR <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Mammals_10k_range_weighted.rds")
Amphibians_SR <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Amphibians_10k_range_weighted.rds")

hist(values(log(Birds_SR)))
hist(values(log(Amphibians_SR)))
hist(values(log(Mammals_SR)))
hist(values(log(Reptiles_SR)))

hist(values(Birds_SR))
max(values(Birds_SR))

hist(values(Mammals_SR))
max(values(Mammals_SR))

hist(values(Birds_SR))
max(values(Birds_SR))

hist(values(Amphibians_SR))
max(values(Amphibians_SR), na.rm=T)

hist(values(Reptiles_SR))
max(values(Reptiles_SR))

values(Amphibians_SR)[values(Amphibians_SR)>0.5] <- NA 
values(Mammals_SR)[values(Mammals_SR)>0.5] <- NA 
values(Reptiles_SR)[values(Reptiles_SR)>0.5] <- NA 
values(Birds_SR)[values(Birds_SR)>0.5] <- NA 

hist(values(log(Amphibians_SR)))
hist(values(log(Mammals_SR)))
hist(values(log(Reptiles_SR)))
hist(values(log(Birds_SR)))

max(values(log(Birds_SR)))
values(log(Birds_SR))[values(log(Birds_SR))>0.5]

max(values(log(Birds_SR))[values(log(Birds_SR))<0.5])
max(values(log(Amphibians_SR))[values(log(Amphibians_SR))<0.5])
max(values(log(Reptiles_SR))[values(log(Reptiles_SR))<0.5])

max(values(log(Amphibians_SR)))
values(log(Amphibians_SR))[values(log(Amphibians_SR))>0.5]

# max(values(log(Mammals_SR)))
# values(log(Mammals_SR))[values(log(Mammals_SR))>0.5]

max(values(log(Reptiles_SR)))
values(log(Reptiles_SR))[values(log(Reptiles_SR))>0.5]



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                                         #
#       Calculate mean trait completeness by class 10 x 10 km resolution                  #      
#                                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

trait.completeness <- readRDS("../../Results/Traits_to_map/traits_tomap_SpeciesAfterCuts.rds")

Birds_mean_comp <- Mean_ApplyToChunks(rangeMapFiles=allFiles$Birds, mask, nCores, ChunkSize=500, trait.completeness[["Birds"]])
saveRDS(Birds_mean_comp, "../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Birds.rds")

Reptiles_mean_comp <- Mean_ApplyToChunks(rangeMapFiles=allFiles$Reptiles, mask, nCores, ChunkSize=500, trait.completeness[["Reptiles"]])
saveRDS(Reptiles_mean_comp, "../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Reptiles.rds")

Amphibians_mean_comp <- Mean_ApplyToChunks(rangeMapFiles=allFiles$Amphibians, mask, nCores, ChunkSize=500, trait.completeness[["Amphibians"]])
saveRDS(Amphibians_mean_comp, "../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Amphibians.rds")

Mammals_mean_comp <- Mean_ApplyToChunks(rangeMapFiles=allFiles$Mammals, mask, nCores, ChunkSize=500, trait.completeness[["Mammals"]])
saveRDS(Mammals_mean_comp, "../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Mammals.rds")
plot(Mammals_mean_comp)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                                         #
#       Generate species richness rasters at 50 km resolution                             #      
#                                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

Birds_SR_50km <- SumSpecies_50km(mask = mask, rangeMapFiles = allFiles[["Birds"]], nCores = 5)
saveRDS(Birds_SR_50km, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Birds_50k.rds")

Reptiles_SR_50km <- SumSpecies_50km(mask = mask, rangeMapFiles = allFiles[["Reptiles"]], nCores = 5)
saveRDS(Reptiles_SR_50km, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Reptiles_50k.rds")

Mammals_SR_50km <- SumSpecies_50km(mask = mask, rangeMapFiles = allFiles[["Mammals"]], nCores = 5)
saveRDS(Mammals_SR_50km, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Mammals_50k.rds")

Amphibians_SR_50km <- SumSpecies_50km(mask = mask, rangeMapFiles = allFiles[["Amphibians"]], nCores = 5)
saveRDS(Amphibians_SR_50km, "../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Amphibians_50k.rds")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                                         #
#  Calculate mean, variance and median trait completeness by class 50 x 50 km resolution  #      
#                                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# NB - for variance and median, 10km resolution is too computationally intensive
# because not possible to work on chunks of rasters as median and variance are not distributive

nCores <- parallel::detectCores()-2

Amphibians_MedianVar_comp <- Median_Var_completeness_50k(rangeMapFiles=allFiles$Amphibians, mask, nCores, trait.completeness[["Amphibians"]])
saveRDS(Amphibians_MedianVar_comp$Variance, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceAmphibians.rds")
saveRDS(Amphibians_MedianVar_comp$Median, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianAmphibians.rds")
saveRDS(Amphibians_MedianVar_comp$Mean, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanAmphibians.rds")
rm(Amphibians_MedianVar_comp)

Mammals_MedianVar_comp <- Median_Var_completeness_50k(rangeMapFiles=allFiles$Mammals, mask, nCores, trait.completeness[["Mammals"]])
saveRDS(Mammals_MedianVar_comp$Variance, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceMammals.rds")
saveRDS(Mammals_MedianVar_comp$Median, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianMammals.rds")
saveRDS(Mammals_MedianVar_comp$Mean, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanMammals.rds")
rm(Mammals_MedianVar_comp)

Birds_MedianVar_comp <- Median_Var_completeness_50k(rangeMapFiles=allFiles$Birds, mask, nCores, trait.completeness[["Birds"]])
saveRDS(Birds_MedianVar_comp$Variance, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceBirds.rds")
saveRDS(Birds_MedianVar_comp$Median, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianBirds.rds")
saveRDS(Birds_MedianVar_comp$Mean, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanBirds.rds")

Reptiles_MedianVar_comp <- Median_Var_completeness_50k(rangeMapFiles=allFiles$Reptiles, mask, nCores, trait.completeness[["Reptiles"]])
saveRDS(Reptiles_MedianVar_comp$Variance, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceReptiles.rds")
saveRDS(Reptiles_MedianVar_comp$Median, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianReptiles.rds")
saveRDS(Reptiles_MedianVar_comp$Mean, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanReptiles.rds")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                                                                         #
#         Calculate standard deviation from variance rasters at 50 x 50 km resolution     #      
#                                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


var_completeness_Birds <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceBirds.rds")
var_completeness_Mammals <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceMammals.rds")
var_completeness_Reptiles <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceReptiles.rds")
var_completeness_Amphibians <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/VarianceAmphibians.rds")

GetSD <- function(r_variance) {
  
  NAVal <- !is.na(values(r_variance))
  Minus1val <- !values(r_variance)==-1
  
  r_variance[NAVal & Minus1val] <- sqrt(r_variance[NAVal & Minus1val])
  return(r_variance)
  
}

sd_completeness_Birds <- GetSD(var_completeness_Birds)
sd_completeness_Mammals <- GetSD(var_completeness_Mammals)
sd_completeness_Reptiles <- GetSD(var_completeness_Reptiles)
sd_completeness_Amphibians <- GetSD(var_completeness_Amphibians)


saveRDS(sd_completeness_Birds, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Birds.rds")
saveRDS(sd_completeness_Mammals, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Mammals.rds")
saveRDS(sd_completeness_Reptiles, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Reptiles.rds")
saveRDS(sd_completeness_Amphibians, "../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Amphibians.rds")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#                                    PLOTS                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

rm(list = ls())


PlotMaps <- function(r_mammals, r_birds, r_amphibians, r_reptiles, filename, metric) {
  
  filename <- paste0("../../Results/Maps_completeness/", filename)
  
  pdf(filename, width=6, height=10, family="Times", pointsize=12)
  
  par(family='serif', tcl=0.2, cex.lab=1.3, mgp=c(1.5,0.2,0), cex.main=1.8, mfrow=c(4,1), oma=c(1,1,1,13), mar=c(1.5,0,1.5,0))

  brks <- c(-1,80,90,91,92,93,94,95,96,97,98,99,100)
  cols <- c("#cccccc",brewer.pal(n = 11,name = "RdYlBu"))
  
  plot(r_mammals, col=cols, breaks=brks, axes=F, legend=FALSE)
  title("(a) Mammals", adj=0, line=0.5, cex=8)
  plot(r_birds, col=cols, breaks=brks, axes=F, legend=FALSE)
  title("(b) Birds", adj=0, line=0.5, cex=8)
  
  par(xpd=NA)
  legend(1.9e7,2e7,c(
    "NA","80 - 90%", "90 - 91%","91 - 92%","92 - 93%","93 - 94%","94 - 95%",
    "95 - 96%","96 - 97%","97 - 98%",
    "98 - 99%","99 - 100%"),fill = cols,cex=1.5, 
    bty="n", title = paste(metric, "completeness \n(birds and mammals):"))
  
  brks <- c(-1,0,10,20,30,40,50,60,70,80,90,100)
  cols <- c("#cccccc",brewer.pal(n = 10,name = "RdYlBu"))
  
  plot(r_amphibians, col=cols, breaks=brks, axes=F, legend=FALSE)
  title("(c) Amphibians", adj=0, line=0.5, cex=8)
  plot(r_reptiles, col=cols, breaks=brks, axes=F, legend=FALSE)
  title("(d) Reptiles", adj=0, line=0.5, cex=8)
  
  par(xpd=NA)
  legend(1.9e7,2e7,c(
    "NA","0 - 10%","10 - 20%","20 - 30%","30 - 40%",
    "40 - 50%","50 - 60%","60 - 70%",
    "70 - 80%","80 - 90%","90 - 100%"),fill = cols,cex=1.5, 
    bty="n", title = paste(metric, "completeness \n(herptiles):"))
  
  dev.off()
  
}

Plot2Maps <- function(r_amphibians, r_reptiles, filename, metric) {
  
  filename <- paste0("../../Results/Maps_completeness/", filename)
  
  pdf(filename, width=10, height=7, family="Times", pointsize=12)
  
  par(family='serif', tcl=0.2, cex.lab=1.3, mgp=c(1.5,0.2,0), cex.main=1.8, mfrow=c(2,1), oma=c(1,1,1,13), mar=c(1.5,0,1.5,0))
  
  brks <- c(-1,0,10,20,30,40,50,60,70,80,90,100)
  cols <- c("#cccccc",brewer.pal(n = 10,name = "RdYlBu"))
  
  plot(r_amphibians, col=cols, breaks=brks, axes=F, legend=FALSE)
  title("(a) Amphibians", adj=0, line=0.5, cex=5)
  plot(r_reptiles, col=cols, breaks=brks, axes=F, legend=FALSE)
  title("(b) Reptiles", adj=0, line=0.5, cex=5)
  
  par(xpd=NA)
  legend(1.9e7,2e7,c(
    "NA","0 - 10%","10 - 20%","20 - 30%","30 - 40%",
    "40 - 50%","50 - 60%","60 - 70%",
    "70 - 80%","80 - 90%","90 - 100%"),fill = cols,cex=1.5,
    bty="n", title = paste(metric, "completeness:"))
  
  dev.off()
  
}



# 4 maps - mean completeness, 10km and 50km resolution

mean_completeness_Birds <- readRDS("../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Birds.rds")
mean_completeness_Mammals <- readRDS("../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Mammals.rds")
mean_completeness_Reptiles <- readRDS("../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Reptiles.rds")
mean_completeness_Amphibians <- readRDS("../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Amphibians.rds")

mean_completeness_Birds50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanBirds.rds")
mean_completeness_Mammals50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanMammals.rds")
mean_completeness_Reptiles50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanReptiles.rds")
mean_completeness_Amphibians50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanAmphibians.rds")

PlotMaps(mean_completeness_Mammals, 
         mean_completeness_Birds, 
         mean_completeness_Amphibians,
         mean_completeness_Reptiles, 
         filename = "Mean_map.pdf", "Mean")

PlotMaps(mean_completeness_Mammals50, 
         mean_completeness_Birds50, 
         mean_completeness_Amphibians50,
         mean_completeness_Reptiles50, 
         filename = "Mean_map_50k.pdf", "Mean")

# 4 maps median completeness

median_completeness_Birds <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianBirds.rds")
median_completeness_Mammals <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianMammals.rds")
median_completeness_Reptiles <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianReptiles.rds")
median_completeness_Amphibians <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianAmphibians.rds")

PlotMaps(median_completeness_Mammals, 
         median_completeness_Birds, 
         median_completeness_Amphibians,
         median_completeness_Reptiles, 
         filename = "Median_map.pdf", "Median")


# 2 maps median completeness (herptiles only)
Plot2Maps(median_completeness_Amphibians,
          median_completeness_Reptiles, 
          filename = "Median_2mapsV2.pdf", "Median")



# variance -- get standard deviation from variance and plot standard deviation

sd_completeness_Birds <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Birds.rds")
sd_completeness_Reptiles <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Reptiles.rds")
sd_completeness_Mammals <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Mammals.rds")
sd_completeness_Amphibians <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/SD_Amphibians.rds")

hist(sd_completeness_Birds)
hist(sd_completeness_Reptiles)
hist(sd_completeness_Mammals)
hist(sd_completeness_Amphibians)

## plotting standard deviation around the mean

pdf("../../Results/Maps_completeness/sd_map.pdf", width=6, height=10, family="Times", pointsize=12)

par(family='serif', tcl=0.2, cex.lab=1.3, mgp=c(1.5,0.2,0), cex.main=1.8, mfrow=c(4,1), oma=c(1,1,1,13), mar=c(1.5,0,1.5,0))

brks <- c(-1, 0, 5, 10, 15 ,20, 25, 30)
cols <- c("#cccccc",brewer.pal(n = 6,name = "RdYlBu"))

plot(sd_completeness_Mammals, col=cols, breaks=brks, axes=F, legend=FALSE)
title("(a) Mammals", adj=0, line=0.5, cex=8)
plot(sd_completeness_Birds, col=cols, breaks=brks, axes=F, legend=FALSE)
title("(b) Birds", adj=0, line=0.5, cex=8)

par(xpd=NA)
legend(1.9e7, 2e7,c(
  "NA","0 - 5%","5 - 10%","10 - 15%","15 - 20%",
  "20 - 25%","25 - 30%"),fill = cols, cex=1.5, 
  bty="n", title = "Standard deviation \nof completeness \n(birds and mammals):")

brks <- c(-1, 0, 10, 15 ,20, 25, 30, 35, 40, 50, 60, 75)
cols <- c("#cccccc", brewer.pal(n = 10,name = "RdYlBu"))

plot(sd_completeness_Amphibians, col=cols, breaks=brks, axes=F, legend=FALSE)
title("(c) Amphibians", adj=0, line=0.5, cex=8)
plot(sd_completeness_Reptiles, col=cols, breaks=brks, axes=F, legend=FALSE)
title("(d) Reptiles", adj=0, line=0.5, cex=8)

par(xpd=NA)
legend(1.9e7,2e7,c(
  "NA","0- 10%","10 - 15%","15 - 20%","20 - 25%",
  "25 - 30%","30 - 35%","35 - 40%",
  "40 - 50%","50 - 60%","60 - 75%"),fill = cols,cex=1.5, 
  bty="n", title = "Standard deviation \nof completeness \n(herptiles):")

dev.off()


## plot correlation between species richness and standard deviation
SRMammals50 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Mammals_50k.rds")
SRMammals10 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Mammals_10k.rds")

SRBirds50 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Birds_50k.rds")
SRBirds10 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Birds_10k.rds")

SRAmphibians50 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Amphibians_50k.rds")
SRAmphibians10 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Amphibians_10k.rds")
values(SRAmphibians10)[values(SRAmphibians10>1000)] <- NA

SReptiles50 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Reptiles_50k.rds")
SReptiles10 <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Reptiles_10k.rds")

Metrics_VS_SR <- function(r_mammals, r_birds, r_amphibians, r_reptiles, filename, metric, res, pdf) {
  
  filename <- paste0("../../Results/Maps_completeness/", filename)
  
  if(pdf){
    pdf(filename, width=8, height=6, family="Times", pointsize=12)
  }
  
  else {
    png(filename, width=500, height=400, family="Times", pointsize=12)
  }
  
  par(family='serif', tcl=0.2, cex.lab=1.3, mgp=c(1.5,0.2,0), cex.main=1.8, mfrow=c(2,2), oma=c(2,2,2,2), mar=c(3,3,3,3))
  
  if(res==50){
  
  plot(values(SRMammals50), values(r_mammals), pch=19, xlab="Species richness", ylab=metric)
  title("(a) Mammals", adj=0, line=0.5, cex=3)
  plot(values(SRBirds50), values(r_birds), pch=19, xlab="Species richness", ylab=metric)
  title("(b) Birds", adj=0, line=0.5, cex=3)
  plot(values(SRAmphibians50), values(r_amphibians), pch=19,  xlab="Species richness", ylab=metric)
  title("(c) Amphibians", adj=0, line=0.5, cex=3)
  plot(values(SReptiles50), values(r_reptiles), pch=19, xlab="Species richness", ylab=metric)
  title("(d) Reptiles", adj=0, line=0.5, cex=3)
  
  dev.off()
  }
  
  if(res==10){
    
    plot(values(SRMammals10), values(r_mammals), pch=19, xlab="Species richness", ylab=metric)
    title("(a) Mammals", adj=0, line=0.5, cex=3)
    plot(values(SRBirds10), values(r_birds), pch=19, xlab="Species richness", ylab=metric)
    title("(b) Birds", adj=0, line=0.5, cex=3)
    plot(values(SRAmphibians10), values(r_amphibians), pch=19,  xlab="Species richness", ylab=metric)
    title("(c) Amphibians", adj=0, line=0.5, cex=3)
    plot(values(SReptiles10), values(r_reptiles), pch=19, xlab="Species richness", ylab=metric)
    title("(d) Reptiles", adj=0, line=0.5, cex=3)
    
    dev.off()
  }
  
}

# plotting SR against standard deviation of the mean, at 50km resolution
Metrics_VS_SR(sd_completeness_Mammals, sd_completeness_Birds, sd_completeness_Amphibians, sd_completeness_Reptiles,
              "VS_SR/sd_SR.png", "SD completeness (%)", res=50, pdf=FALSE)

# plotting SR against median, at 50km resolution
Metrics_VS_SR(median_completeness_Mammals, median_completeness_Birds, median_completeness_Amphibians, median_completeness_Reptiles,
              "VS_SR/median_SR.png", "Median completeness (%)", res=50, pdf=FALSE)

# plotting SR vs mean (at 50x50 km resolution)
mean_completeness_Birds50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanBirds.rds")
mean_completeness_Mammals50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanMammals.rds")
mean_completeness_Reptiles50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanReptiles.rds")
mean_completeness_Amphibians50 <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanAmphibians.rds")

Metrics_VS_SR(mean_completeness_Mammals50, mean_completeness_Birds50, mean_completeness_Amphibians50, mean_completeness_Reptiles50,
              "VS_SR/mean_SR_50km.png", "Mean completeness (%)", res=50, pdf=FALSE)

# plotting SR vs mean (at 10x10 km resolution)
Metrics_VS_SR(mean_completeness_Mammals, mean_completeness_Birds, mean_completeness_Amphibians, mean_completeness_Reptiles,
              "VS_SR/mean_SR_10km.png", "Mean completeness (%)", res=10, pdf=FALSE)


# plottig mean vs median
plot(values(mean_completeness_Mammals50), values(median_completeness_Mammals), pch=19)
abline(a=0, b=1, col="red")
hist(mean_completeness_Mammals50)
hist(median_completeness_Mammals)





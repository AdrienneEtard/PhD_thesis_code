## Function to put the raster files together for a given class, and calculate species richness from stacking the raster files

SumSpecies <- function(mask, rangeMapFiles, nCores){

  # take all cells that have non NA values in "mask" (raster file of the world, in Behrman projection)
  mask.vals <- !is.na(values(mask))
  
  # verify extent and resolution
  verif <- raster::raster(rangeMapFiles[1])
  
  if(extent(mask)!=extent(verif)) {
    print("Different extent between species raster and mask")
    stop()
  }
  
  if((res(mask)[1]!=res(verif)[1]) |(res(mask)[2]!=res(verif)[2])) {
    print("Different resolution between species raster and mask")
    stop()
  }
  
  # else{print("Start processing.")}
  
  rm(verif)
  
  # # # # # # # # 
  
  # parallelise and apply clusterFunc to range files
  cl <- parallel::makeCluster(nCores)
  
  parallel::clusterExport(cl = cl,
                          varlist = c("rangeMapFiles",
                                      "mask.vals"),
                          envir = environment())

 
  
  # Define function "clusterFunc" that takes a raster file, then intersects it with mask where the species is present
  
  clusterFunc <- function(rangeMaps){
    
    # open range file
    r <- raster::raster(rangeMaps)
    
    # from the raster, take all cells that have non NA value (where the species is present)
    sp.pres <- !is.na(raster::values(r))
    
    # create grid with NAs
    sp.pres2 <- rep(NA,length(mask.vals))
    
    # put 1 where the species is present and intersects with mask
    sp.pres2[(mask.vals) & (sp.pres)] <- 1
    
    # return the grid with cells=1 where the species is present
    return(sp.pres2)
  }
  environment(clusterFunc) <- .GlobalEnv
  
  # Run clusterFunc in parallel for all range files
  # (returns a list of numeric elements for species presence in each range file )
  SpPres <- ParallelLogger::clusterApply(cl = cl,
                                         x = as.list(rangeMapFiles),
                                         fun = clusterFunc,
                                         progressBar=TRUE)
  
  stopCluster(cl)
  
  print("Stacking and summing across rasters.")

  # put results together from list to matrix (each column = data for 1 species)
  SpPres <- do.call('cbind', SpPres)
  # calculate species richness in each grid cell by summing on rows
  SpRich <- apply(SpPres, 1, sum, na.rm=TRUE)
  
  # Map species richness onto mask raster layer
  OutRaster <- mask
  values(OutRaster) <- SpRich
  
  return(OutRaster)

}

SumSpecies_range_weighted <- function(mask, rangeMapFiles, nCores, rangeSizes){
  
  #browser()
  
  ## for debugging
  # mask
  # rangeMapFiles <- allFiles[["Birds"]]
  # rangeSizes <- RangeSizes
  # nCores <- 3
  
  # take all cells that have non NA values in "mask" (raster file of the world, in Behrman projection)
  mask.vals <- !is.na(values(mask))
  
  # verify extent and resolution
  verif <- raster::raster(rangeMapFiles[1])
  
  if(extent(mask)!=extent(verif)) {
    print("Different extent between species raster and mask")
    stop()
  }
  
  if((res(mask)[1]!=res(verif)[1]) |(res(mask)[2]!=res(verif)[2])) {
    print("Different resolution between species raster and mask")
    stop()
  }
  
  # else{print("Start processing.")}
  
  rm(verif)
  
  # # # # # # # # 
  
  # parallelise and apply clusterFunc to range files
  cl <- parallel::makeCluster(nCores)
  
  parallel::clusterExport(cl = cl,
                          varlist = c("rangeMapFiles",
                                      "rangeSizes",
                                      "mask.vals"),
                          envir = environment())
  
  
  
  # Define function "clusterFunc" that takes a raster file, then intersects it with mask where the species is present
  
  clusterFunc <- function(rangeMaps, rangeSizes){
    
    # rangeMaps <- rangeMapFiles[1]
    #browser()
    
    # open range file
    r <- raster::raster(rangeMaps)
    
    # get range size for the species
    Range_size <- rangeSizes$Range_area_sq_km[rangeSizes$Index==gsub(names(r), replacement = "", pattern = "sp")]
    
    # from the raster, take all cells that have non NA value (where the species is present)
    sp.pres <- !is.na(raster::values(r))
    
    # create grid with NAs
    sp.pres2 <- rep(NA,length(mask.vals))
    
    # put 1/range size where the species is present and intersects with mask
    sp.pres2[(mask.vals) & (sp.pres)] <- 1/Range_size
    
    # return the grid with cells=1 where the species is present
    return(sp.pres2)
  }
  environment(clusterFunc) <- .GlobalEnv
  
  # Run clusterFunc in parallel for all range files
  # (returns a list of numeric elements for species presence in each range file )
  SpPres <- ParallelLogger::clusterApply(cl = cl,
                                         x = as.list(rangeMapFiles),
                                         rangeSizes=rangeSizes,
                                         fun = clusterFunc,
                                         progressBar=TRUE)
  
  stopCluster(cl)
  
  print("Stacking and summing across rasters.")
  
  # put results together from list to matrix (each column = data for 1 species)
  SpPres <- do.call('cbind', SpPres)
  # calculate species richness in each grid cell by summing on rows
  SpRich <- apply(SpPres, 1, sum, na.rm=TRUE)
  
  # Map species richness onto mask raster layer
  OutRaster <- mask
  values(OutRaster) <- SpRich
  
  return(OutRaster)
  
}


## SumSpecies but with change in resolution

SumSpecies_50km <- function(mask, rangeMapFiles, nCores){
  
  # verify extent and resolution
  verif <- raster::raster(rangeMapFiles[1])
  
  if(extent(mask)!=extent(verif)) {
    print("Different extent between species raster and mask")
    stop()
  }
  
  if((res(mask)[1]!=res(verif)[1]) |(res(mask)[2]!=res(verif)[2])) {
    print("Different resolution between species raster and mask")
    stop()
  }
  
  # else{print("Start processing.")}
  
  rm(verif)
  
  ## aggregate mask at 50km resolution
  mask <- raster::aggregate(mask, fact=5, fun=max, na.rm=TRUE)
  
  # take all cells that have non NA values in "mask" (raster file of the world, in Behrman projection)
  mask.vals <- !is.na(values(mask))
  
  # # # # # # # # 
  
  # parallelise and apply clusterFunc to range files
  cl <- parallel::makeCluster(nCores)
  
  parallel::clusterExport(cl = cl,
                          varlist = c("rangeMapFiles",
                                      "mask.vals"),
                          envir = environment())
  
  
  
  # Define function "clusterFunc" that takes a raster file, then intersects it with mask where the species is present
  
  clusterFunc <- function(rangeMaps){
    
    # open range file and aggregate at 50 km resolution
    r <- raster::raster(rangeMaps)
    r <- raster::aggregate(r, fact=5, fun=max, na.rm=TRUE)
    
    # from the raster, take all cells that have non NA value (where the species is present)
    sp.pres <- !is.na(raster::values(r))
    
    # create grid with NAs
    sp.pres2 <- rep(NA,length(mask.vals))
    
    # put 1 where the species is present and intersects with mask
    sp.pres2[(mask.vals) & (sp.pres)] <- 1
    
    # return the grid with cells=1 where the species is present
    return(sp.pres2)
  }
  environment(clusterFunc) <- .GlobalEnv
  
  # Run clusterFunc in parallel for all range files
  # (returns a list of numeric elements for species presence in each range file )
  SpPres <- ParallelLogger::clusterApply(cl = cl,
                                         x = as.list(rangeMapFiles),
                                         fun = clusterFunc,
                                         progressBar=TRUE)
  
  stopCluster(cl)
  
  print("Stacking and summing across rasters.")
  
  # put results together from list to matrix (each column = data for 1 species)
  SpPres <- do.call('cbind', SpPres)
  # calculate species richness in each grid cell by summing on rows
  SpRich <- apply(SpPres, 1, sum, na.rm=TRUE)
  
  # Map species richness onto mask raster layer
  OutRaster <- mask
  values(OutRaster) <- SpRich
  
  return(OutRaster)
  
}



## function to apply SumSpecies onto chunks of range files (to avoid memory size limitations)

ApplyToChunks <- function(rangeMapFiles, mask, nCores, ChunkSize) {
  
  # divide rangeMapFiles into smaller chunks, otherwise not enough memory size
  MapsChunks <- split(rangeMapFiles, ceiling(seq_along(rangeMapFiles)/ChunkSize))
  
  # initialize a raster layer with 0 everywhere
  Cumul <- mask
  values(Cumul) <- 0
  
  # cumulative species richness
  for (i in 1:length(MapsChunks)) {
    gc()
    print(paste("Processing chunk", i, "on", length(MapsChunks)))
    ResChunk <- SumSpecies(rangeMapFiles = MapsChunks[[i]], mask=mask, nCores = nCores) 
    Cumul <- stack(Cumul, ResChunk)
    rm(ResChunk)
    Cumul <- calc(Cumul, sum)
  }
  
  return(Cumul)
}

ApplyToChunks_rangeweighted <- function(rangeMapFiles, mask, nCores, ChunkSize, rangeSizes) {
  
  # divide rangeMapFiles into smaller chunks, otherwise not enough memory size
  MapsChunks <- split(rangeMapFiles, ceiling(seq_along(rangeMapFiles)/ChunkSize))
  
  # initialize a raster layer with 0 everywhere
  Cumul <- mask
  values(Cumul) <- 0
  
  # cumulative species richness
  for (i in 1:length(MapsChunks)) {
    gc()
    print(paste("Processing chunk", i, "on", length(MapsChunks)))
    #browser()
    ResChunk <- SumSpecies_range_weighted(rangeMapFiles = MapsChunks[[i]], mask=mask, nCores = nCores, rangeSizes = rangeSizes) 
    Cumul <- stack(Cumul, ResChunk)
    rm(ResChunk)
    Cumul <- calc(Cumul, sum)
  }
  
  return(Cumul)
}


# # # # # # # # # # # # # 


## functions to calculate mean trait completeness at 10 km x 10 km resolution

Mean_completeness_spatial <- function(mask, rangeMapFiles, nCores, trait.completeness){
  
  # take all cells that have non NA values in "mask" (raster file of the world)
  mask.vals <- !is.na(values(mask))
  
  # Define function "clusterFunc" that takes a range file, converts it to raster, then intersects it with mask where the species is present
  
  clusterFunc.central_completeness <- function(rangefile, trait.completeness){
    
    # take range file and convert to raster
    r <- raster::raster(rangefile)
    
    # from the raster, take all cells that have non NA value (where the species is present)
    sp.pres <- !is.na(raster::values(r))
    sp.pres2 <- rep(NA, length(mask.vals))
    
    # get species name
    species <-  strsplit(rangefile,'[/]')[[1]]
    species <- species[grepl("sp", species)]
    
    Row <- which(trait.completeness$Index==species)
    # get completeness for that species
    if (length(Row)!=0) {
      compl <- trait.completeness$completeness[Row]
    } else{
      compl <- NA
    }
    
    # Put trait completeness of the corresponding species where the species is present
    sp.pres2[(mask.vals) & (sp.pres)] <- compl
    
    # return numeric vec with trait completeness where the species is present    
    return(sp.pres2)
    
  }
  
  # parallelise
  
  cl <- parallel::makeCluster(nCores)
  
  parallel::clusterExport(cl = cl,varlist = c("rangeMapFiles","mask.vals", "clusterFunc.central_completeness", "trait.completeness"),
                          envir = environment())
  
  
  SpPres <- ParallelLogger::clusterApply(cluster = cl,
                                         fun = clusterFunc.central_completeness,
                                         x = as.list(rangeMapFiles),
                                         trait.completeness=trait.completeness,
                                         progressBar = TRUE)
  
  stopCluster(cl)
  
  # put results together from list to matrix (each column = data for a species)
  SpMap <- do.call('cbind', SpPres)
  rm(SpPres)
  
  
  # calculate mean trait completeness in each grid cell
  gc()
  memory.limit(size=8e6)
  
  print("Computing mean")
  
  Completeness.Mean <- apply(SpMap, 1, mean, na.rm=TRUE)
  rm(SpMap)
  
  # # if there are "holes", map them as -1 (because NAs represent the seas and oceans)
  # Completeness.Mean[is.na(Completeness.Mean) & mask.vals] <- -1
  
  # Map species richness onto main raster layer
  outRaster.mean <- mask
  values(outRaster.mean) <- Completeness.Mean
  
  return(outRaster.mean)
}


Mean_ApplyToChunks <- function(rangeMapFiles, mask, nCores, ChunkSize, trait.completeness) {
  
  # divide rangeMapFiles into smaller chunks, otherwise not enough memory size
  MapsChunks <- split(rangeMapFiles, ceiling(seq_along(rangeMapFiles)/ChunkSize))
  
  Chunks <- list()
  
  # mean for each chunk
  for (i in 1:length(MapsChunks)) {
    gc()
    print(paste("Processing chunk", i, "on", length(MapsChunks)))
    Chunks[[i]] <- Mean_completeness_spatial(rangeMapFiles = MapsChunks[[i]], mask=mask, nCores = nCores, trait.completeness=trait.completeness) 
  }
  
  MeanComp <- stack(Chunks)
  rm(Chunks)
  MeanComp <- calc(MeanComp, mean, na.rm=TRUE)
  
  # # if there are "holes", map them as -1 (because NAs represent the seas and oceans)
  mask.vals <- !is.na(values(mask))
  MeanComp[is.na(values(MeanComp)) & mask.vals] <- -1
  
  return(MeanComp)
}



# # # # # # # # # # # # # 


## Function to put the range files together for a given class, and to calculate median or variance trait completeness
## from the raster files across species in each grid cell AT 50 X 50 km resolution (aggregating the rasters)

## Need to adapt this version for median and mean at lower resolution

Median_Var_completeness_50k <- function(mask, rangeMapFiles, nCores, trait.completeness){
  
  # resample mask at 50x50k resolution
  mask_50 <- raster::aggregate(mask, fact=5, fun=max, na.rm=TRUE)
  rm(mask)
  
  # take all cells that have non NA values in "mask" (raster file of the world)
  mask.vals.50 <- !is.na(values(mask_50))
  
  # Define function "clusterFunc" that takes a range file, converts it to raster, then intersects it with mask where the species is present
  
  clusterFunc.central_completeness50 <- function(rangefile, trait.completeness){
    
    # take range file and convert to raster
    r <- raster::raster(rangefile)
    
    # aggregate raster at 50x50 km
    r <- raster::aggregate(r, fact=5, fun=max, na.rm=TRUE)
    
    # from the raster, take all cells that have non NA value (where the species is present)
    sp.pres <- !is.na(raster::values(r))
    sp.pres2 <- rep(NA, length(mask.vals.50))
    
    # get species name
    species <-  strsplit(rangefile,'[/]')[[1]]
    species <- species[grepl("sp", species)]
    
    Row <- which(trait.completeness$Index==species)
    # get completeness for that species
    if (length(Row)!=0) {
      compl <- trait.completeness$completeness[Row]
    } else{
      compl <- NA
    }
    
    # Put trait completeness of the corresponding species where the species is present
    sp.pres2[(mask.vals.50) & (sp.pres)] <- compl
    
    # return numeric vec with trait completeness where the species is present    
    return(sp.pres2)
    
  }
  
  print("Intersecting rasters")
  
  # parallelise
  cl <- parallel::makeCluster(nCores)
  
  parallel::clusterExport(cl = cl,varlist = c("rangeMapFiles","mask.vals.50", "clusterFunc.central_completeness50", "trait.completeness"),
                          envir = environment())
  
  SpPres <- ParallelLogger::clusterApply(cluster = cl,
                                         fun = clusterFunc.central_completeness50,
                                         x = as.list(rangeMapFiles),
                                         trait.completeness=trait.completeness,
                                         progressBar = TRUE)
  
  
  stopCluster(cl)
  
  # browser()
  
  
  # put results together from list to matrix (each column = data for a species)
  SpMap <- do.call('cbind', SpPres)
  rm(SpPres)
  
  # calculate median or variance trait completeness in each grid cell
  gc()
  memory.limit(size=8e6)
  
  print("Calculating median at 50km resolution")
  Median_Comp <- apply(SpMap, 1, median, na.rm=TRUE)
  
  print("Calculating mean at 50km resolution")
  Mean_Comp <- apply(SpMap, 1, mean, na.rm=TRUE)
  
  print("Calculating variance at 50km resolution")
  Var_Comp <- apply(SpMap, 1, var, na.rm=TRUE)
  
  rm(SpMap)
  gc()
  
  # if there are "holes", map them as -1 (because NAs represent the seas and oceans)
  Median_Comp[is.na(Median_Comp) & mask.vals.50] <- -1
  Var_Comp[is.na(Var_Comp) & mask.vals.50] <- -1
  Mean_Comp[is.na(Mean_Comp) & mask.vals.50] <- -1
  
  # Map species richness onto main raster layer
  outRasterMedian <- mask_50
  outRasterVariance <- mask_50
  outRasterMean <- mask_50
  
  values(outRasterMedian) <- Median_Comp
  values(outRasterVariance) <- Var_Comp
  values(outRasterMean) <- Mean_Comp
  
  return(list(Mean=outRasterMean, Median=outRasterMedian, Variance=outRasterVariance))
}



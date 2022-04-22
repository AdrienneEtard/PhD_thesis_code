#### RASTERISING WWF REALMS - (NB. already run)

realms <- readOGR(dsn="../../Data/Ecoregions", layer="tnc_terr_ecoregions")
crs(realms)
behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
realms <- spTransform(realms, CRSobj = behrCRS)
# reproject realms using Behrman's projection
realms <- spTransform(realms, crs(state_boundary_us))
## at 10km resolution
Raster_completeness <- readRDS("../../Results/SpatialAnalyses/MeanCompletenessRasters_10k/Amphibians.rds")
# convert shapefile to raster with the same number of cells than the other rasters
r <- raster(Raster_completeness)
values(r) <- NA
raster_realms <- rasterize(realms, r, "WWF_REALM2")
plot(raster_realms)
saveRDS(raster_realms, "../../Results/Ecoregions_rasters/realms_10k.rds")
## at 50 km resolution
Raster_completeness <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianAmphibians.rds")
r <- raster(Raster_completeness)
values(r) <- NA
raster_realms <- rasterize(realms, r, "WWF_REALM2")
saveRDS(raster_realms, "../../Results/Ecoregions_rasters/realms_50k.rds")



##### FUNCTION TO FIT SPATIAL MODELS

FitSpatialModel <- function(Raster_completeness, Raster_richness, samplesize, cutoff, realmraster, Realm_Inter, frac) {
  
  # for debugging the function
  # Raster_completeness <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianReptiles.rds")
  # Raster_richness <- Richness_reptiles
  # samplesize <- 2000
  # cutoff <- 3
  # realmraster <- readRDS("../../Results/Ecoregions_rasters/realms_50k.rds")
  #

  ## Model with and without quadratic term, and null model: compared using loglikelihood ratio tests
  
  # get the data from the raster layers, and sample (for a smaller dataset with 10,000 datapoints)
  mod.data <- data.frame(Completeness=raster::values(Raster_completeness),
                         SpeciesRichness=raster::values(Raster_richness),
                         Longitude=raster::coordinates(Raster_richness)[,1],
                         Latitude=raster::coordinates(Raster_richness)[,2])
  
  # if models allow for interaction with realms, add realms
  if(Realm_Inter) {
    
    realm.data <- data.frame(realm=values(realmraster),
                             Longitude=coordinates(realmraster)[,1],
                             Latitude=coordinates(realmraster)[,2])
    
    mod.data$Realm <- realm.data$realm
    mod.data$Realm <- factor(mod.data$Realm, levels=c(1:8), labels=c("Afrotropic",
                                                                     "Antarctic",
                                                                     "Australasia",
                                                                     "Indo-Malay",
                                                                     "Neartic",
                                                                     "Neotropic",
                                                                     "Oceania",
                                                                     "Palearctic"))
    # only consider realms for which sample size > sample size
    Getrid <- mod.data %>% 
      dplyr::group_by(Realm) %>%
      dplyr::summarise(Count=n()) %>% 
      dplyr::filter(!Count>=samplesize) %>% 
      dplyr::select(Realm) %>% 
      as.data.frame()
    Getrid <- Getrid$Realm %>% 
      as.character()
    
    mod.data <- mod.data %>%
      filter(!Realm %in% Getrid) %>% 
      filter(!is.na(Realm))
    
  }
  
  # filter out NAs (non terrestrial areas) and -1 (where completeness is NA)
  mod.data$Completeness[mod.data$Completeness==-1] <- NA
  mod.data <- na.omit(mod.data)
  
  # take all areas where SR>cutoff to avoid sampling biases
  mod.data <- mod.data[mod.data$SpeciesRichness>=cutoff,]
  mod.data$Completeness <- mod.data$Completeness / 100

  # sample randomly for sample size data points in each realm and fit the three models
  if(Realm_Inter) {
    mod.data <- mod.data %>%
      dplyr::group_by(Realm) %>%
      dplyr::sample_frac(frac,replace=FALSE) %>%
      as.data.frame()
    # ggplot(mod.data, aes(log(SpeciesRichness), Completeness)) + geom_point() + facet_wrap(~Realm,scales = "free" )

    # fit a spatial simultaneous autoregressive lag model with realm as interacting facot

    # coordinates as neighbours "The function uses the deldir package to
    # convert a matrix of two-dimensional coordinates into a neighbours list of class nb
    # with a list of integer vectors containing neighbour region number ids."
    nb <- spdep::tri2nb(coords = mod.data[,c('Longitude','Latitude')])
    lw <- spdep::nb2listw(neighbours = nb)

    # # fit the model with quadratic term
    # print("Fitting model 1")
    # Model1 <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~
    #                                  log(mod.data$SpeciesRichness)+
    #                                  I(log(mod.data$SpeciesRichness)^2) +
    #                                  mod.data$Realm +
    #                                  mod.data$Realm:log(mod.data$SpeciesRichness)+
    #                                  mod.data$Realm:I(log(mod.data$SpeciesRichness)^2),
    #                                listw = lw)

    # model without quadratic term
    print("Fitting model 2")
    Model2 <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~
                                     log(mod.data$SpeciesRichness)+
                                     mod.data$Realm +
                                     mod.data$Realm:log(mod.data$SpeciesRichness),
                                     listw = lw)

    # # Null model
    # print("Fitting model 0")
    # Model0 <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~ 1, listw = lw)


    }


  # sample randomly for samplesize data points
  else{

    mod.data <- mod.data[sample(x = 1:nrow(mod.data), size = samplesize, replace = FALSE), ]

  # fit a spatial simultaneous autoregressive lag model

  # coordinates as neighbours "The function uses the deldir package to
  # convert a matrix of two-dimensional coordinates into a neighbours list of class nb
  # with a list of integer vectors containing neighbour region number ids."
  nb <- spdep::tri2nb(coords = mod.data[,c('Longitude','Latitude')])
  lw <- spdep::nb2listw(neighbours = nb)

  # fit the model with quadratic term
  print("Fitting model 1")
  Model1norealm <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~
                                   log(mod.data$SpeciesRichness)+
                                   I(log(mod.data$SpeciesRichness)^2),
                                 listw = lw)

  # model without quadratic term
  print("Fitting model 2")
  Model2norealm <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~
                                   log(mod.data$SpeciesRichness),
                                 listw = lw)

  # Null model
  print("Fitting model 0")
  Model0norealm <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~ 1, listw = lw)


 }


  ## Model comparisons
  # anova and log likelyhood ratio test
  # Anova_results_10 <- anova(Model0,Model1)
  # Anova_results_12 <- anova(Model2,Model1)
  # Anova_results_20 <- anova(Model0,Model2)

  return(list(data=mod.data,
              #model1=Model1,
              model2=Model2
              #NullModel=Model0,
              #Anova10=Anova_results_10,
              #Anova12=Anova_results_12,
              #Anova20=Anova_results_20
              ))
  
}


SampleSize <- function(Raster_completeness, Raster_richness, samplesize, cutoff, realmraster, Realm_Inter, frac) {
  
  ## Model with and without quadratic term, and null model: compared using loglikelihood ratio tests
  
  # get the data from the raster layers, and sample (for a smaller dataset with 10,000 datapoints)
  mod.data <- data.frame(Completeness=raster::values(Raster_completeness),
                         SpeciesRichness=raster::values(Raster_richness),
                         Longitude=raster::coordinates(Raster_richness)[,1],
                         Latitude=raster::coordinates(Raster_richness)[,2])
  
  # if models allow for interaction with realms, add realms
  if(Realm_Inter) {
    
    realm.data <- data.frame(realm=values(realmraster),
                             Longitude=coordinates(realmraster)[,1],
                             Latitude=coordinates(realmraster)[,2])
    
    mod.data$Realm <- realm.data$realm
    mod.data$Realm <- factor(mod.data$Realm, levels=c(1:8), labels=c("Afrotropic",
                                                                     "Antarctic",
                                                                     "Australasia",
                                                                     "Indo-Malay",
                                                                     "Neartic",
                                                                     "Neotropic",
                                                                     "Oceania",
                                                                     "Palearctic"))
    # only consider realms for which sample size > sample size
    Getrid <- mod.data %>% 
      dplyr::group_by(Realm) %>%
      dplyr::summarise(Count=n()) %>% 
      dplyr::filter(!Count>=samplesize) %>% 
      dplyr::select(Realm) %>% 
      as.data.frame()
    Getrid <- Getrid$Realm %>% 
      as.character()
    
    mod.data <- mod.data %>%
      filter(!Realm %in% Getrid) %>% 
      filter(!is.na(Realm))
    
  }
  
  # filter out NAs (non terrestrial areas) and -1 (where completeness is NA)
  mod.data$Completeness[mod.data$Completeness==-1] <- NA
  mod.data <- na.omit(mod.data)
  
  # take all areas where SR>cutoff to avoid sampling biases
  mod.data <- mod.data[mod.data$SpeciesRichness>=cutoff,]
  mod.data$Completeness <- mod.data$Completeness / 100
  
  # sample randomly for sample size data points in each realm and fit the three models
  if(Realm_Inter) {
    mod.data <- mod.data %>%
      dplyr::group_by(Realm) %>%
      dplyr::sample_frac(frac,replace=FALSE) %>%
      as.data.frame()
  }
  
  
  # sample randomly for samplesize data points
  else{
    mod.data <- mod.data[sample(x = 1:nrow(mod.data), size = samplesize, replace = FALSE), ]
  }
  
  
  return(mod.data)
  
}

##### FUNCTIONS TO PLOT THE BEST FITTING SPATIAL MODELS
GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


InSamplePredictions <- function(Model, mod.data, Inter.realms){
  
  browser()
  
  # in-sample predictions for the model
  predictions <- spatialreg::predict.sarlm(object = Model,
                                           newdata=NULL,
                                           listw = NULL,
                                           zero.policy = TRUE)
  
  
  # create dataframe to store results
  
  # trend
  Spatial_predictions <-  attr(predictions, "trend") %>%
    as.data.frame() %>% 
    setNames(., "trend")
  
  # signal
  Spatial_predictions$signal <- attr(predictions, "signal")
  
  # fit = trend + signal
  Spatial_predictions$fit <- Spatial_predictions$signal + Spatial_predictions$trend
  
  # add species richness
  Spatial_predictions$SR <- mod.data$SpeciesRichness
  
  # add realms if interactions are taken into account
  if(Inter.realms){
    Spatial_predictions$Realm <- mod.data$Realm
  }
  
  # back-transforming
  Spatial_predictions$fit_BT <- (sin(Spatial_predictions$fit))^2
  Spatial_predictions$signal_BT <- (sin(Spatial_predictions$signal))^2
  Spatial_predictions$trend_BT <- (sin(Spatial_predictions$trend))^2
  
  
  return(Spatial_predictions)
}

PlotInSamplePreds <- function(Preds, Inter.realms) {
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"),
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  if(Inter.realms) {
    p <- ggplot(Preds, aes(x=SR, y=trend_BT*100, col=Realm)) + geom_line()+ GGPoptions 
  }
  
  else{
    p <- ggplot(Preds, aes(x=SR, y=trend_BT*100)) + geom_line()+ GGPoptions 
  }
  
  return(p)
}

# # # # Plotting with 95% intervals or SE interval constructed manually

Preds_Summary_NoQuad <- function(Model, ModelData, cutoff, Class, SE) {
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"),
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  x <- summary(Model)
  Coef <- x$Coef %>% 
    as.data.frame()
  
  if(!SE){
    Coef$CI <-  Coef$`Std. Error`*1.96
  }
  
  if (SE) {
    Coef$CI <-  Coef$`Std. Error`
  }
  
  
  ## slope estimates - log SR
  Slopes <- Coef[grepl("SpeciesRichness", rownames(Coef)),]
  Slopes$Signif <- ifelse(Slopes$`Pr(>|z|)`<0.05, TRUE, FALSE)
  rownames(Slopes) <- c("Afrotropic", "Australasia", "Indo-Malay", "Neartic", "Neotropic", "Palearctic")
  
  # for(i in 2:nrow(Slopes)) {
  # 
  #   if(Slopes$Signif[i]==TRUE) {
  #     Slopes$Estimate[i] <- Slopes$Estimate[i] + Slopes$Estimate[1]
  #   }
  # 
  #   else{
  #     Slopes$Estimate[i] <- Slopes$Estimate[1]
  #   }
  # 
  # }
  
 Slopes$Estimate[2:nrow(Slopes)] <- Slopes$Estimate[2:nrow(Slopes)] + Slopes$Estimate[1]
  
  
  ## intercept estimates
  Intercept <- Coef[!grepl("SpeciesRichness", rownames(Coef)),]
  Intercept$Signif <- ifelse(Intercept$`Pr(>|z|)`<0.05, TRUE, FALSE)
  rownames(Intercept) <- c("Afrotropic", "Australasia", "Indo-Malay", "Neartic", "Neotropic", "Palearctic")
  
  Intercept$Estimate[2:nrow(Intercept)] <- Intercept$Estimate[2:nrow(Intercept)] + Intercept$Estimate[1]
  
  # for(i in 2:nrow(Intercept)) {
  # 
  #   if(Intercept$Signif[i]==TRUE) {
  #     Intercept$Estimate[i] <- Intercept$Estimate[i] + Intercept$Estimate[1]
  #   }
  # 
  #   else{
  #     Intercept$Estimate[i] <- Intercept$Estimate[1]
  #   }
  # }
  # 
  
  ## predictions for each realm
  SR_preds <- ModelData %>% 
    group_by(Realm) %>% 
    summarise(Max=max(SpeciesRichness)) %>% 
    as.data.frame()
  
  Preds <- matrix(ncol=nrow(SR_preds), nrow=max(SR_preds$Max))
  colnames(Preds) <- SR_preds$Realm
  Preds <- as.data.frame(Preds)
  
  # for CI up
  Preds_CI.up <- Preds 
  Preds_CI.low <- Preds 
  
  
  for (r in 1:length(colnames(Preds))) {
    Preds[cutoff:SR_preds$Max[r],colnames(Preds)[r]] <- c(cutoff:SR_preds$Max[r])
    Preds_CI.up[cutoff:SR_preds$Max[r],colnames(Preds_CI.up)[r]] <- c(cutoff:SR_preds$Max[r])
    Preds_CI.low[cutoff:SR_preds$Max[r],colnames(Preds_CI.low)[r]] <- c(cutoff:SR_preds$Max[r])
    
    ## predictions
    Preds[1:SR_preds$Max[r],colnames(Preds)[r]] <-
      Intercept$Estimate[rownames(Intercept)==colnames(Preds)[r]] +
      Slopes$Estimate[rownames(Slopes)==colnames(Preds)[r]] * log(Preds[1:SR_preds$Max[r],colnames(Preds)[r]])
    
    Preds_CI.up[1:SR_preds$Max[r],colnames(Preds_CI.up)[r]] <-
      Intercept$Estimate[rownames(Intercept)==colnames(Preds_CI.up)[r]] + Intercept$CI[rownames(Intercept)==colnames(Preds_CI.up)[r]] +
      (Slopes$Estimate[rownames(Slopes)==colnames(Preds_CI.up)[r]] + Slopes$CI[rownames(Slopes)==colnames(Preds_CI.up)[r]]) * log(Preds_CI.up[1:SR_preds$Max[r],colnames(Preds_CI.up)[r]])
    
    Preds_CI.low[1:SR_preds$Max[r],colnames(Preds_CI.low)[r]] <-
      Intercept$Estimate[rownames(Intercept)==colnames(Preds_CI.low)[r]] - Intercept$CI[rownames(Intercept)==colnames(Preds_CI.low)[r]] +
      (Slopes$Estimate[rownames(Slopes)==colnames(Preds_CI.low)[r]]- Slopes$CI[rownames(Slopes)==colnames(Preds_CI.low)[r]]) * log(Preds_CI.low[1:SR_preds$Max[r],colnames(Preds_CI.low)[r]])
    
  }
  
  
  Preds_melt <- reshape2::melt(data = Preds) %>% 
    group_by(variable) %>% 
    mutate(SR=seq(1:length(variable))) %>% 
    as.data.frame() %>% 
    setNames(c("Realms", "Trend", "Species richness")) %>% 
    filter(!is.na(Trend))
  
  Preds_melt_CI.up <- reshape2::melt(data = Preds_CI.up) %>% 
    group_by(variable) %>% 
    mutate(SR=seq(1:length(variable))) %>% 
    as.data.frame() %>% 
    setNames(c("Realms", "Trend", "Species richness")) %>% 
    filter(!is.na(Trend))
  colnames(Preds_melt_CI.up)[2] <- "CI_up"
  
  Preds_melt_CI.low <- reshape2::melt(data = Preds_CI.low) %>% 
    group_by(variable) %>% 
    mutate(SR=seq(1:length(variable))) %>% 
    as.data.frame() %>% 
    setNames(c("Realms", "Trend", "Species richness")) %>% 
    filter(!is.na(Trend))
  colnames(Preds_melt_CI.low)[2] <- "CI_low"
  
  Preds_melt$CI_up <- Preds_melt_CI.up$CI_up
  Preds_melt$CI_low <- Preds_melt_CI.low$CI_low
  Preds_melt$Class <- Class
  
  # # Backtransform predictions 
  # Preds_melt_BT <- Preds_melt
  # Preds_melt_BT[,c("Trend", "CI_up", "CI_low")] <- (sin(Preds_melt_BT[,c("Trend", "CI_up", "CI_low")])^2)*100
  
  # adding significance
  
  # if all significant
  if(all(Intercept$Signif)) {
    Preds_melt$Intercept.Significance <- TRUE
  }
  else{
    NotSignifInter <- rownames(Intercept[Intercept$Signif==FALSE,])
    Preds_melt$Intercept.Significance <- TRUE
    Preds_melt$Intercept.Significance[Preds_melt$Realms %in% NotSignifInter] <- FALSE
  }
  
  if(all(Slopes$Signif)) {
    Preds_melt$Slope.Significance <- TRUE
  }
  else{
    NotSignifSlope <- rownames(Slopes[Slopes$Signif==FALSE,])
    Preds_melt$Slope.Significance <- TRUE
    Preds_melt$Slope.Significance[Preds_melt$Realms %in% NotSignifSlope] <- FALSE
  }
  
  return(Preds_melt)
  
  
}

Preds_Summary_Quad <- function(Model, ModelData, cutoff, Class, SE) {
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"),
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  x <- summary(Model)
  Coef <- x$Coef %>% 
    as.data.frame()
  
  if(nrow(Coef)==0) {
    Coef <- x$coefficients %>%  as.data.frame()
  }
  
  if(!SE){
    Coef$CI <-  Coef$`Std. Error`*1.96
  }
  
  if (SE) {
    Coef$CI <-  Coef$`Std. Error`
  }
  
  ## slope estimates - log SR
  Slopes <- Coef[grepl("SpeciesRichness", rownames(Coef)),]
  
  if(!is.null(Slopes$`Pr(>|z|)`)){
    Slopes$Signif <- ifelse(Slopes$`Pr(>|z|)`<0.05, TRUE, FALSE)
  } else{
    Slopes$Signif <- ifelse(Slopes$`Pr(>|t|)`<0.05, TRUE, FALSE)
  }
  
  # slopes for logSR
  Slopes1 <- Slopes[!grepl("2", rownames(Slopes)),]
  
  # slopes for logSR^2
  Slopes2 <- Slopes[grepl("2", rownames(Slopes)),]
  
  rownames(Slopes1) <- c("Afrotropic", "Australasia", "Indo-Malay", "Neartic", "Neotropic", "Palearctic")
  rownames(Slopes2) <- c("Afrotropic", "Australasia", "Indo-Malay", "Neartic", "Neotropic", "Palearctic")
  
  for(i in 1:nrow(Slopes1)) {
    if(Slopes1$Signif[i]==TRUE){
      Slopes1$Estimate[i] <- Slopes1$Estimate[i] + Slopes1$Estimate[1]
    } 
  }
  for(i in 1:nrow(Slopes2)) {
    if(Slopes2$Signif[i]==TRUE){
      Slopes2$Estimate[i] <- Slopes2$Estimate[i] + Slopes2$Estimate[1]
    } 
  }
  
  
  # Slopes1$Estimate[2:nrow(Slopes1)] <- Slopes1$Estimate[2:nrow(Slopes1)] + Slopes1$Estimate[1]
  # Slopes2$Estimate[2:nrow(Slopes2)] <- Slopes2$Estimate[2:nrow(Slopes2)] + Slopes2$Estimate[1]
  
  ## intercept estimates
  Intercept <- Coef[!grepl("SpeciesRichness", rownames(Coef)),]
  
  if(!is.null(Intercept$`Pr(>|z|)`)){
    Intercept$Signif <- ifelse(Intercept$`Pr(>|z|)`<0.05, TRUE, FALSE)
  } else{
    Intercept$Signif <- ifelse(Intercept$`Pr(>|t|)`<0.05, TRUE, FALSE)
  }
  
  
  rownames(Intercept) <- c("Afrotropic", "Australasia", "Indo-Malay", "Neartic", "Neotropic", "Palearctic")
  
  #Intercept$Estimate[2:nrow(Intercept)] <- Intercept$Estimate[2:nrow(Intercept)] + Intercept$Estimate[1]
  
  for(i in 1:nrow(Intercept)) {
    if(Intercept$Signif[i]==TRUE){
      Intercept$Estimate[i] <- Intercept$Estimate[i] + Intercept$Estimate[1]
    } 
  }
  
  ## predictions for each realm
  SR_preds <- ModelData %>% 
    group_by(Realm) %>% 
    summarise(Max=max(SpeciesRichness)) %>% 
    as.data.frame()
  
  Preds <- matrix(ncol=nrow(SR_preds), nrow=max(SR_preds$Max))
  colnames(Preds) <- SR_preds$Realm
  Preds <- as.data.frame(Preds)
  
  # for CI up
  Preds_CI.up <- Preds 
  Preds_CI.low <- Preds 
  
  
  for (r in 1:length(colnames(Preds))) {
    
    Preds[cutoff:SR_preds$Max[r],colnames(Preds)[r]] <- c(cutoff:SR_preds$Max[r])
    Preds_CI.up[cutoff:SR_preds$Max[r],colnames(Preds_CI.up)[r]] <- c(cutoff:SR_preds$Max[r])
    Preds_CI.low[cutoff:SR_preds$Max[r],colnames(Preds_CI.low)[r]] <- c(cutoff:SR_preds$Max[r])
  }
  
  ## predictions
  for (r in 1:length(colnames(Preds))) {
    Preds[1:SR_preds$Max[r],colnames(Preds)[r]] <-
      Intercept$Estimate[rownames(Intercept)==colnames(Preds)[r]] +
      Slopes1$Estimate[rownames(Slopes1)==colnames(Preds)[r]] * log(Preds[1:SR_preds$Max[r],colnames(Preds)[r]]) +
      Slopes2$Estimate[rownames(Slopes2)==colnames(Preds)[r]] * (log(Preds[1:SR_preds$Max[r],colnames(Preds)[r]]))^2
  }
  
  ## confidence interval - upper limit
  
  for (r in 1:length(colnames(Preds_CI.up))){
    Preds_CI.up[1:SR_preds$Max[r],colnames(Preds_CI.up)[r]] <-
      Intercept$Estimate[rownames(Intercept)==colnames(Preds_CI.up)[r]] + Intercept$CI[rownames(Intercept)==colnames(Preds_CI.up)[r]] +
      (Slopes1$Estimate[rownames(Slopes1)==colnames(Preds_CI.up)[r]] + Slopes1$CI[rownames(Slopes1)==colnames(Preds_CI.up)[r]]) * log(Preds_CI.up[1:SR_preds$Max[r],colnames(Preds_CI.up)[r]]) +
      (Slopes2$Estimate[rownames(Slopes2)==colnames(Preds_CI.up)[r]] + Slopes2$CI[rownames(Slopes2)==colnames(Preds_CI.up)[r]]) * (log(Preds_CI.up[1:SR_preds$Max[r],colnames(Preds_CI.up)[r]]))^2
  }
  
  ## confidence interval - lower limit
  for (r in 1:length(colnames(Preds_CI.low))){
    Preds_CI.low[1:SR_preds$Max[r],colnames(Preds_CI.low)[r]] <-
      Intercept$Estimate[rownames(Intercept)==colnames(Preds_CI.low)[r]] - Intercept$CI[rownames(Intercept)==colnames(Preds_CI.low)[r]] +
      (Slopes1$Estimate[rownames(Slopes1)==colnames(Preds_CI.low)[r]]- Slopes1$CI[rownames(Slopes1)==colnames(Preds_CI.low)[r]]) * log(Preds_CI.low[1:SR_preds$Max[r],colnames(Preds_CI.low)[r]]) +
      (Slopes2$Estimate[rownames(Slopes2)==colnames(Preds_CI.low)[r]]- Slopes2$CI[rownames(Slopes2)==colnames(Preds_CI.low)[r]]) * (log(Preds_CI.low[1:SR_preds$Max[r],colnames(Preds_CI.low)[r]]))^2
  }
  
  ## melting 
  
  Preds_melt <- reshape2::melt(data = Preds) %>% 
    group_by(variable) %>% 
    mutate(SR=seq(1:length(variable))) %>% 
    as.data.frame() %>% 
    setNames(c("Realms", "Trend", "Species richness")) %>% 
    filter(!is.na(Trend))
  
  Preds_melt_CI.up <- reshape2::melt(data = Preds_CI.up) %>% 
    group_by(variable) %>% 
    mutate(SR=seq(1:length(variable))) %>% 
    as.data.frame() %>% 
    setNames(c("Realms", "Trend", "Species richness")) %>% 
    filter(!is.na(Trend))
  colnames(Preds_melt_CI.up)[2] <- "CI_up"
  
  Preds_melt_CI.low <- reshape2::melt(data = Preds_CI.low) %>% 
    group_by(variable) %>% 
    mutate(SR=seq(1:length(variable))) %>% 
    as.data.frame() %>% 
    setNames(c("Realms", "Trend", "Species richness")) %>% 
    filter(!is.na(Trend))
  colnames(Preds_melt_CI.low)[2] <- "CI_low"
  
  Preds_melt$CI_up <- Preds_melt_CI.up$CI_up
  Preds_melt$CI_low <- Preds_melt_CI.low$CI_low
  Preds_melt$Class <- Class
  
  # adding significance
  
  # if all significant
  if(all(Intercept$Signif)) {
    Preds_melt$Intercept.Significance <- TRUE
  } else {
    NotSignifInter <- rownames(Intercept[Intercept$Signif==FALSE,])
    Preds_melt$Intercept.Significance <- TRUE
    Preds_melt$Intercept.Significance[Preds_melt$Realms %in% NotSignifInter] <- FALSE
  }
  
  if(all(Slopes1$Signif)) {
    Preds_melt$Slope.Significance <- TRUE
  } else {
    NotSignifSlope1 <- rownames(Slopes1[Slopes1$Signif==FALSE,])
    Preds_melt$Slope.Significance <- TRUE
    Preds_melt$Slope.Significance[Preds_melt$Realms %in% NotSignifSlope1] <- FALSE
  }
  
  if(all(Slopes2$Signif)) {
    Preds_melt$Slope.Significance <- TRUE
  } else {
    NotSignifSlope2 <- rownames(Slopes2[Slopes2$Signif==FALSE,])
    Preds_melt$Slope.Significance <- TRUE
    Preds_melt$Slope.Significance[Preds_melt$Realms %in% NotSignifSlope2] <- FALSE
  }
  
  return(Preds_melt)
  
  
}

# function to plot, facetted by class or by realm
PlotPreds <- function (PredsAmphibians, PredsBirds, PredsMammals, PredsReptiles, facetting, error, backtransform, significance, logt) {
  
  Preds <- rbind(PredsAmphibians, PredsBirds, PredsMammals, PredsReptiles)
  
  Preds$Trend_BT <- sin(Preds$Trend)^2*100
  Preds$CI_up_BT <- sin(Preds$CI_up)^2*100
  Preds$CI_low_BT <- sin(Preds$CI_low)^2*100
  
  if(logt){
    Preds$`Species richness` <- log(Preds$`Species richness`)
  }
  
  Preds <- Preds %>% 
    filter(Class %in% c("Reptiles", "Amphibians"))
  
  if(backtransform){
    
    if(error){
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Realms))+
        geom_line(size=1) +
        geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free")
    }else{
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Realms))+
        geom_line(size=1) +
        #geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free") }
    
  } else{
    
    if(error){
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Realms))+
        geom_line(size=1) +
        geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free")  
    }else{
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Realms))+
        geom_line(size=1) +
        #geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free")
    }
    
  }
  
  
  if(facetting=="realm"){
    
    if(backtransform){
      if(error){
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Class))+
          geom_line(size=1) +
          geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
        
      }else{
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Class))+
          geom_line(size=1) +
          #geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
      }
      
    } else{
      if(error){
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Class))+
          geom_line(size=1) +
          geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
        
      }else{
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Class))+
          geom_line(size=1) +
          #geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
        
      }
      
    }
  }
  
  if(facetting=="realmclass") {
    
    p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Class))+
      geom_line() +
      # geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Class, col=NULL), alpha=0.1) +
      GGPoptions + facet_wrap(~Class+Realms, scales = "free")
  }
  
  p <- p + scale_color_brewer(palette="Dark2")
  
  return(p)
}

PlotPredsv2 <- function (PredsAmphibians, PredsReptiles, PredsMammals,PredsBirds, facetting, error, backtransform, significance, logt, allclasses) {
  
  
  Preds <- rbind(PredsAmphibians, PredsBirds,PredsMammals, PredsReptiles)

  Preds$Realms <- as.character(Preds$Realms)
  Preds$Realms[Preds$Realms=="Neartic"] <-"Nearctic"
    
  Preds$Trend_BT <- sin(Preds$Trend)^2*100
  Preds$CI_up_BT <- sin(Preds$CI_up)^2*100
  Preds$CI_low_BT <- sin(Preds$CI_low)^2*100
  
  if(logt){
    Preds$`Species richness` <- log(Preds$`Species richness`)
  }
  
  if(!allclasses){
    
    Preds <- Preds %>% 
      filter(Class %in% c("Reptiles", "Amphibians"))
    
  }

  if(facetting=="class"){
  
  if(backtransform){
    
    if(error){
      
      if(significance){
        
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Realms))+
          geom_line(size=1, aes(lty=significance)) +
          geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Realms, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Class, scales = "free")
      } else{
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Realms))+
        geom_line(size=1) +
        geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free")}
    }else{
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Realms))+
        geom_line(size=1) +
        #geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free") }
    
  } else{
    
    if(error){
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Realms))+
        geom_line(size=1) +
        geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free")  
    }else{
      p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Realms))+
        geom_line(size=1) +
        #geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Realms, col=NULL), alpha=0.1) +
        GGPoptions + facet_wrap(~Class, scales = "free")
    }
    
  }
  }
  
  
  if(facetting=="realm"){
    
    if(backtransform){
      if(error){
        if(significance){
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Class))+
          geom_line(size=1, aes(lty=Slope.Significance)) +
          scale_linetype_manual(name = "Slope significance", values=c("twodash", "solid")) +
          geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
        
        }else{
          p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Class))+
            geom_line(size=1) +
            geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Class, col=NULL), alpha=0.1) +
            GGPoptions + facet_wrap(~Realms, scales = "free")
        }
      }else{
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Class))+
          geom_line(size=1) +
          #geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
      }
      
    } else{
      if(error){
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Class))+
          geom_line(size=1) +
          geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
        
      }else{
        p <- ggplot(Preds, aes(x=`Species richness`, y=Trend, col=Class))+
          geom_line(size=1) +
          #geom_ribbon(aes(ymin=CI_low, ymax=CI_up, fill=Class, col=NULL), alpha=0.1) +
          GGPoptions + facet_wrap(~Realms, scales = "free")
        
      }
      
    }
  }
  
  if(facetting=="realmclass") {
    
    p <- ggplot(Preds, aes(x=`Species richness`, y=Trend_BT, col=Class))+
      geom_line() +
      # geom_ribbon(aes(ymin=CI_low_BT, ymax=CI_up_BT, fill=Class, col=NULL), alpha=0.1) +
      GGPoptions + facet_wrap(~Class+Realms, scales = "free")
  }
  
  p <- p + scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")
  
  
  return(p)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library(plyr)
library(dplyr)
library(raster)
library(rgdal)
library(ggplot2)
library(ggpubr)
library(spdep)
library(spatialreg)
library(foreach)
library(parallel)


## load data

# median at 50km resolution
# median_completeness_Birds <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianBirds.rds")
# median_completeness_Mammals <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianMammals.rds")
median_completeness_Reptiles <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianReptiles.rds")
median_completeness_Amphibians <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MedianAmphibians.rds")

# mean at 50km resolution
mean_completeness_Birds <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanBirds.rds")
mean_completeness_Mammals <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanMammals.rds")
mean_completeness_Reptiles <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanReptiles.rds")
mean_completeness_Amphibians <- readRDS("../../Results/SpatialAnalyses/MeanMedianVarianceCompletenessRasters_50k/MeanAmphibians.rds")

# SR at 50 km resolution
Richness_amphibians <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Amphibians_50k.rds")
Richness_reptiles <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Reptiles_50k.rds")
Richness_mammals <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Mammals_50k.rds")
Richness_birds <- readRDS("../../Results/SpatialAnalyses/SpeciesRichnessRasters/SR_Birds_50k.rds")


GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

# richness distribution
ValAmphibians <- values(Richness_amphibians) %>%  as.data.frame() %>%  setNames(., "Richness") %>%  filter(Richness!=0)
p1 <- ggplot(ValAmphibians, aes(Richness)) + geom_histogram(binwidth = 1) + 
  GGPoptions + geom_vline(xintercept=10, col="blue", lty="dashed") +
  geom_vline(xintercept=3, col="red", lty="dashed") #+ xlim(c(0, 20)) + ggtitle("Amphibians")

ValReptiles <- values(Richness_reptiles) %>%  as.data.frame() %>%  setNames(., "Richness") %>%  filter(Richness!=0)
p2 <- ggplot(ValReptiles, aes(Richness)) + geom_histogram(binwidth = 1) + 
  GGPoptions + geom_vline(xintercept=10, col="blue", lty="dashed") +
  geom_vline(xintercept=4, col="red", lty="dashed") + xlim(c(0, 20))+ ggtitle("Reptiles")

ValMammals <- values(Richness_mammals) %>%  as.data.frame() %>%  setNames(., "Richness") %>%  filter(Richness!=0)
p3 <- ggplot(ValMammals, aes(Richness)) + geom_histogram(binwidth = 1) + 
  GGPoptions + geom_vline(xintercept=10, col="blue", lty="dashed") +
  geom_vline(xintercept=3, col="red", lty="dashed") #+ xlim(c(0, 20))+ ggtitle("Mammals")

ValBirds <- values(Richness_birds) %>%  as.data.frame() %>%  setNames(., "Richness") %>%  filter(Richness!=0)
p4 <- ggplot(ValBirds, aes(Richness)) + geom_histogram(binwidth = 1) + 
  GGPoptions + geom_vline(xintercept=10, col="blue", lty="dashed") +
  geom_vline(xintercept=3, col="red", lty="dashed") #+ xlim(c(0, 20))+ ggtitle("Birds")

ggarrange(p1, p2, p3, p4)

# raster file for geograpical realms
realmraster <- readRDS("../../Results/Ecoregions_rasters/realms_50k.rds") 

# ## figure out the fraction for sample size -- going with 30% of data points in each realm.
# 
# R <- SampleSize(mean_completeness_Reptiles, Richness_reptiles, 2000, 3, realmraster, TRUE, 0.3)
# R %>% group_by(Realm) %>% 
#   summarise(Count=n()) %>% 
#   as.data.frame() %>% 
#   dplyr::select(Count) %>% 
#   sum()
# 
# A <- SampleSize(mean_completeness_Amphibians, Richness_amphibians, 2000, 3, realmraster, TRUE, 0.3)
# A %>% group_by(Realm) %>% 
#   summarise(Count=n()) %>% 
#   as.data.frame() %>% 
#   dplyr::select(Count) %>% 
#   sum()
# 
# M <- SampleSize(mean_completeness_Mammals, Richness_mammals, 2000, 3, realmraster, TRUE, 0.3)
# M %>% group_by(Realm) %>% 
#   summarise(Count=n()) %>% 
#   as.data.frame() %>% 
#   dplyr::select(Count) %>% 
#   sum()
# 
# B <- SampleSize(mean_completeness_Birds, Richness_birds, 2000, 3, realmraster, TRUE, 0.3)
# B %>% group_by(Realm) %>% 
#   summarise(Count=n()) %>% 
#   as.data.frame() %>% 
#   dplyr::select(Count) %>% 
#   sum()


## fitting the three models with median comp. at 50km resolution, with and interactions with realms 

## With species richness Cut off at 3 for sampling, and sample size of 30%  in each realm -- 

cl <- parallel::makeCluster(2)
clusterEvalQ(cl, {
  library(dplyr)
  library(raster)
  library(spatialreg)
  library(spdep)
})

doParallel::registerDoParallel(cl)

Start <- Sys.time()
Models <- foreach(MedianComp=c(median_completeness_Amphibians, median_completeness_Reptiles),
                  Richness=c(Richness_amphibians, Richness_reptiles),
                  SampleSize=rep(2000,2),
                  CutOff=rep(3,2),
                  realmrasterlist=c(realmraster, realmraster),
                  Interactions=rep(TRUE,2),
                  frac=rep(0.3,2),
                  .verbose = TRUE) %dopar% {
                                        FitSpatialModel(Raster_richness=Richness, 
                                                        Raster_completeness=MedianComp,
                                                        samplesize=SampleSize, 
                                                        cutoff=CutOff,
                                                        realmraster=realmrasterlist,
                                                        Realm_Inter=Interactions,
                                                        frac=frac)
                  }
End <- Sys.time()
stopCluster(cl)
print(Start-End)

saveRDS(Models, "../../Results/SpatialAnalysesModels/Median_Amphibians_Reptiles_cutoff3.rds")

# # # # # # # Models on mean

cl <- parallel::makeCluster(4)
clusterEvalQ(cl, {
  library(dplyr)
  library(raster)
  library(spatialreg)
  library(spdep)
})

doParallel::registerDoParallel(cl)

Start <- Sys.time()
Models <- foreach(MeanComp=c(mean_completeness_Amphibians, mean_completeness_Birds , mean_completeness_Mammals ,mean_completeness_Reptiles),
                  Richness=c(Richness_amphibians, Richness_birds, Richness_mammals, Richness_reptiles),
                  SampleSize=rep(2000,4),
                  CutOff=rep(3,4),
                  realmrasterlist=c(realmraster, realmraster, realmraster, realmraster),
                  Interactions=rep(TRUE,4),
                  frac=rep(0.3,4),
                  .verbose = TRUE) %dopar% {
                    FitSpatialModel(Raster_richness=Richness, 
                                    Raster_completeness=MeanComp,
                                    samplesize=SampleSize, 
                                    cutoff=CutOff,
                                    realmraster=realmrasterlist,
                                    Realm_Inter=Interactions,
                                    frac=frac)
                  }
End <- Sys.time()
stopCluster(cl)
print(Start-End)

saveRDS(Models, "../../Results/SpatialAnalysesModels/Mean_cutoff3.rds")


# # # # # Analysis

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


## median models

Models <- readRDS("../../Results/SpatialAnalysesModels/Median_Amphibians_Reptiles_cutoff3.rds")
ModelReptile <- Models[[2]]
ModelAmphibian <- Models[[1]]

ggplot(ModelReptile$data, aes(SpeciesRichness)) + geom_histogram() + xlim(c(0,10))
ggplot(ModelAmphibian$data, aes(SpeciesRichness)) + geom_histogram() + xlim(c(0,10))

ggplot(ModelAmphibian$data, aes(log(SpeciesRichness), asin(sqrt(Completeness)))) + 
  geom_point() + facet_wrap(~Realm) + theme_bw()

ggplot(ModelReptile$data, aes(log(SpeciesRichness), asin(sqrt(Completeness)))) + 
  geom_point() + theme_bw() #+ facet_wrap(~Realm) + theme_bw()



#ModelAmphibian <- readRDS("../../Results/SpatialAnalysesModels/Amphibians_median.rds")

ModelsMean <- readRDS("../../Results/SpatialAnalysesModels/Mean_cutoff3.rds")


# mean: quadratic term fits better
anova(ModelsMean[[1]]$model1, ModelsMean[[1]]$model2)
anova(ModelsMean[[2]]$model1, ModelsMean[[2]]$model2)
anova(ModelsMean[[3]]$model1, ModelsMean[[3]]$model2)
anova(ModelsMean[[4]]$model1, ModelsMean[[4]]$model2)

summary(ModelsMean[[1]]$model1)
summary(ModelsMean[[2]]$model1)
summary(ModelsMean[[3]]$model1)
summary(ModelsMean[[4]]$model1)

summary(ModelsMean[[1]]$model2)
summary(ModelsMean[[2]]$model2)
summary(ModelsMean[[3]]$model2)
summary(ModelsMean[[4]]$model2)

## plotting models trends
# PredsReptiles <- InSamplePredictions(ModelsMed[[2]]$model2, mod.data = Models[[2]]$data, Inter.realms = TRUE)
# PlotInSamplePreds(PredsReptiles, TRUE)
# 
# PredsAmphibians <- InSamplePredictions(ModelsMed[[1]]$model1, mod.data = ModelsMed[[1]]$data, Inter.realms = TRUE)
# PlotInSamplePreds(PredsAmphibians, TRUE)
 

# plot for amphibians and reptiles 
library(viridis)

PredsAmphibians <- Preds_Summary_NoQuad(ModelAmphibian$model2, ModelAmphibian$data, 3, "Amphibians", TRUE)
PredsReptiles <- Preds_Summary_NoQuad(ModelReptile$model2, ModelReptile$data, 3, "Reptiles", TRUE)

p <- PlotPredsv2(PredsAmphibians,
            PredsReptiles, 
            PredsMammals = NULL,
            PredsBirds = NULL,
            allclasses = FALSE,
            "realm",
            error=TRUE,
            backtransform=TRUE,
            significance = TRUE, 
            logt = FALSE) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  ylab("Trend in completeness (non-spatial smooth)")

ggsave(p, filename="../../Results/SpatialAnalysesModels/Plots/Herptiles2.pdf", width=9, height=5)

x <- summary(ModelAmphibian$model2)
stargazer::stargazer(x$coefficients, digits = 4)

x <- summary(ModelReptile$model2)
stargazer::stargazer(x$Coef, digits = 4)


## means
PredsAmphibiansMean <- Preds_Summary_NoQuad(ModelsMean[[1]]$model2, ModelsMean[[1]]$data, 3, "Amphibians", TRUE)
PredsReptilesMean <- Preds_Summary_NoQuad(ModelsMean[[4]]$model2, ModelsMean[[4]]$data, 3, "Reptiles", TRUE)
PredsMammalsMean <- Preds_Summary_NoQuad(ModelsMean[[3]]$model2, ModelsMean[[3]]$data, 3, "Mammals", TRUE)
PredsBirdsMean <- Preds_Summary_NoQuad(ModelsMean[[2]]$model2, ModelsMean[[2]]$data, 3, "Birds", TRUE)

pmean <- PlotPredsv2(PredsAmphibiansMean,
                 PredsReptilesMean, 
                 PredsMammals = PredsMammalsMean,
                 PredsBirds = PredsBirdsMean,
                 allclasses = TRUE,
                 "realm",
                 error=TRUE,
                 backtransform=TRUE,
                 significance = TRUE, 
                 logt = FALSE) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  ylab("Trend in completeness (non-spatial smooth)")


pmean <- PlotPredsv2(PredsAmphibiansMean,
                     PredsReptilesMean, 
                     PredsMammals = PredsMammalsMean,
                     PredsBirds = PredsBirdsMean,
                     allclasses = FALSE,
                     "realm",
                     error=TRUE,
                     backtransform=TRUE,
                     significance = TRUE, 
                     logt = FALSE) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  ylab("Trend in completeness (non-spatial smooth)")

ggsave(pmean, filename="../../Results/SpatialAnalysesModels/Plots/RegressionLinesMean.pdf", width=9, height=5)



## NB: errors around the lines are standard errors (not 95 ci)


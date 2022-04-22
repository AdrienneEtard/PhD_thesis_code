library(plyr)
library(dplyr)
library(raster)
library(rgdal)
library(ggplot2)
library(ggpubr)
library(spdep)
library(spatialreg)
# library(sppois)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


`%nin%` <- Negate(`%in%`)
## see https://github.com/tidyverse/ggplot2/wiki/plotting-polygon-shapefiles

Amphibians <- read.csv("../../Results/SR_vs_completeness/Amphibians_to_plot_v2.csv")
Reptiles <- read.csv("../../Results/SR_vs_completeness/Reptiles_to_plot_v2.csv")
Birds <- read.csv("../../Results/SR_vs_completeness/Birds_to_plot_v2.csv")
Mammals <- read.csv("../../Results/SR_vs_completeness/Mammals_to_plot_v2.csv")

## load mean completeness
mean_completeness_Amphibians <- readRDS("../../Results/Maps_completeness/mean_completeness/mean_Amphibians.rds")
mean_completeness_Birds <- readRDS("../../Results/Maps_completeness/mean_completeness/mean_Birds.rds")
mean_completeness_Mammals <- readRDS("../../Results/Maps_completeness/mean_completeness/mean_Mammals.rds")
mean_completeness_Reptiles <- readRDS("../../Results/Maps_completeness/mean_completeness/mean_Reptiles.rds")

## species richness
TotalRichness <- readRDS("../../Results/TotalSpeciesRichness.rds")
Richness_class <- readRDS("../../Results/SpeciesRichnessByClass.rds")
Richness_amphibians <- Richness_class[["Amphibians"]]
Richness_reptiles <- Richness_class[["Reptiles"]]
Richness_mammals <- Richness_class[["Mammals"]]
Richness_birds <- Richness_class[["Birds"]]

## to run to generate Amphibians, Birds, Reptiles and Mammals

GetValues <- function(TR, Mean_comp, realms) {

  # replace -1 by NA in the completeness raster (-1 are NAs for terrestrial lands, as opposed to NA in the oceans and seas)
  y <- which(values(Mean_comp)==-1)
  values(Mean_comp)[y] <- NA

  # dataframe to store values of the rasters
  x <- values(TR) %>%
    as.data.frame() %>%
    setNames(., "SR")
  x$Comp <- values(Mean_comp)
  x$Realm <- values(realms)

  # detail realms
  x <- x %>%
    mutate(Realm=ifelse(Realm==1, "Afrotropic",
                        ifelse(Realm==2, "Antarctic",
                               ifelse(Realm==3, "Australasia",
                                      ifelse(Realm==4, "Indo-Malay",
                                             ifelse(Realm==5, "Nearctic",
                                                    ifelse(Realm==6, "Neotropic",
                                                           ifelse(Realm==7, "Oceania",
                                                                  ifelse(Realm==8, "Palearctic", Realm)))))))))

  return(x)
}

## load realm shapefiles
x = readOGR(dsn="../../Data/Ecoregions", layer="tnc_terr_ecoregions")
plot(x)
# # # # # PLOTS

## plotting the realms from the shapefile
x@data$id = rownames(x@data)
x.points = ggplot2::fortify(x, region = "id")
x.df = join(x.points, x@data, by="id")

ggplot(x.df) +
  aes(long,lat,group=group,fill=WWF_REALM2) +
  geom_polygon() +
  coord_equal()


## convert shapefile to raster with the same number of cells than the other rasters
r <- raster(ncol=ncol(mean_completeness_Amphibians), nrow=nrow(mean_completeness_Amphibians))
extent(r) <- extent(x)
ra <- rasterize(x, r, "WWF_REALM2")

## project rasters (TotalRichness and mean completeness rasters) to equal area projection
proj4string(TotalRichness) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
proj4string(Richness_amphibians) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
proj4string(Richness_birds) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
proj4string(Richness_mammals) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
proj4string(Richness_reptiles) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

proj4string(mean_completeness_Amphibians) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
proj4string(mean_completeness_Birds) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
proj4string(mean_completeness_Mammals) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')
proj4string(mean_completeness_Reptiles) <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

## now format the realm raster so that it is in the same format as TotalRichness
realms <- projectRaster(from=ra, to=TotalRichness, method="ngb")
plot(realms)
plot(TotalRichness)

# ## verify that the rasters have the same format (ncell, extent, resolution...)
# compareRaster(realms, TotalRichness)
# compareRaster(realms, mean_completeness_Amphibians)
# compareRaster(TotalRichness, mean_completeness_Amphibians)

## Now, plot TR against mean completeness and color the point by realm
Amphibians <- GetValues(Richness_amphibians, mean_completeness_Amphibians, realms) %>% mutate(Class="Amphibians")
Birds <- GetValues(Richness_birds, mean_completeness_Birds, realms) %>% mutate(Class="Birds")
Mammals <- GetValues(Richness_mammals, mean_completeness_Mammals, realms) %>% mutate(Class="Mammals")
Reptiles <- GetValues(Richness_reptiles, mean_completeness_Reptiles, realms) %>% mutate(Class="Reptiles")

# save these files
write.csv(Amphibians, "../../Results/SR_vs_completeness/Amphibians_to_plot_v2.csv", row.names=FALSE)
write.csv(Reptiles, "../../Results/SR_vs_completeness/Reptiles_to_plot_v2.csv", row.names=FALSE)
write.csv(Birds, "../../Results/SR_vs_completeness/Birds_to_plot_v2.csv", row.names=FALSE)
write.csv(Mammals, "../../Results/SR_vs_completeness/Mammals_to_plot_v2.csv", row.names=FALSE)


## plotting

To_plot <- rbind(Amphibians, Birds, Mammals, Reptiles)
To_plot$Class <- factor(To_plot$Class, levels=c("Mammals", "Birds", "Amphibians", "Reptiles"), labels = c("(a) Mammals", "(b) Birds", "(c) Amphibians", "(d) Reptiles"))

#To_plot$Realm[is.na(To_plot$Realm)] <- "Unknown"

p <-  ggplot(To_plot, aes(SR, Comp, colour=Realm)) + geom_point() + #alpha=0.5
  GGPoptions + facet_wrap(Class ~ ., ncol=2, scales="free") +
  xlab("Species richness") + ylab("Mean trait completeness") 

Colors <- c("#F0E442", "#CC79A7", "#cc0000", "#D55E00",
            "#56B4E9", "#009E73", "#000066","#FF6666")#, "#b3b3b3")
names(Colors) <- c("Afrotropic", "Antarctic", "Australasia", "Indo-Malay", "Nearctic", "Neotropic", "Oceania", "Palearctic")
  
p <- p + scale_colour_manual(name="Realm", values = Colors, na.value="grey")
p
ggsave(p, filename="../Results/SR_vs_completeness/All_classes.pdf", width = 7, height=5)
ggsave(p, filename="../Results/SR_vs_completeness/All_classes.png", dpi=600,  width = 7, height=5)

# mean completeness against species richness
# pdf(file="../../Results/SR_Mean_completeness.pdf", width=8, height=5.5, family="Times", pointsize=15)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,1,2,1), oma=c(1,1,1,1))
# par(mfrow=c(2,2))

## arcsin + square root transformation
ggplot(To_plot, aes(Comp/100)) + geom_histogram() +
  GGPoptions + facet_wrap(Class ~ ., ncol=2, scales="free") 
 
p <-  ggplot(To_plot, aes(asin(sqrt(Comp/100)))) + geom_histogram() +
  GGPoptions + facet_wrap(Class ~ ., ncol=2, scales="free") 


# # # # # # # #  MODELS: spatial simultaneous autoregressive lag models

FitSpatialModel <- function(Raster_mean_completeness, Raster_richness) {

  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"),
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))

  # for debugging only
  # Raster_mean_completeness = mean_completeness_Amphibians
  # Raster_richness = Richness_amphibians

  # TODO model without quadratic term and null model - comparisons. loglikelyhood ratio

  # get the data from the raster layers, and resample (for a smaller dataset with 10000 datapoints)
  mod.data <- data.frame(Completeness=values(Raster_mean_completeness),
                         SpeciesRichness=values(Raster_richness),
                         Longitude=coordinates(Raster_richness)[,1],
                         Latitude=coordinates(Raster_richness)[,2])

  mod.data$Completeness[mod.data$Completeness==-1] <- NA
  mod.data <- na.omit(mod.data)
  mod.data$Completeness <- mod.data$Completeness / 100
  mod.data <- mod.data[sample(x = 1:nrow(mod.data),size = 10000, replace = FALSE), ]

  # fit a spatial simultaneous autoregressive lag model

  # coordinates as neighbours "The function uses the deldir package to
  # convert a matrix of two-dimensional coordinates into a neighbours list of class nb
  # with a list of integer vectors containing neighbour region number ids."
  nb <- spdep::tri2nb(coords = mod.data[,c('Longitude','Latitude')])
  lw <- spdep::nb2listw(neighbours = nb)

  # fit the model
  Model1 <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~
                                  log(mod.data$SpeciesRichness)+
                                  I(log(mod.data$SpeciesRichness)^2),
                                listw = lw)

  # # model without quadratic term
  Model2 <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~
                              log(mod.data$SpeciesRichness),
                              listw = lw)
  #
  # # Null model
  Model0 <- spatialreg::lagsarlm(asin(sqrt(mod.data$Completeness)) ~ 1, listw = lw)

  # predictions for Model1 using predict.sarlm of the spatialreg package
  # fit  = trend + signal
  # trend is the non-spatial smooth
  # signal is the spatial smooth

  predictions <- spatialreg::predict.sarlm(object = Model1,
                                           newdata=NULL,
                                           listw = lw,
                                           zero.policy = TRUE)

  Spatial_preds <-  attr(predictions, "trend") +  attr(predictions, "signal") %>%
    as.data.frame() %>%
    setNames(., "preds")
  Spatial_preds$SR <- mod.data$SpeciesRichness

  # plotting fit and trend (/!\ Need to be back-transformed)

  # back-transforming
  Spatial_preds$preds_BT <- (sin(Spatial_preds$preds))^2

  # the fit = trend + signal
  p1 <- ggplot(Spatial_preds, aes(SR, preds_BT)) +
    geom_point() +
    GGPoptions +
    xlab("Species richness") +
    ylab("Fit (trend + signal)")

  # the trend (non-spatial smooth)
  Spatial_preds$trend <-  attr(predictions, "trend")
  Spatial_preds$trend_BT <- sin(Spatial_preds$trend)^2
  p2 <- ggplot(Spatial_preds, aes(SR, trend_BT)) +
    geom_line() +
    GGPoptions +
    xlab("Species richness") +
    ylab("Non-spatial smooth (trend)")

  ## Model comparison
  # anova and log likelyhood ratio test
  Anova_results_10 <- anova(Model0,Model1)
  Anova_results_12 <- anova(Model2,Model1)
  Anova_results_20 <- anova(Model0,Model2)

  # return
  #return(list(Model=Model, NullModel=Null_Model, Anova_results=Anova_results))

 # browser()
  
  return(list(Plot_fit=p1,
              Plot_trend=p2,
              Model_summary_1=summary(Model1, Nagelkerke=TRUE),
              Model_summary_2=summary(Model2, Nagelkerke=TRUE),
              Model_summary_0=summary(Model0, Nagelkerke=TRUE),
              Anova10=Anova_results_10,
              Anova12=Anova_results_12,
              Anova20=Anova_results_20))

  }
 
 
## Fit models
Start <- Sys.time()

Spatial_model_amphibians <- FitSpatialModel(mean_completeness_Amphibians, Richness_amphibians)
Spatial_model_reptiles <- FitSpatialModel(mean_completeness_Reptiles, Richness_reptiles)
Spatial_model_mammals <- FitSpatialModel(mean_completeness_Mammals, Richness_mammals)
Spatial_model_birds <- FitSpatialModel(mean_completeness_Birds, Richness_birds)

End <- Sys.time()

## save results
saveRDS(Spatial_model_amphibians, "../../Results/SR_vs_completeness/Model_outputs/Amphibians.rds")
saveRDS(Spatial_model_reptiles, "../../Results/SR_vs_completeness/Model_outputs/Reptiles.rds")
saveRDS(Spatial_model_mammals, "../../Results/SR_vs_completeness/Model_outputs/Mammals.rds")
saveRDS(Spatial_model_birds, "../../Results/SR_vs_completeness/Model_outputs/Birds.rds")

p= ggarrange(
Spatial_model_mammals$Plot_trend + xlab("") + ggtitle("(a) Mammals"),
Spatial_model_birds$Plot_trend + xlab("") + ylab("") + ggtitle("(b) Birds"),
Spatial_model_amphibians$Plot_trend + ggtitle("(c) Amphibians"),
Spatial_model_reptiles$Plot_trend + ylab("") + ggtitle("(d) Reptiles")
)
# 
# ggsave(p, filename="../../Results/SR_vs_completeness/Model_outputs/Plots_trends.pdf", width = 7, height=5)

## Load results
Spatial_mod_amphibians <- readRDS("../../Results/SR_vs_completeness/Model_outputs/Amphibians.rds")
Spatial_mod_birds <- readRDS("../../Results/SR_vs_completeness/Model_outputs/Birds.rds")
Spatial_mod_reptiles <- readRDS("../../Results/SR_vs_completeness/Model_outputs/Reptiles.rds")
Spatial_mod_mammals <- readRDS("../../Results/SR_vs_completeness/Model_outputs/Mammals.rds")

## Pseudo R-squared
Spatial_mod_birds[[3]]$NK
Spatial_mod_mammals[[3]]$NK
Spatial_mod_amphibians[[3]]$NK
Spatial_mod_reptiles[[3]]$NK


library(stargazer)

# summaries for the main model
stargazer(Spatial_mod_mammals$Model_summary_1$Coef)
stargazer(Spatial_mod_birds$Model_summary_1$Coef)
stargazer(Spatial_mod_amphibians$Model_summary_1$Coef)
stargazer(Spatial_mod_reptiles$Model_summary_1$Coef)




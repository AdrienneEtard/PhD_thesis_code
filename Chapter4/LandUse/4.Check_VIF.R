## check VIF among traits and range size as well 

#setwd("E:/3.Explanatory_traits/Code/")

#########################################################

library(dplyr)
library(car)
library(stargazer)

stepwise.vif <- function (dataset,
                          metrics,
                          vif.threshold = 5,
                          verbose = F){
  
  #browser()
  
  dataset$dummy <- rnorm(nrow(dataset)) # generates a dummy normal distribution
  output <- metrics # those are all the candidate variables
  step.count <- 1
  output.results <- list() 
  
  repeat {
    
    # VIF scores for the linear model: dummy ~ all variables; returns VIF or GVIF, depending on Degrees of freedom (if any variables has more than one then GVIF is computed rather than VIF)
    vif.scores <- vif(lm(as.formula(paste0(
      "dummy~", paste0(output,
                       collapse = "+")
    )), data = dataset))
    na.coefficients <- Reduce('|', is.nan(vif.scores)) # for scores that are NA -- if so stop the function
    if (na.coefficients) {
      stop("NA coefficient in a regression model.")
    }
    
    # Select VIF scores
    vif.scores <- vif.scores[,1] 
    
    ## output.results stores VIF calculated at each step and output stores variable names that are selected 
    output.results[[step.count]] <-
      sort(vif.scores, decreasing = F) # sort VIF scores
    
    vif.scores <- vif.scores[vif.scores >= vif.threshold]
    
    ## If all VIF scores are under the threshold values then keep all variables and stop the function. Else drop the variable
    # that has a VIF above the threshold and proceed iteratively
    if (length(vif.scores) == 0)
      break
    
    drop.var <-
      names(vif.scores[vif.scores == max(vif.scores)]) # select variable names that have VIF score more than the threshold
    if (verbose) {
      print(paste0(
        "Step ",
        step.count,
        " - Exclude ",
        drop.var,
        " (VIF = ",
        max(vif.scores),
        ")"
      ))
    }
    step.count <- step.count + 1
    output <- output[!output %in% drop.var]
  }
  
  names(output.results) <- paste0("Iteration ", 1:step.count)
  names(output.results)[length(output.results)] <- "Final"
  
  if(length(output.results)==1) {output.results <- as.data.frame(output.results)}
  
  return(list(Selected_vars=output, VIF=output.results))
}

#########################################################

#setwd("D:/PhD/PhD_R_work/3.Explanatory_traits/Code")

## range sizes
Index <- read.csv("D:/PhD/PhD_R_projects/Range_maps_work/Results/7_GenerateSpeciesIndices/SpeciesIndices.csv")
RS <- read.csv("../Data/Range_sizes.csv")
RS$Species <- paste0("sp", RS$Index)
RS <- RS %>% 
  filter(Range_area_sq_km!=0)
RS <- left_join(RS, Index)
colnames(RS)[5] <- "Best_guess_binomial"

## combine traits and range sizes
Traits <- readRDS("../Results/Imputed_traits_transformed.rds")[[8]]
Traits <- left_join(Traits, RS)
colnames(Traits)
Traits$log10_Range_area <- log10(Traits$Range_area_sq_km)

## load PREDICTS and combine with trait data
Predicts <- readRDS("../Results/Predicts_merged_sites.rds")
levels(Predicts$Predominant_land_use)
Predicts$Occurrence <- ifelse(Predicts$Measurement > 0, 1, 0)
Predicts$LandUse <- as.character(Predicts$Predominant_land_use)
Predicts$LandUse[Predicts$LandUse=="Cannot decide"] <- NA
#Predicts$LandUse[Predicts$LandUse=="Secondary vegetation (indeterminate age)"] <- NA
Predicts$LandUseGrouped <- as.character(Predicts$LandUse)

Predicts$LandUse <- factor(Predicts$LandUse,
                           levels=c("Primary vegetation", 
                                    "Mature secondary vegetation",
                                    "Intermediate secondary vegetation",
                                    "Young secondary vegetation",
                                    "Plantation forest" ,
                                    "Pasture", 
                                    "Cropland",
                                    "Urban"))
#Predicts <- subset(Predicts, !is.na(LandUse))

Predicts$Use_intensity <- as.character(Predicts$Use_intensity)
Predicts$Use_intensity[Predicts$Use_intensity=="Cannot decide"] <- NA
Predicts$Use_intensity <- factor(Predicts$Use_intensity, levels=c("Minimal use", "Light use", "Intense use"))

Predicts$LandUseGrouped[Predicts$LandUseGrouped=="Mature secondary vegetation"|
                   Predicts$LandUseGrouped=="Intermediate secondary vegetation"|
                   Predicts$LandUseGrouped=="Young secondary vegetation"|
                   Predicts$LandUseGrouped=="Secondary vegetation (indeterminate age)"] <- "Secondary vegetation"
Predicts$LandUseGrouped[Predicts$LandUseGrouped=="Pasture"|
                   Predicts$LandUseGrouped=="Cropland"] <- "AGR"

Predicts$LandUseGrouped <- factor(Predicts$LandUseGrouped,
                           levels=c("Primary vegetation", 
                                    "Secondary vegetation",
                                    "Plantation forest" ,
                                    "AGR",
                                    "Urban"))

levels(Predicts$LandUseGrouped )

any(is.na(Predicts$Use_intensity))
any(is.na(Predicts$LandUseGrouped))

Predicts$LandUseGrouped %>% table()
Predicts$LandUseGrouped[is.na(Predicts$LandUseGrouped)] %>% nrow()

Predicts <- subset(Predicts, !is.na(Use_intensity))
Predicts <- subset(Predicts, !is.na(LandUseGrouped))

## combine with the trait data
Predicts <- left_join(Predicts, Traits[, -c(1,2,3,4,15,17)], by="Best_guess_binomial")
colnames(Predicts)

## set primary diet levels
Predicts$Primary_diet <- as.character(Predicts$Primary_diet)
Predicts$Primary_diet[Predicts$Primary_diet=="SE"|Predicts$Primary_diet=="PL"] <- "PL|SE"
Predicts$Primary_diet[Predicts$Primary_diet=="FR"|Predicts$Primary_diet=="NE"] <- "FR|NE"


# ## center and scale continuous variables (across classes)
# Predicts_zscores_all <- Predicts
# Predicts_zscores_all$log10_Body_mass_g <- scale(Predicts_zscores_all$log10_Body_mass_g, center = TRUE, scale = TRUE)
# Predicts_zscores_all$log10_Lifespan_proxy <- scale(Predicts_zscores_all$log10_Lifespan_proxy, center = TRUE, scale = TRUE)
# Predicts_zscores_all$log10_Litter_size <- scale(Predicts_zscores_all$log10_Litter_size, center = TRUE, scale = TRUE)
# Predicts_zscores_all$log10_Range_area <- scale(Predicts_zscores_all$log10_Range_area, center = TRUE, scale = TRUE)
# Predicts_zscores_all$sqrt_Habitat_breadth_IUCN <- scale(Predicts_zscores_all$sqrt_Habitat_breadth_IUCN, center = TRUE, scale = TRUE)
# Predicts_zscores_all$sqrt_Diet_breadth <- scale(Predicts_zscores_all$sqrt_Diet_breadth, center = TRUE, scale = TRUE)

# ## check VIF scores across classes
# Predicts_zscores_all$Diel_activity %>%  unique
# Predictors <-  c("log10_Body_mass_g",
#                  "log10_Lifespan_proxy",
#                  "log10_Litter_size",
#                  "log10_Range_area",
#                  "sqrt_Habitat_breadth_IUCN",
#                  "sqrt_Diet_breadth",
#                  "Specialisation", 
#                  "Diel_activity",
#                  "Trophic_level", 
#                  "Primary_diet", 
#                  "LandUse",
#                  "Use_intensity")
# 
# stepwise.vif(dataset = Predicts_zscores_all, metrics = Predictors, vif.threshold = 5, verbose = TRUE)
# 
# saveRDS(Predicts_zscores_all, "../Results/Model_data_allvertebrates.rds")

## center and scale continuous variables (within classes) ## don't scale or center
# Predicts_zscores_classes <- Predicts
# Predicts_zscores_classes <- Predicts_zscores_classes %>% 
#   group_by(Class) %>% 
#   mutate(log10_Body_mass_g=scale(log10_Body_mass_g, center = TRUE, scale = TRUE),
#          log10_Lifespan_proxy=scale(log10_Lifespan_proxy, center = TRUE, scale = TRUE),
#          log10_Litter_size=scale(log10_Litter_size, center = TRUE, scale = TRUE),
#          log10_Range_area=scale(log10_Range_area, center = TRUE, scale = TRUE),
#          sqrt_Habitat_breadth_IUCN=scale(sqrt_Habitat_breadth_IUCN, center = TRUE, scale = TRUE),
#          sqrt_Diet_breadth=scale(sqrt_Diet_breadth, center = TRUE, scale = TRUE))
# var(Predicts_zscores_classes$log10_Body_mass_g[Predicts_zscores_classes$Class=="Aves"])
# var(Predicts_zscores_classes$log10_Body_mass_g[Predicts_zscores_classes$Class=="Reptilia"])


## check VIF scores within classes
Predictors <-  c("log10_Body_mass_g",
                 "log10_Lifespan_proxy",
                 "log10_Litter_size",
                 "log10_Range_area",
                 "sqrt_Habitat_breadth_IUCN",
                 "sqrt_Diet_breadth",
                 "Specialisation", 
                 "Diel_activity",
                 #"Trophic_level", 
                 "Primary_diet", 
                 "LandUseGrouped",
                 "Use_intensity")

#stepwise.vif(dataset = Predicts, metrics = Predictors, vif.threshold = 5, verbose = TRUE)

##### for mammals -- all predictors can be included
stepwise.vif(dataset = Predicts[Predicts$Class=="Mammalia",], 
             metrics = Predictors, vif.threshold = 5, verbose = TRUE)$VIF %>% 
  stargazer(summary = FALSE, digits = 1)

##### for birds -- all predictors can be included
stargazer(stepwise.vif(dataset = Predicts[Predicts$Class=="Aves",], 
                       metrics = Predictors, vif.threshold = 5, verbose = TRUE)$VIF,
          summary = FALSE, digits = 1)

##### for reptiles -- exclude primary diet
Predicts$Primary_diet[Predicts$Class=="Reptilia"] %>%  table()

X <- stepwise.vif(dataset = Predicts[Predicts$Class=="Reptilia",], 
                  metrics = Predictors, vif.threshold = 5, verbose = TRUE)$VIF$`Iteration 1` %>% 
  as.data.frame()
stargazer(X,
          summary = FALSE, digits = 1)

X <- stepwise.vif(dataset = Predicts[Predicts$Class=="Reptilia",], 
                  metrics = Predictors, vif.threshold = 5, verbose = TRUE)$VIF$Final %>% 
  as.data.frame()
stargazer(X,
          summary = FALSE, digits = 1)

##### for amphibians -- exclude primary diet

X <- stepwise.vif(dataset = Predicts[Predicts$Class=="Amphibia",], 
                  metrics = Predictors, vif.threshold = 5, verbose = TRUE)$VIF$`Iteration 1` %>% 
  as.data.frame()
stargazer(X,
          summary = FALSE, digits = 1)

X <- stepwise.vif(dataset = Predicts[Predicts$Class=="Amphibia",], 
                  metrics = Predictors, vif.threshold = 5, verbose = TRUE)$VIF$Final %>% 
  as.data.frame()
stargazer(X,
          summary = FALSE, digits = 1)

saveRDS(Predicts, "../Results/Model_data_within_class.rds")





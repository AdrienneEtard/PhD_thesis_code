## Pre-process traits (choose the 8th dataset because it is the one I used in the FD analyses)
library(dplyr)
`%nin%` <- Negate(`%in%`)

## Function to bind and transform the 8 imputed trait datasets
## Function to transform  results
Transform <- function(TraitDF, Trait, Transf) {
  
  if (Transf=="log10"){
    TraitDF[,paste("log10", Trait, sep="_")] <- as.numeric(log10(TraitDF[,Trait]))
    #TraitDF[, paste("log10", Trait, sep="_")] <- scale(TraitDF[, paste("log10", Trait, sep="_")], center=TRUE, scale=TRUE)
    #TraitDF[, paste("log10", Trait, sep="_")] <- as.numeric(TraitDF[, paste("log10", Trait, sep="_")])
  }
  
  if(Transf=="sqrt") {
    TraitDF[,Trait] <- as.numeric(TraitDF[,Trait])
    TraitDF[, paste("sqrt", Trait, sep="_")]  <- as.numeric(sqrt(TraitDF[,Trait]))
    #TraitDF[, paste("sqrt", Trait, sep="_")] <- scale(TraitDF[,Trait] , center=TRUE, scale=TRUE)
    #TraitDF[, paste("sqrt", Trait, sep="_")] <- as.numeric(TraitDF[, paste("sqrt", Trait, sep="_")])
    
  }
  
  return(TraitDF)
}

Bind_transform_traits <- function(List_imputed_datasets) {
  
  Candidate_traits <-  c("Body_mass_g",
                         "Lifespan_proxy",
                         "Litter_size",
                         "Habitat_breadth_IUCN",
                         "Specialisation",
                         "Diel_activity",
                         "Trophic_level",
                         "Primary_diet",
                         "Diet_breadth")
  
  Candidate_traits_transformed <-  c("log10_Body_mass_g",
                                     "log10_Lifespan_proxy",
                                     "log10_Litter_size",
                                     "sqrt_Habitat_breadth_IUCN",
                                     "Specialisation",
                                     "Diel_activity",
                                     "Trophic_level",
                                     "Primary_diet", 
                                     "sqrt_Diet_breadth")
  
  Taxo <- c("Class", "Order", "Family", "Genus", "Best_guess_binomial")
  
  List_datasets <- list()
  
  for (i in 1:8) {
  
    Imputed <- rbind(List_imputed_datasets[[i]][["A"]]$Imputed.Dataset[, c(Candidate_traits, Taxo)],
                     List_imputed_datasets[[i]][["B"]]$Imputed.Dataset[, c(Candidate_traits, Taxo)],
                     List_imputed_datasets[[i]][["M"]]$Imputed.Dataset[, c(Candidate_traits, Taxo)],
                     List_imputed_datasets[[i]][["R"]]$Imputed.Dataset[, c(Candidate_traits, Taxo)]
                     )
    
    Imputed <- Transform(Imputed, "Body_mass_g", "log10")
    Imputed <- Transform(Imputed, "Lifespan_proxy", "log10")
    Imputed <- Transform(Imputed, "Litter_size", "log10")
    Imputed <- Transform(Imputed, "Habitat_breadth_IUCN", "sqrt")
    Imputed <- Transform(Imputed, "Diet_breadth", "sqrt")
    
    Traits <- Imputed[, c(Taxo, Candidate_traits_transformed)]
    
    List_datasets[[i]] <- Traits
    
  }
  
  return(List_datasets)
}


# # # # # # 

Imputed <- readRDS("../../Results/Imputed_traits/List_of_8_sets.rds")

Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")
Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
             "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")

## Add class
## Create Habitat breadth and diet breadth as sum of habitat and diet related variables (not anymore!)

for(i in 1:8){
  
  Imputed[[i]][["M"]]$Imputed.Dataset$Class <- "Mammals"
  colnames(Imputed[[i]][["M"]]$Imputed.Dataset)[8] <- "Lifespan_proxy"
  # Imputed[[i]][["M"]]$Imputed.Dataset$Diet_breadth <- apply( Imputed[[i]][["M"]]$Imputed.Dataset[, Diet], 1, sum, na.rm=TRUE)
  # Imputed[[i]][["M"]]$Imputed.Dataset[, Habitat] <- apply( Imputed[[i]][["M"]]$Imputed.Dataset[, Habitat], 2, as.numeric)
  # Imputed[[i]][["M"]]$Imputed.Dataset$Habitat_breadth_IUCN <- apply( Imputed[[i]][["M"]]$Imputed.Dataset[, Habitat], 1, sum, na.rm=TRUE)
 
  Imputed[[i]][["B"]]$Imputed.Dataset$Class <- "Birds"
  colnames(Imputed[[i]][["B"]]$Imputed.Dataset)[9] <- "Lifespan_proxy"
  # Imputed[[i]][["B"]]$Imputed.Dataset$Diet_breadth <- apply( Imputed[[i]][["B"]]$Imputed.Dataset[, Diet], 1, sum, na.rm=TRUE)
  # Imputed[[i]][["B"]]$Imputed.Dataset[, Habitat] <- apply( Imputed[[i]][["B"]]$Imputed.Dataset[, Habitat], 2, as.numeric)
  # Imputed[[i]][["B"]]$Imputed.Dataset$Habitat_breadth_IUCN <- apply( Imputed[[i]][["B"]]$Imputed.Dataset[, Habitat], 1, sum, na.rm=TRUE)
  
  Imputed[[i]][["R"]]$Imputed.Dataset$Class <- "Reptiles"
  colnames(Imputed[[i]][["R"]]$Imputed.Dataset)[9] <- "Lifespan_proxy"
  # Imputed[[i]][["R"]]$Imputed.Dataset$Diet_breadth <- apply( Imputed[[i]][["R"]]$Imputed.Dataset[, Diet], 1, sum, na.rm=TRUE)
  # Imputed[[i]][["R"]]$Imputed.Dataset[, Habitat] <- apply( Imputed[[i]][["R"]]$Imputed.Dataset[, Habitat], 2, as.numeric)
  # Imputed[[i]][["R"]]$Imputed.Dataset$Habitat_breadth_IUCN <- apply( Imputed[[i]][["R"]]$Imputed.Dataset[, Habitat], 1, sum, na.rm=TRUE)
  
  Imputed[[i]][["A"]]$Imputed.Dataset$Class <- "Amphibians"
  colnames(Imputed[[i]][["A"]]$Imputed.Dataset)[10] <- "Lifespan_proxy"
  # Imputed[[i]][["A"]]$Imputed.Dataset$Diet_breadth <- apply( Imputed[[i]][["A"]]$Imputed.Dataset[, Diet], 1, sum, na.rm=TRUE)
  # Imputed[[i]][["A"]]$Imputed.Dataset[, Habitat] <- apply( Imputed[[i]][["A"]]$Imputed.Dataset[, Habitat], 2, as.numeric)
  # Imputed[[i]][["A"]]$Imputed.Dataset$Habitat_breadth_IUCN <- apply( Imputed[[i]][["A"]]$Imputed.Dataset[, Habitat], 1, sum, na.rm=TRUE)
}

## check diet breadth -- check those that have 0 diet breadth
Imputed[[8]][["A"]]$Imputed.Dataset$Diet_breadth %>%  unique()
Imputed[[8]][["B"]]$Imputed.Dataset$Diet_breadth %>%  unique()
Imputed[[8]][["M"]]$Imputed.Dataset$Diet_breadth %>%  unique()
Imputed[[8]][["R"]]$Imputed.Dataset$Diet_breadth %>%  unique()

## check habitat breadth -- check those that have 0 habitat breadth
Imputed[[8]][["A"]]$Imputed.Dataset$Habitat_breadth_IUCN %>%  table()
Imputed[[8]][["B"]]$Imputed.Dataset$Habitat_breadth_IUCN %>%  table()
Imputed[[8]][["M"]]$Imputed.Dataset$Habitat_breadth_IUCN %>%  table()
Imputed[[8]][["R"]]$Imputed.Dataset$Habitat_breadth_IUCN %>%  table()

## bind and transform
Imputed_t <- Bind_transform_traits(Imputed)
saveRDS(Imputed_t, "../Results/Imputed_traits_transformed.rds")

## checking distributions
ex <- Imputed_t[[8]]
hist(ex$log10_Body_mass_g)
hist(ex$log10_Lifespan_proxy)
hist(ex$log10_Litter_size)
hist(ex$sqrt_Diet_breadth)
hist(ex$sqrt_Habitat_breadth_IUCN)

table(ex$Trophic_level)
table(ex$Diel_activity)
table(ex$Specialisation)
table(ex$Primary_diet)

## plotting diet before imputations

Mammals <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Mammals.csv") %>% 
  mutate(Class="Mammals")
Birds <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Birds.csv") %>% 
  mutate(Class="Birds")
colnames(Birds)[37] <- "Trophic_level"
Amphibians <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Amphibians.csv") %>% 
  mutate(Class="Amphibians")
Reptiles <- read.csv("../../Data/TraitsWithDiet/Diet_redefined/Reptiles.csv") %>% 
  mutate(Class="Reptiles")

Dietdata <- rbind(Mammals[, c("Class", "Primary_diet", "Best_guess_binomial")], 
                  Birds[, c("Class", "Primary_diet", "Best_guess_binomial")], 
                  Amphibians[, c("Class", "Primary_diet", "Best_guess_binomial")], 
                  Reptiles[, c("Class", "Primary_diet", "Best_guess_binomial")])

Dietdata$Primary_diet <- as.character(Dietdata$Primary_diet)
Dietdata$Primary_diet[Dietdata$Primary_diet=="SE"] <- "PL|SE"
Dietdata$Primary_diet[Dietdata$Primary_diet=="PL"] <- "PL|SE"
Dietdata$Primary_diet[Dietdata$Primary_diet=="NE"] <- "FR|NE"
Dietdata$Primary_diet[Dietdata$Primary_diet=="FR"] <- "FR|NE"
Dietdata$Imputed <- "Collected data"
Dietdata$Class

table(Dietdata$Class, Dietdata$Primary_diet)

## plottong diet after imputations

ModelData <-  readRDS("../../Results/Imputed_traits_transformed.rds")[[8]]
ModelData$Primary_diet %>%  levels
ModelData$Imputed <- "After imputations"

ModelData$Primary_diet <- as.character(ModelData$Primary_diet)
ModelData$Primary_diet[ModelData$Primary_diet=="SE"] <- "PL|SE"
ModelData$Primary_diet[ModelData$Primary_diet=="PL"] <- "PL|SE"
ModelData$Primary_diet[ModelData$Primary_diet=="NE"] <- "FR|NE"
ModelData$Primary_diet[ModelData$Primary_diet=="FR"] <- "FR|NE"


## plot for thesis / manuscript -- on diet: distribution of diet before and after imputations
X1 <- Dietdata[, c("Primary_diet", "Class", "Imputed", "Best_guess_binomial")] 
table(X1)

X2 <-  ModelData[, c("Primary_diet", "Class", "Imputed", "Best_guess_binomial")]
table(X2)

Data <- rbind(X1, X2)
Data$Imputed <- factor(Data$Imputed, levels=c("Collected data", "After imputations"))

table(Data$Class, Data$Primary_diet, Data$Imputed)

pPD2 <- ggplot(Data, aes(x=Primary_diet, group=Imputed, fill=Imputed), stat = "count") + 
  geom_bar() + 
  facet_grid(Imputed~Class) + GGPoptions + xlab("") + ylab("Number of species") +
  theme(axis.text.x = element_text(angle = 60,  hjust =1)) + 
  geom_text( stat='count', aes(label=..count..), vjust=-1) + ylim(0, 11000) +
  ggtitle("Distribution of diets among vertebrate classes") +
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_fill_viridis_d(begin = 0.1, end = 0.7, name="") + ylim(c(0, 12500)) +
  theme(legend.position = "bottom")
  

ggsave(pPD2, filename = "../Results/DietDataCoverage/diet_redefined/Diet_distribution_collected_imputed.pdf", width = 11, height= 6)

## coverage for PREDICTS data

## load PREDICTS 
Predicts <- readRDS("../../Results/Predicts_merged_sites.rds")
Predicts$Occurrence <- ifelse(Predicts$Measurement > 0, 1, 0)
Predicts$LandUse <- as.character(Predicts$Predominant_land_use)
Predicts$LandUse[Predicts$LandUse=="Cannot decide"] <- NA
Predicts$LandUse[Predicts$LandUse=="Secondary vegetation (indeterminate age)"] <- NA
Predicts$LandUse <- factor(Predicts$LandUse,
                           levels=c("Primary vegetation", 
                                    "Mature secondary vegetation",
                                    "Intermediate secondary vegetation",
                                    "Young secondary vegetation",
                                    "Plantation forest" ,
                                    "Pasture", 
                                    "Cropland",
                                    "Urban"))
Predicts <- subset(Predicts, !is.na(LandUse))

Predicts$Use_intensity <- as.character(Predicts$Use_intensity)
Predicts$Use_intensity[Predicts$Use_intensity=="Cannot decide"] <- NA
Predicts$Use_intensity <- factor(Predicts$Use_intensity, levels=c("Minimal use", "Light use", "Intense use"))

## match species with trait

DataPredicts <- subset(Data, Best_guess_binomial %in% unique(Predicts$Best_guess_binomial))

pPD3 <- ggplot(DataPredicts, aes(x=Primary_diet, group=Imputed, fill=Imputed), stat = "count") + 
  geom_bar() + 
  facet_grid(Imputed~Class) + GGPoptions + xlab("") + ylab("Number of species") +
  theme(axis.text.x = element_text(angle = 60,  hjust =1)) + 
  geom_text( stat='count', aes(label=..count..), vjust=-1) + 
  ggtitle("Distribution of diets among vertebrate classes (PREDICTS database)") +
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_fill_viridis_d(begin = 0.1, end = 0.7, name="") + ylim(c(0, 1800)) +
  theme(legend.position = "bottom")

ggsave(pPD3, filename = "../Results/DietDataCoverage/diet_redefined/Diet_distribution_collected_imputed_PREDICTS.pdf", width = 11, height= 6)



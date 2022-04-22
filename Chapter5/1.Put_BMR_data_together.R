library(dplyr)
library(stringr)
library(data.table)

## first source of data is Stark 2020 (https://doi.org/10.1111/geb.13069)
Stark2020 <- read.csv("D:/4.BMR/Data/Stark_2020_BMR_Herptiles.csv")
colnames(Stark2020)[1] <- "Species"
colnames(Stark2020)[12] <- "BMR_ml_O2_per_h"
colnames(Stark2020)[13] <- "FMR_kJ_per_d"
Stark2020 <- Stark2020 %>% 
  select(Species, Class, Order, BMR_ml_O2_per_h, FMR_kJ_per_d, Mean.Body.Mass..grams.) 
Stark2020 <- Stark2020[-nrow(Stark2020),]
Stark2020 %>%  
  group_by(Class) %>% 
  filter(!is.na(BMR_ml_O2_per_h)) %>% 
  dplyr::summarise(C=n())
plot(log(Stark2020$BMR_ml_O2_per_h)~log(Stark2020$Mean.Body.Mass..grams.), col=as.factor(Stark2020$Class))
Stark2020 <- Stark2020 %>% 
  filter(Species!="")

## second source of data is PanTHERIA (for mammals)
Pantheria <- read.csv("D:/4.BMR/Data/Pantheria.csv")
colnames(Pantheria)
Pantheria <- Pantheria %>% 
  select(MSW05_Order, MSW05_Family, MSW05_Binomial, X18.1_BasalMetRate_mLO2hr, X5.1_AdultBodyMass_g)
colnames(Pantheria) <- c("Order", "Family", "Species", "BMR_ml_O2_per_h", "Body_mass_g")
Pantheria[Pantheria==-999] <- NA
Pantheria <- Pantheria %>% 
  filter(!is.na(BMR_ml_O2_per_h))
Pantheria$Class <- "Mammals"
plot(log(Pantheria$BMR_ml_O2_per_h)~log(Pantheria$Body_mass_g))


## third - McNab et al. 2009
McNab2009 <- read.csv("D:/4.BMR/Data/McNab2009.csv")
colnames(McNab2009)[1] <- "Species"
colnames(McNab2009)[2] <- "Mass_g"
McNab2009$Species <- gsub(McNab2009$Species, pattern="Â", replacement="")
McNab2009$Species <- str_trim(McNab2009$Species, side = "left")
McNab2009 <- McNab2009[McNab2009$Mass_g!="",]
McNab2009 <- McNab2009 %>% 
  filter(!is.na(Species))
McNab2009$BMR..kJ.h. <- McNab2009$BMR..kJ.h. %>% as.character() %>%  as.numeric()
plot(log(McNab2009$Mass_g), log(McNab2009$BMR..kJ.h.))
# for McNab, need to transform BMR into mL_O2_per_hour
McNab2009$BMR_ml_O2_per_h <- as.numeric(McNab2009$BMR..kJ.h.)/0.02
plot(log(McNab2009$BMR_ml_O2_per_h)~log(as.numeric(McNab2009$Mass_g)), pch=19)


## Fristoe 2015
Fristoe2015 <- read.csv("D:/4.BMR/Data/Fristoe2015.csv")
Fristoe2015$Mass..g. <- Fristoe2015$Mass..g. %>%  as.character() %>%  as.numeric()
Fristoe2015$BMR_mlO2_per_h <- Fristoe2015$BMR_mlO2_per_h %>%  as.character() %>%  as.numeric()
plot(log(as.numeric(Fristoe2015$Mass..g.)), log(as.numeric(Fristoe2015$BMR_mlO2_per_h)), pch=19)

colnames(Fristoe2015)[1] <- "Species"
Fristoe2015 <- Fristoe2015 %>% 
  filter(Species!="")
Fristoe2015$Species <- str_trim(Fristoe2015$Species, side = "both")
Fristoe2015$Class <- NA
Fristoe2015$Class[1:178] <- "Mammals"
Fristoe2015$Class[179:391] <- "Birds"
Fristoe2015 <- Fristoe2015 %>% 
  filter(Order!="")
Fristoe2015$BMR_mlO2_per_h <- as.numeric(Fristoe2015$BMR_mlO2_per_h)
Fristoe2015$Mass..g. <- as.numeric(Fristoe2015$Mass..g.)
plot(log(Fristoe2015$Mass..g.), log(Fristoe2015$BMR_mlO2_per_h), col=as.factor(Fristoe2015$Class), pch=19)


## Londono 2015 (BMR is in Watts, 1W=20.1J/mlO2) but they were measured per minute and not per hour
Londono2015 <- read.csv("D:/4.BMR/Data/Londono2015.csv")
colnames(Londono2015)[1] <- "Species"
Londono2015$Species <- gsub(Londono2015$Species, pattern="[*]", replacement="")
Londono2015$Species <- gsub(Londono2015$Species, pattern="Â", replacement="")
Londono2015$Species <- str_trim(Londono2015$Species, side = "left")
Londono2015$Species <- str_trim(Londono2015$Species, side = "right")

# convert units
Londono2015$BMR_Watts_per_hour <- Londono2015$BMR_Watts*60*60
Londono2015$BMR_ml_O2_per_h <- Londono2015$BMR_Watts_per_hour/(0.02*1000) # was measured in Joules and not in kJ
#Londono2015$BMR_ml_O2_per_h <- Londono2015$BMR_Watts_per_hour/(0.02) # was measured in Joules and not in kJ
Londono2015$Class <- "Birds"
plot(log(Londono2015$BMR_ml_O2_per_h)~log(Londono2015$Mean.mass..g.), pch=19)

length(unique(Londono2015$Species))
Londono2015 <- Londono2015 %>% 
  group_by(Class, Species) %>% 
  summarise(BMR_ml_O2_per_h=mean(BMR_ml_O2_per_h),
            Mean.mass..g.=mean(BMR_ml_O2_per_h))


#############################################################################
## put the data together

# mammals
Mammals_BMR_data <- Pantheria[, c("Species", "BMR_ml_O2_per_h")]
length(unique(Mammals_BMR_data$Species))
Fristoe_Mammals <- subset(Fristoe2015, Class=="Mammals")
Fristoe_Mammals <- Fristoe_Mammals[, c("Species", "BMR_mlO2_per_h")]
colnames(Fristoe_Mammals)[2] <- "BMR_ml_O2_per_h"
Stark2020_mammals <- subset(Stark2020, Class=="Mammalia")
Stark2020_mammals <- Stark2020_mammals[, c("Species", "BMR_ml_O2_per_h")]
Stark2020_mammals[!is.na(Stark2020_mammals$BMR_ml_O2_per_h),] %>%  nrow
Stark2020_mammals[!is.na(Stark2020_mammals$FMR_kJ_per_d),] %>%  nrow

Mammals <- rbind(Fristoe_Mammals, Mammals_BMR_data, Stark2020_mammals)
length(unique(Mammals$Species))
Mammals <- Mammals %>% 
  group_by(Species) %>% 
  summarise(BMR_ml_O2_per_h=mean(BMR_ml_O2_per_h))
Mammals$Class <- "Mammals"
Mammals <- as.data.frame(Mammals)
length(unique(Mammals$Species))
Mammals <- Mammals[!is.na(Mammals$BMR_ml_O2_per_h),]


# birds
MacNabData <- McNab2009[, c("Species", "BMR_ml_O2_per_h")]
MacNabData$Class <- "Birds"
Fristoe_Birds <- subset(Fristoe2015, Class=="Birds")
Fristoe_Birds <- Fristoe_Birds[, c("Species", "BMR_mlO2_per_h", "Class")]
colnames(Fristoe_Birds)[2] <- "BMR_ml_O2_per_h"
Londono2015Birds <- Londono2015[, -4] %>%  as.data.frame()
Stark2020_birds <- subset(Stark2020, Class=="Aves")
Stark2020_birds <- Stark2020_birds[, c("Species", "BMR_ml_O2_per_h")]
Stark2020_birds$Class <- "Birds"
Birds_BMR_data <- rbind(MacNabData, Fristoe_Birds, Londono2015Birds, Stark2020_birds)

length(unique(Birds_BMR_data$Species))

Birds <- Birds_BMR_data %>% 
  group_by(Species) %>% 
  summarise(BMR_ml_O2_per_h=mean(BMR_ml_O2_per_h)) %>% 
  as.data.frame()
length(unique(Birds$Species))
Birds <- Birds[!is.na(Birds$BMR_ml_O2_per_h),]
Birds$Class <- "Birds"

# amphibians
Amphibians <- Stark2020 %>% 
  filter(Class=="Amphibia") %>%  
  select(Species, BMR_ml_O2_per_h)
Amphibians$Class <- "Amphibians"
Amphibians <- Amphibians %>% 
  filter(!is.na(BMR_ml_O2_per_h))
length(unique(Amphibians$Species))

# reptiles
Reptiles <- Stark2020 %>% 
  filter(Class=="Reptilia") %>% 
  select(Species, BMR_ml_O2_per_h)
Reptiles$Class <- "Reptilia"
Reptiles <- Reptiles %>% 
  filter(!is.na(BMR_ml_O2_per_h))
length(unique(Reptiles$Species))


## Assemble dataset
BMR_litterature_data <- rbind(Birds, Mammals, Amphibians, Reptiles)
write.csv(BMR_litterature_data, "D:/4.BMR/Results/1.BMR_litt_data.csv", row.names=FALSE)






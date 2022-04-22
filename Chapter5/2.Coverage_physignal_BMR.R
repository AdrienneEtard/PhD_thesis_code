## BMR data coverage
library(ggplot2)
library(phytools)
library(dplyr)

BMR_data <- read.csv("../Results/1.BMR_litt_data.csv")
Predicts <- readRDS("../Data/PredictsVertebrates.rds")
Traits <- readRDS("../Data/Imputed_traits.rds")[[8]]

Vars <- c("Order", "Family", "Best_guess_binomial", "Body_mass_g")

Traits <- rbind(Traits$M$Imputed.Dataset[Vars],
                Traits$B$Imputed.Dataset[Vars],
                Traits$R$Imputed.Dataset[Vars],
                Traits$A$Imputed.Dataset[Vars])

colnames(Traits)[3] <- "Species"

BMR_data$Class <- as.character(BMR_data$Class)
BMR_data$Class[BMR_data$Class=="Reptilia"] <- "Reptiles"
 
# join traits and BMR data
BMR_data <- left_join(BMR_data, Traits, by="Species" )
BMR_data$Thermoregulation <- ifelse(BMR_data$Class %in% c("Mammals", "Birds"), "Endotherms", "Ectotherms")

# plot BMR against body mass
ggplot(BMR_data, aes(log(Body_mass_g), log(BMR_ml_O2_per_h))) + geom_point() + facet_wrap(~Class) +
  xlab("Body mass (g, log)") +
  ylab("Resting metabolic rate (mL(O2)/hour, log) ") + theme_bw()  + theme(panel.spacing = unit(0, "lines")) 


## Now check coverage / PREDICTS species
Predicts_Mammals <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Mammalia"])
BMR_Mammals <- BMR_data$Species[BMR_data$Class=="Mammals"]
length(intersect(Predicts_Mammals, BMR_Mammals))/length(Predicts_Mammals)*100

Predicts_Birds <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Aves"])
BMR_Birds <- BMR_data$Species[BMR_data$Class=="Birds"]
length(intersect(Predicts_Birds, BMR_Birds))/length(Predicts_Birds)*100

Predicts_Amphibians <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Amphibia"])
BMR_Amphibians <- BMR_data$Species[BMR_data$Class=="Amphibians"]
length(intersect(Predicts_Amphibians, BMR_Amphibians))/length(Predicts_Amphibians)*100

Predicts_Reptiles <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Reptilia"])
BMR_Reptiles <- BMR_data$Species[BMR_data$Class=="Reptiles"]
length(intersect(Predicts_Reptiles, BMR_Reptiles))/length(Predicts_Reptiles)*100

## measuring phylogenetic signal in BMR
# function to calculate phylogenetic signal for each of a 100 sampled trees
PhySignal <- function(Trees, N, BMR_data, Class) {
  
  Data <- BMR_data %>% 
    filter(Class==Class)
  DataVec <- log(Data$BMR_ml_O2_per_h)
  names(DataVec) <- Data$Species
  
  .Format_tiplabels <- function (Phylogeny) {
    Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
    return(Phylogeny)
  }
  Res_lamba <- c()
  for(i in 1:length(N)) {
    Tree <- .Format_tiplabels(Trees[[N[i]]])
    Signal <- phytools::phylosig(Tree, DataVec, method="lambda", test = FALSE) %>% 
      unlist()
    Res_lamba <- c(Res_lamba, Signal$lambda)
    print(i)
  }
  return(Res_lamba)
}

TreesMammals <- read.nexus("D:/PhD/PhD_R_work/4.BMR/Data/phylogenies/Complete_phylogeny.nex")
TreesBirds <- read.nexus("D:/PhD/PhD_R_work/4.BMR/Data/phylogenies/birds1000trees.nex")
TreesAmphibians <- read.nexus("D:/PhD/PhD_R_work/4.BMR/Data/phylogenies/amphibians1000trees.nex")
TreesReptiles <- read.nexus("D:/PhD/PhD_R_work/4.BMR/Data/phylogenies/reptiles1000trees.nex")

N <- sample(100, x=1:1000)
lambda_mammals <- PhySignal(TreesMammals, N, BMR_data, "Mammals")
median(lambda_mammals)
quantile(lambda_mammals, c(0.025, 0.975))

N <- sample(100, x=1:1000)
lambda_birds <- PhySignal(TreesBirds, N, BMR_data, "Birds")
median(lambda_birds)
quantile(lambda_birds, c(0.025, 0.975))

N <- sample(100, x=1:1000)
lambda_amphibians <- PhySignal(TreesAmphibians, N, BMR_data, "Amphibians")
median(lambda_amphibians)
quantile(lambda_amphibians, c(0.025, 0.975))

N <- sample(100, x=1:1000)
lambda_reptiles <- PhySignal(TreesReptiles, N, BMR_data, "Reptiles")
median(lambda_reptiles)
quantile(lambda_reptiles, c(0.025, 0.975))


## prepare the data for imputations (by adding PREDICTS species that are not represented in the BMR data)
Mammals_to_add <- setdiff(Predicts_Mammals, BMR_data$Species) %>% 
  as.data.frame() %>% 
  setNames("Species") %>% 
  mutate(Class="Mammals")

Birds_to_add <- setdiff(Predicts_Birds, BMR_data$Species)%>% 
  as.data.frame() %>% 
  setNames("Species") %>% 
  mutate(Class="Birds")

Amphibians_to_add <- setdiff(Predicts_Amphibians, BMR_data$Species)%>% 
  as.data.frame() %>% 
  setNames("Species") %>% 
  mutate(Class="Amphibians")

Reptiles_to_add <- setdiff(Predicts_Reptiles, BMR_data$Species)%>% 
  as.data.frame() %>% 
  setNames("Species") %>% 
  mutate(Class="Reptiles")

To_add <- rbind(Mammals_to_add, Birds_to_add, Reptiles_to_add, Amphibians_to_add)
To_add$Imputed <- TRUE
To_add <- left_join(To_add, Traits, by="Species")
To_add$BMR_ml_O2_per_h <- NA
To_add$Thermoregulation <- ifelse(To_add$Class %in% c("Mammals", "Birds"), "Endotherms", "Ectotherms")
BMR_data$Imputed <- FALSE

BMR_data_to_impute <- rbind(BMR_data, To_add)
write.csv(BMR_data_to_impute, "../Results/2.BMR_data_to_impute.csv", row.names=FALSE)

ggplot(BMR_data, aes(log(Body_mass_g), log(BMR_ml_O2_per_h/Body_mass_g))) + geom_point() + facet_wrap(~Class)








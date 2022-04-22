# imputations of missing BMR values within each class using phylogenetic position, body mass and taxonomic orders
X <- c("dplyr",
       "phytools",
       "picante", 
       "stringr", 
       "PVR",
       "missForest")
lapply(X, library, character.only=TRUE); rm(X)

library(ggplot2)
library(ggpubr)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


## Functions to extract phylogenetic eigenvectors from the phylogenies and return as a dataframe - for species that match the trait dataset
Extract_eigenvectors <- function(TraitDF, Phylo, N) {
  
  library(PVR)
  
  ## arguments:
  # N <- number of eigenvectors to extract: 10 is enough to maximise imputation accuracy (Penone et al. 2014)
  # Phylo: phylogeny considered
  # TraitDF: trait dataset with species names in "Species"
  
  # to format tip label names of phylogenies
  .Format_tiplabels <- function (Phylogeny) {
    Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
    return(Phylogeny)
  }
  Phylo <- .Format_tiplabels(Phylo)
  
  ## Prune species from tree that do not intersect with trait dataset
  row.names(TraitDF) <- TraitDF$Species
  Prune_Taxa <- match.phylo.data(Phylo, TraitDF)
  Phylo <- Prune_Taxa$phy
  
  ## Get phylogenetic eigenvectors from the phylogeny and select N first eigenvectors - PVR package (Thiago Santos)
  print("Eigenvector decomposition.")
  EigenV <- PVR::PVRdecomp(Phylo)
  
  Eigenvectors <- EigenV@Eigen$vectors
  Eigenvectors <- as.data.frame(Eigenvectors)
  Eigenvectors <- Eigenvectors[, 1:N]
  
  # rename columns in trait dataset: EV_1, ..., EV_N
  for (i in 1:N) {colnames(Eigenvectors)[i] <- paste0("EV_",i)}
  
  # add species names and reorder

  Eigenvectors$Species <- Prune_Taxa$data$Species
  Eigenvectors <- Eigenvectors[order(Eigenvectors$Species), c(ncol(Eigenvectors), 1:N)]
  
  return(Eigenvectors)
  
}

## Function to add Eigenvectors to trait dataset
Add_eigenvectors <- function (TraitDF, EV) {
  ColN <- vector()
  for (i in 1:5) {ColN <- c(ColN, paste("EV", i, sep="_")) }
  Species <- as.character(EV$Species)
  x <- which(TraitDF$Species %in% Species)
  TraitDF[, ColN] <- NA
  TraitDF[x, ColN] <- EV[, c(2:ncol(EV))]
  return(TraitDF)
}

## Function to impute missing trait values
Impute <- function(data, tree, NEigen) {
  
  EV <- Extract_eigenvectors(data, tree, N=NEigen)
  
  data <- Add_eigenvectors(data, EV)
  
  data$log_BMR <- log(data$BMR_ml_O2_per_h)
  
  data_to_impute <- data %>% 
    dplyr::select(-Family, -Thermoregulation, -Species, -Class, -Imputed)
  
  data_to_impute$Order<- as.factor(data_to_impute$Order)
  data_to_impute$Order <- droplevels(data_to_impute$Order)
  
  data_imputed <- missForest(data_to_impute, variablewise = TRUE)
  
  to_return <- cbind(data_imputed$ximp, data[,c("Family", "Thermoregulation", "Species", "Class", "Imputed")])
  to_return$MSE_logBMR <- data_imputed$OOBerror[9]
  
  return(list(data=to_return, error=data_imputed$OOBerror))
}


###################################################################

# phylogenetic trees
TreesMammals <- read.newick("../Data/phylogenies/consensus/Mammals.nwk")
TreesBirds <- read.newick("../Data/phylogenies/consensus/Birds.nwk")
TreesAmphibians <- read.newick("../Data/phylogenies/consensus/Amphibians.nwk")
TreesReptiles <- read.newick("../Data/phylogenies/consensus/Reptiles.nwk")

# data to impute
BMR_data <- read.csv("../Results/2.BMR_data_to_impute.csv")

BMR_data %>%  filter(!is.na(BMR_ml_O2_per_h)) %>%  group_by(Class) %>%  summarise(N=n())

# plotting
ggplot2::ggplot(BMR_data, aes(log(Body_mass_g),log(BMR_ml_O2_per_h))) + geom_point() + theme_bw() +
  xlab("Body mass (g, log)") + ylab("RMR (mL O2/h, g)") + 
  geom_smooth(method='lm', col="red", fill="red", alpha=0.3)+
  facet_wrap(~Class) +
  stat_cor(label.y = 12, size=3.5)


# imputations
Mammals_imputed <- Impute(BMR_data[BMR_data$Class=="Mammals",], TreesMammals, NEigen = 5)
Birds_imputed <- Impute(BMR_data[BMR_data$Class=="Birds",], TreesBirds, NEigen = 5)
Amphibians_imputed <- Impute(BMR_data[BMR_data$Class=="Amphibians",], TreesAmphibians, NEigen = 5)
Reptiles_imputed <- Impute(BMR_data[BMR_data$Class=="Reptiles",], TreesReptiles, NEigen = 5)

Imputed_data <- rbind(Mammals_imputed$data, Birds_imputed$data, Amphibians_imputed$data, Reptiles_imputed$data)

## sensitivity to number of eigenvectors
res_imp <- list()
for(EV in c(2:10)) {
  M_imp <- Impute(BMR_data[BMR_data$Class=="Mammals",], TreesMammals, NEigen = EV)
  B_imp <- Impute(BMR_data[BMR_data$Class=="Birds",], TreesBirds, NEigen = EV)
  A_imp <- Impute(BMR_data[BMR_data$Class=="Amphibians",], TreesAmphibians, NEigen = EV)
  R_imp <- Impute(BMR_data[BMR_data$Class=="Reptiles",], TreesReptiles, NEigen = EV)
  # res_errors[[EV]] <- data.frame(Mammals=unique(M_imp$data$MSE_logBMR),
  #                                Birds=unique(B_imp$data$MSE_logBMR), 
  #                                Amphibians=unique(A_imp$data$MSE_logBMR),
  #                                Reptiles=unique(R_imp$data$MSE_logBMR), 
  #                                EV_n=EV)
  res_imp[[EV]] <- rbind(M_imp$data, B_imp$data, A_imp$data, R_imp$data)
  print(EV)
}

plot(log(res_imp[[2]]$Body_mass_g), res_imp[[2]]$log_BMR, col=as.factor(res_imp[[2]]$Imputed), pch=19)
plot(log(res_imp[[2]]$Body_mass_g), res_imp[[2]]$log_BMR, col=as.factor(res_imp[[2]]$Imputed), pch=19)

Corrs <- function(res_imp, i, VClass){
  res_imp <- res_imp[[i]]
  res_imp <- subset(res_imp, Class==VClass)
  d1 <- cor(log(res_imp$Body_mass_g[res_imp$Imputed=="FALSE"]), res_imp$log_BMR[res_imp$Imputed=="FALSE"])
  d2 <- cor(log(res_imp$Body_mass_g[res_imp$Imputed=="TRUE"]), res_imp$log_BMR[res_imp$Imputed=="TRUE"])
  return(abs(d1-d2))
}

Corrs(res_imp,2, "Birds") 
Corrs(res_imp,2, "Mammals") 
Corrs(res_imp,3, "Mammals")
Corrs(res_imp,4, "Mammals")
Corrs(res_imp,5, "Mammals")
Corrs(res_imp,6, "Mammals")
Corrs(res_imp,7, "Mammals")
Corrs(res_imp,8, "Mammals")
Corrs(res_imp,9, "Mammals")
Corrs(res_imp,10, "Mammals")

Imputed_data5 <- res_imp[[5]]
Imputed_data6 <- res_imp[[6]]
Imputed_data7 <- res_imp[[7]]
Imputed_data8 <- res_imp[[8]]
Imputed_data9 <- res_imp[[9]]
Imputed_data10 <- res_imp[[10]]

plot(Imputed_data5$log_BMR~Imputed_data10$log_BMR, pch=19)
cor(Imputed_data5$log_BMR,Imputed_data10$log_BMR)
cor(Imputed_data6$log_BMR,Imputed_data10$log_BMR)
cor(Imputed_data5$log_BMR,Imputed_data6$log_BMR)
cor(Imputed_data5$log_BMR,Imputed_data10$log_BMR)

plot(Imputed_data6$log_BMR~Imputed_data10$log_BMR, pch=19)
plot(Imputed_data7$log_BMR~Imputed_data10$log_BMR, pch=19)
plot(Imputed_data8$log_BMR~Imputed_data10$log_BMR, pch=19)
plot(Imputed_data9$log_BMR~Imputed_data10$log_BMR, pch=19)

ggplot(Imputed_data, aes(log(Body_mass_g),log_BMR, col=Imputed)) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~Class) + GGPoptions +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  xlab("Body mass (g, log)") + ylab("RMR (mL O2/hour, log)") +
  scale_colour_manual(values=c( "#0072B2","#E69F00"), name="") + 
  theme(legend.position = "top")


res_errors <- data.table::rbindlist(res_errors)
res_errors$median <- apply(res_errors[,1:4], 1, median)
plot(res_errors$EV_n, res_errors$Mammals)
plot(res_errors$EV_n, res_errors$Birds)
plot(res_errors$EV_n, res_errors$Amphibians)
plot(res_errors$EV_n, res_errors$Reptiles)
plot(res_errors$EV_n, res_errors$median)
plot(res_errors$EV_n, res_errors$mean)

## plotting ot visualise imputed BMR against BM

ggplot(Imputed_data, aes(log(Body_mass_g), log_BMR, col=Imputed)) + 
  geom_point() + 
  facet_wrap(~Class) 

## get mass-standardised BMR
Imputed_data$log_BMR_standardised <- log(exp(Imputed_data$log_BMR)/Imputed_data$Body_mass_g)
write.csv(Imputed_data, "D:/PhD/PhD_R_work/4.BMR/Results/3.BMR_data_imputed.csv", row.names = FALSE)

Imputed_data <- "../Results/3.BMR_data_imputed.csv" %>%  read.csv()
Imputed_data$Imputed <- ifelse(Imputed_data$Imputed==FALSE, "Collected data point", "Imputed data point")

ggplot(Imputed_data, aes(log(Body_mass_g),log_BMR, col=Imputed)) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~Class) + GGPoptions +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  xlab("Body mass (g, log)") + ylab("RMR (mL O2/hour, log)") +
  scale_colour_manual(values=c( "#0072B2","#E69F00"), name="") + 
  theme(legend.position = "top")


cor(log(Imputed_data$Body_mass_g),Imputed_data$log_BMR)

cor(log(Imputed_data$Body_mass_g[Imputed_data$Class=="Birds"]),Imputed_data$log_BMR[Imputed_data$Class=="Birds"])
cor(log(Imputed_data$Body_mass_g[Imputed_data$Class=="Mammals"]),Imputed_data$log_BMR[Imputed_data$Class=="Mammals"])
cor(log(Imputed_data$Body_mass_g[Imputed_data$Class=="Reptiles"]),Imputed_data$log_BMR[Imputed_data$Class=="Reptiles"])
cor(log(Imputed_data$Body_mass_g[Imputed_data$Class=="Amphibians"]),Imputed_data$log_BMR[Imputed_data$Class=="Amphibians"])




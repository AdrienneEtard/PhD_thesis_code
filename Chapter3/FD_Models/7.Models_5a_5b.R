#library(betapart)
library(dplyr)
library(geometry)
library(rcdd)
library(fastmatch)
library(utils)
library(alphahull)
library(plyr)
library(reshape2)
library(ggplot2)
library(viridis)

source("Functions.R")

na.omit.list <- function(y) {
  return(y[!sapply(y, function(x) all(is.na(x)))])
}

## functions for calculation of functional loss and functional gain

functional_volumes <- function (x, traits) {

  ## check conditions
  D <- ncol(traits)
  Si <- rowSums(x)
  if (any(Si <= D)) {
    stop(paste("'community ", row.names(x)[which(Si <= D)],
               " must contain at least ", D + 1, " species", sep = "")) }

  if ( nrow(x) < 2) {stop("Computing dissimilairty requires at least 2 communities", call. = TRUE)}


  ## Volumes calculations (FRi) -- N is the number of sites (2, because pairwise comparisons of 2 land uses)
  N <- nrow(x)
  FRi <- rep(NA, N)
  names(FRi) <- row.names(x)
  coord_vert_i <- list()

  for (i in 1:N) {
    tr_i <- traits[which(x[i, ] == 1), ]
    vert0 <- convhulln(tr_i, "Fx TO 'vert.txt'")
    vert1 <- scan("vert.txt", quiet = TRUE)
    verti <- (vert1 + 1)[-1]
    coord_vert_i[[i]] <- tr_i[verti, ]
    FRi[i] <- convhulln(tr_i[verti, ],"FA")$vol
  }

  ## Intersection
  inter <- function(set1, set2) {
    set1rep <- d2q(cbind(0, cbind(1, set1)))
    set2rep <- d2q(cbind(0, cbind(1, set2)))
    polytope1 <- redundant(set1rep, representation = "V")$output
    polytope2 <- redundant(set2rep, representation = "V")$output
    H_chset1 <- scdd(polytope1, representation = "V")$output
    H_chset2 <- scdd(polytope2, representation = "V")$output
    H_inter <- rbind(H_chset1, H_chset2)
    V_inter <- scdd(H_inter, representation = "H")$output
    vert_1n2 <- q2d(V_inter[, -c(1, 2)])
    coord_vert_inter <- rep(NA, ncol(set1))
    vol_inter <- 0
    if (is.matrix(vert_1n2))
      if (nrow(vert_1n2) > ncol(vert_1n2)) {
        coord_vert_inter <- vert_1n2
        vol_inter <- convhulln(vert_1n2, "FA")$vol
      }
    res <- list(coord_vert_inter = coord_vert_inter, vol_inter = vol_inter)
    return(res)
  }

  comb2 <- combn(N, 2)

  vol_inter2_mat <- matrix(0, N, N, dimnames = list(row.names(x), row.names(x)))
  vol_inter2 <- rep(0, ncol(comb2))
  coord_vert_inter2 <- list()

  for (k in 1:ncol(comb2)) {
    i <- comb2[1, k]
    j <- comb2[2, k]
    seti <- traits[which(x[i, ] == 1), ]
    setj <- traits[which(x[j, ] == 1), ]
    interij <- inter(seti, setj)
    vol_inter2_mat[j, i] <- interij$vol_inter
    vol_inter2[k] <- interij$vol_inter
    coord_vert_inter2[[k]] <- interij$coord_vert_inter
  }

  ## volumes
  FRX <- FRi[1]
  FRO <- FRi[2]
  Inter <- vol_inter2

  return(list(FRX=FRX, FRO=FRO, Inter=Inter))

}

Compute_betafunc <- function(Gower_distances, Community_matrix, N, Predicts.info, Use_intensity) {

  Predicts.info$Predominant_land_use <- as.character(Predicts.info$Predominant_land_use)
  Predicts.info$Predominant_land_use[Predicts.info$Predominant_land_use=="Primary vegetation"] <- "X"

  # get information: sites names and use intensities
  Predicts.info$Use_intensity <- as.character(Predicts.info$Use_intensity)
  Sites <- rownames(Community_matrix)
  UseIntensities <-  Predicts.info$Use_intensity[Predicts.info$SSBS %in% Sites] %>% as.character

  # add this information to the matrix for later use
  Community_matrix <- as.data.frame(Community_matrix)
  Community_matrix$Sites <- Sites
  Community_matrix$Use_intensity <- UseIntensities
  Community_matrix <- as.matrix(Community_matrix)

  # get land uses in Community_matrix & use intensity
  for(i in 1:nrow(Community_matrix)){
    Rn <- rownames(Community_matrix)[i]
    rownames(Community_matrix)[i] <- Predicts.info$Predominant_land_use[Predicts.info$SSBS==Rn] %>% as.character %>% unique()
  }


  # if no Primary vegetation, no need to calculate beta diversity for that study
  RN <- rownames(Community_matrix)

  if(!any(grepl("X", RN))){
    return(NA)
  } else{

    # comparisons with PV only: create pairs
    Community_matrix <- as.data.frame(Community_matrix)
    Pairs <- expand.grid(rownames(Community_matrix), rownames(Community_matrix))
    Pairs <- Pairs %>%  filter(grepl("X", Var1)|grepl("X", Var2))
    Pairs <- unique(t(apply(Pairs, 1, sort))) %>%  as.data.frame() %>%  setNames(., c("Var1", "Var2"))

    # filter out pairs of same site
    Same <- which(as.character(Pairs$Var1)==as.character(Pairs$Var2))
    Pairs <- Pairs[-Same,]

    ## reorder so that PV (renamed X) is always the reference
    if(any(grepl("Young",  Pairs$Var2))){
      Pairs$Var1 <- as.character(Pairs$Var1)
      Pairs$Var2 <- as.character(Pairs$Var2)

      x <- which(grepl("Young.secondary.vegetation", Pairs$Var2))
      Toplace <- Pairs$Var2[x]
      Pairs$Var2[x] <- Pairs$Var1[x]
      Pairs$Var1[x] <- Toplace

      ## consider all pairwise possiblities for PV
      if(length(unique(Pairs$Var2))>1){
        ToAdd <- subset(Pairs, grepl("X", Pairs$Var1) & grepl("X", Pairs$Var2)) %>%  rev()
        colnames(ToAdd) <- c("Var1", "Var2")
        Pairs <- rbind(Pairs, ToAdd)
      }
    }

    Community_matrix_2 <- as.data.frame(Community_matrix)
    Community_matrix <- as.data.frame(Community_matrix) %>%
      dplyr::select(-Sites, -Use_intensity) %>%
      as.matrix()
    Community_matrix <- apply(Community_matrix, 2, as.numeric)
    rownames(Community_matrix) <- rownames(Community_matrix_2)

    PairList <- list()
    for (i in 1:nrow(Pairs)){
      Df <- rbind(Community_matrix[rownames(Community_matrix)==Pairs$Var2[i],], Community_matrix[rownames(Community_matrix)==Pairs$Var1[i],])
      PairList[[i]] <- Df

      # put land uses and use intensities as names
      Name1 <- paste0(Community_matrix_2$Use_intensity[rownames(Community_matrix)==Pairs$Var2[i]], "/",
                    Community_matrix_2$Use_intensity[rownames(Community_matrix)==Pairs$Var1[i]])

      Name2 <- paste0(Community_matrix_2$Sites[rownames(Community_matrix)==Pairs$Var2[i]], "/",
                    Community_matrix_2$Sites[rownames(Community_matrix)==Pairs$Var1[i]])
      Name <- paste0(Name2, "/", Name1)
      names(PairList)[[i]]<- Name
    }

    ## computing beta diversity
    Computebetadiv <- function(comm.matrix, Gowdis){

      # subset distance matrix for species in the community
      Dist <- usedist::dist_subset(Gowdis, colnames(comm.matrix))

      # check if euclidian, and if not apply cailliez correction
      if(!ade4::is.euclid(Dist)){
        Dist <- ade4::cailliez(Dist)
      }

      # dimensionality reduction (n dim <= 4; nspecies > ndim)
      nSp <- apply(comm.matrix, 1, sum)
      nSp <- min(nSp)

      # perform dimensionality reduction and select the first few axes (maximum number of dimensions)
      x.pco <- ade4::dudi.pco(Dist, scannf = FALSE, full = TRUE)

      # number of sepcies in each sites must be striclty > than D

      # get coordinates and retain first axes
      if(nSp > 2){
        Traits <- x.pco$li[,1:2]
      }

      # if too many communities (sites)
      if(nrow(comm.matrix)>N){
        return(NA)
      }

      if(nSp<=2){
        return(NA)
      } else {

        #print(paste("Number of sites:", nrow(Community_matrix)))
        Corebeta <- functional_volumes(x=comm.matrix, traits = as.matrix(Traits))
        if(file.exists("vert.txt")){rm("vert.txt")}
        return(Corebeta)
      }

    }

    BetaDiv <- lapply(PairList, Computebetadiv, Gowdis=Gower_distances)

    return(BetaDiv)
  }
}

## function to process results
ProcessResults <- function(Res, Info){

  List_elements <- function(Y){
    Res_int <- list()
    for(i in 1:length(Y)) {
      X <- Y[[i]]
      Names <- names(Y[i])
      if(is.na(X)){Res_int[[i]] <- NA}
      else{
        Xun <- unlist(X) %>% as.data.frame() %>% t() %>%  as.data.frame()
        Xun$LU1 <- sub(x=colnames(Xun)[1], pattern="FRX.", replacement="")
        Xun$LU2 <- sub(x=colnames(Xun)[2], pattern="FRO.", replacement="")
        colnames(Xun) <- c("Vol.PV", "Vol.LU2", "Vol.Inter", "LU1", "LU2")
        #browser()
        Xun$Use_intensity <- Names
        Names <- str_split(Names, pattern = "/") %>%  unlist()
        Xun$LU1 <- Names[1]
        Xun$LU2 <- Names[2]
        Xun$Use_intensity_1 <- Names[3]
        Xun$Use_intensity_2 <- Names[4]
        Res_int[[i]] <- Xun
      }
    }
    Res_int <- na.omit.list(Res_int)
    Res_int <- data.table::rbindlist(Res_int)
    return(Res_int)
  }

  Results_all <- list()
  for(j in 1:length(Res)){
    Results_all[[j]] <- List_elements(Res[[j]])
    Results_all[[j]]$SS <- names(Res)[[j]]
  }
  to_return <- data.table::rbindlist(Results_all)

  to_return$Land_use_1 <- NA
  to_return$Land_use_1 <- NA

  for(i in 1:nrow(to_return)){
    to_return$Land_use_1[i] <- Info$Predominant_land_use[Info$SSBS==to_return$LU1[i]] %>% as.character()
    to_return$Land_use_2[i] <- Info$Predominant_land_use[Info$SSBS==to_return$LU2[i]] %>% as.character()
  }

  return(to_return)
}

#########################################################################################################################

## Estimating functional loss and gain, within vertebrate classes
Communities_Occu_Mammals <- readRDS("../../Data/Presence_absence_Mammals.rds") %>% na.omit.list()
Communities_Occu_Birds <- readRDS("../../Data/Presence_absence_Birds.rds") %>% na.omit.list()
Communities_Occu_Amphibians <- readRDS("../../Data/Presence_absence_Amphibians.rds") %>% na.omit.list()
Communities_Occu_Reptiles <- readRDS("../../Data/Presence_absence_Reptiles.rds") %>% na.omit.list()

Gower_distance_matrices_Mammals <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Mammals.rds")
Gower_distance_matrices_Birds <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Birds.rds")
Gower_distance_matrices_Amphibians <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Amphibians.rds")
Gower_distance_matrices_Reptiles <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets_Reptiles.rds")

Predictsinfo <- read.csv("../../Results/Predicts_site_info_ER.csv")
Reg <- read.csv("../../Results/dbFD_indices_for_analysis/dbFD_region_std.csv")

## Running the function for each class

## mammals
Betadiv_Mammals <- list()
for(i in 1:length(Communities_Occu_Mammals)) {
  Betadiv_Mammals[[i]] <- Compute_betafunc(Gower_distance_matrices_Mammals[[8]],
                                   Communities_Occu_Mammals[[i]],
                                   N=20,
                                   Predicts.info = Predictsinfo)
  print(paste("Study", i, "on", length(Communities_Occu_Mammals)))
}
names(Betadiv_Mammals) <- names(Communities_Occu_Mammals)
Betadiv_Mammals <- na.omit.list(Betadiv_Mammals)
saveRDS(Betadiv_Mammals, "../../Results/Betadiv/Betadiv_Mammals.rds")
Results_mammals <- ProcessResults(Betadiv_Mammals, Predictsinfo)
write.csv(Results_mammals, "../../Results/Betadiv/Betadiv_Mammals_results.csv", row.names = FALSE)

## amphibians
Betadiv_amphibians <- list()
for(i in 1:13) {
  Betadiv_amphibians[[i]] <- Compute_betafunc(Gower_distance_matrices_Amphibians[[8]],
                                           Communities_Occu_Amphibians[[i]],
                                           N=20,
                                           Predicts.info = Predictsinfo)
  print(paste("Study", i, "on", length(Communities_Occu_Amphibians)))
}
Betadiv_amphibians[[14]] <- NA
for(i in 15:length(Communities_Occu_Amphibians)) {
  Betadiv_amphibians[[i]] <- Compute_betafunc(Gower_distance_matrices_Amphibians[[8]],
                                              Communities_Occu_Amphibians[[i]],
                                              N=20,
                                              Predicts.info = Predictsinfo)
  print(paste("Study", i, "on", length(Communities_Occu_Amphibians)))
}
names(Betadiv_amphibians) <- names(Communities_Occu_Amphibians)
Betadiv_amphibians <- na.omit.list(Betadiv_amphibians)
saveRDS(Betadiv_amphibians, "../../Results/Betadiv/Betadiv_Amphibians.rds")
Results_amphibians <- ProcessResults(Betadiv_amphibians, Predictsinfo)
write.csv(Results_amphibians, "../../Results/Betadiv/Betadiv_Amphibians_results.csv", row.names = FALSE)

## reptiles
Betadiv_reptiles <- list()
for(i in 1:length(Communities_Occu_Reptiles)) {
  Betadiv_reptiles[[i]] <- Compute_betafunc(Gower_distance_matrices_Reptiles[[8]],
                                              Communities_Occu_Reptiles[[i]],
                                              N=20,
                                              Predicts.info = Predictsinfo)
  print(paste("Study", i, "on", length(Communities_Occu_Reptiles)))
}
names(Betadiv_reptiles) <- names(Communities_Occu_Reptiles)
Betadiv_reptiles <- na.omit.list(Betadiv_reptiles)
saveRDS(Betadiv_reptiles, "../../Results/Betadiv/Betadiv_Reptiles.rds")
Results_reptiles <- ProcessResults(Betadiv_reptiles, Predictsinfo)
write.csv(Results_reptiles, "../../Results/Betadiv/Betadiv_reptiles_results.csv", row.names = FALSE)


## birds
Betadiv_birds <- list()
for(i in 1:length(Communities_Occu_Birds)) {
  Betadiv_birds[[i]] <- Compute_betafunc(Gower_distance_matrices_Birds[[8]],
                                              Communities_Occu_Birds[[i]],
                                              N=20,
                                              Predicts.info = Predictsinfo)
  print(paste("Study", i, "on", length(Communities_Occu_Birds)))
}
names(Betadiv_birds) <- names(Communities_Occu_Birds)
Betadiv_birds <- na.omit.list(Betadiv_birds)
saveRDS(Betadiv_birds, "../../Results/Betadiv/Betadiv_Birds.rds")
Results_birds <- ProcessResults(Betadiv_birds, Predictsinfo)
write.csv(Results_birds, "../../Results/Betadiv/Betadiv_birds_results.csv", row.names = FALSE)


#########################################################################################################################
#########################################################################################################################
Results_amphibians$Class <- "Amphibians"
Results_birds$Class <- "Birds"
Results_mammals$Class <- "Mammals"
Results_reptiles$Class <- "Reptiles"

Results_all_classes <- rbind(Results_amphibians, Results_birds, Results_mammals, Results_reptiles)

Results_all_classes$Use_intensity_1[Results_all_classes$Use_intensity_1=="Cannot decide"] <- NA
Results_all_classes$Use_intensity_2[Results_all_classes$Use_intensity_2=="Cannot decide"] <- NA
Results_all_classes$Use_intensity_1 %>%  unique
Results_all_classes$Use_intensity_2 %>%  unique

Results_all_classes$Land_use_2[grepl("econdary", Results_all_classes$Land_use_2)] <- "Secondary vegetation"
Results_all_classes$Land_use_1 %>%  unique
Results_all_classes$Land_use_2 %>%  unique
Results_all_classes$Land_use_2 <- factor(Results_all_classes$Land_use_2, levels=c("Primary vegetation",
                                                                                  "Secondary vegetation",
                                                                                  "Plantation forest",
                                                                                  "Pasture",
                                                                                  "Cropland",
                                                                                  "Urban"))

## sample sizes per use intensities
Results_all_classes %>%
  filter(Use_intensity_1==Use_intensity_2) %>%
  group_by(Class, Use_intensity_1, Land_use_2) %>%
  dplyr::summarise(Count=n())

## add tropical versus temperate divide
Reg <- read.csv("../../Results/dbFD_indices_for_analysis/dbFD_region_std.csv")
Temperate <- Reg$SS[Reg$Biome=="Temperate"] %>%  unique()
Tropical <- Reg$SS[Reg$Biome=="Tropical"] %>%  unique()

Results_all_classes$Biome <- NA
Results_all_classes$Biome[Results_all_classes$SS %in% Temperate] <- "Temperate"
Results_all_classes$Biome[Results_all_classes$SS %in% Tropical] <- "Tropical"

write.csv(Results_all_classes, "../../Results/Betadiv/Betadiv_results_all_classes.csv", row.names = FALSE)

#########################################################################################################################
#########################################################################################################################
library(StatisticalModels)
library(scales)
library(ggpubr)

## running models


Results_all_classes$Nestedness <- (Results_all_classes$Vol.PV - Results_all_classes$Vol.Inter)/Results_all_classes$Vol.PV*100
Results_all_classes$Turnover <- (Results_all_classes$Vol.LU2 - Results_all_classes$Vol.Inter)/Results_all_classes$Vol.LU2*100

## transformation for improving normality and bounding results
hist(Results_all_classes$Nestedness)
hist(Results_all_classes$Turnover)

Results_all_classes$asin_sqrt_loss <- asin(sqrt(Results_all_classes$Nestedness/100))
Results_all_classes$asin_sqrt_gain <- asin(sqrt(Results_all_classes$Turnover/100))
hist(Results_all_classes$asin_sqrt_loss)
hist(Results_all_classes$asin_sqrt_gain)

## fitting model to explain loss and gain by land use, use intensity, class and region

Results <- Results_all_classes %>%
  filter(Use_intensity_1==Use_intensity_2)

Results$Use_intensity_2 <- factor(Results$Use_intensity_2, levels = c("Minimal use", "Light use", "Intense use"))
Results$Land_use_2 <- as.character(Results$Land_use_2)
Results$Land_use_2[Results$Land_use_2=="Pasture"] <- "Agri"
Results$Land_use_2[Results$Land_use_2=="Cropland"] <- "Agri"
Results$Land_use_2 <- factor(Results$Land_use_2, levels=c("Primary vegetation",
                                                          "Secondary vegetation",
                                                          "Plantation forest",
                                                          "Agri",
                                                          "Urban"))

## sample sizes
length(unique(Results$SS[Results$Class=="Amphibians"]))
length(unique(Results$SS[Results$Class=="Birds"]))
length(unique(Results$SS[Results$Class=="Mammals"]))
length(unique(Results$SS[Results$Class=="Reptiles"]))

length(unique(Results$SS[Results$Class=="Amphibians" & Results$Biome=="Temperate"]))
length(unique(Results$SS[Results$Class=="Amphibians" & Results$Biome=="Tropical"]))
length(unique(Results$SS[Results$Class=="Birds" & Results$Biome=="Temperate"]))
length(unique(Results$SS[Results$Class=="Birds" & Results$Biome=="Tropical"]))
length(unique(Results$SS[Results$Class=="Mammals" & Results$Biome=="Temperate"]))
length(unique(Results$SS[Results$Class=="Mammals" & Results$Biome=="Tropical"]))
length(unique(Results$SS[Results$Class=="Reptiles" & Results$Biome=="Temperate"]))
length(unique(Results$SS[Results$Class=="Reptiles" & Results$Biome=="Tropical"]))

print(Results %>%
  filter(Class=="Birds", Biome=="Tropical") %>%
  group_by(Biome, Class, Land_use_2, Use_intensity_2) %>%
  dplyr::summarise(Count=n()))



## model selection (loss)
model4a <- GLMERSelect(modelData = Results,
                       responseVar = "asin_sqrt_loss",
                       fitFamily = "gaussian",
                       fixedFactors = c("Land_use_2", "Biome", "Class", "Use_intensity_2"),
                       fitInteractions=TRUE,
                       randomStruct = "(1|SS)",
                       verbose = TRUE)
model4a <- model4a$model

## model selection (gain)
model4b <- GLMERSelect(modelData = Results,
                       responseVar = "asin_sqrt_gain",
                       fitFamily = "gaussian",
                       fixedFactors = c("Land_use_2", "Biome", "Class", "Use_intensity_2"),
                       fitInteractions=TRUE,
                       randomStruct = "(1|SS)",
                       verbose = TRUE)
model4b <- model4b$model

## plotting predictions

Newdata <- expand.grid(Land_use_2=levels(Results$Land_use_2),
                       Biome=c("Temperate", "Tropical"),
                       Class=c("Amphibians", "Birds", "Mammals", "Reptiles"),
                       Use_intensity_2=c("Minimal use", "Light use", "Intense use"),
                       asin_sqrt_loss=0,
                       asin_sqrt_gain=0)


Predict_effects <- function(Newdata, Model, rescale) {

  preds <- sapply(X=1:1000, FUN=function(i){

    coefs <- mvrnorm(n = 1, mu = fixef(object=Model), Sigma = vcov(object=Model))
    mm <- model.matrix(terms(Model), Newdata)

    # drop coefs that couldn't be estimated
    print(setdiff(colnames(mm), names(coefs)))
    to_drop <- print(setdiff(colnames(mm), names(coefs)))
    mm <- as.data.frame(mm)
    mm <- mm[, -which(colnames(mm) %in% to_drop)]
    mm <- as.matrix(mm)

    y <- mm%*%coefs

    # backtranforming
    y <- sin(y)^2

    # rescaling
    if(rescale){
      # initialisation
      seq <- 1:5
      y[seq] <- y[seq]-y[seq[1]]
      # for loop to rescale all values
      for(i in 1:(nrow(Newdata)/5-1)){
        seq <- seq + 5
        y[seq] <- y[seq]-y[seq[1]]
      }
    }
    return(y)
  })

  preds <- data.frame(Median=apply(X=preds, MARGIN=1, FUN=median),
                      Upper=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.975),
                      Lower=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.025))
  preds <- cbind(preds, Newdata)
  return(preds)

}

## predictions and plotting for functional loss
cols <- scales::show_col(viridis(option = "plasma", n=8))
cbPalette <- c("#000000", viridis(option = "plasma", n=8)[c(2,3,5,6)])
show_col(cbPalette)

preds_loss <- Predict_effects(Newdata, model4a, TRUE)

# [1] "Land_use_2Urban:BiomeTropical"                          "Land_use_2Urban:ClassReptiles"
# [3] "Land_use_2Plantation forest:Use_intensity_2Intense use"

preds_loss$Median[preds_loss$Land_use_2=="Urban" & preds_loss$Biome=="Tropical"] <- NA
preds_loss$Upper[preds_loss$Land_use_2=="Urban" & preds_loss$Biome=="Tropical"] <- NA
preds_loss$Lower[preds_loss$Land_use_2=="Urban" & preds_loss$Biome=="Tropical"] <- NA

preds_loss$Median[preds_loss$Land_use_2=="Urban" & preds_loss$Class=="Reptiles"] <- NA
preds_loss$Upper[preds_loss$Land_use_2=="Urban" & preds_loss$Class=="Reptiles"] <- NA
preds_loss$Lower[preds_loss$Land_use_2=="Urban" & preds_loss$Class=="Reptiles"] <- NA

preds_loss$Median[preds_loss$Land_use_2=="Plantation forest" & preds_loss$Use_intensity_2=="Intense use"] <- NA
preds_loss$Upper[preds_loss$Land_use_2=="Plantation forest" & preds_loss$Use_intensity_2=="Intense use"] <- NA
preds_loss$Lower[preds_loss$Land_use_2=="Plantation forest" & preds_loss$Use_intensity_2=="Intense use"] <- NA

preds_loss$Median <- preds_loss$Median*100
preds_loss$Upper <- preds_loss$Upper*100
preds_loss$Lower <- preds_loss$Lower*100

ploss <- ggplot(preds_loss, aes(Land_use_2, Median, ymin = Lower, ymax = Upper, col=Land_use_2, group=Use_intensity_2, shape=Use_intensity_2)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  GGPoptions +
  facet_grid(Biome~Class, scales="free")+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  guides(colour=FALSE) +
  ylab("Difference to the reference (%)") + xlab("") +
  scale_x_discrete(labels=c("PV\n/PV", "PV\n/SV", "PV\n/PF","PV\n/AGR", "PV\n/UR")) +
  ggtitle("(a) Functional loss")

## predictions and plotting for functional gain

preds_gain <- Predict_effects(Newdata, model4b, TRUE)

# [1] "Land_use_2Urban:ClassReptiles" "Land_use_2Plantation forest:Use_intensity_2Intense use"
preds_gain$Median[preds_gain$Land_use_2=="Urban" & preds_gain$Class=="Reptiles"] <- NA
preds_gain$Upper[preds_gain$Land_use_2=="Urban" & preds_gain$Class=="Reptiles"] <- NA
preds_gain$Lower[preds_gain$Land_use_2=="Urban" & preds_gain$Class=="Reptiles"] <- NA

preds_gain$Median[preds_gain$Land_use_2=="Plantation forest" & preds_gain$Use_intensity_2=="Intense use"] <- NA
preds_gain$Upper[preds_gain$Land_use_2=="Plantation forest" & preds_gain$Use_intensity_2=="Intense use"] <- NA
preds_gain$Lower[preds_gain$Land_use_2=="Plantation forest" & preds_gain$Use_intensity_2=="Intense use"] <- NA

preds_gain$Median <- preds_gain$Median*100
preds_gain$Upper <- preds_gain$Upper*100
preds_gain$Lower <- preds_gain$Lower*100

pgain <- ggplot(preds_gain[preds_gain$Biome=="Temperate",], aes(Land_use_2, Median, ymin = Lower, ymax = Upper, col=Land_use_2, group=Use_intensity_2, shape=Use_intensity_2)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  GGPoptions +
  facet_grid(~Class, scales="free")+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  guides(colour=FALSE) +
  ylab("Difference to the reference (%)") +
  ggtitle("(b) Functional gain (similar effects in tropical and temperate regions)")+
  scale_x_discrete(labels=c("PV\n/PV", "PV\n/SV", "PV\n/PF","PV\n/AGR", "PV\n/UR")) + xlab("")

p_classes <- ggarrange(ploss, pgain, common.legend = TRUE, nrow=2, heights = c(0.6, 0.4))

ggsave(p_classes, filename="../Revisions/PDF_figures/Figure6.pdf",
       height=7, width=10)


## diagnostic plots
library(performance)
check_model(model4a, check = c("qq", "normality"), colors = c("red", "blue"))
check_model(model4b, check = c("qq", "normality"), colors = c("red", "blue"))

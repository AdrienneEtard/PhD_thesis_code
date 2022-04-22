## Script to assess the phylogenetic signal in Traits (continuous and categorical)

## Code is parallelised
library(parallel)
library(dplyr)
library(phytools)
source("Functions_Borges_et_al_2018_delta_statistic.R")

## START CLUSTER
Cluster <- makeCluster(6)

## EXCECUTE ANY PRE PROCESSING CODE NECESSARY
clusterEvalQ(Cluster, {
  library(dplyr)
  library(phytools)
  library(picante)
  library(geiger)
  library(ape)
  library(pbmcapply)
  library(pbapply)
})

## Preamble

# Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  return(Phylogeny)
}

# Function to apply for continuous traits: PhySignal
PhySignal <- function(Traitdata, Names, Phylo) {
  names(Traitdata) <- Names # names are binomial species names
  Signal <- phytools::phylosig(Phylo, Traitdata, method="lambda", test = TRUE) %>%
    unlist()
  return(Signal)
}

# Function to apply for categorical traits: PhySignal_Cat. Based on Borges et al 2018, in Bioinformatics
PhySignal_Cat <- function(Traitdata, Names, Phylo, n) {

 # browser()

  ## The Borges function does not work when branches have 0 length
  ## in that case, add a very small number to these branches
  #Phylo$edge.length[Phylo$edge.length==0] <- 10e-10

  # n = number of simulations
  ## For the current trait, match and prune phylogeny
  Traitdata <- as.character(Traitdata)
  names(Traitdata) <- Names
  Match <- picante::match.phylo.data(Phylo, Traitdata)
  Phylo <- Match$phy
  Trait <- Match$data
  rm(Match)

  ## Priors and parameters
  lambda0 <- 0.1   #rate parameter of the proposal
  se      <- 0.5   #standard deviation of the proposal
  sim     <- 10000 #number of iterations
  thin    <- 10    #we kept only each 10th iterate
  burn    <- 100   #100 iterates are burned-in

  ## Run the function (Borges et al, 2018, Bioinformatics: Measuring phylogenetic signal between categorical traits and phylogenies)
  ## with trycatch to avoid the process crashing when encountering errors.
  Delta <- tryCatch(expr={delta(Trait, Phylo, lambda0, se, sim, thin, burn)}, error = function(e) {NA})

  ## If the signal is not NA, then calculate a null distribution of delta values for the trait
  if(!is.na(Delta)) {

    # Generate randomised trait values - n times to generate a null distribution of delta -- stored in a list
    Func <- function(Trait){
      N <- length(Trait)
      L <- levels(as.factor(Trait))
      return(sample(L, size=N, replace=TRUE))
    }
    ListRandom <- lapply(rep(list(Trait), n), Func)

    Func_delta_toapply <- function(trait, tree, lambda0, se, sim, thin, burn) {
      Result <- tryCatch(expr = {delta(trait, tree, lambda0, se, sim, thin, burn)},
                         error=function(e){NA})
      return(Result)
    }

    Random_Delta <- pbmapply(FUN=Func_delta_toapply,
                             trait=ListRandom,
                             tree=rep(list(Phylo), n),
                             lambda0=rep(list(lambda0), n),
                             se=rep(list(se), n),
                             sim=rep(list(sim), n),
                             thin=rep(list(thin), n),
                             burn=rep(list(burn), n)) %>%
      as.data.frame()
    return(list(Delta=Delta, Delta0=Random_Delta))

  }

  else{return(Delta)}
}

# Read trait data, corrected for taxonomy
Amphibians <- read.csv("../../Data/GlobalGaps_traitdata/Amphibians.csv")
Birds <- read.csv("../../Data/GlobalGaps_traitdata/Birds.csv")
Mammals <- read.csv("../../Data/GlobalGaps_traitdata/Mammals.csv")
Reptiles <- read.csv("../../Data/GlobalGaps_traitdata/Reptiles.csv")

## Load phylogenies, the ones used in the Global Gaps article (consensus trees).
Phylo_Mammals <- read.nexus("../../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Mammals_complete_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Birds <- read.nexus("../../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Birds_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Amphibians <- read.nexus("../../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Amphibians_TreeAnnotator.nex")  %>% .Format_tiplabels()
Phylo_Reptiles <- read.nexus("../../Data/GlobalGaps_phylogenies/Consensus_Trees_TreeAnnotator/Reptiles_TreeAnnotator.nex")  %>% .Format_tiplabels()

# # Are there any polytomies in the trees?
# is.binary(Phylo_Birds)      # true
# is.binary(Phylo_Mammals)    # true
# is.binary(Phylo_Reptiles)   # true
# is.binary(Phylo_Amphibians) # true

# Traits
Continuous.Traits <- c("Body_mass_g",
                       #"Longevity_d",
                       "Litter_size",
                       #"Diet_breadth",
                       "Habitat_breadth_IUCN")

Categorical.Traits <- c("Specialisation",
                        "Trophic_level",
                        "Diel_activity")
                        #"Primary_diet")

# log-transform continuous traits
Transf <- function(Data, Traits_log10) {
  Data[, Traits_log10] <- log10(Data[, Traits_log10])
  Data$Habitat_breadth_IUCN <- sqrt(Data$Habitat_breadth_IUCN)
  return(Data)
}

Mammals <- Transf(Mammals, Traits_log10=c(Continuous.Traits[-which(Continuous.Traits=="Habitat_breadth_IUCN")],"Generation_length_d"))
Birds <- Transf(Birds, Traits_log10=c(Continuous.Traits[-which(Continuous.Traits=="Habitat_breadth_IUCN")],"Generation_length_d"))
Amphibians <- Transf(Amphibians, Traits_log10=c(Continuous.Traits[-which(Continuous.Traits=="Habitat_breadth_IUCN")],"Body_length_mm", "Maturity_d", "Max_longevity_d"))
Reptiles <- Transf(Reptiles, Traits_log10=c(Continuous.Traits[-which(Continuous.Traits=="Habitat_breadth_IUCN")], "Maturity_d", "Max_longevity_d", "Longevity_d"))

# Names
Names.Mammals <- Mammals$Best_guess_binomial
Names.Birds <- Birds$Best_guess_binomial
Names.Reptiles <- Reptiles$Best_guess_binomial
Names.Amphibians <- Amphibians$Best_guess_binomial

## Export variables in all clusters
clusterExport(cl=Cluster, list(".Format_tiplabels","PhySignal", "PhySignal_Cat",
                               "Mammals", "Birds", "Reptiles", "Amphibians",
                               "Phylo_Mammals", "Phylo_Birds", "Phylo_Reptiles", "Phylo_Amphibians",
                               "Continuous.Traits", "Categorical.Traits",
                               "Names.Mammals", "Names.Birds", "Names.Reptiles", "Names.Amphibians",
                               "nentropy", "lpalpha", "lpbeta", "mhalpha", "mhbeta", "emcmc", "ratematrix", "rtrait", "delta"), envir=environment())


## PARALLEL CALCULATIONS WITH parApply

## 1. Phylogenetic signal in continuous traits (takes up to 6 hours for birds)

print("starting estimations on continuous traits")

print(Sys.time())
start_time <- Sys.time()
Lambda_Mammals_continuous <- parApply(Cluster, Mammals[, c(Continuous.Traits,"Generation_length_d")], 2, PhySignal, Names=Names.Mammals, Phylo=Phylo_Mammals)
Lambda_Mammals_continuous <- data.table::rbindlist(lapply(Lambda_Mammals_continuous, function(x){return(as.data.frame(unlist(x)[1:4]))}))
Lambda_Mammals_continuous$Trait <- c("BM", "LCS", "HB", "GL")
end_time <- Sys.time()
timeMcont <- end_time-start_time
print(timeMcont)

write.csv(Lambda_Mammals_continuous, "../../Results/Physignal_traits/ContinuousMammals_log.csv",row.names = FALSE)

print(Sys.time())
start_time <- Sys.time()
TrAm <- c(Continuous.Traits, "Body_length_mm", "Maturity_d", "Max_longevity_d")
Lambda_Amphibians_continuous <- parApply(Cluster, Amphibians[, TrAm], 2, PhySignal, Names=Names.Amphibians, Phylo=Phylo_Amphibians)
Lambda_Amphibians_continuous <- data.table::rbindlist(lapply(Lambda_Amphibians_continuous, function(x){return(as.data.frame(unlist(x)[1:4]))}))
Lambda_Amphibians_continuous$Trait <- c("BM", "LCS", "HB", "BL", "MA", "ML")
end_time <- Sys.time()
timeAcont <- end_time-start_time
print(timeAcont)
write.csv(Lambda_Amphibians_continuous, "../../Results/Physignal_traits/ContinuousAmphibians_log.csv",row.names = FALSE)

print(Sys.time())
start_time <- Sys.time()
Lambda_Reptiles_continuous <- parApply(Cluster, Reptiles[, c(Continuous.Traits, "Maturity_d", "Max_longevity_d", "Longevity_d")], 2, PhySignal, Names=Names.Reptiles, Phylo=Phylo_Reptiles)
Lambda_Reptiles_continuous <- data.table::rbindlist(lapply(Lambda_Reptiles_continuous, function(x){return(as.data.frame(unlist(x)[1:4]))}))
Lambda_Reptiles_continuous$Trait <- c("BM", "LCS", "HB", "MA", "ML", "L")
end_time <- Sys.time()
timeRcont <- end_time-start_time
print(timeRcont)
write.csv(Lambda_Reptiles_continuous, "../../Results/Physignal_traits/ContinuousReptiles_log.csv",row.names = FALSE)


print(Sys.time())
start_time <- Sys.time()
Lambda_Birds_continuous <- parApply(Cluster, Birds[, c(Continuous.Traits, "Generation_length_d")], 2, PhySignal, Names=Names.Birds, Phylo=Phylo_Birds)
end_time <- Sys.time()
timeBcont <- end_time-start_time
print(timeBcont)
Lambda_Birds_continuous <- data.table::rbindlist(lapply(Lambda_Birds_continuous, function(x){return(as.data.frame(unlist(x)[1:4]))}))
Lambda_Birds_continuous$Trait <- c("BM", "LCS", "HB", "GL")
write.csv(Lambda_Birds_continuous, "../../Results/Physignal_traits/ContinuousBirds_log.csv",row.names = FALSE)


# ## 2. Phylogenetic signal in categorical traits
#
# print("starting estimations on categorical traits")
#
print(Sys.time())
start_time <- Sys.time()
delta_Mammals_categorical <- parApply(Cluster, Mammals[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Mammals, Phylo=Phylo_Mammals, n=50)
end_time <- Sys.time()
timeMcat <- end_time-start_time
print(timeMcat)
saveRDS(delta_Mammals_categorical, "../../Results/Physignal_traits/CategoricalMammals.rds")

print(Sys.time())
start_time <- Sys.time()
delta_Birds_categorical <- parApply(Cluster, Birds[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Birds, Phylo=Phylo_Birds, n=50)
end_time <- Sys.time()
timeBcat <- end_time-start_time
print(timeBcat)
saveRDS(delta_Birds_categorical, "../../Results/Physignal_traits/CategoricalBirds.rds")

print(Sys.time())
start_time <- Sys.time()
delta_Reptiles_categorical <- parApply(Cluster, Reptiles[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Reptiles, Phylo=Phylo_Reptiles, n=50)
end_time <- Sys.time()
timeRcat <- end_time-start_time
print(timeRcat)
saveRDS(delta_Reptiles_categorical, "../../Results/Physignal_traits/CategoricalReptiles.rds")

print(Sys.time())
start_time <- Sys.time()
delta_Amphibians_categorical <- parApply(Cluster, Amphibians[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Amphibians, Phylo=Phylo_Amphibians, n=50)
end_time <- Sys.time()
timeAcat <- end_time-start_time
print(timeAcat)
saveRDS(delta_Amphibians_categorical, "../../Results/Physignal_traits/CategoricalAmphibians.rds")

## DESTROY CLUSTER
stopCluster(Cluster)

## Testing whether phylogenetic signal is significant

## categorcial traits
wilcox.test(mu=delta_Mammals_categorical$Specialisation$Delta, delta_Mammals_categorical$Specialisation$Delta0$., alternative = "less")
wilcox.test(mu=delta_Mammals_categorical$Trophic_level$Delta, delta_Mammals_categorical$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=delta_Mammals_categorical$Diel_activity$Delta, delta_Mammals_categorical$Diel_activity$Delta0$., alternative = "less")

# shapiro.test(delta_Mammals_categorical$Specialisation$Delta0$.)
# shapiro.test(delta_Mammals_categorical$Trophic_level$Delta0$.)
# shapiro.test(delta_Mammals_categorical$Diel_activity$Delta0$.)

wilcox.test(mu=delta_Birds_categorical$Specialisation$Delta, delta_Birds_categorical$Specialisation$Delta0$., alternative = "less")
wilcox.test(mu=delta_Birds_categorical$Trophic_level$Delta, delta_Birds_categorical$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=delta_Birds_categorical$Diel_activity$Delta, delta_Birds_categorical$Diel_activity$Delta0$., alternative = "less")

hist(delta_Birds_categorical$Diel_activity$Delta0$.)
delta_Birds_categorical$Diel_activity$Delta

wilcox.test(mu=delta_Amphibians_categorical$Specialisation$Delta, delta_Amphibians_categorical$Specialisation$Delta0$., alternative = "less")
wilcox.test(mu=delta_Amphibians_categorical$Trophic_level$Delta, delta_Amphibians_categorical$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=delta_Amphibians_categorical$Diel_activity$Delta, delta_Amphibians_categorical$Diel_activity$Delta0$., alternative = "less")

wilcox.test(mu=delta_Reptiles_categorical$Specialisation$Delta, delta_Reptiles_categorical$Specialisation$Delta0$., alternative = "less")
wilcox.test(mu=delta_Reptiles_categorical$Trophic_level$Delta, delta_Reptiles_categorical$Trophic_level$Delta0$., alternative = "less")
wilcox.test(mu=delta_Reptiles_categorical$Diel_activity$Delta, delta_Reptiles_categorical$Diel_activity$Delta0$., alternative = "less")

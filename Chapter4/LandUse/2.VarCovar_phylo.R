## Prepare inverse of variance - covariance matrix from phylogeny, to add as random effects in mcmcglmm models

library(geiger)
library(phytools)
library(dplyr)
library(MCMCglmm)

# # # # #  f u n c t i o n s  # # # # #

# format tip labels of phylogenies

.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) stringr::word(x, 1, 2))
  return(Phylogeny)}

# keep phylogenetic tips for species in PREDICTS only

.DropTipsPredicts <- function(Phylo, Predicts, Class) {

  ## tips to drop
  if(!is.null(Class)) {
    x <- unique(Predicts$Best_guess_binomial[Predicts$Class==Class])
  }
  
  else{x <- unique(Predicts$Best_guess_binomial)}
  to_drop <- dplyr::setdiff(Phylo$tip.label, x) %>% unlist()
  
  ## drop tips to drop
  Phylo <- ape::drop.tip(Phylo, to_drop)
  return(Phylo)
}


# # # # # # # # # # # # # # # # # # # # #

## PREDICTS (already preprocessed with mergeSites)
Predicts <- readRDS("../Results/Predicts_merged_sites.rds")

## vertebrate classes separately
Mammals <- read.tree("../Data/Phylogenies_corrected/Mammals.nwk") %>% .Format_tiplabels()
Birds <- read.tree("../Data/Phylogenies_corrected/Birds.nwk") %>% .Format_tiplabels()
Reptiles <- read.tree("../Data/Phylogenies_corrected/Reptiles.nwk") %>% .Format_tiplabels()
Amphibians <- read.tree("../Data/Phylogenies_corrected/Amphibians.nwk") %>% .Format_tiplabels()
All <- read.tree("../Data/Phylogeny_all_vertebrates/Vertebrates.nwk") %>% .Format_tiplabels()


length(All$tip.label)
length(setdiff(unique(Predicts$Best_guess_binomial), All$tip.label)) # 102 species 


## prune phylogenies to keep only species in PREDICTS
Mammals <- .DropTipsPredicts(Mammals, Predicts, "Mammalia") 
Birds <- .DropTipsPredicts(Birds, Predicts, "Aves") 
Amphibians <- .DropTipsPredicts(Amphibians, Predicts, "Amphibia") 
Reptiles <- .DropTipsPredicts(Reptiles, Predicts, "Reptilia") 
AllV <- .DropTipsPredicts(All, Predicts, Class=NULL) 

## species in PREDICTS that are not represented in the phylogenies
length(setdiff(unique(Predicts$Best_guess_binomial[Predicts$Class=="Mammalia"]), Mammals$tip.label)) 
length(setdiff(unique(Predicts$Best_guess_binomial[Predicts$Class=="Aves"]), Birds$tip.label)) 
length(setdiff(unique(Predicts$Best_guess_binomial[Predicts$Class=="Amphibia"]), Amphibians$tip.label)) 
length(setdiff(unique(Predicts$Best_guess_binomial[Predicts$Class=="Reptilia"]), Reptiles$tip.label))

## overall 102 species in PREDICTS are not represented in the phylogenies
length(setdiff(unique(Predicts$Best_guess_binomial), AllV$tip.label))

## where BL=0 add a very small number
Mammals$edge.length[Mammals$edge.length==0] <- 10e-10
Birds$edge.length[Birds$edge.length==0] <- 10e-10
Reptiles$edge.length[Reptiles$edge.length==0] <- 10e-10
Amphibians$edge.length[Amphibians$edge.length==0] <- 10e-10
AllV$edge.length[AllV$edge.length==0] <- 10e-10


# # # # # # Generate variance - covariance matrix and take the inverse

AinvMammals <-inverseA(Mammals,nodes="TIPS",scale=TRUE)$Ainv
AinvBirds <-inverseA(Birds,nodes="TIPS",scale=TRUE)$Ainv
AinvReptiles <-inverseA(Reptiles,nodes="TIPS",scale=TRUE)$Ainv
AinvAmphibians <-inverseA(Amphibians,nodes="TIPS",scale=TRUE)$Ainv
AinvAll <-inverseA(AllV,nodes="TIPS",scale=TRUE)$Ainv

## save
saveRDS(AinvMammals, "../Results/Inv_var_covar_phylo/Mammals.rds")
saveRDS(AinvAmphibians, "../Results/Inv_var_covar_phylo/Amphibians.rds")
saveRDS(AinvBirds, "../Results/Inv_var_covar_phylo/Birds.rds")
saveRDS(AinvReptiles, "../Results/Inv_var_covar_phylo/Reptiles.rds")
saveRDS(AinvAll, "../Results/Inv_var_covar_phylo/Vertebrates.rds")



# ## variance-covariance matrix: expected trait variances and trait covariances among species under Brownian-motion evolution
# vcv.Mammals <- ape::vcv.phylo(Mammals)
# vcv.Birds <- ape::vcv.phylo(Birds)
# vcv.Amphibians <- ape::vcv.phylo(Amphibians)
# vcv.Reptiles <- ape::vcv.phylo(Reptiles)
# 
# ## solve for the inverse of the var-covar matrix
# vcvInv.Mammals <- solve(vcv.Mammals)
# vcvInv.Birds <- solve(vcv.Birds)
# vcvInv.Amphibians <- solve(vcv.Amphibians)
# vcvInv.Reptiles <- solve(vcv.Reptiles)






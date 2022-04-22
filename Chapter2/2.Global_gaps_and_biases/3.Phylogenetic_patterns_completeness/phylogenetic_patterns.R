## Patterns in missing trait values with regards to the phylogenetic position of the species

# plot per traits (percent and number of species with missing values within families for each trait)
# plot completeness (within family median)

# setwd("../../Gaps_biases_trait_data/Code/3.Phylogenetic_patterns_completeness/")

library(phytools)
library(dplyr)
library(ape)
library(picante)
library(geiger)
library(ggplot2)
library(RColorBrewer)
library(ggtree)
library(grid)
library(ggpubr)
`%nin%` <- Negate(`%in%`)
source("Functions_patterns_missing_values.R")

# ## Load phylogenies -- these are from the Time Tree of Life
# PhyloAmphibians <- read.newick("../../Data/Phylogenies_processed/Amphibians.nwk") %>% .Format_tiplabels()
# PhyloMammals <- read.newick("../../Data/Phylogenies_processed/Mammals.nwk")  %>% .Format_tiplabels()
# PhyloBirds <- read.newick("../../Data/Phylogenies_processed/Birds.nwk")  %>% .Format_tiplabels() 
# PhyloReptiles <- read.newick("../../Data/Phylogenies_processed/Reptiles.nwk")  %>% .Format_tiplabels() 

## class specific phylogenies from VertLife, processed with TreeAnnotator to get a consensus tree

# for TreeAnnotator - to format into nexus file, run the following:
# PhyloAmphibians2 <-  read.tree("../../Data/Phylogenies_class_specific/Amphibians/Amphibians_unzipped/amph_shl_new_Posterior_7238.1000.trees")
# writeNexus(PhyloAmphibians2, "../../Data/Phylogenies_class_specific/Amphibians/amphibians1000trees.nex")
# PhyloReptiles2 <-  read.tree("../../Data/Phylogenies_class_specific/Reptiles/Squamates_unzipped/squam_shl_new_Posterior_9755.1000.trees")
# writeNexus(PhyloReptiles2, "../../Data/Phylogenies_class_specific/Reptiles/reptiles1000trees.nex")
# PhyloBirds2 <- read.tree("../../Data/Phylogenies_class_specific/Birds/Stage2_MayrAll_Ericson_set10_decisive/Stage2_MayrAll_Ericson_set10_decisive.tre")
# writeNexus(PhyloBirds2, "../../Data/Phylogenies_class_specific/Birds/birds1000trees.nex")

## Mammals (Phylacine_1.2)
Phy_MammalsComplete <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Mammals_complete_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_MammalsSmall <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Mammals_small_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_Amphibians <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Amphibians_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_Reptiles <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Reptiles_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_Birds <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Birds_TreeAnnotator.nex") %>%
 .Format_tiplabels()

## Load traits with completeness

Tr <- readRDS("../../Results/Traits_to_map/traits_completeness_V2.rds")
Tr_amphibians <- Tr[["Amphibians"]]
Tr_birds <- Tr[["Birds"]]
Tr_mammals <- Tr[["Mammals"]]
Tr_reptiles <- Tr[["Reptiles"]]

# ## Phylogenetic signal in trait completeness?
PhySignal <- function(Traitdata, Names, Phylo) {
  names(Traitdata) <- Names # names are binomial species names
  Signal <- phytools::phylosig(Phylo, Traitdata, method="lambda", test = TRUE) %>%
    unlist()
  
  # # http://blog.phytools.org/2012/11/testing-for-pagels-10.html
  # # comparing to Brownian motion
  # # get likelihood for bm model
  # browser()
  # 
  # tree <- Phylo
  # tree$mapped.edge <- matrix(tree$edge.length, nrow(tree$edge),1,dimnames=list(NULL,"1"))
  # bm.logL<-brownie.lite(tree,Traitdata)$logL1
  # # conduct hypothesis test using chi-square
  # LR <- 2*(bm.logL-Signal$logL)
  # P <- pchisq(LR,df=1,lower.tail=F)
  
 # return(list(Signal=Signal, p_BM=P))
  return(Signal)
  
  }

# lambda_amphibians <- PhySignal(Tr_amphibians$completeness, Tr_amphibians$Best_guess_binomial, Phy_Amphibians)
# lambda_reptiles <- PhySignal(Tr_reptiles$completeness, Tr_reptiles$Best_guess_binomial, Phy_Reptiles)
# 
# lambda_mammals_C <- PhySignal(Tr_mammals$completeness, Tr_mammals$Best_guess_binomial, Phy_MammalsComplete)
# lambda_mammals_S <- PhySignal(Tr_mammals$completeness, Tr_mammals$Best_guess_binomial, Phy_MammalsSmall)
# 
# lambda_birds <- PhySignal(Tr_birds$completeness, Tr_birds$Best_guess_binomial, Phy_Birds)
# 
# L <- cbind(lambda_amphibians, lambda_mammals_C, lambda_mammals_S, lambda_reptiles, lambda_birds)
# write.csv(L, "../../Results/lambda_completeness_V2.csv", row.names = TRUE)

# with Tree of Life phylogeny
Signal_TTOL <- read.csv("../../Results/lambda_completeness.csv")

# with VertLife
Signal_VertLife <- read.csv("../../Results/lambda_completeness_V2.csv")


## Phylogenetic signal using a distribution of trees

Bootstrap_lambda <- function(Phy, TraitData, nCores){
  
  PhySignal <- function(Traitdata, Phylo) {
    Signal <- phytools::phylosig(Phylo, Traitdata, method="lambda", test = TRUE)
    return(unlist(Signal))
  }
  
  Phy <- sample(Phy, 50, replace = FALSE)
  Phy <- lapply(Phy, .Format_tiplabels)
  
  # Signal <- list()
  # for(i in 1:100){
  #   Signal[[i]] <- PhySignal(Traitdata=TraitData$completeness, TraitData$Best_guess_binomial, Phylo=Phy[[i]])
  #   print(i)
  # }
  
  ## parallelise
  # browser()
  # CompletenessList <- rep(list(TraitData$completeness),100)
  # NamesList <- rep(list(TraitData$Best_guess_binomial), 100)
  
  CompletenessList <- TraitData$completeness
  NamesList <- TraitData$Best_guess_binomial
  names(CompletenessList) <- NamesList # names are binomial species names
  #CompletenessList <- rep(list(CompletenessList),100)
  
   cl <- parallel::makeCluster(nCores)
  
  parallel::clusterExport(cl = cl,
                          varlist = 
                            c("PhySignal",
                              "Phy", 
                              "CompletenessList"),
                              #"NamesList"),
                          envir = environment())
  
  
  Signal <- ParallelLogger::clusterApply(cluster = cl,
                                         fun = PhySignal,
                                         x = Phy,
                                         Traitdata=CompletenessList,
                                         #Names=NamesList,
                                         progressBar = TRUE)
  
  #browser()
  
  ParallelLogger::stopCluster(cl)
  
  Signal.df <- lapply(Signal, function(x){
    y <- as.data.frame(unlist(x[1:4])) %>%  t() %>%  as.data.frame()
    colnames(y) <- c("Lambda", "LogL", "Log0", "pvalue")# %>%  t() 
    return(y)})
  Signal.df <- data.table::rbindlist(Signal.df)
  return(Signal.df)
}

gc()
memory.limit(size=500000)

DistPhyAmphibians <- read.nexus("../../Data/Phylogenies_class_specific/Amphibians/amphibians1000trees.nex")
SignalDistAmphibians <- Bootstrap_lambda(Phy=DistPhyAmphibians, Tr_amphibians, nCores=6)
write.csv(SignalDistAmphibians, "../../Results/Lambda_dist_amphibians.csv", row.names = FALSE)

DistPhyMammals <- read.nexus("../../Data/Phylogenies_class_specific/Mammals/doi_10.5061_dryad.bp26v20__v1/Complete_phylogeny.nex")
SignalDistMammals <- Bootstrap_lambda(Phy=DistPhyMammals, Tr_mammals, nCores = 6)
write.csv(SignalDistMammals, "../../Results/Lambda_dist_mammals.csv", row.names = FALSE)

DistPhyBirds <- read.nexus("../../Data/Phylogenies_class_specific/Birds/birds1000trees.nex")
SignalDistBirds1 <- Bootstrap_lambda(Phy=DistPhyBirds, Tr_birds, nCores = 6)
write.csv(SignalDistBirds1, "../../Results/Lambda_dist_birds1.csv", row.names = FALSE)

SignalDistBirds2 <- Bootstrap_lambda(Phy=DistPhyBirds, Tr_birds, nCores = 6)
SignalDistBirds <- rbind(SignalDistBirds1, SignalDistBirds2)
write.csv(SignalDistBirds, "../../Results/Lambda_dist_birds.csv", row.names = FALSE)

DistPhyReptiles <- read.nexus("../../Data/Phylogenies_class_specific/Reptiles/reptiles1000trees.nex")
PhyR1 <- sample(DistPhyReptiles, 100, replace = FALSE)
PhyR11 <- PhyR1[1:50]
PhyR12 <- PhyR1[51:100]

SignalDistReptiles1 <- Bootstrap_lambda(Phy=PhyR11, Tr_reptiles, nCores = 6)
write.csv(SignalDistReptiles1, "../../Results/Lambda_dist_reptiles1.csv", row.names = FALSE)

SignalDistReptiles2 <- Bootstrap_lambda(Phy=PhyR12, Tr_reptiles, nCores = 6)
SignalDistReptiles <- rbind(SignalDistReptiles1, SignalDistReptiles2)
write.csv(SignalDistReptiles, "../../Results/Lambda_dist_reptiles.csv", row.names = FALSE)


SignalDistAmphibians <- read.csv("../../Results/Lambda_dist_amphibians.csv")[51:100,]
SignalDistMammals <- read.csv("../../Results/Lambda_dist_mammals.csv")[51:100,]
SignalDistBirds <- read.csv("../../Results/Lambda_dist_birds.csv")[51:100,]
SignalDistReptiles <- read.csv("../../Results/Lambda_dist_reptiles.csv")

hist(SignalDistAmphibians$Lambda, xlab="")
hist(SignalDistMammals$Lambda)
hist(SignalDistBirds$Lambda)
hist(SignalDistReptiles$Lambda)

shapiro.test(SignalDistAmphibians$Lambda)
shapiro.test(SignalDistMammals$Lambda)
shapiro.test(SignalDistBirds$Lambda)
shapiro.test(SignalDistReptiles$Lambda)

gmodels::ci(SignalDistAmphibians$Lambda)
gmodels::ci(SignalDistMammals$Lambda)
gmodels::ci(SignalDistBirds$Lambda)
gmodels::ci(SignalDistReptiles$Lambda)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Phy_MammalsComplete <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Mammals_complete_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_MammalsSmall <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Mammals_small_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_Amphibians <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Amphibians_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_Reptiles <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Reptiles_TreeAnnotator.nex") %>% 
  .Format_tiplabels()
Phy_Birds <- read.nexus("../../Results/Consensus_Trees_TreeAnnotator/Birds_TreeAnnotator.nex") %>%
  .Format_tiplabels()


## 1. Calculate within-family median and mean completeness
Med_comp_amphibians <- Completeness_families(Tr_amphibians, "median")
Med_comp_reptiles <- Completeness_families(Tr_reptiles, "median")
Med_comp_birds <- Completeness_families(Tr_birds,  "median")
Med_comp_mammals <- Completeness_families(Tr_mammals,  "median")

# Mea_comp_amphibians <- Completeness_families(Tr_amphibians, "mean")
# Mea_comp_reptiles <- Completeness_families(Tr_reptiles, "mean")
# Mea_comp_birds <- Completeness_families(Tr_birds,  "mean")
# Mea_comp_mammals <- Completeness_families(Tr_mammals,  "mean")


## variance 
VarcompAmphibians <-  Completeness_families(Tr_amphibians, "variance")
VarcompReptiles <-  Completeness_families(Tr_reptiles, "variance")
VarcompMammals <-  Completeness_families(Tr_mammals, "variance")
VarcompBirds <-  Completeness_families(Tr_birds, "variance")

## SD
VarcompAmphibians$SD <- sqrt(VarcompAmphibians$Variance)
VarcompReptiles$SD <- sqrt(VarcompReptiles$Variance)
VarcompMammals$SD <- sqrt(VarcompMammals$Variance)
VarcompBirds$SD <- sqrt(VarcompBirds$Variance)


pdf(file = "../../Results/Phylogenetic_plots/SD_VS_SR.pdf", width=8, height=7, pointsize = 13)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,1), oma=c(2,1,2,1))
par(mfrow=c(2,2))
plot(VarcompAmphibians$SD~ log(VarcompAmphibians$SR), pch=19, xlab="log(Within-family species richness)", ylab="SD(Within-family Completeness) (%)", main="Amphibians")
plot(VarcompReptiles$SD~ log(VarcompReptiles$SR), pch=19, xlab="log(Within-family species richness)", ylab="SD(Within-family Completeness) (%)", main="Reptiles")
plot(VarcompMammals$SD~ log(VarcompMammals$SR), pch=19, xlab="log(Within-family species richness)", ylab="SD(Within-family Completeness) (%)", main="Mammals")
plot(VarcompBirds$SD~ log(VarcompBirds$SR), pch=19, xlab="log(Within-family species richness)", ylab="SD(Within-family Completeness) (%)", main="Birds")
dev.off()

summary(lm(VarcompAmphibians$Variance~ log(VarcompAmphibians$SR)))
summary(lm(VarcompReptiles$Variance~ log(VarcompReptiles$SR)))
summary(lm(VarcompMammals$Variance~ log(VarcompMammals$SR)))
summary(lm(VarcompBirds$Variance~ log(VarcompBirds$SR)))


## sample sizes
intersect(Tr_amphibians$Best_guess_binomial, Phy_Amphibians$tip.label) %>%  length()
intersect(Tr_reptiles$Best_guess_binomial, Phy_Reptiles$tip.label) %>%  length()
intersect(Tr_mammals$Best_guess_binomial, Phy_MammalsComplete$tip.label) %>%  length()
intersect(Tr_birds$Best_guess_binomial, Phy_Birds$tip.label) %>%  length()



## 2. Plot within-family median completeness against phylogenetic tree


# # Amphibians
# pdf(file="../../Results/Phylogenetic_plots/MEDIAN_Amphibians_completeness.pdf", width=5, height=7, pointsize=9)
# p <- Plot_NA_patterns(PhyloAmphibians, Med_comp_amphibians, Tr_amphibians, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 75, cols=p$cols, title="",  digits=2, lims = c(0,100), x=45, y=10, subtitle="\nwithin family completeness (%)")
# box("outer")
# dev.off()
# 
# pdf(file="../../Results/Phylogenetic_plots/MEAN_Amphibians_completeness.pdf", width=5, height=7, pointsize=9)
# p <- Plot_NA_patterns(PhyloAmphibians, Mea_comp_amphibians, Tr_amphibians, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 75, cols=p$cols, title="",  digits=2, lims = c(0,100), x=45, y=10, subtitle="\nwithin family completeness (%)")
# box("outer")
# dev.off()
# 
# 
# # Reptiles
# pdf(file="../../Results/Phylogenetic_plots/MEDIAN_Reptiles_completeness.pdf", width=5, height=7, pointsize=9)
# Plot_NA_patterns(PhyloReptiles, Med_comp_reptiles, Tr_reptiles, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 75, cols=p$cols, title="",  digits=2, lims = c(0,100), x=30, y=10, subtitle="\nwithin family completeness (%)")
# box("outer")
# dev.off()
# 
# pdf(file="../../Results/Phylogenetic_plots/MEAN_Reptiles_completeness.pdf", width=5, height=7, pointsize=9)
# Plot_NA_patterns(PhyloReptiles, Mea_comp_reptiles, Tr_reptiles, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 75, cols=p$cols, title="",  digits=2, lims = c(0,100), x=30, y=10, subtitle="\nwithin family completeness (%)")
# box("outer")
# dev.off()
# 
# # Birds
# pdf(file="../../Results/Phylogenetic_plots/MEDIAN_Birds_completeness.pdf", width=10, height=20, pointsize=9)
# Plot_NA_patterns(PhyloBirds, Med_comp_birds, Tr_birds, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 50, cols=p$cols, title="",  digits=2, lims = c(0,100), x=10, y=3, subtitle="\nwithin family completeness (%)", fsize=2)
# box("outer")
# dev.off()
# 
# pdf(file="../../Results/Phylogenetic_plots/MEAN_Birds_completeness.pdf", width=10, height=20, pointsize=9)
# Plot_NA_patterns(PhyloBirds, Mea_comp_birds, Tr_birds, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 50, cols=p$cols, title="",  digits=2, lims = c(0,100), x=10, y=3, subtitle="\nwithin family completeness (%)", fsize=2)
# box("outer")
# dev.off()
# 
# # Mammals
# pdf(file="../../Results/Phylogenetic_plots/MEDIAN_Mammals_completeness.pdf", width=10, height=15, pointsize=9)
# Plot_NA_patterns(PhyloMammals, Med_comp_mammals, Tr_mammals, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 50, cols=p$cols, title="",  digits=2, lims = c(0,100), x=20, y=5, subtitle="\nwithin family completeness (%)", fsize=1.5)
# box("outer")
# dev.off()
# 
# pdf(file="../../Results/Phylogenetic_plots/MEAN_Mammals_completeness.pdf", width=10, height=15, pointsize=9)
# Plot_NA_patterns(PhyloMammals, Mea_comp_mammals, Tr_mammals, TRUE, TRUE)
# add.color.bar(prompt=FALSE, 50, cols=p$cols, title="",  digits=2, lims = c(0,100), x=20, y=5, subtitle="\nwithin family completeness (%)", fsize=1.5)
# box("outer")
# dev.off()

# # # # # # # # # #

# ## All together without the tip labels -- median
# pdf(file="../../Results/Phylogenetic_plots/Median_completeness_all.pdf", width=7, height=5, pointsize=12)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(6,1,2,1))
# par(xpd=NA)
# par(mfrow=c(2,2))
# Plot_NA_patterns(PhyloMammals, Med_comp_mammals, Tr_mammals, FALSE, TRUE); title("(a) Mammals", line=0.5, adj=0, font=2)
# Plot_NA_patterns(PhyloBirds, Med_comp_birds, Tr_birds, FALSE, TRUE); title("(b) Birds", line=0.5, adj=0, font=2)
# p <- Plot_NA_patterns(PhyloAmphibians, Med_comp_amphibians, Tr_amphibians, FALSE, TRUE); title("(c) Amphibians", line=0.5, adj=0, font=2)
# Plot_NA_patterns(PhyloReptiles, Med_comp_reptiles, Tr_reptiles, FALSE, TRUE); title("(d) Reptiles", line=0.5, adj=0, font=2)
# 
# add.color.bar(prompt=FALSE, 150, cols=p$cols, title="",  digits=2, lims = c(0,100), x=-50, y=-20, subtitle="\nWithin-family median completeness (%)", fsize=1.2)
# box("outer")
# dev.off()
# 
# ## All together without the tip labels -- mean
# #pdf(file="../../Results/Phylogenetic_plots/Mean_completeness_all.pdf", width=7, height=5, pointsize=12)
# png(filename="../../Results/Phylogenetic_plots/Mean_completeness_all.png", width=7, height=5, pointsize=12, res=400, units="in")
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(6,1,2,1))
# par(xpd=NA)
# par(mfrow=c(2,2))
# Plot_NA_patterns(PhyloMammals, Mea_comp_mammals, Tr_mammals, FALSE, Which = NULL, TRUE); title("(a) Mammals", line=0.5, adj=0, font=2)
# Plot_NA_patterns(PhyloBirds, Mea_comp_birds, Tr_birds, FALSE, Which = NULL, TRUE); title("(b) Birds", line=0.5, adj=0, font=2)
# p <- Plot_NA_patterns(PhyloAmphibians, Mea_comp_amphibians, Tr_amphibians, FALSE, Which = NULL, TRUE); title("(c) Amphibians", line=0.5, adj=0, font=2)
# Plot_NA_patterns(PhyloReptiles, Mea_comp_reptiles, Tr_reptiles, FALSE, Which = NULL, TRUE); title("(d) Reptiles", line=0.5, adj=0, font=2)
# add.color.bar(prompt=FALSE, 150, cols=p$cols, title="",  digits=2, lims = c(0,100), x=-50, y=-20, subtitle="\nWithin-family mean completeness (%)", fsize=1.2)
# box("outer")
# dev.off()
# 
# 
# # with certain clades highlighted
# 
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(6,1,2,1))
# par(xpd=NA)
# par(mfrow=c(2,2))
# Plot_NA_patterns(PhyloMammals, Mea_comp_mammals, Tr_mammals, TRUE, Which="Hippopotamidae",TRUE); title("(a) Mammals", line=0.5, adj=0, font=2)
# Plot_NA_patterns(PhyloBirds, Mea_comp_birds, Tr_birds, TRUE,Which="Falcunculidae", TRUE); title("(b) Birds", line=0.5, adj=0, font=2)
# p <- Plot_NA_patterns(PhyloAmphibians, Mea_comp_amphibians, Tr_amphibians, TRUE,Which=c("Ceratobatrachidae", "Ascaphidae"), TRUE); title("(c) Amphibians", line=0.5, adj=0, font=2)
# Plot_NA_patterns(PhyloReptiles, Mea_comp_reptiles, Tr_reptiles, TRUE, Which=c("Anomalepididae","Calamariidae"),TRUE); title("(d) Reptiles", line=0.5, adj=0, font=2)
# add.color.bar(prompt=FALSE, 150, cols=p$cols, title="",  digits=2, lims = c(0,100), x=-50, y=-20, subtitle="\nWithin-family mean completeness (%)", fsize=1.2)
# box("outer")
# dev.off()


# #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #   #  
# 
# ## 2. Within family trait coverage - both as the percentages of missing trait value within the value and the number of species (log10?)
# 
# ## Amphibians
# 
# X <- c("Body_size", "Life_span_proxy", "Litter_size","Trophic_level","Diel_activity",'Habitat_breadth_IUCN', "Specialisation")
# Titles <- c("A. BS", "B. LS", "C. LCS", "D. TL", "E. DA", "F. HB", "G. Spe")
# 
# pdf(file="../../Results/Plots/Phylogeny_missing_values/Amphibians_coverage.pdf", width=6, height=4, pointsize=9)
# par(mfrow=c(2,4)); par(xpd=NA)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(1,1,2,2))
# for (i in 1:7) {
#   Results <- PercentNA_families(Tr_amphibians, X[i])
#   p <- Plot_NA_patterns(PhyloAmphibians, Results, Tr_amphibians, FALSE, FALSE)
#   title(main=Titles[i], adj=0, line=0.5)
# }
# 
# Results <- Family_rep(Tr_amphibians, TRUE)
# p2 <- Plot_NA_patterns(PhyloAmphibians, Results, Tr_amphibians, FALSE, FALSE)
# title(main="H. % rep", adj=0, line=0.5)
# 
# add.color.bar(prompt=FALSE, 300, cols=p$cols, title="",  digits=1, lims = c(0,100), x=600, y=100, subtitle="\nwithin family missingness (%): \nplots A - K", fsize=1.5)
# add.color.bar(prompt=FALSE, 300, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=600, y=50, subtitle="\nrepresentation (log-10): \nplot L", fsize=1.5)
# box("outer")
# dev.off()
# 
# 
# ## Mammals
# X <- c("Body_mass_g", "Generation_length_d", "Range_size_m2","Trophic_level",'Diel_activity',
#        "Primary_diet", "Habitat_breadth_IUCN", "Specialisation", "Diet_breadth", "Adult_svl_cm", "Litter_size", "Longevity_d")
# Titles <- c("A. BM", "B. GL", "C. RS", "D. TL", "E. DA", "F. PD", "G. HB", "H. Sp", "I. DB", "J. BL", "K. LCS", "L. Longevity")
# 
# pdf(file="../../Results/Plots/Phylogeny_missing_values/Mammals_coverage.pdf", width=6, height=4, pointsize=9)
# par(mfrow=c(3,5)); par(xpd=NA)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(1,1,2,1))
# for (i in 1:12) {
#   Results <- PercentNA_families(Mammals, X[i])
#   p <- Plot_NA_patterns(PhyloMammals, Results, Mammals, FALSE, FALSE)
#   title(main=Titles[i], adj=0, line=0.5)
# }
# Results <- Family_rep(Mammals, TRUE)
# p2 <- Plot_NA_patterns(PhyloMammals, Results, Mammals, FALSE, FALSE)
# title(main="L. % rep", adj=0, line=0.5)
# 
# add.color.bar(prompt=FALSE, 200, cols=p$cols, title="",  digits=1, lims = c(0,100), x=300, y=130, subtitle="\nwithin family missingness (%):\nplots A - K", fsize=1.5)
# add.color.bar(prompt=FALSE, 200, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=300, y=50, subtitle="\nrepresentation (log-10):\nplot L", fsize=1.5)
# box("outer")
# dev.off()
# 
# 
# ## Birds
# X <- c("Habitat_breadth_IUCN", "Specialisation", "Body_mass_g","Range_size_m2","Primary_diet",'Diel_activity',
#        "Generation_length_d", "Trophic_level", "Diet_breadth", "Litter_size", "Longevity_d")
# Titles <- c("A. HB", "B. Sp", "C. BM", "D. RS", "E. PD", "F. DA", "G. GL", "H. TL", "I. DB", "J. LS", "K. Longevity")
# 
# pdf(file="../../Results/Plots/Phylogeny_missing_values/Birds_coverage.pdf", width=6, height=4, pointsize=9)
# par(mfrow=c(3,4)); par(xpd=NA)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(1,1,2,20))
# for (i in 1:11) {
#   Results <- PercentNA_families(Birds, X[i])
#   p <- Plot_NA_patterns(PhyloBirds, Results, Birds, FALSE, FALSE)
#   title(main=Titles[i], adj=0, line=0.5)
# }
# Results <- Family_rep(Birds, TRUE)
# p2 <- Plot_NA_patterns(PhyloBirds, Results, Birds, FALSE, FALSE)
# title(main="L. % rep", adj=0, line=0.5)
# 
# add.color.bar(prompt=FALSE, 100, cols=p$cols, title="",  digits=1, lims = c(0,100), x=200, y=500, subtitle="\n within family missingness (%): \n plots A - K", fsize=1.5)
# add.color.bar(prompt=FALSE, 100, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=200, y=200, subtitle="\n representation (log-10): \n plot L", fsize=1.5)
# box("outer")
# dev.off()
# 
# 
# ## Reptiles
# X <- c("Body_mass_g", "Range_size_m2", "Litter_size","Habitat_breadth_IUCN","Specialisation",'Diel_activity',
#        "Longevity_d", "Trophic_level", "Adult_svl_cm", "Maturity_d")
# Titles <- c("A. BM", "B. RS", "C. LCS", "D. HB", "E. Sp", "F. DA", "G. Longevity", "H. TL", "J. BL", "K. SM")
# 
# pdf(file="../../Results/Plots/Phylogeny_missing_values/Reptiles_coverage.pdf", width=6, height=4, pointsize=9)
# par(mfrow=c(3,4)); par(xpd=NA)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(3,1,2,5))
# for (i in 1:10) {
#   Results <- PercentNA_families(Reptiles, X[i])
#   p <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, FALSE, FALSE)
#   title(main=Titles[i], adj=0, line=0.5)
# }
# Results <- Family_rep(Reptiles, TRUE)
# p2 <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, FALSE, FALSE)
# title(main="L. % rep", adj=0, line=0.5)
# 
# add.color.bar(prompt=FALSE, 150, cols=p$cols, title="",  digits=1, lims = c(0,100), x=330, y=60, subtitle="\nwithin family missingness (%): \nplots A - K", fsize=1.5)
# add.color.bar(prompt=FALSE, 150, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=330, y=20, subtitle="\nrepresentation (log-10): \nplot L", fsize=1.5)
# box("outer")
# dev.off()


# # # # # Alternative
# 
# # for Life span 
# par(mfrow=c(2,2), mar = c(5,1,5,1), oma=c(1,1,2,2)); par(xpd=NA)
# Plot_NA_patterns(PhyloMammals, PercentNA_families(Tr_mammals, "Life_span_proxy"), Tr_mammals, FALSE, FALSE)
# title("Mammals",adj=0, line=0.5)
# Plot_NA_patterns(PhyloBirds, PercentNA_families(Tr_birds, "Life_span_proxy"), Tr_birds, FALSE, FALSE)
# title("Birds",adj=0, line=0.5)
# Plot_NA_patterns(PhyloAmphibians, PercentNA_families(Tr_amphibians, "Life_span_proxy"), Tr_amphibians, FALSE, FALSE)
# title("Amphibians",adj=0, line=0.5)
# Plot_NA_patterns(PhyloReptiles, PercentNA_families(Tr_reptiles, "Life_span_proxy"), Tr_reptiles, FALSE, FALSE)
# title("Reptiles",adj=0, line=0.5)

## Amphibians

# split.screen(c(1,2))
# screen(1)
# split.screen(c(3,4))
# screen(3)
# par(oma=c(4,0,2,4), las=1)
# for (i in 1:11) {
#   screen(i+2)
#   Results <- PercentNA_families(Amphibians, X[i])
#   p <- Plot_NA_patterns(PhyloAmphibians, Results, Amphibians, FALSE, FALSE)
#   mtext(text = Titles[i], line =-12, at=100, font = 2)
# }
# 
# screen(2)
# par(oma=c(4,2,1.5,0), las=1)
# Results <- Family_rep(Amphibians)
# p2 <- Plot_NA_patterns(PhyloAmphibians, Results, Amphibians, TRUE, FALSE)
# 
# par(xpd=NA)
# add.color.bar(prompt=TRUE, 400, cols=p$cols, title="",  digits=2, lims = c(0,100), x=20, y=5, subtitle="within family coverage (%)")
# add.color.bar(prompt=TRUE, 500, cols=p2$cols, title="",  digits=0, lims = range(Results$Percent), x=20, y=5, subtitle="% representation")
# 
# close.screen()
# box("outer")
# mtext("B", side=3, font=2, cex=1.5)
# mtext("A", side=3, font=2, cex=1.5, at=-650)
# 
# dev.off()


## Reptiles

# X <- c("Body_mass_g", "Range_size_m2", "Litter_size","Habitat_breadth_IUCN","Specialisation",'Diel_activity',
#        "Longevity_d", "Trophic_level", "Adult_svl_cm", "Maturity_d")
# Titles <- c("BM", "RS", "LCS", "HB", "Sp", "DA", "L", "TL", "BL", "M")
# 
# pdf(file="../../Results/Plots/Phylogeny_missing_values/Reptiles_coverage.pdf", width=7, height=8, pointsize=9)
# 
# split.screen(c(1,2))
# screen(1)
# split.screen(c(3,4))
# screen(3)
# par(oma=c(4,0,2,4), las=1)
# for (i in 1:11) {
#   screen(i+2)
#   Results <- PercentNA_families(Reptiles, X[i])
#   p <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, FALSE, FALSE)
#   mtext(text = Titles[i], line =-12, at=80, font = 2)
# }
# 
# screen(2)
# par(oma=c(4,2,1.5,0), las=1)
# Results <- Family_rep(Reptiles)
# p2 <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, TRUE, FALSE)
# 
# par(xpd=NA)
# add.color.bar(prompt=FALSE, 300, cols=p$cols, title="",  digits=2, lims = c(0,100), x=-500, y=-5, subtitle="within family missing values (%)")
# add.color.bar(prompt=FALSE, 300, cols=p2$cols, title="",  digits=0, lims = range(Results$Percent), x=150, y=-5, subtitle="% representation")
# 
# close.screen()
# box("outer")
# mtext("B", side=3, font=2, cex=1.5)
# mtext("A", side=3, font=2, cex=1.5, at=-400)
# 
# dev.off()



##  ##  ##  ##  ##  ##  ##  ##
library(rcartocolor)
library(RColorBrewer)
display.brewer.all(colorblindFriendly = T)

# # # # circular phylogenies

# amphibians and reptiles

cols <- brewer.pal(n = 11,name = "RdYlBu")
cols <- cols[c(1,2,3,4,9,10,11)]
#cols <- c(cols, "#000000")
scales::show_col(cols)
breaks=c(0,25,50,60,70,80,90,100)

# amphibians
# X = Plot_NA_patterns_cir_nofam(PhyloAmphibians, Med_comp_amphibians, Tr_amphibians, TRUE, TRUE, Size = 5.3)
X2 = Plot_NA_patterns_cir_nofam(Phy_Amphibians, Med_comp_amphibians, Tr_amphibians, TRUE, TRUE, Size = 5.3)
# x <- X$x
x2 <- X2$x
# xf <- merge(x, x2, by="Taxa", all = TRUE)

## plot with phyloheny from Tree of Life
# Xplot <- X$Tree + #theme_tree(bgcolor = "black", fgcolor = "black") +
#   scale_colour_gradientn(colours=cols, 
#                          breaks=breaks, 
#                          labels = format(breaks))+ ggtree::geom_rootpoint()
# p1 <- Xplot + ggtitle("(a) Amphibians") + xlim(c(0,700))
# rm(p1)

## plot
X2plot <- X2$Tree + #theme_tree(bgcolor = "black", fgcolor = "black") +
  scale_colour_gradientn(colours=cols, 
                         #breaks=breaks, 
                         values=breaks/100)+ ggtree::geom_rootpoint() + theme(legend.position = "right") 

p1 <- X2plot + ggtitle("(a) Amphibians") + xlim(c(0,700))

rm(X, X1, X2, x, x2, xf)

# reptiles

#X = Plot_NA_patterns_cir_nofam(PhyloReptiles, Med_comp_reptiles, Tr_reptiles, TRUE, TRUE, Size = 5.3)
X2 = Plot_NA_patterns_cir_nofam(Phy_Reptiles, Med_comp_reptiles, Tr_reptiles, TRUE, TRUE, Size = 5.3)

# x <- X$x
x2 <- X2$x
# xf <- merge(x, x2, by="Taxa", all = TRUE)

X2 <- X2$Tree + #theme_tree(bgcolor = "black", fgcolor = "black") +
  scale_colour_gradientn(colours=cols, 
                         values=breaks/100)+ ggtree::geom_rootpoint()+ theme(legend.position = "right")
p2 <- X2 + ggtitle("(b) Reptiles") + xlim(c(0,500))

## Figure for manuscript

fig1 <- ggpubr::ggarrange(p1+ 
                            theme(plot.title = element_text(size = 30, face = "bold")),
                          p2+ 
                            theme(plot.title = element_text(size = 30, face = "bold")), 
                          nrow=2, common.legend = TRUE, legend = "right")
ggsave(fig1, filename="../../Results/Phylogenetic_plots/Circular_herptiles_final.pdf", width = 18, height = 24)

# #pdf("../../Results/Phylogenetic_plots/Circular_herptiles.pdf", width = 18, height = 24)
# pdf("../../Results/Phylogenetic_plots/Circular_herptiles_V2.pdf", width = 18, height = 24)
# 
# # plot the legend with base graphics
# plot.new(); par(xpd=NA)
# 
# legend(0.8,0.7,c(
#   "0 - 25%","25 - 50%","50 - 60%","60 - 70%", "70 - 80%",
#   "80 - 90%","90 - 100%"), fill = cols, cex=2, 
#   bty="n", title = "Median completeness:")
# 
# # create an apporpriate viewport
# vp <- viewport(height=unit(1, "npc"), width=unit(0.7, "npc"), just=c("left","top"), y=1, x=0)
# 
# # plot the ggplot using the print command
# print(fig1, vp=vp)
# 
# dev.off()



# # # # # # # # # # # # # #
# mammals

breaks <- c(0,25,50,75,80,90,100)
cols <- brewer.pal(n = 11,name = "RdYlBu")
cols <- cols[c(1,8,9,10,11)]
cols <- c(cols, "#000000")
scales::show_col(cols)

X = Plot_NA_patterns_cir_nofam(Phy_MammalsComplete, Med_comp_mammals, Tr_mammals, TRUE, TRUE, Size = 3)
X$x
X3 <- X$Tree + #theme_tree(bgcolor = "black", fgcolor = "black") +
  scale_colour_gradientn(colours=cols, 
                         values=breaks/100)+ ggtree::geom_rootpoint() + theme(legend.position = "right") +
  labs(col = "Median completeness (%)")
p3 <- X3 + xlim(c(0,300)) 
ggsave(p3, filename="../../Results/Phylogenetic_plots/Circular_mammals_final.pdf", width = 15 ,height = 15)

# plot the legend with base graphics
pdf("../../Results/Phylogenetic_plots/Circular_mammals_final.pdf", width = 15 ,height = 15)
plot.new(); par(xpd=NA)
legend(0.8,0.6,c(
  "0 - 25%","25 - 50%","50 - 75%","75 - 80%", "80 - 90%",
  "90 - 100%"),fill = cols,cex=1,
  bty="n", title = "Median completeness:")
# create an apporpriate viewport
vp <- viewport(height=unit(1, "npc"), width=unit(0.7, "npc"), just=c("left","top"), y=1, x=0)
# plot the ggplot using the print command
print(p3, vp=vp)
dev.off()



# birds
X = Plot_NA_patterns_cir_nofam(Phy_Birds, Med_comp_birds, Tr_birds, TRUE, TRUE, Size = 1)
X4 <- X$Tree +# theme_tree(bgcolor = "black", fgcolor = "black") +
  scale_colour_gradientn(colours=cols, 
                         values=breaks/100)+ ggtree::geom_rootpoint()+ theme(legend.position = "right")+
  labs(col = "Median completeness (%)")

p4 <- X4 + xlim(c(0,120))
ggsave(p4, filename="../../Results/Phylogenetic_plots/Circular_birds_final.pdf", width = 17 ,height = 17)

# plot the legend with b`ase graphics
pdf("../../Results/Phylogenetic_plots/Circular_birds_V2.pdf", width = 17, height = 17)
plot.new(); par(xpd=NA)
legend(0.8,0.6,c(
  "0 - 25%","25 - 90%","90 - 92%","92 - 94%", "94 - 96%",
  "96 - 98%","98 - 100%"),fill = cols,cex=0.3, 
  bty="n", title = "Median completeness:")
# create an apporpriate viewport
vp <- viewport(height=unit(1, "npc"), width=unit(0.7, "npc"), just=c("left","top"), y=1, x=0)
# plot the ggplot using the print command
print(p4, vp=vp)
dev.off()

###### Variance plot for amphibians


# hist(VarcompAmphibians$SD, breaks=50)
# hist(VarcompReptiles$SD, breaks=30)
# 
# safe_colorblind_palette <- c("#CC6677", "#117733", "#332288", "#AA4499",   "#DDCC77", "#88CCEE",
#                              "#999933", "#661100", "#6699CC", "#888888")
# scales::show_col(safe_colorblind_palette)
# cols <- c("#661100","#882255","#AA4499", "#117733", "#999933", "#DDCC77",  "#44AA99", "#6699CC","#888888","#332288", "#000000") 
# scales::show_col(cols)


breaks=c(0,10,20,30,40,50,60,100)
cols <- brewer.pal(n = 11,name = "RdYlBu")
cols <- cols[c(1,2,3,4,9,10,11)]
#cols <- c(cols, "#000000")
scales::show_col(cols)
# amphibians
VarcompAmphibians$Variance <- VarcompAmphibians$SD
X2 = Plot_NA_patterns_cir_nofam(Phy_Amphibians, VarcompAmphibians, Tr_amphibians, TRUE, TRUE, Size = 5.3)
x2 <- X2$x
X2plot <- X2$Tree + #theme_tree(bgcolor = "black", fgcolor = "black") +
  scale_colour_gradientn(colours=cols, 
                         values=breaks/100)+ ggtree::geom_rootpoint() + theme(legend.position = "right")
p1 <- X2plot + ggtitle("(a) Amphibians") + xlim(c(0,650))
rm(x2, X2)
# reptiles
VarcompReptiles$Variance <- VarcompReptiles$SD
X2 = Plot_NA_patterns_cir_nofam(Phy_Reptiles, VarcompReptiles, Tr_reptiles, TRUE, TRUE, Size = 5.3)
x2 <- X2$x
X2 <- X2$Tree + #theme_tree(bgcolor = "black", fgcolor = "black") +
  scale_colour_gradientn(colours=cols, 
                         values=breaks/100)+ ggtree::geom_rootpoint()
p2 <- X2 + ggtitle("(b) Reptiles") + xlim(c(0,400))+ theme(legend.position = "right")

## Figure for manuscript

fig1 <- ggpubr::ggarrange(p1+ 
                            theme(plot.title = element_text(size = 30, face = "bold")),
                          p2+ 
                            theme(plot.title = element_text(size = 30, face = "bold")), 
                          nrow=2, common.legend = TRUE, legend="right")

ggsave(fig1, filename="../../Results/Phylogenetic_plots/Circular_herptiles_SD.pdf", width = 18, height = 24)
# #pdf("../../Results/Phylogenetic_plots/Circular_herptiles.pdf", width = 18, height = 24)
# pdf("../../Results/Phylogenetic_plots/Circular_herptiles_SD.pdf", width = 18, height = 24)
# 
# # plot the legend with base graphics
# plot.new(); par(xpd=NA)
# breaks=c(0,5,10,15,20,25,30,35,40,45,50,60)
# legend(0.7,0.7,c(
#   "0 - 5","5 - 10","10 - 15", "15 - 20","20 - 25","25 - 30", "30 - 35", "35 - 40",
#   "40 - 45","45 - 50", "50 - 60"),fill = cols,cex=2, 
#   bty="n", title = "Standard deviation\n in median completeness (%):")
# 
# 
# # create an apporpriate viewport
# vp <- viewport(height=unit(1, "npc"), width=unit(0.7, "npc"), just=c("left","top"), y=1, x=0)
# 
# # plot the ggplot using the print command
# print(fig1, vp=vp)
# 
# dev.off()

median(VarcompAmphibians$Variance, na.rm=TRUE)
median(VarcompReptiles$Variance, na.rm=TRUE)
mean(VarcompAmphibians$Variance, na.rm=TRUE)
mean(VarcompReptiles$Variance, na.rm=TRUE)














###### all species
safe_colorblind_palette <- c("#CC6677", "#117733", "#332288", "#AA4499",   "#DDCC77", "#88CCEE",
                             "#999933", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)

cols <- c("#661100","#882255","#AA4499", "#117733", "#44AA99", "#6699CC","#332288") 

scales::show_col(cols)


breaks=c(0,25,50,60,70,80,90,100)

pA <- Plot_phylo_species(Phy_Amphibians, Tr_amphibians)+ 
  scale_colour_gradientn(colours=cols,
                         breaks=breaks,
                         labels = format(breaks))+ ggtree::geom_rootpoint() +
  theme(legend.position="right")


pR <- Plot_phylo_species(Phy_Reptiles, Tr_reptiles)+ 
  scale_colour_gradientn(colours=cols,
                         breaks=breaks,
                         labels = format(breaks))+ ggtree::geom_rootpoint() +
  theme(legend.position="right")

pM <- Plot_phylo_species(Phy_MammalsComplete, Tr_mammals)+ 
  scale_colour_gradientn(colours=cols,
                         breaks=breaks,
                         labels = format(breaks))+ ggtree::geom_rootpoint() +
  theme(legend.position="right")

pB <- Plot_phylo_species(Phy_Birds, Tr_birds)+ 
  scale_colour_gradientn(colours=cols,
                         breaks=breaks,
                         labels = format(breaks))+ ggtree::geom_rootpoint() +
  theme(legend.position="right")

pall <- ggpubr::ggarrange(pA,pR,pM,pB, common.legend = TRUE)

ggsave(pall, filename="../../Results/Phylogenetic_plots/plot_splevel.pdf", width = 18, height = 24)




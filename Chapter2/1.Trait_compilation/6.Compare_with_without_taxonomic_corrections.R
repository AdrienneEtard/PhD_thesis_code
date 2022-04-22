## Comparison of basic trait coverage with VS without taxonomic correction
## and of number of species for which traits were compiled

# # #   P  R  E  A  M  B  L  E

X <- c("dplyr", "ggplot2", "ggpubr", "grid", "cowplot", "phytools", "picante")
lapply(X, library, character.only=TRUE); rm(X)
source("Comparison_with_without_taxonomic_corrections_functions.R")
source("Functions_for_phylogenies.R")
`%nin%` <- Negate(`%in%`)

## Load trait files

# No taxonomic correction
UN_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
UN_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
UN_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
UN_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

# With taxonomic correction
C_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
C_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
C_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
C_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

## Load PREDICTS
UN_Predicts <- readRDS("../../Data/PREDICTS_database.rds") %>% 
  filter(Class %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"))

C_Predicts <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

## Phylogenies after taxonomic corrections, without random additions -- to modify EV1
Phylo_Mammals <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloMammals.nwk")  %>% .Format_tiplabels()
Phylo_Birds <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloBirds.nwk")  %>% .Format_tiplabels() 
Phylo_Amphibians <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloAmphibians.nwk")  %>% .Format_tiplabels() 
Phylo_Reptiles <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloReptiles.nwk")  %>% .Format_tiplabels() 

EV_corrected <- function(TraitDF, Phylo) {
  
  EV <- TraitDF$EV_1
  Names=TraitDF$Best_guess_binomial
  names(EV) <- Names
  Match = match.phylo.data(Phylo, EV)
  Phylo <- Match$phy
  EV <- Match$data
  TraitDF$EV_1 <- NA
  
  for (i in 1:length(EV)) {
    TraitDF$EV_1[TraitDF$Best_guess_binomial==names(EV)[i]] <- 1
  }
   return(TraitDF)
}

C_Amphibians_NoAdd <- EV_corrected(C_Amphibians, Phylo_Amphibians)
C_Mammals_NoAdd <- EV_corrected(C_Mammals, Phylo_Mammals)
C_Reptiles_NoAdd <- EV_corrected(C_Reptiles, Phylo_Reptiles)
C_Birds_NoAdd <- EV_corrected(C_Birds, Phylo_Birds)


## 1. Differences in species number

# Amphibians
Delta_Number(UN_Amphibians, C_Amphibians, "Amphibia", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Amphibians, C_Amphibians, "Amphibia", "Predicts", C_Predicts, UN_Predicts)

# Mammals
Delta_Number(UN_Mammals, C_Mammals, "Mammalia", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Mammals, C_Mammals, "Mammalia", "Predicts", C_Predicts, UN_Predicts)

# Birds
Delta_Number(UN_Birds, C_Birds, "Aves", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Birds, C_Birds, "Aves", "Predicts", C_Predicts, UN_Predicts)

# Reptiles
Delta_Number(UN_Reptiles, C_Reptiles, "Reptilia", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Reptiles, C_Reptiles, "Reptilia", "Predicts", C_Predicts, UN_Predicts)


## 2. Plotting coverage

# 2.1. % Species representation in phylogenies
CovPhyloAll <- Phylo_Delta(C_Mammals_NoAdd, C_Birds_NoAdd, C_Reptiles_NoAdd, C_Amphibians_NoAdd, C_Predicts, FALSE,
                           UN_Mammals, UN_Birds, UN_Reptiles, UN_Amphibians, UN_Predicts)

CovPhyloPredicts <- Phylo_Delta(C_Mammals_NoAdd, C_Birds_NoAdd, C_Reptiles_NoAdd, C_Amphibians_NoAdd, C_Predicts, TRUE,
                                UN_Mammals, UN_Birds, UN_Reptiles, UN_Amphibians, UN_Predicts)

p <- PlotPhyloCov(CovPhyloAll,CovPhyloPredicts, 15)
ggsave(p, file="../../Results/Plots/Coverage/Phylogenies/Species_representation.pdf", width=8, height=4.2)
ggsave(p, file="../../Results/Plots_CBER_talk_20119/phylo_species_representation.png", width=8, height=4.2, dpi=1000)


# 2.2. % Trait coverage

# TraitsBirds <-  c("Body_mass_g", "Adult_svl_cm", "Maturity_d", "Longevity_d",
#                   "Litter_size", "Range_size_m2", "Diel_activity", "Primary_diet","Trophic_level","Specialisation",
#                   "Habitat_breadth_IUCN")
# 
# TraitsMammals <-  c("Body_mass_g", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm",
#                     "Generation_length_d", "Maturity_d", "Longevity_d", "AFR_d",
#                     "Litter_size", "Range_size_m2", "Diel_activity", "Primary_diet","Trophic_level","Specialisation",
#                     "Habitat_breadth_IUCN")
# 
# TraitsAmphibians <-  c("Body_mass_g", "Body_length_mm","Svl_length_mm", "Maturity_d", "Longevity_d",
#                        "Litter_size", "Range_size_m2", "Diel_activity","Primary_diet", "Trophic_level","Specialisation",
#                        "Habitat_breadth_IUCN")
# 
# TraitsReptiles <-  c("Body_mass_g", "Adult_svl_cm", "Maturity_d", "Longevity_d",
#                      "Litter_size", "Range_size_m2", "Diel_activity", "Primary_diet","Trophic_level", "Specialisation",
#                      "Habitat_breadth_IUCN")
#
#
#
# ## Delta in trait coverage before and after taxonomic correction -- FOR ALL TRAITS
#
# # For all species
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(2,2))
# Plot.Delta.Cov(UN_Mammals, C_Mammals, TraitsMammals, FALSE, UN_Predicts, C_Predicts, "Mammals")
# Plot.Delta.Cov(UN_Birds, C_Birds, TraitsBirds, FALSE, UN_Predicts, C_Predicts, "Birds")
# Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TraitsReptiles, FALSE, UN_Predicts, C_Predicts, "Reptiles")
# Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TraitsAmphibians, FALSE, UN_Predicts, C_Predicts, "Amphibians")
# 
# 
# ## For PREDICTS species
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(2,2))
# Plot.Delta.Cov(UN_Mammals, C_Mammals, TraitsMammals, TRUE, UN_Predicts, C_Predicts, "Mammals")
# Plot.Delta.Cov(UN_Birds, C_Birds, TraitsBirds, TRUE, UN_Predicts, C_Predicts, "Birds")
# Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TraitsReptiles, TRUE, UN_Predicts, C_Predicts, "Reptiles")
# Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TraitsAmphibians, TRUE, UN_Predicts, C_Predicts, "Amphibians")



## Delta in trait coverage before and after taxonomic correction -- FOR TARGET TRAITS

TargetTraits <- c("Body_mass_g",
                  "Longevity_d",
                  "Litter_size", 
                  "Diel_activity",
                  "Trophic_level",
                  "Diet_breadth",
                  "Specialisation",
                  "Habitat_breadth_IUCN",
                  "Primary_diet")

# For all species

pdf(file="../../Results/Plots/Coverage/Target_traits/All_species.pdf", width=7, height=5.5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(7,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TargetTraits, FALSE, UN_Predicts, C_Predicts, "A")
Plot.Delta.Cov(UN_Birds, C_Birds, TargetTraits, FALSE, UN_Predicts, C_Predicts, "B")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TargetTraits, FALSE, UN_Predicts, C_Predicts, "C")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TargetTraits, FALSE, UN_Predicts, C_Predicts, "D")

mtext(at=50, line=-11, "% coverage", cex=0.8)
mtext(at=-140, line=-11, "% coverage", cex=0.8)

par(xpd=NA)
legend(x=-110, y=-6, title="Trait coverage",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("#F8766D", "#00BFC4"),  bty="n")

dev.off()


# For PREDICTS species
pdf(file="../../Results/Plots/Coverage/Target_traits/Predicts.pdf", width=7, height=5.5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TargetTraits, TRUE, UN_Predicts, C_Predicts, "A")
Plot.Delta.Cov(UN_Birds, C_Birds, TargetTraits, TRUE, UN_Predicts, C_Predicts, "B")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TargetTraits, TRUE, UN_Predicts, C_Predicts, "C")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TargetTraits, TRUE, UN_Predicts, C_Predicts, "D")

mtext(at=50, line=-11, "% coverage", cex=0.8)
mtext(at=-140, line=-11, "% coverage", cex=0.8)

par(xpd=NA)
legend(x=-110, y=-6, title="Trait coverage",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("#F8766D", "#00BFC4"),  bty="n")

dev.off()



## Delta in trait coverage before and after taxonomic correction -- FOR TRAITS USED IN MISSFOREST IMPUTATIONS only (predictor traits)

Traits_cont <-  c("Body_mass_g", "Longevity_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Diet_breadth")
Traits_cat <- c("Specialisation", "Diel_activity","Trophic_level", "Primary_diet")

TMammalsI <- c(Traits_cont, Traits_cat, "Generation_length_d", "Adult_svl_cm")
TBirdsI <- c(Traits_cont, Traits_cat, "Generation_length_d")
TReptilesI <- c(Traits_cont, Traits_cat, "Adult_svl_cm", "Maturity_d")
TAmphibiansI <- c(Traits_cont, Traits_cat, "Body_length_mm")

# For all species

pdf(file="../../Results/Plots/Coverage/Predictor_traits/All_species.pdf", width=7, height=6, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TMammalsI, FALSE, UN_Predicts, C_Predicts, "A")
Plot.Delta.Cov(UN_Birds, C_Birds, TBirdsI, FALSE, UN_Predicts, C_Predicts, "B")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TReptilesI, FALSE, UN_Predicts, C_Predicts, "C")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TAmphibiansI, FALSE, UN_Predicts, C_Predicts, "D")

mtext(at=50, line=-12.5, "% coverage", cex=0.8)
mtext(at=-140, line=-12.5, "% coverage", cex=0.8)

par(xpd=NA)
legend(x=-110, y=-6, title="Trait coverage",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("#F8766D", "#00BFC4"), bty="n")
dev.off()


# For PREDICTS species
pdf(file="../../Results/Plots/Coverage/Predictor_traits/Predicts.pdf", width=7, height=6, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TMammalsI, TRUE, UN_Predicts, C_Predicts, "A")
Plot.Delta.Cov(UN_Birds, C_Birds, TBirdsI, TRUE, UN_Predicts, C_Predicts, "B")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TReptilesI, TRUE, UN_Predicts, C_Predicts, "C")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TAmphibiansI, TRUE, UN_Predicts, C_Predicts, "D")

mtext(at=50, line=-12.5, "% coverage", cex=0.8)
mtext(at=-140, line=-12.5, "% coverage", cex=0.8)

par(xpd=NA)
legend(x=-110, y=-6, title="Trait coverage",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("#F8766D", "#00BFC4"),  bty="n")
dev.off()




# png(file="../../Results/Plots_CBER_talk_20119/Coverage.png", width =7, height=6, units="in", res=600)
# par(tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
# par(mfrow=c(2,2))
# Plot.Delta.Cov(UN_Mammals, C_Mammals, TMammalsI, FALSE, UN_Predicts, C_Predicts, "Mammals")
# Plot.Delta.Cov(UN_Birds, C_Birds, TBirdsI, FALSE, UN_Predicts, C_Predicts, "Birds")
# Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TReptilesI, FALSE, UN_Predicts, C_Predicts, "Reptiles")
# Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TAmphibiansI, FALSE, UN_Predicts, C_Predicts, "Amphibians")
# 
# mtext(at=50, line=-12.5, "% coverage", cex=0.8)
# mtext(at=-140, line=-12.5, "% coverage", cex=0.8)
# 
# par(xpd=NA)
# legend(x=-110, y=-6, title="Trait coverage",
#        legend = c("Without taxonomic correction", "With taxonomic correction"),
#        fill = c("#F8766D", "#00BFC4"), bty="n")
# dev.off()



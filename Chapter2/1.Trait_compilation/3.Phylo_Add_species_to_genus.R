## Phylogenies: adding species to their genus, when possible, at the root (phytools::add.species.to.genus function)
## Relevant for species in the trait dataset that are not represented in the phylogeny initially

# # # # P r e a m b l e 

X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ape", "treeio", "ngram","phylobase", "PDcalc")
lapply(X, library, character.only=TRUE); rm(X)
source("Functions_for_phylogenies.R")


# ## Load phylogenies, with added species (code commented below if these datasets have to be generated again)
# 
## These phylogenies are resolved and ultrametric but contain branches of length zero that need to be corrected for

C_Mammals <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Mammals.nwk")
C_Birds <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Birds.nwk")
C_Amphibians <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Amphibians.nwk")
C_Reptiles <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Reptiles.nwk")

BL_0_n(C_Mammals)
BL_0_n(C_Birds)
BL_0_n(C_Amphibians)
BL_0_n(C_Reptiles)

UN_Mammals <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Mammals.nwk")
UN_Birds <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Birds.nwk")
UN_Amphibians <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Amphibians.nwk")
UN_Reptiles <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Reptiles.nwk")

BL_0_n(UN_Mammals)
BL_0_n(UN_Birds)
BL_0_n(UN_Amphibians)
BL_0_n(UN_Reptiles)


# Original phylogenies: nor all resolved, all ultrametric, contain branches of length zero
Mammals_original <- read.newick("../../Data/Mammals/Phylogenies/TTOL_mammals_smoothed_interpolated.nwk")
Amphibians_original <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk")
Birds_original <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk")
Reptiles_original <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk")

BL_0_n(Mammals_original)
BL_0_n(Birds_original)
BL_0_n(Amphibians_original)
BL_0_n(Reptiles_original)

## how many 0 BL have been added due to these species addition?
BL_0_n(C_Mammals) - BL_0_n(Mammals_original)
BL_0_n(C_Birds)-BL_0_n(Birds_original)
BL_0_n(C_Amphibians)-BL_0_n(Amphibians_original)
BL_0_n(C_Reptiles)-BL_0_n(Reptiles_original)


## BUT there are some branch length of 0 still because of addition of species to the root that generates polytomies
## Plus there are some branches of length 0 in the original trees.


# # # # # SCRIPT TO RUN TO RANDOMLY ADD SPECIES # # # #
# 
# ## Load trait data
# 
# # No taxonomic correction
# UN_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Amphibians.csv")
# UN_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Birds.csv")
# UN_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Mammals.csv")
# UN_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Reptiles.csv")
# 
# # With taxonomic correction
# C_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Amphibians.csv")
# C_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Birds.csv")
# C_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Mammals.csv")
# C_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Reptiles.csv")
# 
# ## Load phylogenies
# 
# # uncorrected - all ultrametric, with BL, all binary (resolved) except Amphibians
# PhyloMammal_UN <- read.newick("../../Data/Mammals/Phylogenies/TTOL_mammals_smoothed_interpolated.nwk")
# PhyloAmphibian_UN <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk")
# PhyloBird_UN <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk")
# PhyloReptile_UN <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk")
# 
# # corrected - all ultrametric, with BL, all binary (resolved) except Amphibians
# Phylo_Mammals_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloMammals.nwk")
# Phylo_Amphibians_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloAmphibians.nwk")
# Phylo_Reptiles_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloReptiles.nwk")
# Phylo_Birds_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloBirds.nwk")
# 
# 
# ## 2. Attach species to their genus if possible, using the function phytools::add.species.to.genus
# ## parameter set to "root" rather than "random" because with "random" branch length information is lost for all phylogenies except birds.
# 
# # Test <- Attach_Species_Phylo(UN_Mammals, PhyloMammal_UN, Where="root")
# # Test2 <- Attach_Species_Phylo(C_Mammals, Phylo_Mammals_C, Where="root")
# 
# # for uncorrected datasets
# PhyloMammal_UN <- Attach_Species_Phylo(UN_Mammals, PhyloMammal_UN,  Where="root")
# PhyloAmphibian_UN <- Attach_Species_Phylo(UN_Amphibians, PhyloAmphibian_UN,  Where="root")
# PhyloBird_UN <- Attach_Species_Phylo(UN_Birds, PhyloBird_UN,  Where="root")
# PhyloReptile_UN <- Attach_Species_Phylo(UN_Reptiles, PhyloReptile_UN,  Where="root")
# 
# # for corrected datasets
# Phylo_Mammals_C <- Attach_Species_Phylo(C_Mammals, Phylo_Mammals_C,  Where="root")
# Phylo_Amphibians_C <- Attach_Species_Phylo(C_Amphibians, Phylo_Amphibians_C,  Where="root")
# Phylo_Birds_C <- Attach_Species_Phylo(C_Birds, Phylo_Birds_C,  Where="root")
# Phylo_Reptiles_C <- Attach_Species_Phylo(C_Reptiles, Phylo_Reptiles_C,  Where="root")
# 
# 
# ## 3. Resolve polytomies with bifurcatr (https://rdrr.io/github/davidnipperess/PDcalc/src/R/bifurcatr.R), when trees are not binary
# ## Function to resolve polytomies while adding branch length information to both edges that have 0 BL and their descendants
# ## Adapted from: Rangel TF, Colwell RK, Graves GR, Fucíková K, Rahbek C, & Diniz-Filho JAF (2015). Phylogenetic uncertainty revisited:
# #  Implications for ecological analyses. Evolution 69:1301-1312.
# 
# 
# # uncorrected trees (ultrametric at this stage)
# if(!is.binary(PhyloMammal_UN)) {PhyloMammal_UN <- bifurcatr(PhyloMammal_UN, runs=1)}
# if(!is.binary(PhyloAmphibian_UN)) {PhyloAmphibian_UN <- bifurcatr(PhyloAmphibian_UN, runs=1)}
# if(!is.binary(PhyloBird_UN)) {PhyloBird_UN <- bifurcatr(PhyloBird_UN, runs=1)}
# if(!is.binary(PhyloReptile_UN)) {PhyloReptile_UN <- bifurcatr(PhyloReptile_UN, runs=1)}
# 
# # corrected trees (ultrametric at this stage)
# if(!is.binary(Phylo_Mammals_C)) {Phylo_Mammals_C <- bifurcatr(Phylo_Mammals_C, runs=1)}
# if(!is.binary(Phylo_Amphibians_C)) {Phylo_Amphibians_C <- bifurcatr(Phylo_Amphibians_C, runs=1)}
# if(!is.binary(Phylo_Birds_C)) {Phylo_Birds_C <- bifurcatr(Phylo_Birds_C, runs=1)}
# if(!is.binary(Phylo_Reptiles_C)) {Phylo_Reptiles_C <- bifurcatr(Phylo_Reptiles_C, runs=1)}
# 
# # save these phylogenies with attached species and resolved politomies
# write.tree(Phylo_Amphibians_C, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Amphibians.nwk")
# write.tree(Phylo_Reptiles_C, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Reptiles.nwk")
# write.tree(Phylo_Birds_C, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Birds.nwk")
# write.tree(Phylo_Mammals_C, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Mammals.nwk")
# 
# write.tree(PhyloAmphibian_UN, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Amphibians.nwk")
# write.tree(PhyloReptile_UN, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Reptiles.nwk")
# write.tree(PhyloBird_UN, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Birds.nwk")
# write.tree(PhyloMammal_UN, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Mammals.nwk")
# 
# # # # #

# ## for all vertebrates: 17262 species in the trait dataset are not represented in the phylogeny
# TraitsAll <- rbind(C_Amphibians[, c("Best_guess_binomial", "Body_mass_g")], C_Mammals[, c("Best_guess_binomial", "Body_mass_g")], C_Birds[, c("Best_guess_binomial", "Body_mass_g")], C_Reptiles[, c("Best_guess_binomial", "Body_mass_g")])
# Vert <- read.tree("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloVertebrates.nwk")
# VertAtt <- Attach_Species_Phylo(TraitsAll, Vert,  Where="root")
# # # # # #
# 
# ## resolve polytomies
# if(!is.binary(VertAtt)) {VertAtt <- bifurcatr(VertAtt, runs=1)}
# write.tree(VertAtt, "../../Results/1.Phylogenies/Corrected/1.Random_additions/PhyloVertebrates.nwk")
# 
# ## after randomly attaching species to the root of their genus
# PredictsCor <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")
# VertAttCor <- .Format_tiplabels(VertAtt)
# length(setdiff(unique(PredictsCor$Best_guess_binomial), VertAttCor$tip.label))
# # 863 species in PREDICTS are not represented in the verterbate phylogeny (previous number) -- 76 should be the actual number -- 102 species.
# # from 863 species to 102 species in PREDICTS not being represented in the phylogeny.





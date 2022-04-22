## Comparing original and corrected phylogenies, with random additions and dropped tips:
## plotting phylogenies with colors for orders

library(ggtree)
library(ggplot2)
library(ggpubr)
library(picante)
library(phytools)

# Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  return(Phylogeny)
}

## Load original phylogenies
Phylo_Mammals_original <- read.newick("../../Data/Phylogenies/TTOL_mammals_smoothed_interpolated_Hedges2015.nwk") %>% .Format_tiplabels() 
Phylo_Amphibians_original <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()
Phylo_Reptiles_original <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()
Phylo_Birds_original <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk") %>% .Format_tiplabels()  

## Load corrected phylogenies
Phylo_Mammals_corrected <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk") %>% .Format_tiplabels() 
Phylo_Amphibians_corrected <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk") %>% .Format_tiplabels()
Phylo_Reptiles_corrected <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk") %>% .Format_tiplabels()
Phylo_Birds_corrected <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk") %>% .Format_tiplabels()  

## how many species added?
length(Phylo_Mammals_corrected$tip.label)-length(Phylo_Mammals_original$tip.label)
length(Phylo_Birds_corrected$tip.label)-length(Phylo_Birds_original$tip.label)
length(Phylo_Amphibians_corrected$tip.label)-length(Phylo_Amphibians_original$tip.label)
length(Phylo_Reptiles_corrected$tip.label)-length(Phylo_Reptiles_original$tip.label)


## Load trait dataframes
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")


### 


Plot_phylogeny_orders <- function(Phylo, TraitDF) {
  
  TraitDF$Order <- tolower(TraitDF$Order)
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  TraitDF$Order <- firstup(TraitDF$Order)
  
  # Prune the tree to keep species that appear in the trait dataset
  Taxa=TraitDF$Order %>% as.vector()
  Names=TraitDF$Best_guess_binomial
  names(Taxa) <- Names
  Match = match.phylo.data(Phylo, Taxa)
  Phylo <- Match$phy
  Taxa <- Match$data
  
  ## get the nodes for each orders
  
  Groups <- as.data.frame(Taxa) %>%
    setNames("Order") %>%
    mutate(Species=names(Taxa))
  Groups <- split(Groups, f=Groups$Order)
  Groups <- lapply(Groups, function(x){return(x[[2]])} )
  
  tree <- groupOTU(Phylo, Groups)
  
  browser()
  
  #  get ancestral nodes or each order to plot the order name
  Func <- function(x, Phylo) {
    return(getMRCA(Phylo, x))
  }
  
  NodeOrders <- mapply(FUN=Func, x=Groups, Phylo=rep(list(Phylo), length(Groups)))
  X <- names(NodeOrders) %>%
    as.data.frame()%>%
    setNames("Order")
  for(i in 1:length(NodeOrders)) {
    if(length(NodeOrders[[i]])!=0){
      X$Node[i] <- NodeOrders[[i]]  
    } else{X$Node[i] <- NA}
  }
  
  X <- X %>% filter(!is.na(Node))
  
  p <- ggtree(tree, layout = "circular", aes(colour=group)) +
    theme(legend.position = "bottom") 
  
  for (i in 1:nrow(X)){
    p <- p + geom_cladelabel(X$Node[i], X$Order[i], )
  }
  
  return(p)
}

po = Plot_phylogeny_orders(Phylo_Mammals_original, Mammals)
pc = Plot_phylogeny_orders(Phylo_Mammals_corrected, Mammals)

pmammals <- ggarrange(po, pc, common.legend = TRUE)

po = Plot_phylogeny_orders(Phylo_Reptiles_original, Reptiles)
pc = Plot_phylogeny_orders(Phylo_Reptiles_corrected, Reptiles)

preptiles <- ggarrange(po, pc, common.legend = TRUE)

po = Plot_phylogeny_orders(Phylo_Amphibians_original, Amphibians)
pc = Plot_phylogeny_orders(Phylo_Amphibians_corrected, Amphibians)

pamphibians <- ggarrange(po, pc, common.legend = TRUE)

po = Plot_phylogeny_orders(Phylo_Birds_corrected, Birds)
pc = Plot_phylogeny_orders(Phylo_Birds_original, Birds)

pbirds <- ggarrange(po, pc, common.legend = TRUE)


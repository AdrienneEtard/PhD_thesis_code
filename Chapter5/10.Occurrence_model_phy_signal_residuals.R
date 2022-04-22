## checking phylogenetic signal in models residuals

library(dplyr)
library(phytools)

Model <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")

### models residuals
Residuals <- residuals(Model)
Names <- Model@frame$Best_guess_binomial

### phylogeny across vertebrates
Phylo <- read.newick("../Data/phylogenies/TTOL_animals_unsmoothed.nwk")

PhySignal <- function(residuals, Names, Tree, simulations, N) {
  
  .Format_tiplabels <- function (Phylogeny) {
    Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
    return(Phylogeny)
  }
  
  Residuals <- residuals
  names(Residuals) <- Names
  
  Res_lamba <- c()
  Tree <- .Format_tiplabels(Tree)
  Signal <- phytools::phylosig(Tree, Residuals, method="lambda", test = TRUE) %>% 
    unlist()
  if(!simulations){return(Signal)} 
  
  if(simulations){
    Res_lambda <- Signal$lambda
  for(i in 1:N) {
    # randomise residuals
    Residuals_sim <- sample(residuals, size = length(residuals), replace = FALSE)
    names(Residuals_sim) <- Names
    Signalsim <- phytools::phylosig(Tree, Residuals_sim, method="lambda", test = FALSE) %>% 
      unlist()
    Res_lamba <- c(Res_lamba, Signalsim$lambda)
    print(i)
  }
    return(Res_lamba)
    }
}

# without simulations 
S <- Sys.time()
lambda_residuals <- PhySignal(Residuals, Names,Tree= Phylo, simulations = FALSE, N=NULL)
E <- Sys.time()
print(E-S) ##  takes less than 10 minutes


# > lambda_residuals
# $lambda
# [1] 0.004270293
# 
# $logL
# [1] -3897.374
# 
# $logL0
# [1] -3898.198
# 
# $P
# [1] 0.1994117


# with simulations 
#lambda_residuals_sim <- PhySignal(Residuals, Names, Tree=Phylo, simulations = TRUE, N=3)

#lambda_residuals_sim[1]
# median(lambda_residuals_sim)
# quantile(lambda_residuals_sim, c(0.025, 0.975))



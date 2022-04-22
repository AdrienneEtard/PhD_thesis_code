na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

## function to run to calculate functional diversity indices
To_run <- function(List_of_communities,
                   Gowerdist,
                   Abundance,
                   SimTRUE,
                   nsim,
                   sim_pool,
                   global_pool_character,
                   Std.FRic,
                   Predicts_site_info)
{

  List <- list()

  for (i in 1:length(List_of_communities)) {

    # print(Sys.time())
    # print(paste("Starting study", i, "of", length(List_of_communities)))

    List[[i]] <- dbFD_to_apply(Distmatrix=Gowerdist,
                               community_matrix=List_of_communities[[i]],
                               Abundance_weighted=Abundance,
                               Simulations=SimTRUE,
                               nsim=nsim,
                               sim_pool=sim_pool,
                               global_pool_character=global_pool_character,
                               Std.FRic=Std.FRic,
                               Predicts_site_info=Predicts_site_info)

    print(paste("Finished study", i, "of", length(List_of_communities)))
    #print(Sys.time())
  }
  return(List)
}


## Function to transform and z-score results
Transform_zscore <- function(TraitDF, Trait, Transf) {

  if (Transf=="log10"){
    TraitDF[,paste("log10", Trait, sep="_")] <- as.numeric(log10(TraitDF[,Trait]))
    TraitDF[, paste("log10", Trait, sep="_")] <- scale(TraitDF[, paste("log10", Trait, sep="_")], center=TRUE, scale=TRUE)
    TraitDF[, paste("log10", Trait, sep="_")] <- as.numeric(TraitDF[, paste("log10", Trait, sep="_")])
  }

  if(Transf=="sqrt") {
    TraitDF[,Trait]  <- as.numeric(sqrt(TraitDF[,Trait]))
    TraitDF[, paste("sqrt", Trait, sep="_")] <- scale(TraitDF[,Trait] , center=TRUE, scale=TRUE)
    TraitDF[, paste("sqrt", Trait, sep="_")] <- as.numeric(TraitDF[, paste("sqrt", Trait, sep="_")])

  }

  return(TraitDF)
}


## Function to bind and transform and z-score the 8 imputed trait datasets
Bind_transform_traits <- function(List_imputed_datasets) {

  Candidate_traits <-  c("Body_mass_g",
                         "Lifespan_proxy",
                         "Litter_size",
                         "Habitat_breadth_IUCN",
                         "Specialisation",
                         "Diel_activity",
                         "Trophic_level")

  Candidate_traits_transformed <-  c("log10_Body_mass_g",
                                     "log10_Lifespan_proxy",
                                     "log10_Litter_size",
                                     "sqrt_Habitat_breadth_IUCN",
                                     "Specialisation",
                                     "Diel_activity",
                                     "Trophic_level")

  Taxo <- c("Class", "Order", "Family", "Genus", "Best_guess_binomial")

  List_datasets <- list()

  for (i in 1:8) {

    Imputed <- List_imputed_datasets[[i]]

    Mammals <- Imputed[["M"]][[1]]
    Mammals <- Mammals[, c(Taxo, Candidate_traits)]

    Birds <- Imputed[["B"]][[1]]
    Birds <- Birds[, c(Taxo, Candidate_traits)]

    Reptiles <- Imputed[["R"]][[1]]
    Reptiles <- Reptiles[, c(Taxo, Candidate_traits)]

    Amphibians <- Imputed[["A"]][[1]]
    Amphibians <- Amphibians[, c(Taxo, Candidate_traits)]

    Traits <- rbind(Mammals, Birds, Reptiles, Amphibians)

    Traits <- Transform_zscore(Traits, "Body_mass_g", "log10")
    Traits <- Transform_zscore(Traits, "Lifespan_proxy", "log10")
    Traits <- Transform_zscore(Traits, "Litter_size", "log10")
    Traits <- Transform_zscore(Traits, "Habitat_breadth_IUCN", "sqrt")

    Traits <- Traits[, c(Taxo, Candidate_traits_transformed)]

    List_datasets[[i]] <- Traits

  }

  return(List_datasets)
}


## Function to bind and transform and z-score the 8 imputed trait datasets - WITHIN class
Bind_transform_traits_BYCLASS <- function(List_imputed_datasets) {


  Candidate_traits <-  c("Body_mass_g",
                         "Lifespan_proxy",
                         "Litter_size",
                         "Habitat_breadth_IUCN",
                         "Specialisation",
                         "Diel_activity",
                         "Trophic_level")

  Candidate_traits_transformed <-  c("log10_Body_mass_g",
                                     "log10_Lifespan_proxy",
                                     "log10_Litter_size",
                                     "sqrt_Habitat_breadth_IUCN",
                                     "Specialisation",
                                     "Diel_activity",
                                     "Trophic_level")

  Taxo <- c("Class", "Order", "Family", "Genus", "Best_guess_binomial")


  List_datasets <- list()

  for (i in 1:8) {

    Imputed <- List_imputed_datasets[[i]]

    Mammals <- Imputed[["M"]][[1]]
    Mammals <- Mammals[, c(Taxo, Candidate_traits)]

    Mammals <- Transform_zscore(Mammals, "Body_mass_g", "log10")
    Mammals <- Transform_zscore(Mammals, "Lifespan_proxy", "log10")
    Mammals <- Transform_zscore(Mammals, "Litter_size", "log10")
    Mammals <- Transform_zscore(Mammals, "Habitat_breadth_IUCN", "sqrt")

    Birds <- Imputed[["B"]][[1]]
    Birds <- Birds[, c(Taxo, Candidate_traits)]

    Birds <- Transform_zscore(Birds, "Body_mass_g", "log10")
    Birds <- Transform_zscore(Birds, "Lifespan_proxy", "log10")
    Birds <- Transform_zscore(Birds, "Litter_size", "log10")
    Birds <- Transform_zscore(Birds, "Habitat_breadth_IUCN", "sqrt")

    Reptiles <- Imputed[["R"]][[1]]
    Reptiles <- Reptiles[, c(Taxo, Candidate_traits)]

    Reptiles <- Transform_zscore(Reptiles, "Body_mass_g", "log10")
    Reptiles <- Transform_zscore(Reptiles, "Lifespan_proxy", "log10")
    Reptiles <- Transform_zscore(Reptiles, "Litter_size", "log10")
    Reptiles <- Transform_zscore(Reptiles, "Habitat_breadth_IUCN", "sqrt")

    Amphibians <- Imputed[["A"]][[1]]
    Amphibians <- Amphibians[, c(Taxo, Candidate_traits)]

    Amphibians <- Transform_zscore(Amphibians, "Body_mass_g", "log10")
    Amphibians <- Transform_zscore(Amphibians, "Lifespan_proxy", "log10")
    Amphibians <- Transform_zscore(Amphibians, "Litter_size", "log10")
    Amphibians <- Transform_zscore(Amphibians, "Habitat_breadth_IUCN", "sqrt")

    Traits <- rbind(Mammals, Birds, Reptiles, Amphibians)

    Traits <- Traits[, c(Taxo, Candidate_traits_transformed)]

    List_datasets[[i]] <- Traits

  }

  return(List_datasets)
}

## Function to select traits based in variance inflation factor (stepwise selection)
# https://rdrr.io/github/software-analytics/Rnalytica/src/R/stepwise.vif.R
# function here just for categorical + continuous data: they are fitted against a dummy variable
stepwise.vif <- function (dataset,
                          metrics,
                          vif.threshold = 5,
                          verbose = F){

  #browser()

  dataset$dummy <- rnorm(nrow(dataset)) # generates a dummy normal distribution
  output <- metrics # those are all the candidate variables
  step.count <- 1
  output.results <- list()

  repeat {

    # VIF scores for the linear model: dummy ~ all variables; returns VIF or GVIF, depending on Degrees of freedom (if any variables has more than one then GVIF is computed rather than VIF)
    vif.scores <- vif(lm(as.formula(paste0(
      "dummy~", paste0(output,
                       collapse = "+")
    )), data = dataset))
    na.coefficients <- Reduce('|', is.nan(vif.scores)) # for scores that are NA -- if so stop the function
    if (na.coefficients) {
      stop("NA coefficient in a regression model.")
    }

    # Select VIF scores
    vif.scores <- vif.scores[,1]

    ## output.results stores VIF calculated at each step and output stores variable names that are selected
    output.results[[step.count]] <-
      sort(vif.scores, decreasing = F) # sort VIF scores

    vif.scores <- vif.scores[vif.scores >= vif.threshold]

    ## If all VIF scores are under the threshold values then keep all variables and stop the function. Else drop the variable
    # that has a VIF above the threshold and proceed iteratively
    if (length(vif.scores) == 0)
      break

    drop.var <-
      names(vif.scores[vif.scores == max(vif.scores)]) # select variable names that have VIF score more than the threshold
    if (verbose) {
      print(paste0(
        "Step ",
        step.count,
        " - Exclude ",
        drop.var,
        " (VIF = ",
        max(vif.scores),
        ")"
      ))
    }
    step.count <- step.count + 1
    output <- output[!output %in% drop.var]
  }

  names(output.results) <- paste0("Iteration ", 1:step.count)
  names(output.results)[length(output.results)] <- "Final"

  if(length(output.results)==1) {output.results <- as.data.frame(output.results)}

  return(list(Selected_vars=output, VIF=output.results))
}

## Function to compute Gower distance matrix from trait datasets, across all species
GowerDist <- function(TraitsDB, Traits) {

  rownames(TraitsDB) <- TraitsDB$Best_guess_binomial
  TraitsDB <- TraitsDB %>%
    dplyr::select(Traits[-which(Traits=="Best_guess_binomial")])

  Dist <- FD::gowdis(TraitsDB)

  # Dist <- gowdis(TraitsDB) %>%
  #   as.matrix() %>%
  #   as.data.frame()

  return(Dist)

}

########################################################################################
## Functions for functional diversity indices

# RandomSimCom a site community nsim times for nsim simulations
Replicate <- function(X, nsim){
  X <- as.data.frame(X)
  Nrow <- nsim
  X[nrow(X)+1:Nrow,] <- X[1,]
  rownames(X) <- c(1:nrow(X))
  return(X)
}

# 1. randomise replicated communities (either just presence-absence or abundance)
RandomSimCom <- function(sim_community, SR) {

  # function to randomise each row of a simulated community
  Func <- function(row) {
    row[row!=0] <- 0
    row[sample(c(1:length(row)), SR, replace = FALSE)] <- 1
    return(row)}

  # Func_ab <- function(row) {
  #   toplace <- row[row!=0]
  #   where <- which(row!=0)
  #   row[where] <- sample(toplace)
  #
  #   return(row)}

  sim_community <- apply(sim_community, 1, Func) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()

  rownames(sim_community) <- c(1:nrow(sim_community))

  # check that after randomisation, all species occur at least at one site; drop species that do not (in the case of species presence-absence data)
  sim_community[nrow(sim_community) + 1, ] <- apply(sim_community, 2, sum, na.rm=T)

  todrop <- which(sim_community[nrow(sim_community),]==0)
  if (length(todrop)!=0){
    sim_community <- sim_community[-nrow(sim_community), -todrop]
  }
  else {sim_community <- sim_community[-nrow(sim_community),]}

  # if there is only one species left in the species pool after randomising, return NA
  if (class(sim_community)!="data.frame") {return(NA)}
  else{return(as.matrix(sim_community))}
}

# 2. for dbFD
Compute_dbFD <- function(Distmatrix, comm_matrix, Abundance_weighted, Std.FRic) {

  # Subset the distance matrix for the species in the "abundance" matrix
  Dist <- usedist::dist_subset(Distmatrix, colnames(comm_matrix))

  FD <- tryCatch(expr={dbFD(x=Dist,
                            a=comm_matrix,
                            m="min",
                            corr = "cailliez",
                            w.abun = Abundance_weighted,
                            calc.FRic = TRUE,
                            stand.FRic = Std.FRic,
                            calc.FDiv=TRUE,
                            messages=FALSE)},
                 error = function(e) {dbFD(Dist,
                                           comm_matrix,
                                           m="min",
                                           corr = "cailliez",
                                           w.abun = Abundance_weighted,
                                           calc.FRic = FALSE,
                                           calc.FDiv=FALSE,
                                           messages=FALSE)})

  return(FD)
}


# 3. other functions that are called

# 95% confidence interval for the medians
Boot_med <- function(X, nsim) {
  bootmed <- apply(matrix(sample(X, rep=TRUE, nsim), nsim), 1, median)
  return(quantile(bootmed, c(.025, 0.975), na.rm = TRUE))
}

# 4. functions to extract results
extract.results <- function(list_results, to_extract) {
  lapply(list_results, function(y) y[to_extract])}

Extract_main_results <- function(List, Predicts_site_info){
  x <- extract.results(List, "results")
  results <- list()
  for (i in 1:length(x)) {
    results[[i]] <- x[[i]]$results
  }
  results <- data.table::rbindlist(results) %>% as.data.frame()
  Y <- match(results$SSBS, Predicts_site_info$SSBS)
  Predicts_site_info <- Predicts_site_info[Y,]
  results <- merge(Predicts_site_info, results, by="SSBS")
  return(results)
}

Extract_simulations <- function(List, to_extract, Sites){

  Simulations <- extract.results(List, to_extract)
  x <- list()
  for (i in 1:length(Simulations)) {
    if(to_extract=="simulations_F_distinctiveness"){
      x[[i]] <- Simulations[[i]]$simulations_F_distinctiveness
    }
    if(to_extract=="simulations_F_rarity"){
      x[[i]] <- Simulations[[i]]$simulations_F_rarity
    }
  }
  simulations <- data.table::rbindlist(x) %>% as.data.frame()

  rownames(simulations) <- Sites
  return(simulations)
}

## for dbFD
Extract <- function(Results, Predicts_site_info){

  Res <- na.omit.list(Results)
  Res <- data.table::rbindlist(Res, fill=TRUE) %>%
    as.data.frame()

  Y <- match(Res$SSBS, Predicts_site_info$SSBS)
  Predicts_site_info <- Predicts_site_info[Y,]
  Final <- merge(Predicts_site_info, Res, by="SSBS")

  return(Final)

}

# # # # MAIN FUNCTIONS  # # # #

## OTHER INDICES USING THE dBFD PACKAGE: FUNCTIONAL DIVERGENCE, DISPERSION + RAO'S QUADRATIC ENTROPY, EVENNESS, RICHNESS
# Ref: Villeger 2008, and, Legendre and Laliberte

## this function applies to one PREDICTS study (equivalent to one community matrix where columns are species and rows are sites)

dbFD_to_apply <- function(Distmatrix, community_matrix, Abundance_weighted, nsim, Simulations, sim_pool, global_pool_character, Std.FRic, Predicts_site_info) {


  # gc()
  # memory.limit(size=50000)

  # # # # # # # #

  # arguments:
  # Distmatrix: distance matrix for all the PREDICTS vertebrates
  # community_matrix: for the considered PREDICTS study. Columns are species and rows are sites
  # Abundance_weighted: whether the calculation of FD takes abundance data (in that case, no simulation). TRUE or FALSE
  # nsim number of simulations - that is, for one site, number of simulated communities
  # Simulations: TRUE or FALSE (if Abundance_weighted is TRUE, then Simulations is FALSE)
  # sim_pool: "global", "regional" or "ecoregional"
  # global_pool_character: if sim_pool=="global"
  # Std.FRic: whether FRic value should be standardised or not (TRUE/FALSE)
  # Predicts_site_info file (containing site-level information including Ecoregions)

  # # # # # # # #


  # # # # # EMPIRICAL FD

  ## Compute dbFD indices for the site -> empirical values of FD
  print("Empirical FD")
  FD <- Compute_dbFD(Dist = Distmatrix, comm_matrix = community_matrix, Abundance_weighted = Abundance_weighted, Std.FRic = Std.FRic)

  ## Abundances or species richness in each community
  if(Abundance_weighted){
    X <- apply(community_matrix, 1, function(x) {return(length(x[x!=0]))}) %>%
      as.data.frame() %>%
      setNames("SR")
  }

  else{ X <- apply(community_matrix, 1, sum) %>%
    as.data.frame() %>%
    setNames("SR")}

  ## Store dbFD results
  X$SSBS <- rownames(community_matrix)
  Sites <- names(FD$nbsp)
  if (length(FD$FRic)!=0) { X$FRic <- FD$FRic }
  if (length(FD$qual.FRic)!=0) { X$qual.FRic <- FD$qual.FRic }
  if (length(FD$FEve)!=0) { X$FEve <- FD$FEve }
  if (length(FD$FDiv)!=0) { X$FDiv <- FD$FDiv }
  if (length(FD$FDis)!=0) { X$FDis <- FD$FDis }
  if (length(FD$RaoQ)!=0) { X$RaoQ <- FD$RaoQ }

  file.remove("vert.txt")

  # # # # # SIMULATED FD


  ## Generate null expectations for each site, via simulations (nsim randomisation of communities at each site). Only if Abundance_weighted is FALSE

  # if(Simulations) {
  #
  #   print("Simulations")
  #
  #   ## prepare dataset to store future results
  #   X$median_sim_FRic <- NA
  #   X$median_sim_FEve <- NA
  #   X$median_sim_FDiv <- NA
  #   X$median_sim_FDis <- NA
  #   X$median_sim_RaoQ <- NA
  #
  #
  #   ## Species randomisation: depends on the species pool. 3 possibilities.
  #   # Here, create a matrix where ncol=number of species in the species pool, nrow=number of sites
  #   # Then, takes each row of the matrix, and for each row creates a new matrix with nsim rows
  #
  #   if(sim_pool=="global") {
  #
  #     path_to_save <- "../../Results/dbFD_indices/Sim_Global/Sim_results/"
  #
  #     # randomise species over global_pool_character
  #     new_comm <- as.data.frame(matrix(nrow=nrow(community_matrix), ncol=length(global_pool_character)))
  #     colnames(new_comm) <- global_pool_character
  #     rownames(new_comm) <- rownames(community_matrix)
  #     new_comm[,] <- 0
  #
  #     List_sim <- split(new_comm, f=rownames(new_comm))
  #     List_sim <- mapply(FUN = Replicate, List_sim, nsim=rep(list(nsim), length(List_sim)), SIMPLIFY = FALSE)
  #
  #     # Now, randomise presence/absence in simulated communities according to species richness in the site
  #     List_simulations <- mapply(FUN=RandomSimCom, List_sim, SR=as.list(X$SR), SIMPLIFY = FALSE)
  #
  #
  #     # function to verify that species are represented at least in one simulation - filter out otherwise
  #     # if SR=1, keep one other species
  #     Filter_0 <- function(x) {
  #
  #       # sum over columns
  #       y <- apply(x, 2, sum)
  #       # are there any species that are not represented in any simulated community?
  #       v <- which(y==0)
  #
  #       if(length(v)>1) {
  #
  #         # remove corresponding columns if SR>1
  #         SR <- unique(apply(x, 1, sum))
  #         if(SR>1) {return(x[,-v])}
  #
  #         # if SR=1, keep one species that is not represented anywhere (so that later, the distance matrix can be subsetted correctly)
  #         if(SR==1) {
  #           keep1 <- sample(v, 1)
  #           keep2 <- which(y!=0)
  #           return(x[,c(keep1, keep2)])
  #         }
  #
  #       }
  #
  #       else{return(x)}
  #
  #     }
  #
  #     List_simulations <- lapply(X=List_simulations, FUN=Filter_0)
  #
  #     # Calculate FD on each simulation set
  #     FD_Sim <- mapply(FUN=Compute_dbFD,
  #                      comm_matrix=List_simulations,
  #                      Dist=rep(list(Distmatrix), length(List_simulations)),
  #                      Abundance_weighted=rep(list(FALSE), length(List_simulations)),
  #                      Std.FRic=rep(list(Std.FRic), length(List_simulations)),
  #                      SIMPLIFY = FALSE)
  #
  #
  #
  #     # Store results: for each community, calculate median metrics for the simulation
  #
  #     for (i in 1:length(FD_Sim)) {
  #
  #       x1 <- median(FD_Sim[[i]]$FRic, na.rm=TRUE)
  #       x2 <- median(FD_Sim[[i]]$FEve, na.rm=TRUE)
  #       x3 <- median(FD_Sim[[i]]$FDiv, na.rm=TRUE)
  #       x4 <- median(FD_Sim[[i]]$FDis, na.rm=TRUE)
  #       x5 <- median(FD_Sim[[i]]$RaoQ, na.rm=TRUE)
  #
  #       if(length(x1)!=0){X$median_sim_FRic[i] <- x1}
  #       if(length(x2)!=0){X$median_sim_FEve[i] <- x2}
  #       if(length(x3)!=0){X$median_sim_FDiv[i] <- x3}
  #       if(length(x4)!=0){X$median_sim_FDis[i] <- x4}
  #       if(length(x5)!=0){X$median_sim_RaoQ[i] <- x5}
  #
  #     }
  #
  #   }
  #
  #   if(sim_pool=="regional") {
  #
  #     path_to_save <- "../../Results/dbFD_indices/Sim_Region/Sim_results/"
  #
  #     # create a list of nrow(community_matrix) length with nsim simulated communities inside each element, as matrices
  #     List_sim <- split(as.data.frame(community_matrix), f=rownames(community_matrix))
  #     List_sim <- mapply(FUN = Replicate, List_sim, nsim=rep(list(nsim), length(List_sim)), SIMPLIFY = FALSE)
  #
  #     # Now, randomise presence/absence in simulated communities according to species richness in the site
  #     List_simulations <- mapply(FUN=RandomSimCom, List_sim, SR=as.list(X$SR), SIMPLIFY = FALSE)
  #
  #     # For each simulated community calculate functional diversity indices (Try to compute functional diversity indices)
  #     FD_Sim <- mapply(FUN=Compute_dbFD,
  #                      comm_matrix=List_simulations,
  #                      Dist=rep(list(Distmatrix), length(List_simulations)),
  #                      Abundance_weighted=rep(list(FALSE), length(List_simulations)),
  #                      Std.FRic=rep(list(Std.FRic), length(List_simulations)),
  #                      SIMPLIFY = FALSE)
  #
  #     # Store results: for each community, calculate median metrics for the simulation
  #
  #     for (i in 1:length(FD_Sim)) {
  #
  #       x1 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FRic
  #       x2 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FEve
  #       x3 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDiv
  #       x4 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDis
  #       x5 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$RaoQ
  #
  #
  #       if(length(x1)!=0){X$median_sim_FRic[i] <- x1}
  #       if(length(x2)!=0){X$median_sim_FEve[i] <- x2}
  #       if(length(x3)!=0){X$median_sim_FDiv[i] <- x3}
  #       if(length(x4)!=0){X$median_sim_FDis[i] <- x4}
  #       if(length(x5)!=0){X$median_sim_RaoQ[i] <- x5}
  #
  #     }
  #   }
  #
  #   if(sim_pool=="ecoregional") {
  #
  #     path_to_save <- "../../Results/dbFD_indices/Sim_Ecoregion/Sim_results/"
  #
  #     # get the ecoregion names from the PREDICTS_site_info file. Load corresponding distance matrices.
  #     # Beware: some studies overlap over several ecoregions. The community matrix is going to need to be split for the simulations.
  #     # So the for loop here executes the whole simulation process over several ecoregions within one study.
  #     # ex: "HP1_2007__Wells 1" is a study that overlaps over several ecoregions.
  #
  #
  #     EcoRegions <- subset(Predicts_site_info, SSBS %in% Sites)
  #
  #     UniqueER <- unique(EcoRegions$EcoRegion_corrected)
  #
  #     for (i in UniqueER) { # do the whole process in the loop! easier than lists of lists of lists of lists....
  #
  #       # load distance matrix file
  #       path_to_file <- paste0("../../Results/Ecoregion_distance_matrices/", i, ".rds")
  #       DistER <- readRDS(path_to_file)
  #       Pool <- attr(DistER,"Labels")
  #
  #       # subset the community matrix for sites in Ecoregion i
  #       SubSites <- EcoRegions$SSBS[EcoRegions$EcoRegion_corrected==i] %>%
  #         unique() %>%
  #         as.character()
  #
  #
  #       community_matrix_sub <- community_matrix[SubSites,] # this should work on matrix class
  #
  #       if(length(SubSites)==1) {
  #         community_matrix_sub <- as.matrix(community_matrix_sub) %>% t()
  #         rownames(community_matrix_sub) <- SubSites
  #       }
  #
  #       # randomise species over the species pool that corresponds to the ecoregion
  #       new_comm <- as.data.frame(matrix(nrow=nrow(community_matrix_sub), ncol=length(Pool)))
  #       colnames(new_comm) <- Pool
  #       rownames(new_comm) <- rownames(community_matrix_sub)
  #       new_comm[,] <- 0
  #
  #       List_sim <- split(new_comm, f=rownames(new_comm))
  #       List_sim <- mapply(FUN = Replicate, List_sim, nsim=rep(list(nsim), length(List_sim)), SIMPLIFY = FALSE)
  #
  #       #  Now, randomise presence/absence in simulated communities according to species richness in the sites
  #
  #       SubSR <- subset(X, SSBS %in% SubSites)
  #
  #       SR_for_func <- as.list(SubSR$SR)
  #       names(SR_for_func) <- rownames(SubSR)
  #
  #       # order the two lists in the same way
  #       SR_for_func <- SR_for_func[names(List_sim)]
  #
  #       List_simulations <- mapply(FUN=RandomSimCom, List_sim, SR=SR_for_func, SIMPLIFY = FALSE)
  #
  #
  #       # function to verify that species are represented at least in one simulation - filter out otherwise
  #       # if SR=1, keep one other species
  #       Filter_0 <- function(x) {
  #
  #         # sum over columns
  #         y <- apply(x, 2, sum)
  #         # are there any species that are not represented in any simulated community?
  #         v <- which(y==0)
  #
  #         if(length(v)>1) {
  #
  #           # remove corresponding columns if SR>1
  #           SR <- unique(apply(x, 1, sum))
  #           if(SR>1) {return(x[,-v])}
  #
  #           # if SR=1, keep one species that is not represented anywhere (so that later, the distance matrix can be subsetted correctly)
  #           if(SR==1) {
  #             keep1 <- sample(v, 1)
  #             keep2 <- which(y!=0)
  #             return(x[,c(keep1, keep2)])
  #           }
  #
  #         }
  #
  #         else{return(x)}
  #
  #       }
  #
  #       List_simulations <- lapply(X=List_simulations, FUN=Filter_0)
  #
  #       # Calculate FD on each simulation set
  #       FD_Sim <- mapply(FUN=Compute_dbFD,
  #                        comm_matrix=List_simulations,
  #                        Dist=rep(list(DistER), length(List_simulations)),
  #                        Abundance_weighted=rep(list(FALSE), length(List_simulations)),
  #                        Std.FRic=rep(list(Std.FRic), length(List_simulations)),
  #                        SIMPLIFY = FALSE)
  #
  #       # Store results: for each community, calculate median metrics for the simulation
  #       for (i in 1:length(FD_Sim)) {
  #         x1 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FRic
  #         x2 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FEve
  #         x3 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDiv
  #         x4 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDis
  #         x5 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$RaoQ
  #
  #         N <- names(FD_Sim)[i]
  #
  #         if(length(x1)!=0){X$median_sim_FRic[which(rownames(X)==N)] <- x1}
  #         if(length(x2)!=0){X$median_sim_FEve[which(rownames(X)==N)] <- x2}
  #         if(length(x3)!=0){X$median_sim_FDiv[which(rownames(X)==N)] <- x3}
  #         if(length(x4)!=0){X$median_sim_FDis[which(rownames(X)==N)] <- x4}
  #         if(length(x5)!=0){X$median_sim_RaoQ[which(rownames(X)==N)] <- x5}
  #
  #       }
  #
  #     }
  #
  #   }
  #
  #   if (file.exists("vert.txt")) {file.remove("vert.txt")}
  #
  # } # end of "if(Simulations)"

  ## save the simulation results if possible and return the main results

  FDTrue <- exists("FD_Sim")
  if(Simulations && FDTrue) {

    Name_SS <- unique(Predicts_site_info$SS[Predicts_site_info$SSBS %in% rownames(community_matrix)])
    path_to_save <- paste0(path_to_save, Name_SS, ".rds")
    saveRDS(FD_Sim, path_to_save)
  }

  ## also save the main results, to have back up if the function crashes....

  return(X)
}

dbFD_to_apply_Ecoregional_sim <- function(community_matrix, nsim, Simulations, sim_pool, Std.FRic, Predicts_site_info) {

  gc()
  memory.limit(size=50000)

  # # # # # SIMULATED FD


  ## Generate null expectations for each site, via simulations (nsim randomisation of communities at each site). Only if Abundance_weighted is FALSE

  if(Simulations) {

    print("Simulations")

    ## prepare dataset to store future results

    X <- apply(community_matrix, 1, sum) %>%
      as.data.frame() %>%
      setNames("SR")
    X$SSBS <- rownames(community_matrix)
    X$FRic <- NA
    X$FDis <- NA
    X$RaoQ <- NA
    X$FEve <- NA
    X$FDiv <- NA
    X$QualFRic <- NA

    X$median_sim_FRic <- NA
    X$median_sim_FEve <- NA
    X$median_sim_FDiv <- NA
    X$median_sim_FDis <- NA
    X$median_sim_RaoQ <- NA


    ## Species randomisation: depends on the species pool. 3 possibilities.
    # Here, create a matrix where ncol=number of species in the species pool, nrow=number of sites
    # Then, takes each row of the matrix, and for each row creates a new matrix with nsim rows

    if(sim_pool=="ecoregional") {

      if(Std.FRic){
      path_to_save <- "../../Results/dbFD_indices/Sim_Ecoregion/Sim_results/Standardised/"
      } else{
        path_to_save <- "../../Results/dbFD_indices/Sim_Ecoregion/Sim_results/Not_standardised/"
      }

      # get the ecoregion names from the PREDICTS_site_info file. Load corresponding distance matrices.
      # Beware: some studies overlap over several ecoregions. The community matrix is going to need to be split for the simulations.
      # So the for loop here executes the whole simulation process over several ecoregions within one study.
      # ex: "HP1_2007__Wells 1" is a study that overlaps over several ecoregions.

      Sites <- rownames(community_matrix)

      EcoRegions <- subset(Predicts_site_info, SSBS %in% Sites)

      UniqueER <- unique(EcoRegions$EcoRegion_corrected)

      for (i in UniqueER) { # do the whole process in the loop! easier than lists of lists of lists of lists....

        # load distance matrix file
        path_to_file <- paste0("../../Results/Ecoregion_distance_matrices/", i, ".rds")
        DistER <- readRDS(path_to_file)
        Pool <- attr(DistER,"Labels")

        # subset the community matrix for sites in Ecoregion i
        SubSites <- EcoRegions$SSBS[EcoRegions$EcoRegion_corrected==i] %>%
          unique() %>%
          as.character()

        community_matrix_sub <- community_matrix[SubSites,] # this should work on matrix class

        if(length(SubSites)==1) {
          community_matrix_sub <- as.matrix(community_matrix_sub) %>% t()
          rownames(community_matrix_sub) <- SubSites
        }


        # randomise species over the species pool that corresponds to the ecoregion
        new_comm <- as.data.frame(matrix(nrow=nrow(community_matrix_sub), ncol=length(Pool)))
        colnames(new_comm) <- Pool
        rownames(new_comm) <- rownames(community_matrix_sub)
        new_comm[,] <- 0

        List_sim <- split(new_comm, f=rownames(new_comm))
        List_sim <- mapply(FUN = Replicate, List_sim, nsim=rep(list(nsim), length(List_sim)), SIMPLIFY = FALSE)

        #  Now, randomise presence/absence in simulated communities according to species richness in the sites
        SubSR <- subset(X, SSBS %in% SubSites)
        SR_for_func <- as.list(SubSR$SR)
        names(SR_for_func) <- rownames(SubSR)

        # order the two lists in the same way
        SR_for_func <- SR_for_func[names(List_sim)]

        # species for the actual community
        species <- apply(community_matrix, 1, FUN = function(x){return(names(which(x==1)))})

        # randomise replicated communities (either just presence-absence or abundance)
        RandomSimComm <- function(sim_community, SR, species) {

          # function to randomise each row of a simulated community
          Func <- function(row) {
            row[row!=0] <- 0
            row[sample(c(1:length(row)), SR, replace = FALSE)] <- 1
            return(row)}

          sim_community <- apply(sim_community, 1, Func) %>%
            as.data.frame() %>%
            t() %>%
            as.data.frame()

          rownames(sim_community) <- c(1:nrow(sim_community))

          # for the first row of sim community, add 1 for "species" - this represents the actual community for that site
          sim_community[1,] <- 0
          sim_community[1, species] <- 1

          # check that after randomisation, all species occur at least at one site;
          # drop species that do not (in the case of species presence-absence data)
          sim_community[nrow(sim_community) + 1, ] <- apply(sim_community, 2, sum, na.rm=T)

          todrop <- which(sim_community[nrow(sim_community),]==0)
          if (length(todrop)!=0){
            sim_community <- sim_community[-nrow(sim_community), -todrop]
          }
          else {sim_community <- sim_community[-nrow(sim_community),]}

          # if there is only one species left in the species pool after randomising, return NA
          if (class(sim_community)!="data.frame") {return(NA)}
          else{return(as.matrix(sim_community))}
        }

        List_simulations <- mapply(FUN=RandomSimComm, List_sim, SR=SR_for_func, species, SIMPLIFY = FALSE)

        # function to verify that species are represented at least in one simulation - filter out otherwise
        # if SR=1, keep one other species
        Filter_0 <- function(x) {

          # sum over columns
          y <- apply(x, 2, sum)
          # are there any species that are not represented in any simulated community?
          v <- which(y==0)

          if(length(v)>1) {

            # remove corresponding columns if SR>1
            SR <- unique(apply(x, 1, sum))
            if(SR>1) {return(x[,-v])}

            # if SR=1, keep one species that is not represented anywhere (so that later, the distance matrix can be subsetted correctly)
            if(SR==1) {
              keep1 <- sample(v, 1)
              keep2 <- which(y!=0)
              return(x[,c(keep1, keep2)])
            }

          }

          else{return(x)}

        }

        List_simulations <- lapply(X=List_simulations, FUN=Filter_0)

        # Calculate FD on each simulation set
        FD_Sim <- mapply(FUN=Compute_dbFD,
                         comm_matrix=List_simulations,
                         Dist=rep(list(DistER), length(List_simulations)),
                         Abundance_weighted=rep(list(FALSE), length(List_simulations)),
                         Std.FRic=rep(list(Std.FRic), length(List_simulations)),
                         SIMPLIFY = FALSE)

        # Store results: for each community, calculate median metrics for the simulation
        for (i in 1:length(FD_Sim)) {
          x1 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FRic
          x2 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FEve
          x3 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDiv
          x4 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDis
          x5 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$RaoQ

          N <- names(FD_Sim)[i]

          if(length(x1)!=0){X$median_sim_FRic[which(rownames(X)==N)] <- x1}
          if(length(x2)!=0){X$median_sim_FEve[which(rownames(X)==N)] <- x2}
          if(length(x3)!=0){X$median_sim_FDiv[which(rownames(X)==N)] <- x3}
          if(length(x4)!=0){X$median_sim_FDis[which(rownames(X)==N)] <- x4}
          if(length(x5)!=0){X$median_sim_RaoQ[which(rownames(X)==N)] <- x5}

          # add main results (empirical results at each site - the first element)
          xFRic <- FD_Sim[[i]]$FRic[1]
          xFDis <- FD_Sim[[i]]$FDis[1]
          xRaoQ <- FD_Sim[[i]]$RaoQ[1]
          xFEve <- FD_Sim[[i]]$FEve[1]
          xFDiv <- FD_Sim[[i]]$FDiv[1]
          xQualFRic <- FD_Sim[[i]]$qual.FRic[1]

          if(length(xFRic)!=0){X$FRic[which(rownames(X)==N)] <- xFRic}
          if(length(xFDis)!=0){X$FDis[which(rownames(X)==N)] <- xFDis}
          if(length(xRaoQ)!=0){X$RaoQ[which(rownames(X)==N)] <- xRaoQ}
          if(length(xFEve)!=0){X$FEve[which(rownames(X)==N)] <- xFEve}
          if(length(xFDiv)!=0){X$FDiv[which(rownames(X)==N)] <- xFDiv}
          if(length(xQualFRic)!=0){X$QualFRic[which(rownames(X)==N)] <- xQualFRic}

        }

      }

    }

    else{print("This function is for randomisations using ecoregional species pool.")}

    if (file.exists("vert.txt")) {file.remove("vert.txt")}

  } # end of "if(Simulations)"


  ## save the simulation results if possible and return the main results

  FDTrue <- exists("FD_Sim")
  if(Simulations && FDTrue) {

    Name_SS <- unique(Predicts_site_info$SS[Predicts_site_info$SSBS %in% rownames(community_matrix)])
    path_to_save <- paste0(path_to_save, Name_SS, ".rds")
    saveRDS(FD_Sim, path_to_save)
  }

  ## also save the main results, to have back up if the function crashes....

  return(X)
}

dbFD_to_apply_Regional <- function(Distmatrix, community_matrix, nsim, Simulations, sim_pool, Std.FRic, Predicts_site_info) {

  # gc()
  # memory.limit(size=50000)

  X <- apply(community_matrix, 1, sum) %>%
    as.data.frame() %>%
    setNames("SR")

  X$SSBS <- rownames(community_matrix)
  X$FRic <- NA
  X$qual.FRic <- NA
  X$FEve <- NA
  X$FDiv <- NA
  X$FDis <- NA
  X$RaoQ <- NA

  # # # # # SIMULATED FD

  ## Generate null expectations for each site, via simulations (nsim randomisation of communities at each site). Only if Abundance_weighted is FALSE

  if(Simulations) {

    print("Simulations")

    ## prepare dataset to store future results
    X$median_sim_FRic <- NA
    X$median_sim_FEve <- NA
    X$median_sim_FDiv <- NA
    X$median_sim_FDis <- NA
    X$median_sim_RaoQ <- NA


    ## Species randomisation: depends on the species pool. 3 possibilities.
    # Here, create a matrix where ncol=number of species in the species pool, nrow=number of sites
    # Then, takes each row of the matrix, and for each row creates a new matrix with nsim rows

    if(sim_pool=="regional") {

      if(Std.FRic){
      path_to_save <- "../../Results/dbFD_indices/Sim_Region/Sim_results/Standardised/"
      } else{
        path_to_save <- "../../Results/dbFD_indices/Sim_Region/Sim_results/Not_standardised/"
      }

      # create a list of nrow(community_matrix) length with nsim simulated communities inside each element, as matrices
      List_sim <- split(as.data.frame(community_matrix), f=rownames(community_matrix))
      List_sim <- mapply(FUN = Replicate, List_sim, nsim=rep(list(nsim), length(List_sim)), SIMPLIFY = FALSE)

      # species for the actual community
      species <- apply(community_matrix, 1, FUN = function(x){return(names(which(x==1)))})

      # study 45
      if(class(species)=="matrix") {

        print(class(species))
        List <- list()
        List[[1]] <- species[,1]
        List <- rep(List, ncol(species))
        names(List) <- colnames(species)
        species <- List
      }

      # randomise replicated communities (either just presence-absence or abundance)
      RandomSimComm <- function(sim_community, SR, species) {

        # function to randomise each row of a simulated community
        Func <- function(row) {
          row[row!=0] <- 0
          row[sample(c(1:length(row)), SR, replace = FALSE)] <- 1
          return(row)}

        sim_community <- apply(sim_community, 1, Func) %>%
          as.data.frame() %>%
          t() %>%
          as.data.frame()

        rownames(sim_community) <- c(1:nrow(sim_community))

        # for the first row of sim community, add 1 for "species" - this represents the actual community for that site
        sim_community[1,] <- 0
        sim_community[1, species] <- 1

        # check that after randomisation, all species occur at least at one site;
        # drop species that do not (in the case of species presence-absence data)
        sim_community[nrow(sim_community) + 1, ] <- apply(sim_community, 2, sum, na.rm=T)

        todrop <- which(sim_community[nrow(sim_community),]==0)
        if (length(todrop)!=0){
          sim_community <- sim_community[-nrow(sim_community), -todrop]
        }
        else {sim_community <- sim_community[-nrow(sim_community),]}

        # if there is only one species left in the species pool after randomising, return NA
        if (class(sim_community)!="data.frame") {return(NA)}
        else{return(as.matrix(sim_community))}
      }

      # Now, randomise presence/absence in simulated communities according to species richness in the site
      List_simulations <- mapply(FUN=RandomSimComm, List_sim, SR=as.list(X$SR), species, SIMPLIFY = FALSE)

      # function to verify that species are represented at least in one simulation - filter out otherwise
      # if SR=1, keep one other species
      Filter_0 <- function(x) {

        # sum over columns
        y <- apply(x, 2, sum)
        # are there any species that are not represented in any simulated community?
        v <- which(y==0)

        if(length(v)>1) {

          # remove corresponding columns if SR>1
          SR <- unique(apply(x, 1, sum))
          if(SR>1) {return(x[,-v])}

          # if SR=1, keep one species that is not represented anywhere (so that later, the distance matrix can be subsetted correctly)
          if(SR==1) {
            keep1 <- sample(v, 1)
            keep2 <- which(y!=0)
            return(x[,c(keep1, keep2)])
          }

        }

        else{return(x)}

      }

      List_simulations <- lapply(X=List_simulations, FUN=Filter_0)


      # For each simulated community calculate functional diversity indices (Try to compute functional diversity indices)
      FD_Sim <- mapply(FUN=Compute_dbFD,
                       comm_matrix=List_simulations,
                       Dist=rep(list(Distmatrix), length(List_simulations)),
                       Abundance_weighted=rep(list(FALSE), length(List_simulations)),
                       Std.FRic=rep(list(Std.FRic), length(List_simulations)),
                       SIMPLIFY = FALSE)

      # Store results: for each community, calculate median metrics for the simulation
      for (i in 1:length(FD_Sim)) {

        x1 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FRic
        x2 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FEve
        x3 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDiv
        x4 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$FDis
        x5 <- lapply(FD_Sim[[i]], median, na.rm=TRUE)$RaoQ


        if(length(x1)!=0){X$median_sim_FRic[i] <- x1}
        if(length(x2)!=0){X$median_sim_FEve[i] <- x2}
        if(length(x3)!=0){X$median_sim_FDiv[i] <- x3}
        if(length(x4)!=0){X$median_sim_FDis[i] <- x4}
        if(length(x5)!=0){X$median_sim_RaoQ[i] <- x5}

        # add main results (empirical results at each site = the first element)
        xFRic <- FD_Sim[[i]]$FRic[1]
        xFDis <- FD_Sim[[i]]$FDis[1]
        xRaoQ <- FD_Sim[[i]]$RaoQ[1]
        xFEve <- FD_Sim[[i]]$FEve[1]
        xFDiv <- FD_Sim[[i]]$FDiv[1]
        xQualFRic <- FD_Sim[[i]]$qual.FRic[1]

        if(length(xFRic)!=0){X$FRic[i] <- xFRic}
        if(length(xFDis)!=0){X$FDis[i] <- xFDis}
        if(length(xRaoQ)!=0){X$RaoQ[i] <- xRaoQ}
        if(length(xFEve)!=0){X$FEve[i] <- xFEve}
        if(length(xFDiv)!=0){X$FDiv[i] <- xFDiv}
        if(length(xQualFRic)!=0){X$qual.FRic[i] <- xQualFRic}

      }
    }

    if (file.exists("vert.txt")) {file.remove("vert.txt")}

  } # end of "if(Simulations)"

  else{print("This function is for simulations using regional (study) species pool")}
  ## save the simulation results if possible and return the main results

  FDTrue <- exists("FD_Sim")
  if(Simulations && FDTrue) {

    Name_SS <- unique(Predicts_site_info$SS[Predicts_site_info$SSBS %in% rownames(community_matrix)])
    path_to_save <- paste0(path_to_save, Name_SS, ".rds")
    saveRDS(FD_Sim, path_to_save)
  }

  return(X)
}


GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


# Function to reorder factor levels and get rid of Predominant LU NAs
Clean <- function(DF) {

  DF$Predominant_land_use <- as.character(DF$Predominant_land_use)
  DF$Predominant_land_use[DF$Predominant_land_use=="Cannot decide"] <- NA
  DF$Predominant_land_use[DF$Predominant_land_use=="Secondary vegetation (indeterminate age)"] <- NA
  DF$Predominant_land_use <- factor(DF$Predominant_land_use,
                                    levels=c( "Primary vegetation",
                                              "Mature secondary vegetation",
                                              "Intermediate secondary vegetation",
                                              "Young secondary vegetation",
                                              "Plantation forest",
                                              "Pasture",
                                              "Cropland",
                                              "Urban"))

  # DF$Use_intensity <- as.character(DF$Use_intensity)
  # DF$Use_intensity[DF$Use_intensity=="Cannot decide"] <- NA
  DF$Use_intensity <- factor(DF$Use_intensity, levels=c("Minimal use",
                                                        "Light use",
                                                        "Intense use"))


  # DF <- DF %>%
  #   dplyr::filter(!is.na(Predominant_land_use))

  DF <- subset(DF, !is.na(Predominant_land_use))

  return(DF)

}

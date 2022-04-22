## FUNCTIONS USED FOR IMPUTATIONS OF MISSING TRAIT VALUES

## Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  # Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2)) %>% unlist()
  return(Phylogeny)
}

## Functions to extract phylogenetic eigenvectors from the phylogenies and return as a dataframe - for species that match the trait dataset
Extract_eigenvectors <- function(TraitDF, Phylo, N) {

  ## arguments:
  # N <- number of eigenvectors to extract: 10 is enough to maximise imputation accuracy (Penone et al. 2014)
  # Phylo: phylogeny considered
  # TraitDF: trait dataset with species names in "Best_guess_binomial"

  ## Prune species from tree that do not intersect with trait dataset
  row.names(TraitDF) <- TraitDF$Best_guess_binomial
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
  Eigenvectors$Best_guess_binomial <- Prune_Taxa$data$Best_guess_binomial
  Eigenvectors <- Eigenvectors[order(Eigenvectors$Best_guess_binomial), c(11, 1:10)]

  return(Eigenvectors)

}


## Function to add EV to trait dataset
Add_eigenvectors <- function (TraitDF, EV) {

  ColN <- vector()
  for (i in 1:10) {ColN <- c(ColN, paste("EV", i, sep="_")) }

  Species <- as.character(EV$Best_guess_binomial)
  x <- which(TraitDF$Best_guess_binomial %in% Species)
  TraitDF[, ColN] <- NA
  TraitDF[x, ColN] <- EV[, c(2:ncol(EV))]

  return(TraitDF)

}


## Function to impute missing trait values using missForest

## ARG:
# Formatted phylogeny (class phylo)
# Formatted trait dataset (class dataframe)
# Taxinfo: genus, family, order
# Traits to impute (character string): continuous and categorical + eigenvectors to retain
# ErrorTrue: error format returned by missForest (variablewise: TRUE or FALSE)
# std: whether reprocessed diet breadth or not ("raw" diet breadth)

## RETURNS:
# Imputed trait values (dataframe)
# OOB errors from the imputations

Imputations_missForest <- function (TraitDF, Taxinfo, Traits_cont, Traits_cat, EV, ErrorTrue, DietTRUE, std) {

  ## order by alphabetical order
  TraitDF <-  TraitDF[order(TraitDF$Best_guess_binomial),]

  ## Select traits of interest, to impute, and phylogenetic eigenvectors, and taxinfo
  To_impute <- TraitDF[, colnames(TraitDF) %in% c(Taxinfo, Traits_cont, Traits_cat, EV)]
  rownames(To_impute) <- TraitDF$Best_guess_binomial

  ## Set habitat and diet variables as binary (TRUE for 1 and FALSE for 0)

  ## Set class as numeric for all continuous traits
  To_impute[, Traits_cont] <- apply(To_impute[, Traits_cont], 2, as.numeric)

  ## Set habitat and diet variables as factors. Otherwise, when they are numeric or logical, outputs are numeric

  if(DietTRUE){

    Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")
    for (i in Diet) {
      To_impute[,i] <- factor(To_impute[,i], levels=c(0,1))

    }
  }


  Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
               "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
               "Artificial","Introduced.vegetation","Other.Unknown")

  for (i in Habitat) {
    To_impute[,i] <- factor(To_impute[,i], levels=c(0,1))
    }

  ## Impute missing values
  print("Imputing missing values.")

  if (ErrorTrue) {
    R.Imputed <- missForest(To_impute, variablewise = TRUE)
    Imputed <- R.Imputed$ximp
    Errors <- R.Imputed$OOBerror %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()
    colnames(Errors) <- paste(colnames(Imputed), names(R.Imputed$OOBerror))
    row.names(Errors) <- 1
  }
  else {
    R.Imputed <- missForest(To_impute, variablewise = FALSE)
    Errors <- R.Imputed$OOBerror
    Imputed <- R.Imputed$ximp
  }

  ## Select traits and variables of interest after imputations
  Imputed$Best_guess_binomial <- rownames(Imputed)
  Imputed <- Imputed[order(Imputed$Best_guess_binomial),
                     c(Taxinfo, "Best_guess_binomial",
                       Traits_cont, Traits_cat, "EV_1")]

  ## Add a column for phylogenetic information (yes or no) to know if it was available during the imputations
  Imputed$Phylo_info <- TraitDF$EV_1
  Imputed$Phylo_info[!is.na(Imputed$Phylo_info)] <- "YES"
  Imputed$Phylo_info[is.na(Imputed$Phylo_info)] <- "NO"

  ## Add taxonomic information
  Imputed$Order <- TraitDF$Order
  Imputed$Family <- TraitDF$Family
  Imputed$Genus <- TraitDF$Genus

  ## Reprocess Primary diet and diet breadth (sqrt + normalise) if diet is included
  ## Then reorder columns.

  if(DietTRUE){

    Func <- function(X) {
      names(X) <- Diet
      ToPaste <- names(X)[which(X==1)]
      return(paste(ToPaste, collapse = "|"))
    }

    # Reprocess primary diet, for comparison with imputed values
    Imputed$Primary_diet_reprocessed <- apply(Imputed[,Diet], 1, Func)
    # Imputed$Primary_diet_reprocessed[Imputed$Primary_diet_reprocessed==""] <- "OM"

    # Reprocess diet breadth, for comparison with imputed values
    Imputed[, Diet] <- apply(Imputed[, Diet], 2, as.numeric)
    Imputed$Diet_breadth_reprocessed <- apply(Imputed[, Diet], 1, sum, na.rm=T)

    if(std) {
      Imputed$Diet_breadth_reprocessed <- sqrt(Imputed$Diet_breadth_reprocessed)
      Imputed$Diet_breadth_reprocessed <- scale(Imputed$Diet_breadth_reprocessed, center = TRUE, scale = TRUE)
      colnames(Imputed)[colnames(Imputed)=="Diet_breadth_reprocessed"] <- "sqrt_Diet_breadth_reprocessed"
    }

    ## Rearrange columns
    Imputed <- Imputed[, unique(c("Order", "Family", "Genus", "Best_guess_binomial",
                                  colnames(Imputed[grepl("Body_mass_g", colnames(Imputed))]),
                                  colnames(Imputed[grepl("Longevity_d", colnames(Imputed))]),
                                  colnames(Imputed[grepl("Litter_size", colnames(Imputed))]),
                                  colnames(Imputed[grepl("Range_size_m2", colnames(Imputed))]),
                                  "Diel_activity",
                                  "Trophic_level",
                                  colnames(Imputed[grepl("Diet_breadth", colnames(Imputed))]),
                                  colnames(Imputed[grepl("Diet_breadth_reprocessed", colnames(Imputed))]),
                                  "Primary_diet",
                                  colnames(Imputed[grepl("Primary_diet_reprocessed", colnames(Imputed))]),
                                  Diet,
                                  "Specialisation",
                                  colnames(Imputed[grepl("Habitat_breadth_IUCN", colnames(Imputed))]),
                                  Habitat,
                                  "Phylo_info"))]

  }

  else {
    Imputed <- Imputed[, c("Order", "Family", "Genus", "Best_guess_binomial", Traits_cont, Traits_cat,"Phylo_info")]
  }

  rownames(Imputed) <- c(1:nrow(Imputed))
  ToReturn <- list(Imputed.Dataset=Imputed, Imputation.errors=Errors)
  ToReturn <- list(ToReturn)
  return(ToReturn)
}


## Function to apply in parallel (runs the above function Imputations_missForest)
To_apply_parallel_imputations <- function (List_of_arguments) {

  Imputations_results <- pbmapply (FUN=Imputations_missForest,
                                   TraitDF=List_of_arguments[["TraitDF"]],
                                   Taxinfo=List_of_arguments[["Taxinfo"]],
                                   Traits_cont=List_of_arguments[["Traits_cont"]],
                                   Traits_cat=List_of_arguments[["Traits_cat"]],
                                   EV=List_of_arguments[["EV"]],
                                   ErrorTrue=List_of_arguments[["ErrorTrue"]],
                                   DietTRUE=List_of_arguments[["DietTRUE"]],
                                   std=List_of_arguments[["std"]])

  return (Imputations_results)
}

## Functions for resolving taxonomic matching (typos, synonyny, and all other operations)

# ------------------------------------------------------------------------------- DATA CLEANING

## Format tip labels in phylogenies
.Format_tiplabels <- function (Phylogeny) {
  
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2))
  
  return(Phylogeny)
  
}


## Species original names: check that all have two words (binomial), and return original species list

SubClass <- function(Predictsdata, Class) {
  
  Data <- unique(Predictsdata$Best_guess_binomial[Predicts$Class==Class]) %>%
    as.character() %>% as.data.frame() %>% setNames(., "Original")
  
  Data$Original <- Data$Original %>% as.character()
  
  Data$NChar <- sapply(strsplit(Data$Original, " "), length)
  
  if (length(unique(Data$NChar))==1) {
    print("All names are binomial. Returning original names.")
    Data <- Data %>% select(Original)
    return(Data)
  } else {
    print("Not all names are binomial. Manual verification needed. Returning original names and NChar.")
    return(Data)
  }
  
}


## Species original names for phylogenies and trait datasets

SpNames <- function(SpList) {
  
  Data <- unique(SpList) %>%
    as.character() %>% as.data.frame() %>% setNames(., "Original")
  
  Data$Original <- Data$Original %>% as.character()
  
  Data$NChar <- sapply(strsplit(Data$Original, " "), length)
  
  if (length(unique(Data$NChar))==1) {
    print("All names are binomial. Returning original names.")
    Data <- Data %>% select(Original)
    return(Data)
  } else {
    print("Not all names are binomial. Manual verification needed. Returning original names and NChar.")
    return(Data)
  }
  
}


## Checking for typos.

CheckTypos <- function(SpDF) {
  
  print("checking for typos...")
  
  SpDF$CorrectedTypos <- gnr_resolve(SpDF$Original, best_match_only = TRUE, canonical = TRUE)$matched_name2
  
  # When only one character in corrected names - replace by the original name
  SpDF$NCharCorrected <- sapply(strsplit(SpDF$CorrectedTypos, " "), length)
  SpDF$CorrectedTypos[SpDF$NCharCorrected!=2] <- SpDF$Original[SpDF$NCharCorrected!=2]
  
  SpDF <- SpDF %>% dplyr::select(Original, CorrectedTypos)
  
  # Order by alphabetical order
  SpDF <- SpDF[order(SpDF$Original),]
  rownames(SpDF) <- c(1:nrow(SpDF))
  
  SpDF$IsCorrected[SpDF$Original!=SpDF$CorrectedTypos] <- TRUE
  SpDF$IsCorrected[SpDF$Original==SpDF$CorrectedTypos] <- FALSE
  
  return(SpDF)
}


ForBirds <- function(DF) {
  
  for (i in 1:nrow(DF)) {
    
    X <- gnr_resolve(DF$Original[i], best_match_only = TRUE, canonical = TRUE)$matched_name2
    
    if (is.null(X)){ DF$CorrectedTypos[i] <- "NULL"
    } else {DF$CorrectedTypos[i] <- X}
    
    print(i)
  }
  
  DF$IsCorrected[DF$Original!=DF$CorrectedTypos] <- TRUE
  DF$IsCorrected[DF$Original==DF$CorrectedTypos] <- FALSE
  
  return(DF)
}


# ------------------------------------------------------------------------------- EXTRACT SYNONYMS

## Function to extract synonyms # WITH RED LIST

Redlist_synonyms <- function(SpeciesDF) {
  
 FunctionToApply_RL = function(x) {  
   
# browser()
    
    print(paste("Processing", x, ": species", which(SpeciesDF$CorrectedTypos==x),"on", nrow(SpeciesDF)))  
    
    # Several attempts in case the connection is not established straightaway
   
   FuzzyMatch <- FALSE
   
    X <- NULL
    
    while(is.null(X)) {  
      
      try(X <- rl_synonyms(name=x))
      
    }
    
    ## 1. No synonym found OR no Red List record
    
    if (X$count==0) {
      
      # Check if there is no synonym OR if the species does not appear in the RedList database, by querying the IUCN website
      # several attempts in case the connection is not established straightaway
      
      Search <- NULL
      attempt <- 1
      while( is.null(Search) && attempt <= 3 ) {
        
        try(Search <- rl_search(name = x)$result)
        attempt <- attempt + 1
        
      } 
      
      id <- Search$taxonid
      
        # 1.1. If id is null: no record in the Red List -- 
        
      if (is.null(id)) {
        
          InRedList <- FALSE
          IsAccepted <- NA
          IsSynonym <- NA
          Accepted <- NA
          Synonyms <- NA
        
          Family <- NA
          Order <- NA
       }
      
      
        # 1.2. If id is not null: the species does not have recorded synonyms in the Red List
        
      else {
        
          InRedList <- TRUE
          IsAccepted <- TRUE
          IsSynonym <- FALSE
          Accepted <- X$name
          Synonyms <- "No recorded synonym in the Red List."
        
          if (is.null(Search$family)) { Family <- NA }
          else {Family <- unique(Search$family)}
          
          if (is.null(Search$order)) { Order <- NA }
          else {Order <- unique(Search$order) }
      }
      
      
    } 
    
    
    ## 2. When synonyms are found
    
    else {
      
      Name <- X$name
      
      ## Typo in the red list
      # if (Name=="aotus lemurinus") {Name=="Aotus lemurinus"}
      
      InRedList <- TRUE
      
        # Find out status of the name (accepted or synonym)
      
        # 2.1. If the name is the accepted scientific name
      
        # with "any" because in some cases, names differ sightly (with Acomys cahirinus for example)
        if (any(grepl(tolower(Name), tolower(X$result$accepted_name)))) {
        
          IsAccepted <- TRUE
          IsSynonym <- FALSE
        
          Accepted <- Name
          
          Synonyms <- X$result$synonym[X$result$synonym != Name]
          Synonyms <- paste(Synonyms, collapse = " | ")
          
          }
      
      
        # 2.2 If the name is a synonym that is not accepted  
        else { 
        
          IsAccepted <- FALSE
          IsSynonym <- TRUE
        
          # Look for the accepted name corresponding to that synonym -- a few special cases treated separatly
          
          ## Special cases -- when there are several accepted names returned: Find the best match between the name and the accepted names
          
          ToLook <- X$result$accepted_name %>% unique
          
          if (length(ToLook)!=1) {
            
            ToLook <- word(ToLook, 1,2) %>% unique # when genus, species + subspecies 
            
            if (length(ToLook)!=1) { # when different species: find the best match by computing the distances between the character string
              
              MinStringDist <- stringdist(Name, ToLook) %>% which.min()
              ToLook <- X$result$accepted_name[MinStringDist]
              
              FuzzyMatch <- TRUE
              
            }
            
            }
          
          # if (Name=="Bunopithecus hoolock") {
          #   ToLook <- X$result$accepted_name[2] # Hoolock hoolock
          # } 
          # 
          # if (Name=="Callicebus cupreus") {
          #   ToLook <- X$result$accepted_name[4] # Plecturocebus cupreus
          # } 
          # 
          # if (Name=="Callicebus torquatus") {
          #   ToLook <- X$result$accepted_name[4] # Cheracebus torquatus
          # }
          
          ## end of special cases
          
          else {
            
            ToLook <- X$result$accepted_name %>% unique
          }
        
          Y <- try(rl_synonyms(name=ToLook))
        
          if (Y$count==0) {
          
              Accepted <- Y$`name`
              Synonyms <- X$`name`
          
          } 
          
          else {
            
              # if there are more than one accepted names appearing because of subspecies: choose ToLook as the accepted name
            
              if (length(unique(Y$result$accepted_name))!=1) {
                Accepted <- ToLook
                Synonyms <- Y$result$synonym[Y$result$synonym != ToLook]
                Synonyms <- paste(Synonyms, collapse = " | ")
              }
              
              else {
                Accepted <- unique(Y$result$accepted_name)
                Synonyms <- paste(Y$result$synonym, collapse = " | ")
              }
            
          }
        
        }
      
      # Family and order information for the accepted name
      
      Search2 <- NULL
      attempt <- 1

      while( is.null(Search2) && attempt <= 3 ) {

        try(Search2 <- rl_search(name = Accepted)$result)
        attempt <- attempt + 1
      }

      if (is.null(Search2$family)) { Family <- NA }
      else {Family <- unique(Search2$family)}

      if (is.null(Search2$order)) { Order <- NA }
      else {Order <- unique(Search2$order)}
        
      } 
    
    return(list(InRedList=InRedList, 
                IsAccepted=IsAccepted, 
                IsSynonym=IsSynonym, 
                Accepted=Accepted, 
                Synonyms=Synonyms, 
                Order=Order, 
                Family=Family, 
                FuzzyMatch=FuzzyMatch))
    
  } # End of function to apply
  

  List <- pblapply(SpeciesDF$CorrectedTypos, FunctionToApply_RL)
 
  
  SpeciesDF$InRedList <- unlist(List)[attr(unlist(List),"names")=="InRedList"]
  SpeciesDF$IsAccepted <- unlist(List)[attr(unlist(List),"names")=="IsAccepted"]
  SpeciesDF$IsSynonym <- unlist(List)[attr(unlist(List),"names")=="IsSynonym"]
  SpeciesDF$Accepted <- unlist(List)[attr(unlist(List),"names")=="Accepted"]
  SpeciesDF$Synonyms <- unlist(List)[attr(unlist(List),"names")=="Synonyms"]
  SpeciesDF$Order <- unlist(List)[attr(unlist(List),"names")=="Order"]
  SpeciesDF$Family <- unlist(List)[attr(unlist(List),"names")=="Family"]
  SpeciesDF$FuzzyMatch <- unlist(List)[attr(unlist(List),"names")=="FuzzyMatch"]
  
  .firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  if (!is.na(SpeciesDF$Accepted)) {
    SpeciesDF$Accepted <- .firstup(SpeciesDF$Accepted)
  }
  return(SpeciesDF)
  
  }




## Function to run: runs the Redlist_synonym function on subsets (splits) of the original data

RunSyn <- function (SpeciesDF) {
  
  # browser()
  
  # Split the dataframe into smaller DF
  n <- 100
  nr <- nrow(SpeciesDF)
  DF <- split(SpeciesDF, rep(1:ceiling(nr/n), each=n, length.out=nr))
  
  # Run the function
  DF_results <- pblapply(DF, Redlist_synonyms)
  DF_results <- rbind.fill(DF_results)
  
  return(DF_results)
}

## Modification of this function (RunSyn) to save at each iteration the split datasets with synonyms: useful when 
## imputations times are very long and when errors frequently occur (type error in curl...)

RunSynSave <- function (SpeciesDF, Path) { # path is either "Birds_splits/Birds_" or "Reptiles_splits/Reptiles_"
  
  # browser()
  
  # Split the dataframe into smaller DFs
  n <- 100
  nr <- nrow(SpeciesDF)
  DF <- split(SpeciesDF, rep(1:ceiling(nr/n), each=n, length.out=nr))
  
  # Run the function
  # DF_results <- pblapply(DF, Redlist_synonyms)
  # DF_results <- rbind.fill(DF_results)
  
  for (i in 1:length(DF)) {
      
    Results <- Redlist_synonyms(DF[[i]])
    write.csv(Results,
              paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/", Path, i, ".csv", sep=""), row.names=FALSE)
    
    print(paste("Saved dataset", i, "of", length(DF)))
    
  }
}




## Function to extract synonyms # WITH ITIS

## Function to look for ITIS synonyms for remaining taxa: species that have not been found in the Red List db

FunctionToApply_ITIS <- function(x) {
  
  #browser()
  
  # 1. Get the Taxonomic Serial Number
  TSN <- tryCatch(expr={get_tsn(x)[1]},
                  error=function(e) {NA})
  
  
  # 1.1. TSN is NA: taxon not in ITIS
  if (is.na(TSN)) { 
    
    InITIS <- FALSE
    IsAccepted <- NA
    IsSynonym <- NA
    Accepted <- NA
    Synonyms <- NA
    # Family <- NA
    # Order <- NA
    
  } 
  
  # 1.2. TSN not NA: taxon is in ITIS
  
  else {
    
    InITIS <- TRUE
    
    # A. Look for the accepted name
    Search <- tryCatch(expr={itis_acceptname(TSN)},
                       error=function(e) {NA})
    
    ## A.1. The name is already the accepted name: returns NA in acceptedname -- x is the accepted name
    if (is.na(Search$acceptedname)) {
      
      # A.1.1. Look for the synonyms of this accepted name
      Syn <- synonyms(x, db="itis")
      Syn <- Syn[[1]]$syn_name
      
      # A.1.2. Store results
      IsAccepted <- TRUE
      IsSynonym <- FALSE
      Accepted <- x
      Synonyms <- paste(Syn, collapse = " | ")
      
    }
    
    ##  A.2. The name is a synonym: get the accepted name and look for the synonyms of that accepted name      
    else {
      
      # A.2.1. Look for the synonyms of this accepted name
      Accepted <- Search$acceptedname[1]
      Syn <- synonyms(Accepted, db="itis")
      Syn <- Syn[[1]]$syn_name
      
      # A.2.2. Store results
      IsAccepted <- FALSE
      IsSynonym <- TRUE
      Synonyms <- paste(Syn, collapse = " | ")
    }
    
  }
  
  return(list(InITIS=InITIS, 
              IsAccepted=IsAccepted, 
              IsSynonym=IsSynonym, 
              Accepted=Accepted, 
              Synonyms=Synonyms ))
  
}


## The whole function (TO RUN):

ITIS_synonyms <- function(SpeciesDF) {
  
 #browser()
  
  SpeciesDF$CorrectedTypos <- as.character(SpeciesDF$CorrectedTypos)
  
  List <- pblapply(SpeciesDF$CorrectedTypos, FunctionToApply_ITIS)
  
  SpeciesDF$InITIS <- unlist(List)[attr(unlist(List),"names")=="InITIS"]
  SpeciesDF$IsAccepted <- unlist(List)[attr(unlist(List),"names")=="IsAccepted"]
  SpeciesDF$IsSynonym <- unlist(List)[attr(unlist(List),"names")=="IsSynonym"]
  SpeciesDF$Accepted <- unlist(List)[attr(unlist(List),"names")=="Accepted"]
  SpeciesDF$Synonyms <- unlist(List)[attr(unlist(List),"names")=="Synonyms"]
  
  return(SpeciesDF)
  
}


## Function to split the synonym dataset and run the ITIS_synonym function on the subset where RedList==FALSE

Complement_ITIS <- function (SynDF, Split, N) {
  
  SplitRedList <- split(SynDF , f=SynDF$InRedList )
  RedListFALSE <- SplitRedList[[1]]
  RedListTRUE <- SplitRedList[[2]]
  
  print(paste("There are", nrow(RedListFALSE), "species that were not found on the Red List."))
  
  # Retrieve synonyms for RedListFALSE: divide the dataset in smaller portions to save at each iterations
  
  if (Split){
    
    # Split the dataframe into smaller DFs (for birds) => 16 smaller datasets
    n <- 100
    nr <- nrow(RedListFALSE)
    DF <- split(RedListFALSE, rep(1:ceiling(nr/n), each=n, length.out=nr))
    
    for (i in N:length(DF)) {
      
      Results <- ITIS_synonyms(DF[[i]])
      write.csv(Results,
                paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Birds_splits/BirdITIS_", i, ".csv", sep=""), row.names=FALSE)
      
      print(paste("Saved dataset", i, "of", length(DF)))
      
    }
    
    files <- list.files(path = "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Birds_splits/", pattern = ".csv")
    files <- paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL+ITIS/Birds_splits/", files, sep="")
    RedListFALSE <- do.call(rbind,lapply(files,read.csv))
  }
  
  else (RedListFALSE <- ITIS_synonyms(RedListFALSE))
 
  
  RedListTRUE$InITIS <- "Assessed_RL"
  SynDF <- rbind(RedListTRUE, RedListFALSE)
  SynDF <- SynDF[order(SynDF$Original),]
  
  return(SynDF)
  
}

## Transform to characters

ToChar <- function(DF) {
  
  DF$Original <- as.character(DF$Original)
  DF$CorrectedTypos <- as.character(DF$CorrectedTypos)
  DF$Accepted <- as.character(DF$Accepted)
  DF$Synonyms <- as.character(DF$Synonyms)
  DF$Order <- as.character(DF$Order)
  DF$Family <- as.character(DF$Family)
  
  return(DF)
}


## Function to check that all species names have 2 words (binomial)

Check_binomial <- function(DF) {
  
 N_Words_Original <- sapply(strsplit(DF$Original, " "), length)
 N_Words_Corrected <- sapply(strsplit(DF$CorrectedTypos, " "), length) 
 N_Words_Accepted <- sapply(strsplit(DF$Accepted, " "), length) 
 
 if(length(unique(N_Words_Original))==1) { print("All original names have 2 words") }
 else { print("Not all original names have two words") }
 
 if(length(unique(N_Words_Corrected))==1) { print("All corrected names have 2 words") }
 else { print("Not all corrected names have two words") }
 
 if(length(unique(N_Words_Accepted))==1) { print("All accepted names have 2 words") }
 else { print("Not all accepted names have two words") }
 
}


## Function for species entered under the form genus + sp. / + cf.

ToGenus <- function(DF) {
  
  DF$Genus_level <- FALSE
  
  DF$Accepted[grep("cf.", DF$Original, fixed=TRUE)] <- ""
  DF$Accepted[grep("sp.", DF$Original, fixed=TRUE)] <- ""
  DF$Synonyms[grep("cf.", DF$Original, fixed=TRUE)] <- ""
  DF$Synonyms[grep("sp.", DF$Original, fixed=TRUE)] <- ""
  DF$Genus_level[grep("cf.", DF$Original, fixed=TRUE)] <- TRUE
  DF$Genus_level[grep("sp.", DF$Original, fixed=TRUE)] <- TRUE
  
  return(DF)
}

# ## Function to subset species that have not been assessed by either RL or ITIS, to be processed by GBIF backbone
# 
# NotAssessed_RL_ITIS <- function(DF, Taxon) {
#   
#   DF <- subset(DF, InITIS==FALSE)
#   
#   DF <- DF %>% select(CorrectedTypos) %>%
#     setNames("scientificName")
#   
#   write.csv(DF, paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Not_assessed_RL_ITIS/", 
#                       Taxon, ".csv", sep=""), row.names=FALSE)
#   
# }
# 
# ## Function to subset and save the names that need manual checks
# 
# Require_manual_check <- function(DF, Taxon) {
#   
#   `%nin%` <- Negate(`%in%`)
#   DF <- subset(DF, status %nin% "ACCEPTED")
#   print(paste("There are", nrow(DF), "species that require manual checks for", Taxon))
#   write.csv(DF, paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Not_assessed_RL_ITIS/GBIF_outputs/Manual_checks/", 
#                       Taxon, ".csv", sep=""), row.names=FALSE)
# 
# }
# 
# 
# 
# ## Function that pass original names as accepted names for species not found in the Red list or in ITIS 
# 
# Remaining_species <- function(DF) {
#   
#   DF <- split(DF, f=DF$InITIS)
#   
#   NotITIS <- DF[[2]]
#   Assessed <- rbind(DF[[1]], DF[[3]])
#   
#   NotITIS$Accepted <- NotITIS$CorrectedTypos
#   
#   DF <- rbind(NotITIS, Assessed)
#   DF <- DF[order(DF$Original),]
#   
#   return(DF)
# }

## Function to add a genus column

GenusCol <- function(DF) {
  
  DF$Genus[DF$Genus_level==TRUE] <- word(DF$CorrectedTypos[DF$Genus_level==TRUE], 1)
  DF$Genus[DF$Genus_level==FALSE] <- word(DF$Accepted[DF$Genus_level==FALSE], 1)
  
  return(DF)
}

## Function to get a subset of the data for which order and family information is missing (prepare for GBIF)
# 
# MissingOrder <- function(DF, Taxon) {
#   
#   DF <- subset(DF, is.na(Order))
#   DF <- DF %>% select(CorrectedTypos, Accepted) %>%
#     setNames(c("CorrectedTypos", "scientificName"))
#   
#   write.csv(DF, paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Missing_order_family/", 
#                   Taxon, ".csv", sep=""), row.names=FALSE)
#   
# }


# ## Function to add order and family information 
# 
# AddOrderFam <- function(SynDF, Info) {
#   
#   SynDF_missing <- SynDF %>% filter(is.na(Order))
#   SynDF_complete <- SynDF %>% filter(!is.na(Order))
#   
#   Names <- SynDF_missing$CorrectedTypos %>% unique
#   
#   for (i in 1:length(Names)) {
#     
#     SynDF_missing$Order[SynDF_missing$CorrectedTypos==Names[i]] <- Info$order[Info$verbatimScientificName==Names[i]]
#     SynDF_missing$Family[SynDF_missing$CorrectedTypos==Names[i]] <- Info$family[Info$verbatimScientificName==Names[i]]
#     
#     print(i)
#     
#   }
#   
#   SynDF <- rbind(SynDF_missing, SynDF_complete)
#   SynDF <- SynDF[order(SynDF$Original),]
# }


## Fuction to prepare a Genus-Family-Order information dataset and return a list of genus to check against gbif backbone

Order_family_info <- function(DF) {
  
  Family_order <- DF %>% select(Genus, Family, Order, CorrectedTypos) %>%
  group_by(Genus) %>%
  as.data.frame()
Family_order <- unique(Family_order[c("Genus", "Family", "Order")])

# Genus that appear only once: do not filter for NAs
Counts <- Family_order %>% group_by(Genus) %>%
  summarise(Count=count(Genus)) %>%
  as.data.frame() 

ToFilter <- Counts$Count$x[Counts$Count$freq==2]
NotToFilter <- Counts$Count$x[Counts$Count$freq==1]

Filtered <- subset(Family_order, Genus %in% ToFilter) %>% filter(!is.na(Family))

Not_Filtered <- subset(Family_order, Genus %nin% ToFilter)

GBIFInfo <- subset(Not_Filtered, is.na(Family))
TaxInfo <- rbind(Filtered, Not_Filtered)
GBIFInfo <- GBIFInfo[order(GBIFInfo$Genus),] %>%
  select(Genus)

FunctionX <- function(x) {
  Y <- DF$CorrectedTypos[DF$Genus==x][1]
  return(Y)
}

# TaxInfo$Name <- sapply(TaxInfo$Genus, FunctionX)

ForGBIF=GBIFInfo
colnames(ForGBIF) <- "scientificName"
ForGBIF$scientificName <- sapply(ForGBIF$scientificName, FunctionX)

print(length(unique(TaxInfo$Genus)))
print(length(unique(Family_order$Genus)))
print(length(unique(GBIFInfo$Genus)))

return(list(TaxInfo=TaxInfo, ForGBIF=ForGBIF, GBIFInfo=GBIFInfo))
  
}


## Complement with GBIF data

Complement_with_GBIF <- function(DF, GBIF) {
  
  FunctionY <- function(y) {
    return(GBIF$family[GBIF$genus==y] %>% toupper())
  }
  
  FunctionZ <- function(z) {
    return(GBIF$order[GBIF$genus==z] %>% toupper())
  }
  
  DF$Family <- apply(DF, 1, FunctionY)
  DF$Order <- apply(DF, 1, FunctionZ)
  
  return(DF)  
}


## Function to add additional taxonomic information (family and order) to synonym datasets where lacking

AddFamilyOrder <- function(Syn, ToAdd) {

  Missing <- c()
  Lacking <- subset(Syn, is.na(Order)) 
  
  for (i in 1:nrow(Lacking)){
    
    if(any(grepl(Lacking$Genus[i], ToAdd$Genus))) {
      
      Lacking$Family[i] <- ToAdd$Family[ToAdd$Genus==Lacking$Genus[i]] %>% unique
      Lacking$Order[i] <- ToAdd$Order[ToAdd$Genus==Lacking$Genus[i]]%>% unique
    }
    
    else(Missing <- c(Missing, Lacking$CorrectedTypos[i]))
    
  }
  
  SynToreturn <- rbind(Lacking, subset(Syn, !is.na(Syn$Order)))
  SynToreturn <- SynToreturn[order(SynToreturn$Original),]
  
  return(list(Syn=SynToreturn, Missing=Missing))
}


## Second round of GBIF addition
GBIF2 <- function(Syn, GBIF) {
  
  # browser()
  
  for (i in 1:nrow(GBIF)) {
    
    Syn$Family[Syn$CorrectedTypos==GBIF$verbatimScientificName[i]] <- GBIF$family[i] %>% toupper
    Syn$Order[Syn$CorrectedTypos==GBIF$verbatimScientificName[i]] <- GBIF$order[i] %>% toupper
  }
  
  return(Syn)
  
}



## Functions to replace names in original trait datasets / PREDICTS / phylogenies

Replace_by_accepted_name <- function(Syn, TargetDF, colnamespecies) {
  
  No_match <- c()
  
  Syn$Original <- as.character(Syn$Original)
  Syn$Accepted <- as.character(Syn$Accepted)
  Syn$CorrectedTypos <- as.character(Syn$CorrectedTypos)
  
  if (class(TargetDF)=="phylo") {
    
    Names <- TargetDF$tip.label %>% unlist
    
    for (i in 1:length(Names)) {
      
      if (!grepl("cf.|sp.", Names[i])) {
        
        if(any(grepl(Names[i], Syn$Original))) {
          TargetDF$tip.label[TargetDF$tip.label==Names[i]] <- Syn$Accepted[Syn$Original==Names[i]]
        }
        
      }
      
    # print(i)
      
    }
    
    print(paste("Delta in species number:", length(unique(Names))- length(unique(TargetDF$tip.label))))
    
  }
  
  else {
    
    TargetDF[ , colnamespecies] <- as.character(as.factor(TargetDF[ , colnamespecies]))
  
    Names <- TargetDF[ , colnamespecies] %>% unique
  
      for (i in 1:length(Names)) {
        
       if(Names[i]=="Desmodus draculae") {
         
         TargetDF$Best_guess_binomial[TargetDF[ , colnamespecies]==Names[i]] <-"Desmodus draculae"; next()}
        
        if(Names[i]=="Pholidoscelis maynardi") {
          
          TargetDF$Best_guess_binomial[TargetDF[ , colnamespecies]==Names[i]] <-"Pholidoscelis maynardi"; next()}
       
        if(any(grepl(Names[i], Syn$Original, fixed=TRUE))){
          TargetDF$Best_guess_binomial[TargetDF[ , colnamespecies]==Names[i]] <- Syn$Accepted[Syn$Original==Names[i]]
        }
        
        else {
          
          if(any(grepl(Names[i], Syn$CorrectedTypos))){
            TargetDF$Best_guess_binomial[TargetDF[ , colnamespecies]==Names[i]] <- Syn$Accepted[Syn$CorrectedTypos==Names[i]]
          }
          
          else {
            TargetDF$Best_guess_binomial[TargetDF[ , colnamespecies]==Names[i]] <- Names[i]
            No_match <- c(Names[i], No_match)
          }
        }
       
        print(i)
      }
    
    print(No_match)
  print(paste("Delta in species number:", length(unique(Names))- length(unique(TargetDF$Best_guess_binomial))))
  
  }
  
  if(length(No_match!=0)) {return(list(No_match=No_match, Data=TargetDF))}
  else{(TargetDF)}
  
}

## Replace by accepted names for PREDICTS

ReplacePredicts <- function (Predictssub, syndf) {
  
  syndf$Accepted <- as.character(syndf$Accepted)
  Names <- unique(Predictssub$Best_guess_binomial)
  
  for (i in 1:length(Names)) {
    
    Predictssub$Best_guess_binomial[Predictssub[ , "Best_guess_binomial"]==Names[i]] <- syndf$Accepted[syndf$Original==Names[i]]
    # print(paste(i, "on", length(Names)))
  }
  
  return(Predictssub)
}




# ## Function to check which Predicts species are not represented in the phylogenies
# 
# Not_in_phylogeny <- function(Predictsdata, VClass, Phylogeny) {
#   
#   SubsetClass <- subset(Predictsdata, Class==VClass)
#   SubsetClass <- subset(SubsetClass, Best_guess_binomial!="")
#   Species <- unique(SubsetClass$Best_guess_binomial)
#   
#   print(paste("There are", length(Species), "species known by binomial name in Predicts for", VClass))
# 
#   Diff <- setdiff(Species, Phylogeny$tip.label)
# 
#   print(paste("Among these, there are", length(Diff), "species that do not have a match in the phylogeny"))
#   
#   return(Diff)  
# }
# 








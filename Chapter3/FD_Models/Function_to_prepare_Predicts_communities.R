## FUNCTIONS TO PREPARE PREDICTS COMMUNITY MATRICES, IN A LIST

# output = a list of PREDICTS studies; in each element, a list of PREDICTS community matrix at each site

## Function to apply to a list of PREDICTS studies, to prepare community matrices -- with an abundance option
# to return only community matrices with abundances (excluding sites for which only presence/absence was measured)
Create_assemblage_matrix <- function(Predicts_SS_data, Abundance) {
  
  # from predicts site level data, prepare a community matrix: 
  # each row is a site, each column is a species.
  # species presence/absence as a binary variable
  # if abundance option, then abundance rather than presence/absence
  
  Community <- Predicts_SS_data %>%
    dplyr::select("Best_guess_binomial", "SSBS", "Effort_corrected_measurement")
  
  ## If abundance only is wanted (arg Abundance=TRUE), then return NA if metrics is presence/absence only
  if(Abundance) {
    
    if(unique(Predicts_SS_data$Diversity_metric_type)!="Abundance") {
      return(NA); 
      next
      }
    
  } 
    
    ## verify that the species pool has more than one species; otherwise not possible to compute functional richness / possibly other functional diversity indices
    n <- length(unique(Community$Best_guess_binomial))
    
    if(n==1){
      
      print("Species has only one species")
      return(NA)
      }
    
    else{
      
      ## verify that the same set of species is sampled across all the sites of the study; otherwise not possible to compute functional richness
      check <- Community %>%
        dplyr::group_by(SSBS, Best_guess_binomial) %>%
        dplyr::summarise(Count=n()) %>%
        dplyr::group_by(Best_guess_binomial) %>%
        dplyr::summarise(Count=n())
      
      if (length(unique(check$Count))!=1) {
        print("different set of species sampled across the site")
        return(NA)
        }
      
      ## if all checks are good, then create assemblage matrix for the study (site by species)
      ## ie: for studies that have more than one species in the species pool and for which the same sets of species have been sampled, prepare assemblage matrix
      else {
        
        # order within sites by alphabetical order
        Community <- Community %>% 
          group_by(SSBS) %>%
          arrange(Best_guess_binomial, .by_group = TRUE)
        
        # put 1 for presence of the species (rather than the metric provided, which can be abundance)
        if(!Abundance){
          Community <- Community %>%
            mutate(Effort_corrected_measurement=ifelse(Effort_corrected_measurement!=0, 1, 0))
        } 
        
        # select unique sites and sum occurrence over similar sites -- in case there are replicates (should not happen)
        if(length(unique(Community$SSBS))!=length(Community$SSBS)) {
          print("replicated sites!")
          }
        
        Community <- unique(Community[c("Best_guess_binomial", "SSBS", "Effort_corrected_measurement")])
        
        Community <- Community %>% 
          dplyr::group_by(Best_guess_binomial, SSBS) %>%
          dplyr::summarise(Sum=sum(Effort_corrected_measurement))
        
        ## Now prepare community matrix: rows are sites, columns are species
        Comm.matrix <- matrix(ncol=length(unique(Community$Best_guess_binomial)),
                              nrow=length(unique(Community$SSBS))) %>%
          as.data.frame()
        
        row.names(Comm.matrix) <- unique(Community$SSBS)
        colnames(Comm.matrix) <- unique(Community$Best_guess_binomial)
        
        ## fill in the matrix with abundance / presence-absence data
        for (i in 1:ncol(Comm.matrix)){
          Comm.matrix[,i] <- Community$Sum[Community$Best_guess_binomial==colnames(Comm.matrix)[i]]
        }
        
        ## Additional checks
        
        # check that all communities (sites) have non-zero total abundance/presence
        # (at least one species occurs); 
        # select sites that have non-zero total abundance if it happens that total abundance is 0
        Comm.matrix$TotAb <- apply(Comm.matrix, 1, sum, na.rm=TRUE)
        if(any(Comm.matrix$TotAb)==0){
          print("0 total abundance!")
        }
        
        Comm.matrix <- subset(Comm.matrix, TotAb != 0)
        Comm.matrix <- Comm.matrix %>%
          dplyr::select(-TotAb)
        
        # Check that species occur at least at one site;
        # select species that occur at least at one site
        Comm.matrix[nrow(Comm.matrix) + 1, ] <- apply(Comm.matrix, 2, sum, na.rm=TRUE)
        Todrop <- which(Comm.matrix[nrow(Comm.matrix),]==0)
        
        if (length(Todrop)!=0){
          print("")
          Comm.matrix <- Comm.matrix[-nrow(Comm.matrix), -Todrop]
        }  
        
        else {Comm.matrix <- Comm.matrix[-nrow(Comm.matrix),]}
        
        # if after removing species that did not occur across all sites, there is only
        # one species left in the community, return NA (when there is only one species left in the species pool)
        if (class(Comm.matrix)!="data.frame") 
        {
          return(NA)
        }
        
        else{
          return(as.matrix(Comm.matrix))
        }
        
      }
      
    }
}

## Function to run to generate the list of community matrices (applies the above function)
# Outputs in a list: community dataframe (by studies: each element contains a matrix of communities
# for each site of a study)
Generate_communities <- function(Predictsdata, AbundanceTrue) {
  
  # each study as a list element
  SS_List <- split(Predictsdata, f=Predictsdata$SS, drop=TRUE)
  
  # apply function Community matrix
  x <- lapply(FUN =Community_matrix, X = SS_List, Abundance=AbundanceTrue)
  
  return(x)
}



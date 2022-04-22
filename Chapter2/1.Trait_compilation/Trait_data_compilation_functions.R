## Functions to prepare the trait datasets and the data extraction:
## (necessary runs of these functions, sourced to the main script, are at the end atm)


#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
# Normalise trait datasets -------------------------------------------------------------------------

## This function changes column names, and format trait datasets
## so that column names are unique and identifiable across datasets

.Normalise_TDB <- function(Kissling, Elton_MD, Pantheria, Pacifici, Myhrvold, 
                           Amphibio, Cooper, Senior, Bickford, 
                           Butchart_BM, Butchart_GL, Elton_BD, # Sekercoiglu_Diet, 
                           Scharf, Vidan, Stark, Schwarz, Novosolov, Novosolov_2, Slavenko, Meiri, Meiri2018, Feldman,
                           #Range_mammal, Range_amphibian, Range_bird,
                           Format_range) {
  
  
  # Diet variables
  Diet <- c("FR", "NE", "SE", "PL", "IN", "VE")
  
  # Kissling
  colnames(Kissling)[22] <- "Trophic_level"
  Kissling$Trophic_level %<>% as.character()

  # Elton (Kissling replaced by Elton trait for mammalian diet)
  colnames(Elton_MD)[14] <- "Body_mass_g"
  Elton_MD$Primary_diet %<>% as.character()
  Elton_MD <- Elton_MD %>%
    mutate(Diel_activity=ifelse(Activity.Nocturnal==1,"Nocturnal","Other"))
  Mammal_traits.Elton <- c("Diel_activity","Body_mass_g", Diet, "Primary_diet", "Diet_breadth")
  
  # Pantheria
  Pantheria[Pantheria==-999] <- NA
  Mammal_traits.Pantheria <- c("Diel_activity","Body_mass_g", "Forearm_length_mm", "Head_length_mm", "Habitat_breadth",
                               "Home_range_group_km2", "Home_range_ind_km2", "Litter_size", "Max_longevity_d", "Maturity_d",
                               "Terrestriality")
  colnames(Pantheria)[c(1, 2, 3, 6, 7, 8, 9, 17:19, 21, 23, 28, 31)] <- c("Order", "Family", "Genus", Mammal_traits.Pantheria)
  Pantheria$Max_longevity_d <- Pantheria$Max_longevity_d*30.44
  # Pantheria <- Pantheria %>%
  #   mutate(Diel_activity=ifelse(Diel_activity==1, "Nocturnal", Diel_activity)) %>%
  #   mutate(Diel_activity=ifelse(Diel_activity==2, "Else", Diel_activity)) %>%
  #   mutate(Diel_activity=ifelse(Diel_activity==3, "Diurnal", Diel_activity))
  
  Pantheria <- Pantheria %>%
    mutate(Diel_activity=ifelse(Diel_activity==1, "Nocturnal", "Other"))
  
  # Pacifici
  Mammal_traits.Pacifici <- c("Body_mass_g", "Max_longevity_d","AFR_d", "Generation_length_d")
  colnames(Pacifici)[c(6, 8, 11, 14)] <- Mammal_traits.Pacifici
  Pacifici[Pacifici=="no information"] <- NA
  Pacifici$Max_longevity_d  %<>% as.character %>% as.numeric()  
  Pacifici$Generation_length_d  %<>% as.character %>% as.numeric()  
  Pacifici$AFR_d  %<>% as.character %>% as.numeric()  
  
  # Myhrvold
  Myhrvold[Myhrvold==-999] <- NA
  Myhrvold$female_maturity_d[Myhrvold$female_maturity_d<0] <- NA
  Myhrvold$Maturity_d <- apply(Myhrvold[, c("female_maturity_d", "male_maturity_d")], 1, mean, na.rm=T)
  Myhrvold$Maturity_d[is.na(Myhrvold$Maturity_d)] <- Myhrvold$no_sex_maturity_d[is.na(Myhrvold$Maturity_d)]
  Traits.Myhrvold <- c("Litter_size", "Body_mass_g", "Max_longevity_d", "Longevity_d", "Adult_svl_cm")
  colnames(Myhrvold)[c(9, 11, 12, 20, 29)] <- Traits.Myhrvold
  colnames(Myhrvold)[c(2:4)] <- c("Order", "Family", "Genus")
  Myhrvold$Longevity_d <- Myhrvold$Longevity_d * 365.25
  Myhrvold$Max_longevity_d <- Myhrvold$Max_longevity_d * 365.25
  Traits.Myhrvold <- c(Traits.Myhrvold, "Maturity_d")
  
  # Amphibio
  
  # Transforming activity times into one variable
  # Amphibio <- Amphibio %>% 
  #   mutate(Diel_activity="Else") %>%
  #   mutate(Diel_activity=ifelse((Diu==1 & is.na(Noc) & is.na(Crepu)), "Diurnal", Diel_activity)) %>%
  #   mutate(Diel_activity=ifelse((Noc==1 & is.na(Diu) & is.na(Crepu)), "Nocturnal", Diel_activity)) 
  Amphibio <- Amphibio %>% 
    mutate(Diel_activity="Other") %>%
    mutate(Diel_activity=ifelse((Noc==1 & is.na(Diu) & is.na(Crepu)), "Nocturnal", Diel_activity))
  Amphibio$Primary_diet %<>% as.character()
  
  colnames(Amphibio)[c(26, 29)] <- c("Body_length_mm", "Max_longevity_d")
  Amphibio %<>% select(-id)
  Amphibio$Maturity_d <- apply(Amphibio[, c("Age_at_maturity_min_y", "Age_at_maturity_min_y")], 1, mean, na.rm=T) * 365.25
  Amphibio$Max_longevity_d <- Amphibio$Max_longevity_d * 365.25
  Amphibio$Litter_size <- apply(Amphibio[, c("Litter_size_min_n", "Litter_size_max_n")], 1, mean, na.rm=T) 
  Amphibian_traits.Amphibio <- c("Body_mass_g", "Litter_size", "Body_length_mm", "Max_longevity_d", "Maturity_d",
                                 "Fos", "Ter", "Aqu", "Arb", 
                                 Diet, "Diet_breadth", "Trophic_level", "Primary_diet", "Diel_activity")
  
  # Cooper
  Amphibian_traits.Cooper <- c("Habitat_breadth", "Svl_length_mm", "Litter_size")
  colnames(Cooper)[c( 4, 7, 9)] <- Amphibian_traits.Cooper
  # Cooper$Range_size_m2 <- Cooper$Range_size_m2 * 1000000
  Cooper$Genus <- word(Cooper$Best_guess_binomial, 1)
  
  # Senior
  Senior %<>% filter(Rank %in% "Species") %>% select(-Rank)
  Amphibian_traits.Senior <- c("Svl_length_mm")
  colnames(Senior)[11] <- Amphibian_traits.Senior
  
  # Bickford
  Amphibian_traits.Bickford <- c("Body_length_mm", "Terrestriality")
  colnames(Bickford)[c(9, 10)] <- Amphibian_traits.Bickford
  # Bickford$Range_size_m2 <- as.numeric(Bickford$Range_size_m2) * 1000000
  
  # Scharf reptile
  Reptile_traits.Scharf <- c("Max_longevity_d", "Body_mass_g", "Diel_activity", "Trophic_level", "Litter_size", "Maturity_d")
  colnames(Scharf)[c(5,6,8,9,11,14)] <- Reptile_traits.Scharf
  Scharf$Trophic_level %<>% as.character()
  Scharf <- Scharf %>% mutate(Trophic_level=ifelse(Trophic_level=="Carnivorous", "Carnivore", 
                                                   ifelse(Trophic_level=="Herbivorous", "Herbivore", 
                                                          ifelse(Trophic_level=="Omnivorous", "Omnivore", NA))))
  Scharf$Max_longevity_d <- Scharf$Max_longevity_d * 365.25
  Scharf$Body_mass_g <- exp(Scharf$Body_mass_g)
  Scharf$Maturity_d <- Scharf$Maturity_d * 30.44 
  Scharf$Diel_activity <- as.character(Scharf$Diel_activity)
  Scharf$Diel_activity[Scharf$Diel_activity=="Diurnal"] <- "Other"
  Scharf$Diel_activity[Scharf$Diel_activity=="Cathemeral"] <- "Other"
  
  # # Meiri reptile -- the non published data
  # Meiri %<>% filter(Rank=="Species")
  # colnames(Meiri)[8] <- "Body_mass_g"
  # Meiri$Body_mass_g <- 10^Meiri$Body_mass_g
  
  # Vidan reptiles (diel activity)
  colnames(Vidan)[2] <- "Diel_activity"
  Vidan$Diel_activity <- as.character(Vidan$Diel_activity)
  Vidan$Diel_activity[Vidan$Diel_activity=="Diurnal"] <- "Other"
  Vidan$Diel_activity[Vidan$Diel_activity=="Cathemeral"] <- "Other"
  
  # Stark reptiles: diel activity; max longevity; mean BM; 
  Traits.Stark <- c("Max_longevity_d", "Diel_activity", "Body_mass_g", "Litter_size")
  colnames(Stark)[c(4, 9, 16, 17)] <- Traits.Stark
  Stark$Max_longevity_d <- Stark$Max_longevity_d * 365.25
  Stark$Diel_activity <- as.character(Stark$Diel_activity)
  Stark$Diel_activity [Stark$Diel_activity=="Diurnal"] <- "Other"
  Stark$Diel_activity [Stark$Diel_activity=="Crepuscular "] <- "Other"
  Stark$Diel_activity [Stark$Diel_activity=="Crepuscular"] <- "Other"
  Stark$Diel_activity [Stark$Diel_activity=="Cathemeral"] <- "Other"
  Stark$Diel_activity [Stark$Diel_activity==""] <- NA
  
  # Schwarz: clutch size
  colnames(Schwarz)[5] <- "Litter_size"
  Schwarz$Litter_size <- 10^Schwarz$Litter_size
  
  # Novosolov (1): body mass, range size, trophic level
  Traits.Novosolov <- c("Body_mass_g", "Trophic_level")
  colnames(Novosolov)[c(4,16)] <- Traits.Novosolov
 # Novosolov$Range_size_m2 <- Novosolov$Range_size_m2 * 1000000
  Novosolov$Trophic_level %<>% as.character()
  Novosolov$Trophic_level[Novosolov$Trophic_level=="Carnivorous"] <- "Carnivore"
  Novosolov$Trophic_level[Novosolov$Trophic_level=="Omnivorous"] <- "Omnivore"
  Novosolov$Trophic_level[Novosolov$Trophic_level=="Herbivorous"] <- "Herbivore"

  # Novosolov (2): clutch size
  colnames(Novosolov_2)[6] <- "Litter_size"
  
  # Slavenko: BM (it is log10)
  colnames(Slavenko)[5] <- "Body_mass_g"
  Slavenko$Body_mass_g <- 10^(Slavenko$Body_mass_g) # back transforming
  
  # Meiri 2015 Evolutionary Biology
  Traits.Meiri <- c("Diel_activity", "Trophic_level")
  colnames(Meiri)[c(13,15)] <- Traits.Meiri
  Meiri$Trophic_level %<>% as.character()
  Meiri$Trophic_level[Meiri$Trophic_level=="Carnivorous"] <- "Carnivore"
  Meiri$Trophic_level[Meiri$Trophic_level=="Omnivorous"] <- "Omnivore"
  Meiri$Trophic_level[Meiri$Trophic_level=="Herbivorous"] <- "Herbivore"
  Meiri$Diel_activity <- as.character(Meiri$Diel_activity)
  Meiri$Diel_activity [Meiri$Diel_activity=="Cathemeral"] <- "Other"
  Meiri$Diel_activity [Meiri$Diel_activity=="Diurnal"] <- "Other"
  
  # Meiri 2018 GEB - added after reviewers suggestion
  # activity time, extinct/extant/EW, TL, clutch size, age at 1st breeding, (max SVL)
  Traits.Meiri.GEB <- c("Diel_activity", "Trophic_level", "Litter_size", "Maturity_d")
  
  Meiri2018$Litter_size <- apply(Meiri2018[, c("smallest.clutch", "largest.clutch")], 1, mean, na.rm=T)
  Meiri2018$Diel_activity <- ifelse(Meiri2018$Activity.time=="Nocturnal", "Nocturnal", 
                                    ifelse(is.na(Meiri2018$Activity.time), NA, "Other"))
  Meiri2018$Trophic_level <- ifelse(Meiri2018$diet=="Carnivorous", "Carnivore", 
                                    ifelse(Meiri2018$diet=="Omnivorous", "Omnivore",
                                           ifelse(Meiri2018$diet=="Herbivorous", "Herbivore", NA)))
  
  Meiri2018$Maturity_d <- apply(Meiri2018[, c("youngest.age.at.first.breeding..months.",
                                              "oldest.age.at.first.breeding..months.")], 1, mean, na.rm=T)* 30.44 # converting into days

  
  # Feldman 2016 - added after reviewers suggestion
  Traits.Feldman <- "Body_mass_g"
  colnames(Feldman)[12] <- Traits.Feldman
  
  # Butchart Avian body mass
  colnames(Butchart_BM)[29] <- "Body_mass_g"
  
  # Butchart Avian generation length
  colnames(Butchart_GL)[2] <- "Generation_length_d"
  Butchart_GL$Generation_length_d <- Butchart_GL$Generation_length_d * 365.25
  
 # # Sekercioglu - we chose to use diet data from Elton traits instead
 # Traits.Sekercioglu <- c(Diet, "Diet_breadth", "Trophic_level")
  
  # Elton bird diet
  colnames(Elton_BD)[26] <- "Body_mass_g"
  Elton_BD$Primary_diet %<>% as.character()
  Elton_BD <- Elton_BD %>%
    mutate(Diel_activity=ifelse(Nocturnal==1,"Nocturnal","Other"))
  Bird_traits.Elton <- c("Body_mass_g", Diet, "Primary_diet", "Diet_breadth", "Trophic_level", "Diel_activity")
  
  
  # # AnAge -- we dedided not to include AnAge
  # AnAge[AnAge==-9999] <- NA
  # AnAge$Maturity_d <- apply(AnAge[, c("Male.maturity..days.", "Female.maturity..days.")], 1, mean, na.rm=T)
  # AnAge$Species <-  paste(AnAge$Genus, AnAge$Species, sep= " ")
  # AnAge.Traits <- c("Litter_size", "Max_longevity_d", "Body_mass_g")
  # colnames(AnAge)[c(8, 14, 21, 29)] <- c("Binomial_name", AnAge.Traits)
  # AnAge <- AnAge %>% filter(Class %in% c("Mammalia", "Reptilia", "Amphibia", "Aves"))
  # AnAge$Max_longevity_d <- AnAge$Max_longevity_d * 365.25
  # AnAge.Traits <- c(AnAge.Traits, "Maturity_d")
  
  # # Format the Range size datasets
  # if (Format_range==TRUE) {
  #     Format_range <- function (Range){
  #     colnames(Range) <- c("Species", "Range_size_m2", "Best_guess_binomial")
  #     return(Range)
  #     }
  #     
  #     Range_mammal <- Format_range(Range_mammal)
  #     Range_amphibian <- Format_range(Range_amphibian)
  #     Range_bird <- Format_range(Range_bird)
  # }
  

  # OUTPUTS
  return(list(Traits.Myhrvold=Traits.Myhrvold,
              
              Kissling=Kissling, 
              Elton_MD=Elton_MD,
              Pantheria=Pantheria, 
              Pacifici=Pacifici, 
              Myhrvold=Myhrvold,
              
              Amphibio=Amphibio,
              Cooper=Cooper, 
              Senior=Senior,
              Bickford=Bickford,
              
              Scharf=Scharf,
              Meiri=Meiri,
              Vidan=Vidan,
              Schwarz=Schwarz,
              Stark=Stark,
              Traits.Stark=Traits.Stark,
              Novosolov=Novosolov, 
              Traits.Novosolov=Traits.Novosolov,
              Novosolov_2=Novosolov_2,
              Slavenko=Slavenko,
              Meiri=Meiri,
              Traits.Meiri=Traits.Meiri,
              Feldman=Feldman,
              Traits.Feldman=Traits.Feldman,
              Meiri2018=Meiri2018, 
              Traits.Meiri.GEB=Traits.Meiri.GEB,
              
              Butchart_BM=Butchart_BM,
              Butchart_GL=Butchart_GL,
              Elton_BD=Elton_BD,
              Bird_traits.Elton=Bird_traits.Elton,
              #Sekercioglu_Diet=Sekercioglu_Diet,
              #Traits.Sekercioglu=Traits.Sekercioglu,
              
              # Range_mammal=Range_mammal,
              # Range_amphibian=Range_amphibian,
              # Range_bird=Range_bird,
              
              Mammal_traits.Pantheria=Mammal_traits.Pantheria,
              Mammal_traits.Pacifici=Mammal_traits.Pacifici,
              Mammal_traits.Elton=Mammal_traits.Elton,
              
              Amphibian_traits.Amphibio=Amphibian_traits.Amphibio, 
              Amphibian_traits.Cooper=Amphibian_traits.Cooper,
              Amphibian_traits.Senior=Amphibian_traits.Senior,
              Amphibian_traits.Bickford=Amphibian_traits.Bickford,
              
              Reptile_traits.Scharf=Reptile_traits.Scharf))
              #AnAge=AnAge))
  
}


## Function to remove duplicate values for continuous traits, then take mean or median

.Remove_duplicates_average <- function(df, mean_or_median) {
  
  # Create dataframe with unique "Best_guess_binomial" to store results
  Out <- unique(df$Best_guess_binomial) %>% 
    as.data.frame() %>% 
    setNames(., "Best_guess_binomial")
  
  # for each continuous traits, remove duplicated values, and then take mean or median
  Traits <- colnames(df)[colnames(df)!="Best_guess_binomial"]
  
  for (i in Traits) {
    
    x <- unique(df[c("Best_guess_binomial", i)])
    x <- x %>%  
      dplyr::group_by(Best_guess_binomial) 
    x <- x[!is.na(x[, i]),] %>% 
      ungroup %>% 
      as.data.frame()

    # x only contains species for which at least there is one estimate for the trait
    dup_unique <- x %>% 
      dplyr::group_by(Best_guess_binomial) %>% 
      dplyr::summarise(Count=n())
    dup_unique <- unique(dup_unique$Count) 
    dup_unique <- dup_unique[order(dup_unique)]
    print(dup_unique)
    
    # duplicates 
    print(paste(i, nrow(df[!is.na(df[,i]),]) - nrow(x)))
    
    if(mean_or_median=="mean"){
       xf <- setDT(x)[, lapply(.SD, mean, na.rm=TRUE), by = Best_guess_binomial]  %>%
         as.data.frame() 
    }
    if(mean_or_median=="median"){
      xf <- setDT(x)[, lapply(.SD, median, na.rm=TRUE), by = Best_guess_binomial]  %>%
        as.data.frame() 
    }
    # join with NAs where the species is not present in xf
    Out <- dplyr::left_join(Out, xf, by="Best_guess_binomial")
  }
  
    return(Out)
}


#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 


Plot_df1_versus_df2 <- function(df1, df2, Traits, Xlab, Ylab) {
  
  df1 <- df1 %>% dplyr::select(Traits)
  df2 <- df2 %>% dplyr::select(Traits)
  
  par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2.5,2,1), oma=c(1,1,2,1))
  par(mfrow=c(2,2))
  
  for (i in 1:ncol(df1)) {
    plot(log(df1[, i])~log(df2[, i]), 
         main=paste(colnames(df1)[i], "(log)"), 
         xlab=Xlab, ylab=Ylab,
         pch=19)
    abline(a=0, b=1, col="red", lwd=2)
  }
  box("outer")
} 



#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

# Reduce IUCN Habitat traits and diet redundancy -----------------------------------

.MutHabIUCN <- function(CatTrait) {
  
  CatTrait$Specialisation %<>% as.character()
  
  CatTrait  %<>% group_by(Best_guess_binomial)%>%
    # mutate(Habitat_breadth_IUCN=
    #          ifelse(length(unique(Habitat_breadth_IUCN))==1, unique(Habitat_breadth_IUCN), unique(Habitat_breadth_IUCN[!is.na(Habitat_breadth_IUCN)]))) %>%
    dplyr::mutate(Specialisation=
             ifelse(length(unique(Specialisation))==1, unique(Specialisation), unique(Specialisation[!is.na(Specialisation)]))) %>%
    dplyr::mutate(Forest=
             ifelse(length(unique(Forest))==1, unique(Forest), unique(Forest[!is.na(Forest)]))) %>%
    dplyr::mutate(Savanna=
             ifelse(length(unique(Savanna))==1, unique(Savanna), unique(Savanna[!is.na(Savanna)]))) %>%
    dplyr::mutate(Shrubland=
             ifelse(length(unique(Shrubland))==1, unique(Shrubland), unique(Shrubland[!is.na(Shrubland)]))) %>%
    dplyr::mutate(Grassland=
             ifelse(length(unique(Grassland))==1, unique(Grassland), unique(Grassland[!is.na(Grassland)]))) %>%
    dplyr::mutate(Wetland=
             ifelse(length(unique(Wetland))==1, unique(Wetland), unique(Wetland[!is.na(Wetland)]))) %>%
    dplyr::mutate(Rocky.areas=
             ifelse(length(unique(Rocky.areas))==1, unique(Rocky.areas), unique(Rocky.areas[!is.na(Rocky.areas)]))) %>%
    dplyr::mutate(Caves.and.subterranean=
             ifelse(length(unique(Caves.and.subterranean))==1, unique(Caves.and.subterranean), unique(Caves.and.subterranean[!is.na(Caves.and.subterranean)]))) %>%
    dplyr::mutate(Desert=
             ifelse(length(unique(Desert))==1, unique(Desert), unique(Desert[!is.na(Desert)]))) %>%
    dplyr::mutate(Marine=
             ifelse(length(unique(Marine))==1, unique(Marine), unique(Marine[!is.na(Marine)]))) %>%
    dplyr::mutate(Marine.intertidal.or.coastal.supratidal=
             ifelse(length(unique(Marine.intertidal.or.coastal.supratidal))==1, 
                    unique(Marine.intertidal.or.coastal.supratidal), 
                    unique(Marine.intertidal.or.coastal.supratidal[!is.na(Marine.intertidal.or.coastal.supratidal)]))) %>%
    dplyr::mutate(Artificial=
             ifelse(length(unique(Artificial))==1, unique(Artificial), unique(Artificial[!is.na(Artificial)]))) %>%
    dplyr::mutate(Introduced.vegetation=
             ifelse(length(unique(Introduced.vegetation))==1, unique(Introduced.vegetation), unique(Introduced.vegetation[!is.na(Introduced.vegetation)]))) %>%
    dplyr::mutate(Other.Unknown=
             ifelse(length(unique(Other.Unknown))==1, unique(Other.Unknown), unique(Other.Unknown[!is.na(Other.Unknown)])))
  
  return(CatTrait)
}



.MutDiet <- function(CatTrait) {
  
  CatTrait$Trophic_level %<>% as.character()
  
  CatTrait  %<>% group_by(Best_guess_binomial)%>%
    dplyr::mutate(Trophic_level=
             ifelse(length(unique(Trophic_level))==1, unique(Trophic_level), unique(Trophic_level[!is.na(Trophic_level)]))) %>%
    dplyr::mutate(Primary_diet=
             ifelse(length(unique(Primary_diet))==1, unique(Primary_diet), unique(Primary_diet[!is.na(Primary_diet)]))) %>%
    # mutate(Diet_breadth=
    #          ifelse(length(unique(Diet_breadth))==1, unique(Diet_breadth), unique(Diet_breadth[!is.na(Diet_breadth)]))) %>%
    dplyr::mutate(IN=
             ifelse(length(unique(IN))==1, unique(IN), unique(IN[!is.na(IN)]))) %>%
    dplyr::mutate(VE=
             ifelse(length(unique(VE))==1, unique(VE), unique(VE[!is.na(VE)]))) %>%
    dplyr::mutate(FR=
             ifelse(length(unique(FR))==1, unique(FR), unique(FR[!is.na(FR)]))) %>%
    dplyr::mutate(NE=
             ifelse(length(unique(NE))==1, unique(NE), unique(NE[!is.na(NE)]))) %>%
    dplyr::mutate(SE=
             ifelse(length(unique(SE))==1, unique(SE), unique(SE[!is.na(SE)]))) %>%
    dplyr::mutate(PL=
             ifelse(length(unique(PL))==1, unique(PL), unique(PL[!is.na(PL)]))) 
    
  return(CatTrait)
}





#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
# Match by binomial species names ------------------------------------------------------------------

## Function to get species that do not match 

NoMatch <- function(Predicts, TraitDF, VClass) {
  
  Y <- Predicts %>% subset(Class==VClass) %>%
    subset(Best_guess_binomial != "")
  
  Y <- unique(Y$Best_guess_binomial)
  
  X <- intersect(Y, TraitDF$Best_guess_binomial)
  
  cat("Predicts contains", length(Y), " species known by binomial name of which",
      length(X), "species intersect with the Trait dataset for", VClass)
  
  Diff <- setdiff(Y, X)
  
  return(Diff)
}


## Function to add the species that appear only in Predicts
.Add_species <- function(Predicts, TraitsDB, VClass) {

  
  Y <- Predicts %>% subset(Class==VClass) %>% subset(Best_guess_binomial != "") %>%
  distinct(Best_guess_binomial)
  
  X <- intersect(Y$Best_guess_binomial, TraitsDB$Best_guess_binomial)
  
  cat("Predicts contains", nrow(Y), " species known by binomial name of which",
      length(X), "species intersect with the Trait dataset for ", VClass)
  
  rm(X)
  
  Y <- as.data.frame(setdiff(Y$Best_guess_binomial, TraitsDB$Best_guess_binomial))
  colnames(Y) <- "Best_guess_binomial"
  
  cat("\nAdding ", nrow(Y), "species")
  
  if (nrow(Y)>=1) {
  
  X.match <- match(Y$Best_guess_binomial, Predicts$Best_guess_binomial)
  Y[, colnames(TraitsDB)[!colnames(TraitsDB) %in% c("Best_guess_binomial")]] <- NA
  
  TraitsDB <- rbind(TraitsDB, Y)
  TraitsDB <- TraitsDB[order(TraitsDB$Best_guess_binomial),]
  row.names(TraitsDB) <- c(1:nrow(TraitsDB))
  TraitsDB %<>% as.data.frame()
  
  } else {cat("\n - - Retuning Traits dataset with no addition")}
  
  return(list(TraitDB=TraitsDB,Y=Y))

}
  


## This function matches species in Predicts against species 
## in the trait datasets by binomial names and extracts the relevant trait
## information
## ARG: Veterbrate class under consideration, Predicts DB, Trait DB, Trait names
## RETURNS: a dataset with all species known by binomial name in Predicts for the 
## considered class and the relevant trait information extracted from the trait DB

.Match_Species <- function (VClass, Predicts, TraitDB) {
  
  #browser()
  
  # Select species in Predicts, that are known by binomial name
  Sp <- Predicts %>% filter(Class==VClass) %>%  
    distinct(Best_guess_binomial) %>% 
    filter(Best_guess_binomial!="") %>%
    droplevels()
  
  # Match by species binomial name
  Y <- intersect(Sp$Best_guess_binomial, TraitDB$Best_guess_binomial)
  Trait_Predicts <- subset(TraitDB, Best_guess_binomial %in% Y)
  
  
  cat("Matching by best binomial guess provides", length(Y), "matches for", VClass, 
      "on a total of", nrow(Sp)," binomial guesses for all traits with diverse % of coverage")
  
  return (Trait_Predicts)
} 



#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
# Get taxa not known by binomial names in the PREDICTS DB and add to the Traits dataset ------------

.GetOther <-function (Predicts, VClass, TraitDataset) {
  
  Vert <- Predicts %>% filter(Class == VClass )
  
  # Get the other species, that are unknown by binomial name, extract Genus, Family and Order
  # Add them to the Traits dataset
  
  # Species known by genus only
  Other_genus <- Vert %>%
    filter(Best_guess_binomial=="") %>%
    filter(Genus != "") %>% distinct(Order, Family, Genus) %>% droplevels()
  
  # Species known by family only
  Other_family <- Vert %>% 
    filter(Best_guess_binomial=="") %>%
    filter(Genus == "") %>%
    filter(Family != "") %>%
    distinct(Order, Family, Genus) %>% droplevels()
  
  # Species known by Order only
  Other_order <- Vert %>% 
    filter(Best_guess_binomial=="") %>%
    filter(Genus == "") %>%
    filter(Family == "") %>%
    distinct(Order, Family, Genus) %>% droplevels()
  
  Other <- rbind(Other_genus, Other_family, Other_order) 
  rm(Other_family, Other_genus, Other_order)
  Other$Binomial_name <- ""
  # Other[, (ncol(Other)+1):(ncol(TraitDataset))] <- NA
  Other[, colnames(TraitDataset)[!colnames(TraitDataset) %in% c("Order", "Family", "Genus", "Binomial_name")]] <- NA
  
  Other <- rbind(TraitDataset, Other)
  
  return(Other)
  
}






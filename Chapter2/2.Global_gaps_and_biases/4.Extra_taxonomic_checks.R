library(dplyr)

ExtraChecks <- function(UncorData, Taxa) {
  
  require(parallel)
  require(rangeBuilder)
  
  x <- synonymMatch(sp=as.character(UncorData$Best_guess_binomial), 
                    db=Taxa, 
                    returnMultiple=FALSE,
                    nthreads = 8, 
                    fuzzyDist = 4)
  
  UncorData$Accepted <- x
  UncorData$Accepted[is.na(UncorData$Accepted)] <- UncorData$Best_guess_binomial[is.na(UncorData$Accepted)]
  UncorData$Accepted <- gsub(pattern="_", x=UncorData$Accepted, replacement = " ")
  
  UncorDataCont <- UncorData %>% 
    dplyr::select(Accepted,
                  Body_size,
                  Life_span_proxy,
                  Litter_size, 
                  Habitat_breadth_IUCN)
  
  UncorDataCont  <- data.table::setDT(UncorDataCont)[, lapply(.SD, mean, na.rm=TRUE),
                                                                     by = Accepted]  %>% as.data.frame()
  
  UncorDataCat <- UncorData %>% 
    dplyr::select(Accepted, Trophic_level, Diel_activity, Specialisation)
  
  UncorDataCat$Trophic_level <-  as.character(UncorDataCat$Trophic_level)
  UncorDataCat$Diel_activity <-  as.character(UncorDataCat$Diel_activity)
  UncorDataCat$Specialisation <-  as.character(UncorDataCat$Specialisation)
  
  UncorDataCat <- UncorDataCat %>% 
    group_by(Accepted) %>% 
    dplyr::mutate(Trophic_level=
                    ifelse(length(unique(Trophic_level))==1, unique(Trophic_level), unique(Trophic_level[!is.na(Trophic_level)]))) %>% 
    dplyr::mutate(Diel_activity=
                    ifelse(length(unique(Diel_activity))==1, unique(Diel_activity), unique(Diel_activity[!is.na(Diel_activity)]))) %>% 
    dplyr::mutate(Specialisation=
                    ifelse(length(unique(Specialisation))==1, unique(Specialisation), unique(Specialisation[!is.na(Specialisation)])))
  
  UncorDataCat <- UncorDataCat %>% distinct() %>% as.data.frame()
  
  UncorData <- cbind(UncorDataCat, UncorDataCont)
  
  return(UncorData)
}

# datasets corrected for taxonomy with IUCN and ITIS checklists
traits.completeness <- readRDS("../../Results/Traits_to_map/traits_completeness_V2.rds")
Amphibians.corrected1 <- traits.completeness$Amphibians
Birds.corrected1 <- traits.completeness$Birds
Mammals.corrected1 <- traits.completeness$Mammals
Reptiles.corrected1 <- traits.completeness$Reptiles

# datasets uncorrected for taxonomy
ReptilesU <- read.csv("../../Data/Trait_data/DataForSharing/UN_Reptiles_final.csv") %>%
  dplyr::select("Body_mass_g", 
                "Longevity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                "Family",
                "Best_guess_binomial")
colnames(ReptilesU)[c(1,2)] <- c("Body_size", "Life_span_proxy")

AmphibiansU <- read.csv("../../Data/Trait_data/DataForSharing/UN_Amphibians_final.csv") %>%
  dplyr::select("Body_length_mm", 
                "Maturity_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation", 
                "Family",
                "Best_guess_binomial")
colnames(AmphibiansU)[c(1,2)] <- c("Body_size", "Life_span_proxy")

BirdsU <- read.csv("../../Data/Trait_data/DataForSharing/UN_Birds_final.csv")%>%
  dplyr::select("Body_mass_g", 
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation", 
                "Family",
                "Best_guess_binomial") 
colnames(BirdsU)[c(1,2)] <- c("Body_size", "Life_span_proxy")

MammalsU <- read.csv("../../Data/Trait_data/DataForSharing/UN_Mammals_final.csv")%>%
  dplyr::select("Body_mass_g", 
                "Generation_length_d",
                "Litter_size",
                "Trophic_level",
                "Diel_activity",
                "Habitat_breadth_IUCN",
                "Specialisation",
                "Family",
                "Best_guess_binomial") 

colnames(MammalsU)[c(1,2)] <- c("Body_size", "Life_span_proxy")

## applying taxonomic corrections using class specific sources from rangeBuilder package

Reptiles.corrected2 <- ExtraChecks(ReptilesU, "squamates")

# exact match | accepted name		9910
# fuzzy match | accepted name		110
# fuzzy match | synonym		4
# multiple hits | advanced search		1
# multiple hits | synonym		41
# not found		602
# strict match | advanced search		9
# strict match | synonym		341

Amphibians.corrected2 <- ExtraChecks(AmphibiansU, "amphibians")

# exact match | accepted name		6318
# fuzzy match | accepted name		207
# fuzzy match | synonym		13
# multiple hits | advanced search		1
# multiple hits | synonym		11
# not found		251
# strict match | advanced search		9
# strict match | synonym		1778

Mammals.corrected2 <- ExtraChecks(MammalsU, "mammals")

# exact match | accepted name		5291
# fuzzy match | accepted name		33
# fuzzy match | advanced search		1
# fuzzy match | synonym		2
# multiple hits | synonym		5
# not found		319
# strict match | synonym		142

Birds.corrected2 <- ExtraChecks(BirdsU, "birds")

# exact match | accepted name		10757
# fuzzy match | accepted name		204
# fuzzy match | synonym		23
# multiple fuzzy hits | accepted name		1
# multiple hits | synonym		1
# not found		1928
# strict match | synonym		624


# # saving

write.csv(Amphibians.correcte2, "../../Results/Additional_taxonomic_checks_rangeBuilder/Amphibians.csv", row.names = FALSE)
write.csv(Reptiles.correcte2, "../../Results/Additional_taxonomic_checks_rangeBuilder/Reptiles.csv", row.names = FALSE)
write.csv(Mammals.correcte2, "../../Results/Additional_taxonomic_checks_rangeBuilder/Mammals.csv", row.names = FALSE)
write.csv(Birds.correcte2, "../../Results/Additional_taxonomic_checks_rangeBuilder/Birds.csv", row.names = FALSE)

# # plotting difference in coverage
Plot.Delta.Cov <- function(TraitData1, TraitData2, TraitData3, Main ){
  
  GGPoptions <- theme_classic()+ theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=13), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=13),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13)) 
  
  Traits_name  <- c("Body_size",
                    "Life_span_proxy",
                    "Litter_size",
                    "Trophic_level",
                    "Diel_activity",
                    "Habitat_breadth_IUCN",
                    "Specialisation")
  
  Completeness1 <- apply(TraitData1[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness1 <- as.data.frame(Completeness1/nrow(TraitData1)*100)
  colnames(Completeness1) <- "Completeness"
  Completeness1$which <- "No correction"
  Completeness1$Trait <- rownames(Completeness1)
  
  Completeness2 <- apply(TraitData2[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness2 <- as.data.frame(Completeness2/nrow(TraitData2)*100)
  colnames(Completeness2) <- "Completeness"
  Completeness2$which <- "IUCN/ITIS corrections"
  Completeness2$Trait <- rownames(Completeness2)
  
  Completeness3 <- apply(TraitData3[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness3 <- as.data.frame(Completeness3/nrow(TraitData3)*100)
  colnames(Completeness3) <- "Completeness"
  Completeness3$which <- "rangeBuilder corrections"
  Completeness3$Trait <- rownames(Completeness3)
  
  Coverage <- rbind(Completeness1, Completeness2, Completeness3)
  
  # order by highest to lowest mean coverage
  Order <- Coverage %>%  group_by(which) %>% 
    summarise(Mean=mean(Completeness)) %>% 
    as.data.frame()
  Order <- Order[order(Order$Mean, decreasing = FALSE),]
  Order <- Order$which
  
  Plotdata <- Coverage %>%
    group_by(which) %>%
    mutate(position=rank(Completeness))
  
  Plotdata <- Plotdata %>% 
    group_by(Trait) %>% 
    arrange(Completeness) %>% 
    mutate(Group=c(1:3))
  
  Lim <- unique(Plotdata$Trait)
  Lim <- textclean::mgsub(Lim, "_", " ")
  
  # plot for all classes together
  p <- ggplot(Plotdata, aes(position, Completeness, fill=which, group=Group)) + 
    geom_bar(stat="identity", position="dodge", width=0.7) +
    GGPoptions +  
    scale_x_discrete(limits=Lim)+
    coord_flip() + ylab("Coverage (%)") + xlab("") +
    scale_fill_viridis_d() + 
    geom_hline(yintercept=50, lty="dashed") +
    theme(legend.title = element_blank()) + ggtitle(Main)
  
  
  return(p)
  
}

p <- ggpubr::ggarrange(
Plot.Delta.Cov(AmphibiansU, Amphibians.corrected1, Amphibians.corrected2, Main="Amphibians"),
Plot.Delta.Cov(ReptilesU, Reptiles.corrected1, Reptiles.corrected2, Main="Reptiles"),
Plot.Delta.Cov(MammalsU, Mammals.corrected1, Mammals.corrected2, Main="Mammals"),
Plot.Delta.Cov(BirdsU, Birds.corrected1, Birds.corrected2, Main="Birds"),
nrow=2, ncol=2, common.legend = TRUE, legend="right")
p
ggsave(p, filename="../../Results/Additional_taxonomic_checks_rangeBuilder/DeltaCov_F3.pdf", width = 12, height=8)


## Taxonomic corrections: overview

Syn_final <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/SingleFinalDataset.csv")

#### Preamble
library(dplyr)
library(reshape)
library(ggplot2)
library(ggpubr)
`%nin%` <- Negate(`%in%`)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=20, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=20), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=20),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=20))

#### Functions

## Overall differences in species number
Delta <- function(SynM, SynA, SynB, SynR) {
  DF <-as.data.frame(matrix(ncol=3, nrow=4)) %>%
    setNames(., c("original","corrected","accepted"))
  rownames(DF) <- c("mammals","birds","reptiles","amphibians")
  
  DF[1,] <- c(length(unique(SynM$Original)),
              length(unique(SynM$CorrectedTypos)),
              length(unique(SynM$Accepted)))
  
  DF[2,] <- c(length(unique(SynB$Original)),
              length(unique(SynB$CorrectedTypos)),
              length(unique(SynB$Accepted)))
  
  DF[3,] <- c(length(unique(SynR$Original)),
              length(unique(SynR$CorrectedTypos)),
              length(unique(SynR$Accepted)))
  
  DF[4,] <- c(length(unique(SynA$Original)),
              length(unique(SynA$CorrectedTypos)),
              length(unique(SynA$Accepted)))
  DF$Class <- rownames(DF)
  
  toplot <- DF %>% select(-corrected)
  toplot <- reshape::melt(toplot)
  toplot$Class <- factor(toplot$Class, levels=c("birds","reptiles","amphibians","mammals"), labels=c("Birds", "Reptiles",
                                                                                                     "Amphibians", "Mammals"))
  
  p <- ggplot(toplot, aes(Class,value, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") + ylab("Names \nacross all datasets") +
    scale_fill_manual(values=c("cornflowerblue","coral"), name="", labels=c("original", "extracted")) + GGPoptions +
    theme(axis.text.x = element_text(angle = 13, hjust = 1)) +
    geom_text(
      aes(label = value, y = value + 0.05),
      position = position_dodge(0.9),
      vjust = -0.5, size=5
    ) + ylim(0, 14200) +
    theme(legend.background = element_rect(color="black"), legend.position = c(0.86505,0.90505))
  
  return(list(p=p, n=toplot))
  
}

## Species with the most number of replicates
MaxSyn <- function(Syn) {
  browser()
  Table <- Syn %>% 
    filter(Accepted!="") %>%
    group_by(Accepted) %>% 
    summarise(Count=n())
  Max <- Table$Accepted[Table$Count==max(Table$Count)]
  if(length(Max)==1) {  All <- Syn$CorrectedTypos[Syn$Accepted==Max]}
  else{All1 <- Syn$CorrectedTypos[Syn$Accepted==Max[1]] %>% as.character()
  All2 <- Syn$CorrectedTypos[Syn$Accepted==Max[2]]%>% as.character()
  return(list(Accepted=Max, Replicates1=All1, Replicates2=All2))
  }
  return(list(Accepted=Max, Replicates=All))
}

## Function to plot distribution of number of names for each accepted name
Dist_NR <- function(Syn) {
  
  GGPoptions <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=20, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=20), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=20),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=20))
  
  Table <- Syn %>% 
    filter(Accepted!="") %>%
    group_by(Class, Accepted) %>%
    summarise(Count=n())
  
  Table <- table(Table$Count, Table$Class) %>%
    as.data.frame() 
  Table_top <- Table %>% 
    #filter(Var1!=1) %>%
    filter(Freq!=0)

  
 p <- ggplot(Table_top, aes(Var1, Freq, fill=Var2)) + geom_point(size=4,pch=21,alpha=0.5) + GGPoptions +
    scale_y_continuous(trans="log10", breaks=c(1,2,5, 10,25, 50, 100,250, 500, 1000,2000,5000)) +
    scale_fill_discrete(name="Class", labels=c("Amphibians","Birds","Mammals","Reptiles")) +
    ylab("Number of species") + xlab("Number of synonyms") +
    theme(legend.background = element_rect(color="black"), legend.position = c(0.84,0.86))
  
    return(list(p=p, Freq=Table_top))
  
}

#### Data

## Load synonym datasets
Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals.csv") %>% mutate(Class="Mammalia")
colnames(Syn_Mammals)[1] <- "Original"
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds.csv") %>% mutate(Class="Aves")
colnames(Syn_Birds)[15] <- "Manual_edits"
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles2.csv") %>% select(-X) %>% mutate(Class="Reptilia")
Syn_Reptiles$Manual_edits[is.na(Syn_Reptiles$Manual_edits)] <- FALSE
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians.csv")%>% mutate(Manual_edits=FALSE, Notes=NA, Class="Amphibia")

## correct a few mistakes in genera
Syn_Mammals$Genus <- stringr::word(Syn_Mammals$Accepted, 1)
Syn_Birds$Genus <- stringr::word(Syn_Birds$Accepted, 1)
Syn_Reptiles$Genus <- stringr::word(Syn_Reptiles$Accepted, 1)
Syn_Amphibians$Genus <- stringr::word(Syn_Amphibians$Accepted, 1)

Syn <- rbind(Syn_Birds, Syn_Mammals, Syn_Reptiles, Syn_Amphibians)
write.csv(Syn, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/SingleFinalDataset.csv", row.names = FALSE)
write.csv(Syn_Birds, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds_final.csv", row.names = FALSE)
write.csv(Syn_Mammals, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals_final.csv", row.names = FALSE)
write.csv(Syn_Amphibians, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians_final.csv", row.names = FALSE)
write.csv(Syn_Reptiles, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles_final.csv", row.names = FALSE)

length(Syn_Mammals$Manual_edits[Syn_Mammals$Manual_edits==TRUE])
length(Syn_Birds$Manual_edits[Syn_Birds$Manual_edits==TRUE])
length(Syn_Reptiles$Manual_edits[Syn_Reptiles$Manual_edits==TRUE])
length(Syn_Amphibians$Manual_edits[Syn_Amphibians$Manual_edits==TRUE])

# ## Load species in compiled trait datasets
# M <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
# B <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
# R <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")
# A <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
# 
# ## mammals: filter non-terrestrial species out
# M <- M %>%
#   filter(Order %nin% c("SIRENIA")) %>%
#   filter(Family %nin% c("OTARIIDAE", "PHOCIDAE", "ODOBENIDAE", 
#                         "BALAENIDAE", "BALAENOPTERIDAE", "ZIPHIIDAE", 
#                         "NEOBALAENIDAE", "DELPHINIDAE", "MONODONTIDAE",
#                         "ESCHRICHTIIDAE", "INIIDAE", "PHYSETERIDAE","LIPOTIDAE",
#                         "PHOCOENIDAE", "PLATANISTIDAE", "PONTOPORIIDAE"))
# 

M <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Mammals_final.csv")
B <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Birds_final.csv")
A <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Amphibians_final.csv")
R <- read.csv("../../Results/Trait_data_for_ms_GLOBALGAPS/Reptiles_final.csv")

SynM <- Syn_final %>% filter(Accepted %in% M$Best_guess_binomial)
SynB <- Syn_final %>% filter(Accepted %in% B$Best_guess_binomial)
SynR <- Syn_final %>% filter(Accepted %in% R$Best_guess_binomial)
SynA <- Syn_final %>% filter(Accepted %in% A$Best_guess_binomial)

Syn_Mammals <- Syn_final %>% 
  filter(Class=="Mammalia")

Syn_Birds <- Syn_final %>% 
  filter(Class=="Aves")

Syn_Reptiles <- Syn_final %>% 
  filter(Class=="Reptilia")

Syn_Amphibians <- Syn_final %>% 
  filter(Class=="Amphibia")



# # # # Script

# 1. Differences in species number before and after taxonomic corrections
ggarrange(
  Delta(Syn_Mammals, Syn_Amphibians, Syn_Birds, Syn_Reptiles)$p,
  Delta(SynM, SynA, SynB, SynR)$p, common.legend = TRUE
)

# 2. Identified synonyms: which species was the most replicated across datasets? 

# For mammals: Tachyoryctes splendens appeared under 12 different names
# across all datasets (and just for species in the trait datasets: same result)
Mammals <- MaxSyn(Syn_Mammals)
Birds <- MaxSyn(Syn_Birds)
Reptiles <- MaxSyn(Syn_Reptiles)
Amphibians <- MaxSyn(Syn_Amphibians)


# 3. Distribution of number of replicates across accepted names

## Plot
p1 <- Delta(SynM, SynA, SynB, SynR)$p + ylab("Binomial names") + xlab(NULL) + labs(tag = "A") + theme(plot.tag.position = c(0.97, 0.96))
p2 <- Dist_NR(Syn_final)$p + labs(tag = "B") + theme(plot.tag.position = c(0.97, 0.96))

X <- Dist_NR(Syn_final)$Freq

p <- ggarrange(p1, p2, widths = c(1/2, 1/2))
ggsave(p, file="../../Results/Plots/Taxonomic_corrections/tax_corrections.pdf", width = 14, height = 6)


## cber talk
synplot <- Dist_NR(Syn)$p + theme(legend.background = element_rect(color="white"), legend.position = c(0.89,0.88))
ggsave(synplot, filename="../../Results/Plots_CBER_talk_20119/synonyms_dist.png", width = 10, height = 7, dpi =1000)

species_numbers = Delta(Syn_Mammals, Syn_Amphibians, Syn_Birds, Syn_Reptiles)$p + ylab("Binomial names")+ xlab("")+ theme(legend.background = element_rect(color="white"), legend.position = c(0.86505,0.90505))
ggsave(species_numbers, filename="../../Results/Plots_CBER_talk_20119/synonymy_delta_species.png", width = 8, height = 7, dpi =1000)


Syn_final %>%  filter(Manual_edits==TRUE) %>% group_by(Class) %>% summarise(C=n())

Syn_final$Notes[Syn_final$Notes==""] <- NA
write.csv(Syn_final, "../../Results/Trait_data_for_ms_GLOBALGAPS/Synonyms.csv", row.names = FALSE)

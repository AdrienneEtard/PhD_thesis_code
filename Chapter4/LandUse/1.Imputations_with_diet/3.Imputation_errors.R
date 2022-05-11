## Imputation errors across 8 imputed datasets for each class

library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape)
library(ggthemes)
library(ggcorrplot)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

Cattraits <- c("Diel_activity", "Specialisation", "Primary_diet", "Habitat_breadth_IUCN", "Diet_breadth")
Conttraits <- c("Body_mass_g", "Lifespan_proxy", "Litter_size")

## Function to get all imputed datasets for a class
Get_all_results <- function(List, Class=c("A", "B", "M", "R")) {
  
  Results <- list()
  Errors <- list()
  
  #browser()
  
  for (i in 1:length(List)) {
    Results[[i]] <- List[[i]][[Class]]$Imputed.Dataset  
    Errors[[i]] <- List[[i]][[Class]]$Imputation.errors
  }
  
  return(list(Results=Results, Errors=Errors))
} 


############################################################################################################################################

## Load imputed datasets (8); retrieve errors[[8]]
Imputed <- readRDS("../../Results/Imputed_traits/List_of_8_sets.rds")

for(i in 1:8){
  
  Imputed[[i]][["M"]]$Imputed.Dataset$Class <- "Mammals"
  colnames(Imputed[[i]][["M"]]$Imputed.Dataset)[8] <- "Lifespan_proxy"

  Imputed[[i]][["B"]]$Imputed.Dataset$Class <- "Birds"
  colnames(Imputed[[i]][["B"]]$Imputed.Dataset)[9] <- "Lifespan_proxy"
 
  Imputed[[i]][["R"]]$Imputed.Dataset$Class <- "Reptiles"
  colnames(Imputed[[i]][["R"]]$Imputed.Dataset)[9] <- "Lifespan_proxy"
 
  Imputed[[i]][["A"]]$Imputed.Dataset$Class <- "Amphibians"
  colnames(Imputed[[i]][["A"]]$Imputed.Dataset)[10] <- "Lifespan_proxy"
}

## Amphibians
Amphibians <- Get_all_results(Imputed, "A")$Errors
Errors_amp <- data.table::rbindlist(Amphibians)%>% as.data.frame()
Errors_amp8 <- Errors_amp[8,]; names(Errors_amp8) <- colnames(Errors_amp)
Errors_amp8$Class <- "Amphibians"
colnames(Errors_amp8)[7] <- "Lifespan_proxy MSE"

## Reptiles
Reptiles <- Get_all_results(Imputed, "R")$Errors
Errors_rep <- data.table::rbindlist(Reptiles)%>% as.data.frame()
Errors_rep8 <- Errors_rep[8,]; names(Errors_rep8) <- colnames(Errors_rep)
Errors_rep8$Class <- "Reptiles"
colnames(Errors_rep8)[5] <- "Lifespan_proxy MSE"

## Mammals
Mammals <- Get_all_results(Imputed, "M")$Errors
Errors_mam <- data.table::rbindlist(Mammals)%>% as.data.frame()
Errors_mam8 <- Errors_mam[8,]; names(Errors_mam8) <- colnames(Errors_mam)
Errors_mam8$Class <- "Mammals"
colnames(Errors_mam8)[5] <- "Lifespan_proxy MSE"

## Birds
Birds <- Get_all_results(Imputed, "B")$Errors
Errors_bir <- data.table::rbindlist(Birds)%>% as.data.frame()
Errors_bir8 <- Errors_bir[8,]; names(Errors_bir8) <- colnames(Errors_bir)
Errors_bir8$Class <- "Birds"
colnames(Errors_bir8)[6] <- "Lifespan_proxy MSE"

## All PFC errors for the 8th imputed set -- plot
Cattraits <- c(paste(Cattraits, "PFC"), "Class")
PFC_Errors8 <- rbind(Errors_amp8[,which(colnames(Errors_amp8) %in% Cattraits)],
                     Errors_bir8[,which(colnames(Errors_bir8) %in% Cattraits)], 
                     Errors_mam8[,which(colnames(Errors_mam8) %in% Cattraits)],
                     Errors_rep8[,which(colnames(Errors_rep8) %in% Cattraits)])

PFC <- melt(PFC_Errors8)
PFC$value <- PFC$value*100

pPFC <- 
ggplot(PFC, aes(variable, value, fill = Class, group=Class)) +
  geom_bar(stat="identity", position = "dodge")+
  scale_x_discrete(labels=c("Diet breadth", "Habitat breadth", "Primary diet","Use of artificial habitats", "Diel activity")) +
  ylim(0,max(PFC$value)+2) +
  xlab("") + ylab("OOB % falsely classified") +
  GGPoptions + scale_fill_colorblind() +
  geom_text(aes(label = round(value, digits = 1), x = variable),
             position = position_dodge(width = 1), size=3, angle=0, vjust=-1) + ggtitle("(a) Proportion of Falsely Classified (PFC)")
  #geom_text(aes(label = round(value, digits = 1), x=variable, group=Class), vjust = 0,  position = "dodge")

## getting NRMSE from MSE
Conttraits <- c(paste(Conttraits, "MSE"), "Class")
MSE_Errors8 <- rbind(Errors_amp8[,which(colnames(Errors_amp8) %in% Conttraits)],
                     Errors_bir8[,which(colnames(Errors_bir8) %in% Conttraits)], 
                     Errors_mam8[,which(colnames(Errors_mam8) %in% Conttraits)],
                     Errors_rep8[,which(colnames(Errors_rep8) %in% Conttraits)])


## square-rooting to obtain RMSE -- need to sort out longevity
RMSE_Errors8 <- apply(MSE_Errors8[,c(1:3)], 2, sqrt) %>%  
  as.data.frame()
RMSE_Errors8$Class <- MSE_Errors8$Class

## divide by variance of the known trait distribution 
Mammals <- read.csv("../../Results/Traits_with_phy_eigenvectors/Mammals.csv") %>% 
  dplyr::select(-Trophic_level.Elton)
Birds <- read.csv("../../Results/Traits_with_phy_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/Traits_with_phy_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/Traits_with_phy_eigenvectors/Reptiles.csv")

var_LS <- c(var(Amphibians$Litter_size, na.rm = TRUE),
            var(Birds$Litter_size, na.rm = TRUE),
            var(Mammals$Litter_size, na.rm = TRUE),
            var(Reptiles$Litter_size, na.rm = TRUE))

var_BM <- c(var(Amphibians$Body_mass_g, na.rm = TRUE),
            var(Birds$Body_mass_g, na.rm = TRUE),
            var(Mammals$Body_mass_g, na.rm = TRUE),
            var(Reptiles$Body_mass_g, na.rm = TRUE))

var_LP <- c(var(Amphibians$Maturity_d, na.rm = TRUE),
            var(Birds$Generation_length_d, na.rm = TRUE),
            var(Mammals$Generation_length_d, na.rm = TRUE),
            var(Reptiles$Max_longevity_d, na.rm = TRUE))

RMSE_Errors8$NRMSE_Litter_size <- RMSE_Errors8$`Litter_size MSE`/sqrt(var_LS)
RMSE_Errors8$NRMSE_Body_mass_g <- RMSE_Errors8$`Body_mass_g MSE`/sqrt(var_BM)
RMSE_Errors8$NRMSE_Lifespan_proxy <- RMSE_Errors8$`Lifespan_proxy MSE`/sqrt(var_LP)

NRMSE <- melt(RMSE_Errors8[, c(4:7)])

pMSE <- 
ggplot(NRMSE, aes(variable, value, fill = Class)) +
  geom_bar(stat="identity", position = "dodge")+
  scale_x_discrete(labels=c("Litter/clutch size", "Body mass", "Lifespan")) +
  xlab("") + ylab("NRMSE") +
  GGPoptions + scale_fill_colorblind() + ylim(c(0,1))+
  geom_text(aes(label = round(value, digits = 2), x = variable),
            position = position_dodge(width = 1), size=3, angle=0, vjust=-1) +ggtitle("(b) Normalised Root Mean Squared Error (NRMSE)")


# library(patchwork)
# pPFC/pMSE

plotSIerrors <- ggarrange(pPFC, pMSE, common.legend=TRUE, legend = "right", nrow=2)
ggsave(plotSIerrors, filename="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/Imputation_errors.pdf",
       width=9, height=7)

############################################################################################################################################

## imputation congruence among sets of imputed datasets -- continous traits


## get lists of imputed datasets
GetListDatasets <- function(Class, Imputed) {
  ResList <- list()
  for(i in 1:8){
    ResList[[i]] <- Imputed[[i]][[Class]]$Imputed.Dataset
  }
  return(ResList)
} 

MammalsList <- GetListDatasets("M", Imputed)
BirdsList <- GetListDatasets("B", Imputed)
ReptilesList <- GetListDatasets("R", Imputed)
AmphibiansList <- GetListDatasets("A", Imputed)

## to extract imputed values only from imputed datasets
GetImputedValues <- function(TraitsNotImp, TraitsImp, Trait) {
  colnames(TraitsNotImp)[colnames(TraitsNotImp)==Trait] <- "The_trait"
  colnames(TraitsImp)[colnames(TraitsImp)==Trait] <- "The_trait"
  SpeciesImputed <- TraitsNotImp$Best_guess_binomial[is.na(TraitsNotImp$The_trait)]
  Imputed_subset <- TraitsImp[TraitsImp$Best_guess_binomial %in% SpeciesImputed,]
  colnames(Imputed_subset)[colnames(Imputed_subset)=="The_trait"] <- Trait
  return(Imputed_subset)
}

## function to bind together all 8 imputed sets
BindVals <- function(List8Imp, Trait){
  
  Init <- vector(length = nrow(List8Imp[[1]]))
  
  List8Imp <- lapply(List8Imp,
                     FUN = function(x){ 
                       colnames(x)[colnames(x)==Trait] <- "The_trait"
                       return(x)})  
  
  for (i in 1:length(List8Imp)) {
    Init <- cbind(Init, List8Imp[[i]]$"The_trait") 
  }
  return(Init[, c(2:8)])
}


############################################################################################################################################

## Congruence for continous traits

GetCorrelationMatrix <- function(ListImp, TraitsNotImp, Trait){
  ListTrait <- lapply(FUN = GetImputedValues, X=ListImp, TraitsNotImp=TraitsNotImp, Trait=Trait)
  ListTrait <- BindVals(ListTrait, Trait)
  return(cor(ListTrait))
}

## for body mass
bm_mammals <- GetCorrelationMatrix(MammalsList, Mammals, "Body_mass_g")
bm_birds <- GetCorrelationMatrix(BirdsList, Birds, "Body_mass_g")
bm_amphibians <- GetCorrelationMatrix(AmphibiansList, Amphibians, "Body_mass_g")
bm_reptiles <- GetCorrelationMatrix(ReptilesList, Reptiles, "Body_mass_g")

## for litter/clutch size
lcs_mammals <- GetCorrelationMatrix(MammalsList, Mammals, "Litter_size")
lcs_birds <- GetCorrelationMatrix(BirdsList, Birds, "Litter_size")
lcs_amphibians <- GetCorrelationMatrix(AmphibiansList, Amphibians, "Litter_size")
lcs_reptiles <- GetCorrelationMatrix(ReptilesList, Reptiles, "Litter_size")

## for lifespan proxy
colnames(Mammals)[10] <- "Lifespan_proxy"
colnames(Birds)[11] <- "Lifespan_proxy"
colnames(Amphibians)[10] <- "Lifespan_proxy"
colnames(Reptiles)[9] <- "Lifespan_proxy"

lp_mammals <- GetCorrelationMatrix(MammalsList, Mammals, "Lifespan_proxy")
lp_birds <- GetCorrelationMatrix(BirdsList, Birds, "Lifespan_proxy")
lp_amphibians <- GetCorrelationMatrix(AmphibiansList, Amphibians, "Lifespan_proxy")
lp_reptiles <- GetCorrelationMatrix(ReptilesList, Reptiles, "Lifespan_proxy")


## getting min, max and median of correlation coefficients across sets of imputed values for each trait

ExtractCors <- function(A, B, M, R, Trait) {
  
  ReturnCor <- function(X, Class, Trait) {
    Cors <- c(range(X[lower.tri(X, diag = FALSE)]), median(X[lower.tri(X, diag = FALSE)]))
    names(Cors) <- c("Min", "Max", "Median")
    Cors$Class <- Class
    Cors$Trait <- Trait
    Cors <- as.data.frame(Cors)
   return(Cors)
  }
  
  Amphibians <- ReturnCor(A, "Amphibians", Trait = Trait)
  Birds <- ReturnCor(B, "Birds", Trait = Trait)
  Mammals <- ReturnCor(M, "Mammals", Trait = Trait)
  Reptiles <- ReturnCor(R, "Reptiles", Trait = Trait)
  Res <- rbind(Amphibians, Birds, Mammals, Reptiles)
  return(Res)
}

bm_cors <- ExtractCors(bm_amphibians, bm_birds, bm_mammals, bm_reptiles, "Body_mass_g")
lcs_cors <- ExtractCors(lcs_amphibians, lcs_birds, lcs_mammals, lcs_reptiles, "Litter_size")
lp_cors <- ExtractCors(lp_amphibians, lp_birds, lp_mammals, lp_reptiles, "Lifespan_proxy")

AllCors <- rbind(bm_cors, lcs_cors, lp_cors)

## plotting

p_congruence_continuous <- 
  ggplot(AllCors, aes(Trait, Median, col = Class, ymin=Min, ymax=Max)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#EBEBEB", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#EBEBEB", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#EBEBEB", col=NA) +
  # geom_hline(yintercept=0, col="grey") +
  geom_hline(yintercept=0.25, col="grey") +
  geom_hline(yintercept=0.5, col="grey") +
  geom_hline(yintercept=0.75, col="grey") +
  geom_hline(yintercept=1, col="grey") +
  geom_point(position = position_dodge(width=0.4))+
  geom_errorbar(width=0.1, position = position_dodge(width=0.4)) +
  scale_x_discrete(labels=c("Body mass", "Litter/clutch size", "Lifespan proxy")) +
  xlab("") + ylab("Correlation coefficients among pairs \nof imputed value sets (median +/- min and max)") +
  GGPoptions + scale_colour_colorblind()

############################################################################################################################################

## imputation congruence among sets of imputed datasets -- categorical traits

## for each combination of imputed datsets, check if values are similar or not and get proportion of values that differ

GetComparison_Cat <- function(ImputedTraitsList, Trait){
  
  ## create all pairwise comparisons between 1 and 8
  Combinations <- as.data.frame(t(combn(1:8, m=2)))
  Combinations$Results_comparison_same <- NA
  Combinations$Results_comparison_different <- NA
  
  ImputedTraitsList <- lapply(ImputedTraitsList, 
                              FUN = function(x){ 
                                colnames(x)[colnames(x)==Trait] <- "The_trait"
                                return(x)
                                })
 
  for (i in 1:nrow(Combinations)) {
    v1 <- Combinations$V1[i]
    v2 <- Combinations$V2[i]
    Comp <- ImputedTraitsList[[v2]]$The_trait==ImputedTraitsList[[v1]]$The_trait
    Table <- table(Comp)
    if(length(Table)==2) {
      Combinations$Results_comparison_same[i] <- table(Comp)[2]
      Combinations$Results_comparison_different[i] <- table(Comp)[1]  
    }
    
    if (length(Table)==1) {
      if(names(Table)=="TRUE"){
        Combinations$Results_comparison_same[i] <- table(Comp)[1]
      } else{
        Combinations$Results_comparison_different[i] <- table(Comp)[1]  
      }
    }
  }
  
  Combinations$Results_comparison_different[is.na(Combinations$Results_comparison_different)] <- 0
  Combinations$Proportion_classified_different <- Combinations$Results_comparison_different/
    Combinations$Results_comparison_same*100
  return(Combinations)
}

Run_GetComparison <- function(ListImp, TraitsNotImp, Trait, Class){
  ListTrait <- lapply(FUN = GetImputedValues, X=ListImp, TraitsNotImp=TraitsNotImp, Trait=Trait)
  Res <- GetComparison_Cat(ListTrait, Trait)
  Res$Class <- Class
  Res$Trait <- Trait
 return(Res)
}


# ## trophic levels 
# TL_M <- Run_GetComparison(ListImp=MammalsList, TraitsNotImp = Mammals, Trait = "Trophic_level", Class = "Mammals")
# TL_A <- Run_GetComparison(ListImp=AmphibiansList, TraitsNotImp = Amphibians, Trait = "Trophic_level", Class = "Amphibians")
# TL_B <- Run_GetComparison(ListImp=BirdsList, TraitsNotImp = Birds, Trait = "Trophic_level", Class = "Birds")
# TL_R <- Run_GetComparison(ListImp=ReptilesList, TraitsNotImp = Reptiles, Trait = "Trophic_level", Class = "Reptiles")


## habitat breadth
hb_m <- Run_GetComparison(ListImp=MammalsList, TraitsNotImp = Mammals, Trait = "Habitat_breadth_IUCN", Class = "Mammals")
hb_a <- Run_GetComparison(ListImp=AmphibiansList, TraitsNotImp = Amphibians, Trait = "Habitat_breadth_IUCN", Class = "Amphibians")
hb_b <- Run_GetComparison(ListImp=BirdsList, TraitsNotImp = Birds, Trait = "Habitat_breadth_IUCN", Class = "Birds")
hb_r <- Run_GetComparison(ListImp=ReptilesList, TraitsNotImp = Reptiles, Trait = "Habitat_breadth_IUCN", Class = "Reptiles")
hb <- rbind(hb_a, hb_b, hb_m, hb_r)

## specialisation
sp_m <- Run_GetComparison(ListImp=MammalsList, TraitsNotImp = Mammals, Trait = "Specialisation", Class = "Mammals")
sp_a <- Run_GetComparison(ListImp=AmphibiansList, TraitsNotImp = Amphibians, Trait = "Specialisation", Class = "Amphibians")
sp_b <- Run_GetComparison(ListImp=BirdsList, TraitsNotImp = Birds, Trait = "Specialisation", Class = "Birds")
sp_r <- Run_GetComparison(ListImp=ReptilesList, TraitsNotImp = Reptiles, Trait = "Specialisation", Class = "Reptiles")
sp <- rbind(sp_a, sp_b, sp_m, sp_r)

## diel activity
da_m <- Run_GetComparison(ListImp=MammalsList, TraitsNotImp = Mammals, Trait = "Diel_activity", Class = "Mammals")
da_a <- Run_GetComparison(ListImp=AmphibiansList, TraitsNotImp = Amphibians, Trait = "Diel_activity", Class = "Amphibians")
da_b <- Run_GetComparison(ListImp=BirdsList, TraitsNotImp = Birds, Trait = "Diel_activity", Class = "Birds")
da_r <- Run_GetComparison(ListImp=ReptilesList, TraitsNotImp = Reptiles, Trait = "Diel_activity", Class = "Reptiles")
da <- rbind(da_a, da_b, da_m, da_r)

## primary diet
pd_m <- Run_GetComparison(ListImp=MammalsList, TraitsNotImp = Mammals, Trait = "Primary_diet", Class = "Mammals")
pd_a <- Run_GetComparison(ListImp=AmphibiansList, TraitsNotImp = Amphibians, Trait = "Primary_diet", Class = "Amphibians")
pd_b <- Run_GetComparison(ListImp=BirdsList, TraitsNotImp = Birds, Trait = "Primary_diet", Class = "Birds")
pd_r <- Run_GetComparison(ListImp=ReptilesList, TraitsNotImp = Reptiles, Trait = "Primary_diet", Class = "Reptiles")
pd <- rbind(pd_a, pd_b, pd_m, pd_r)

## diet breadth
db_m <- Run_GetComparison(ListImp=MammalsList, TraitsNotImp = Mammals, Trait = "Diet_breadth", Class = "Mammals")
db_a <- Run_GetComparison(ListImp=AmphibiansList, TraitsNotImp = Amphibians, Trait = "Diet_breadth", Class = "Amphibians")
db_b <- Run_GetComparison(ListImp=BirdsList, TraitsNotImp = Birds, Trait = "Diet_breadth", Class = "Birds")
db_r <- Run_GetComparison(ListImp=ReptilesList, TraitsNotImp = Reptiles, Trait = "Diet_breadth", Class = "Reptiles")
db <- rbind(db_a, db_b, db_m, db_r)


Proportion_classified_diff <- expand.grid(Class=c("Amphibians", "Birds", "Mammals", "Reptiles"),
                                           Trait=c("Primary_diet", "Diet_breadth","Diel_activity",
                                                   "Habitat_breadth_IUCN", "Specialisation"))

Proportion_classified_diff$Max <- NA
Proportion_classified_diff$Min <- NA 
Proportion_classified_diff$Median <- NA 

## fill in table
FillDF <- function(Trait, M, B, R, A) {
  
  Proportion_classified_diff$Max[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Mammals"] <- max(M$Proportion_classified_different)
  Proportion_classified_diff$Min[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Mammals"] <- min(M$Proportion_classified_different)
  Proportion_classified_diff$Median[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Mammals"] <- median(M$Proportion_classified_different)
  
  Proportion_classified_diff$Max[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Reptiles"] <- max(R$Proportion_classified_different)
  Proportion_classified_diff$Min[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Reptiles"] <- min(R$Proportion_classified_different)
  Proportion_classified_diff$Median[Proportion_classified_diff$Trait==Trait &
                                      Proportion_classified_diff$Class=="Reptiles"] <- median(R$Proportion_classified_different)
  
  Proportion_classified_diff$Max[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Birds"] <- max(B$Proportion_classified_different)
  Proportion_classified_diff$Min[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Birds"] <- min(B$Proportion_classified_different)
  Proportion_classified_diff$Median[Proportion_classified_diff$Trait==Trait &
                                      Proportion_classified_diff$Class=="Birds"] <- median(B$Proportion_classified_different)
  
  Proportion_classified_diff$Max[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Amphibians"] <- max(A$Proportion_classified_different)
  Proportion_classified_diff$Min[Proportion_classified_diff$Trait==Trait &
                                   Proportion_classified_diff$Class=="Amphibians"] <- min(A$Proportion_classified_different)
  Proportion_classified_diff$Median[Proportion_classified_diff$Trait==Trait &
                                      Proportion_classified_diff$Class=="Amphibians"] <- median(A$Proportion_classified_different)
  
  return(Proportion_classified_diff)
  
}

Proportion_classified_diff <- FillDF(Trait="Primary_diet", M = pd_m, B = pd_b, R = pd_r, A = pd_a)
Proportion_classified_diff <- FillDF(Trait="Diet_breadth", M = db_m, B = db_b, R = db_r, A = db_a)
Proportion_classified_diff <- FillDF(Trait="Diel_activity", M = da_m, B = da_b, R = da_r, A = da_a)
Proportion_classified_diff <- FillDF(Trait="Habitat_breadth_IUCN", M = hb_m, B = hb_b, R = hb_r, A = hb_a)
Proportion_classified_diff <- FillDF(Trait="Specialisation", M = sp_m, B = sp_b, R = sp_r, A = sp_a)

## plotting 

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

p_congruence_categories <- 
ggplot(Proportion_classified_diff, aes(Trait,Median, col = Class, ymin=Min, ymax=Max)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#EBEBEB", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_hline(yintercept=0, col="grey") +
  #geom_hline(yintercept=5, col="grey") +
  geom_hline(yintercept=10, col="grey") +
  #geom_hline(yintercept=15, col="grey") +
  geom_hline(yintercept=20, col="grey") +
  #geom_hline(yintercept=25, col="grey") +
  geom_hline(yintercept=30, col="grey") +
  #geom_hline(yintercept=35, col="grey") +
  geom_hline(yintercept=40, col="grey") +
  #geom_hline(yintercept=45, col="grey") +
  geom_hline(yintercept=50, col="grey") +
  #geom_hline(yintercept=55, col="grey") +
  geom_point(position = position_dodge(width=0.4))+
  geom_errorbar(width=0.1, position = position_dodge(width=0.4)) +
  scale_x_discrete(labels=c("Primary diet", "Diet breadth", "Diel activity", "Habitat breadth", "Use of artifical habitats")) +
  xlab("") + ylab("Proportion of species classified differently (%) \n(median +/- min and max)") +
  GGPoptions + scale_colour_colorblind()


p_congruence_continuous <- 
  ggplot(AllCors, aes(Trait, Median, col = Class, ymin=Min, ymax=Max)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#EBEBEB", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#EBEBEB", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#EBEBEB", col=NA) +
  # geom_hline(yintercept=0, col="grey") +
  geom_hline(yintercept=0.25, col="grey") +
  geom_hline(yintercept=0.5, col="grey") +
  geom_hline(yintercept=0.75, col="grey") +
  geom_hline(yintercept=1, col="grey") +
  geom_point(position = position_dodge(width=0.4))+
  geom_errorbar(width=0.1, position = position_dodge(width=0.4)) +
  scale_x_discrete(labels=c("Body mass", "Litter/clutch size", "Lifespan proxy")) +
  xlab("") + ylab("Correlation coefficients among pairs \nof imputed value sets (median +/- min and max)") +
  GGPoptions + scale_colour_colorblind()

## assembling congruence plots together for SI

plots_congruence <- 
ggarrange( p_congruence_categories + ggtitle("(a) Categorical traits"),
  p_congruence_continuous + ggtitle("(b) Continuous traits"),
          common.legend = TRUE, nrow=2, legend="right")

ggsave(plots_congruence, filename ="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/Imputation_congruence.pdf",
       width=9, height = 8)

## Function for differences in species number
Delta_Number <- function(NotCorrected, Corrected, Class, All_or_Predicts, C_Predictsdata, UN_Predictsdata) {
  
  if(All_or_Predicts=="all species") {
    
    print(paste("Not corrected:", length(unique(NotCorrected$Best_guess_binomial)), "species for", All_or_Predicts, Class))
    print(paste("Corrected:", length(unique(Corrected$Best_guess_binomial)), "species for", All_or_Predicts, Class))
    print(paste("Correcting reduced species redundancy by:", length(unique(NotCorrected$Best_guess_binomial))- length(unique(Corrected$Best_guess_binomial)), "species" ))
    
  }
  
  if(All_or_Predicts=="Predicts") {
    
    C_Sp <- unique(C_Predictsdata$Best_guess_binomial[C_Predictsdata$Class==Class])
    Corrected <- Corrected %>% filter(Best_guess_binomial %in% C_Sp)
    Corrected <- Corrected$Best_guess_binomial
    
    UN_Sp <- unique(UN_Predictsdata$Best_guess_binomial[UN_Predictsdata$Class==Class])
    NotCorrected <- NotCorrected %>% filter(Best_guess_binomial %in% UN_Sp)
    NotCorrected <- NotCorrected$Best_guess_binomial
    
    print(paste("Not corrected:", length(UN_Sp), "species for", All_or_Predicts, Class))
    print(paste("Corrected:", length(C_Sp), "species for", All_or_Predicts, Class))
    print(paste("Correcting reduced species redundancy by:", length(UN_Sp)- length(C_Sp), "species" ))
    print(paste("Corrected:", length(Corrected), "species intersect"))
    print(paste("Uncorrected:", length(NotCorrected), "species intersect"))
    
  }
  
}


## Functions to plot coverage

## 1. Species representation in the phylogenies
Phylo_cov <- function(Mammals, Birds, Reptiles, Amphibians) {
  
  CovM <- sum(!is.na(Mammals$EV_1))/nrow(Mammals)*100 %>%
    as.data.frame() %>%
    setNames(., "Mammals")
  
  CovR <- sum(!is.na(Reptiles$EV_1))/nrow(Reptiles)*100 %>%
    as.data.frame() %>%
    setNames(., "Reptiles")
  
  CovB <- sum(!is.na(Birds$EV_1))/nrow(Birds)*100 %>%
    as.data.frame() %>%
    setNames(., "Birds")
  
  CovA <- sum(!is.na(Amphibians$EV_1))/nrow(Amphibians)*100 %>%
    as.data.frame() %>%
    setNames(., "Amphibians")
  
  Cov <- cbind(CovM, CovR, CovB, CovA)
  Cov <- t(Cov) %>% as.data.frame()
  colnames(Cov) <- "Coverage"
  Cov$Taxon[1:4] <- rownames(Cov)
  rownames(Cov) <- c(1:4)
  
  Cov <- Cov[order(Cov$Coverage, decreasing=FALSE), , drop=FALSE] %>%
    as.data.frame()
  
  return(Cov)
  
  
}

Phylo_Delta <- function(Mammals, Birds, Reptiles, Amphibians, Predicts, PredictsTRUE, UN.Mammals, UN.Birds, UN.Reptiles, UN.Amphibians, UN.Predicts ) {
  
  if(PredictsTRUE) {
    
    # Corrected 
    PCMammals <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Mammalia"])
    Mammals <- subset(Mammals, Best_guess_binomial %in% PCMammals)
    
    PCBirds <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Aves"])
    Birds <- subset(Birds, Best_guess_binomial %in% PCBirds)
    
    PCReptiles <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Reptilia"])
    Reptiles <- subset(Reptiles, Best_guess_binomial %in% PCReptiles)
    
    PCAmphibians <- unique(Predicts$Best_guess_binomial[Predicts$Class=="Amphibia"])
    Amphibians <- subset(Amphibians, Best_guess_binomial %in% PCAmphibians)
    
    # Uncorrected
    PUMammals <- unique(UN.Predicts$Best_guess_binomial[UN.Predicts$Class=="Mammalia"])
    UN.Mammals <- subset(UN.Mammals, Best_guess_binomial %in% PUMammals)
    
    PUBirds <- unique(UN.Predicts$Best_guess_binomial[UN.Predicts$Class=="Aves"])
    UN.Birds <- subset(UN.Birds, Best_guess_binomial %in% PUBirds)
    
    PUReptiles <- unique(UN.Predicts$Best_guess_binomial[UN.Predicts$Class=="Reptilia"])
    UN.Reptiles <- subset(UN.Reptiles, Best_guess_binomial %in% PUReptiles)
    
    PUAmphibians <- unique(UN.Predicts$Best_guess_binomial[UN.Predicts$Class=="Amphibia"])
    UN.Amphibians <- subset(UN.Amphibians, Best_guess_binomial %in% PUAmphibians)
    
  }
  
  
  # Coverage corrected
  Cov_corrected <-  Phylo_cov(Mammals, Birds, Reptiles, Amphibians)
  Cov_corrected$IsCorrected <- "Corrected"
  
  # Coverage uncorrected
  Cov_uncorrected <- Phylo_cov(UN.Mammals, UN.Birds, UN.Reptiles, UN.Amphibians)
  Cov_uncorrected$IsCorrected <- "Uncorrected"
  
  Cov <- rbind(Cov_uncorrected,Cov_corrected)
  Cov <- Cov[order(Cov$Coverage, decreasing=FALSE), ]
  rownames(Cov) <- NULL
  
  return(Cov)
  
  
}

PlotPhyloCov <- function(Cor, Uncor, FontSize) {
  
  # browser()
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=FontSize, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,5,0,"pt"), size=FontSize), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=FontSize),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=FontSize))
  
  p1 <- ggplot(Cor, aes(y = Coverage, x = Taxon, fill=IsCorrected)) +
    geom_hline(yintercept = 50, linetype="dashed", col="darkgrey") +
    geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE, width=0.2), width=0.5) +
    scale_x_discrete(limits = c("Amphibians","Reptiles","Birds","Mammals")) +
    GGPoptions  + xlab("") + ylab("% species represented in phylogenies") +
    geom_hline(yintercept = 100, linetype="dashed") +
    coord_flip() +
    scale_fill_discrete(name = "Taxonomic correction") +
    scale_fill_hue(name="Taxonomy", labels = c("Corrected", "Uncorrected")) #, values = c("firebrick", "royalblue")
    #labs(tag = "A") + theme(plot.tag.position = c(0.13,0.97))
    
  
 p2 <- ggplot(Uncor, aes(y = Coverage, x = Taxon, fill=IsCorrected)) +
    geom_hline(yintercept = 50, linetype="dashed", col="darkgrey") +
    geom_bar(stat = "identity", position = position_dodge2(reverse = TRUE, width=0.2), width=0.5) +
    scale_x_discrete(limits = c("Amphibians","Reptiles","Birds","Mammals"), labels=c("A", "R", "B", "M")) +
    theme_classic() + GGPoptions + 
    geom_hline(yintercept = 100, linetype="dashed") +
    coord_flip() +
    scale_fill_hue(name="Taxonomy", labels = c("Corrected", "Uncorrected"))+ #values = c("firebrick", "royalblue")
    guides(fill=FALSE) + xlab("") + ylab("") 
    #labs(tag = "B")  + 
    #theme(plot.tag.position = c(0.13,0.97))
   

 p <- ggdraw() +
   draw_plot(p1) +
   draw_plot(p2, x = 0.49, y = .16, width = .24, height = .40)
 
  return(p)
}



## 2. Plotting trait coverage  - to select predictors in the imputations: here, no comparison of taxonomic corrections

Plot.Cov <- function(TraitData, Traits, Main, PredictsTrue, Predicts) {
  
  if(PredictsTrue) {
    
    Sp <- unique(Predicts$Best_guess_binomial)
    Sp <- intersect(TraitData$Best_guess_binomial, Sp)
    TraitData <- TraitData %>%
      filter(Best_guess_binomial %in% Sp)
    
    print(paste("Predicts:", length(Sp), "species"))
    
  }
  
  Names <- as.data.frame(c("Body_mass_g", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm","Body_length_mm",
                           "Svl_length_mm","Generation_length_d", "Longevity_d", "Maturity_d", "AFR_d",
                           "Litter_size", "Diel_activity", "Trophic_level","Diet_breadth", "Specialisation","Habitat_breadth_IUCN", "EV_1"))
  colnames(Names) <- "Original"
  Names$FP <- c("Body mass", "Svl length", "Forearm length", "Head length", "Body length", "Svl length", "Generation length", "Longevity", "Sexual maturity age",
                "Age 1st reproduction","Litter/clutch size",
                "Diel activity", "Trophic level", "Diet breadth","Specialisation", "Habitat breadth", "Phylogenetic position")
  
  # Traits <- Traits[-which(Traits=="EV_1")]
  
  Completeness <- apply(TraitData[, Traits], 2,  function(y) sum(!is.na(y))) 
  Completeness <- as.data.frame(Completeness/nrow(TraitData)*100)
  colnames(Completeness) <- "Completeness"
  Completeness <- Completeness[order(Completeness, decreasing=FALSE), , drop=FALSE]
  
  # # Phylo information coverage
  # browser()
  # ToBind <- sum(!is.na(TraitData$EV_1))/nrow(TraitData)*100 %>%
  #   as.data.frame() %>%
  #   setNames(., "Completeness")
  # Completeness <- rbind(Completeness, ToBind)
  # row.names(Completeness)[nrow(Completeness)] <- "EV_1"
  
  Names_plot <- as.data.frame(row.names(Completeness))
  colnames(Names_plot) <- "Or"
  Names$Original <- as.character(Names$Original)
  Names_plot$Or <- as.character(Names_plot$Or)
  for (i in 1:nrow(Names_plot)) {Names_plot$TP[i] <- Names$FP[Names$Original==Names_plot$Or[i]]}
  
  barplot(Completeness$Completeness, horiz = TRUE, 
          xlim = c(0,100), las=1,col="#00BFC4", main = Main, names.arg = Names_plot$TP)
  
  
  abline(v=100, lty="dotted")
  # abline(v=50, lty="dotted", col="grey")
}

## 3. Plotting trait coverage, before VS after taxonomic corrections
Plot.Delta.Cov <- function(TraitData1, TraitData2, Traits_name, PredictsTrue, Predicts_UN, Predicts_C, Main ){#Order=c("By1", "By2")) {
  
  if(PredictsTrue) {
    
    SpUN <- unique(Predicts_UN$Best_guess_binomial)
    SpUN <- intersect(TraitData1$Best_guess_binomial, SpUN)
    TraitData1 <- TraitData1 %>%
      filter(Best_guess_binomial %in% SpUN)
    
    SpC <- unique(Predicts_C$Best_guess_binomial)
    SpC <- intersect(TraitData2$Best_guess_binomial, SpC)
    TraitData2 <- TraitData2 %>%
      filter(Best_guess_binomial %in% SpC)
    
    print(paste("Predicts corrected:", length(SpC), "species"))
    print(paste("Predicts uncorrected:", length(SpUN), "species"))
  }
  
  Names <- c("Body_mass_g", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm","Body_length_mm",
             "Svl_length_mm","Generation_length_d", "Longevity_d", "Maturity_d", "AFR_d",
             "Litter_size", "Diel_activity", "Trophic_level","Diet_breadth","Primary_diet","Specialisation","Habitat_breadth_IUCN") %>%
    as.data.frame() %>%
    setNames(., "Original")
  
  Names$FP <- c("Body mass", "Svl length", "Forearm length", "Head length", "Body length", "Svl length",
                "Generation length", "Longevity", "Sexual maturity age", "Age 1st reproduction",
                "Litter/clutch size","Diel activity", "Trophic level","Diet breadth","Primary diet", "Specialisation", "Habitat breadth")

  
  Completeness1 <- apply(TraitData1[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness1 <- as.data.frame(Completeness1/nrow(TraitData1)*100)
  colnames(Completeness1) <- "Completeness"
  
  Completeness2 <- apply(TraitData2[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness2 <- as.data.frame(Completeness2/nrow(TraitData2)*100)
  colnames(Completeness2) <- "Completeness"
  
  Completeness2 <- Completeness2[order(Completeness2, decreasing=FALSE), , drop=FALSE]
  X <- row.names(Completeness2)
  Completeness1 <- Completeness1[X,] %>% as.data.frame()
  row.names(Completeness1) <- X
  colnames(Completeness1) <- "Completeness"
  
  BarplotArg2 <- vector("character")
  for (i in rownames(Completeness2)){
    BarplotArg2 <- c(BarplotArg2, as.character(as.factor(Names$FP[Names$Original==i])))
  }
  
  barplot(Completeness2$Completeness, horiz = TRUE, 
          xlim = c(0,100), names.arg=BarplotArg2, las=1, col="#00BFC4")
  
  barplot(Completeness1$Completeness, horiz=TRUE, add=T, las=1, col="#F8766D")
  title(main=Main, adj=0)
  
  abline(v=100, lty="dotted")
  abline(v=50, lty="dotted")
  
}



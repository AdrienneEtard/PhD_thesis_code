# # # # Plot trait completeness
library(dplyr)
library(grid)
library(gridBase)
library(moments)
library(RColorBrewer)
library(patchwork)
library(ggplot2)
display.brewer.all(colorblindFriendly = TRUE)


# # # #


# function to plot trait completeness
Plot_completeness <- function(traits, FontSize, BW, Hist) {
  
  traits[[1]]$Class <- "Amphibians"
  traits[[2]]$Class <- "Birds"
  traits[[3]]$Class <- "Mammals"
  traits[[4]]$Class <- "Reptiles"
  
  traits <- data.table::rbindlist(traits)
  
  suppressWarnings(require(ggplot2))
  
  GGPoptions <- theme_classic()+ theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=FontSize, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=FontSize), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=FontSize),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=FontSize)) 
  
  
  if(!Hist){
    
    Med <- traits %>% 
      dplyr::group_by(Class) %>%
      dplyr::summarise(x=median(completeness, na.rm=TRUE))
    print(paste("Median trait completeness", Med))
    
    Mean <- traits %>% 
      dplyr::group_by(Class) %>%
      dplyr::summarise(x=mean(completeness, na.rm=TRUE))
    print(paste("Mean trait completeness", Mean))
    
    p <- ggplot(traits, aes(completeness, col=Class, fill=Class)) +
      geom_density(alpha=0.5, adjust=BW) + GGPoptions +
      geom_vline(data=Mean, aes(xintercept=x, col=Class), linetype="solid", alpha=0.7) +
      geom_vline(data=Med, aes(xintercept=x, col=Class), linetype="dashed", alpha=0.7) +
      xlab("Completeness (%)") + ylab("Density of species")
    
    
  } 
  else {
    p <- ggplot(traits, aes(completeness, col=Class, fill=Class)) +
      geom_histogram(alpha=0.5, bins=BW, position = "dodge") + GGPoptions +
      xlab("Completeness (%)") + ylab("Density of species")
  }
  
  return(p)
  
}

# functions to plot trait coverage

## Plotting trait coverage, before VS after taxonomic corrections
Plot_coverage <- function(traits){
  #browser()
  
  traits <- traits %>%
    dplyr::select(-completeness, -Family, -Best_guess_binomial, -Genus, -Order)
  
  Names <- colnames(traits) %>%
    as.data.frame() %>%
    setNames(., "Original")
  
  Names$FP <- c("Body size", "Life span", "Litter/clutch size", "Trophic level", "Diel activity", "Habitat breadth", "Specialisation")
  
  Coverage <- apply(traits, 2,  function(y) sum(!is.na(y))) 
  Coverage <- as.data.frame(Coverage/nrow(traits)*100)
  colnames(Coverage) <- "Completeness"
  
  Coverage <- Coverage[order(Coverage, decreasing=FALSE), , drop=FALSE]
  
  Names_plot <- as.data.frame(row.names(Coverage))
  colnames(Names_plot) <- "Or"
  Names$Original <- as.character(Names$Original)
  Names_plot$Or <- as.character(Names_plot$Or)
  for (i in 1:nrow(Names_plot)) {Names_plot$TP[i] <- Names$FP[Names$Original==Names_plot$Or[i]]}
  
  
  barplot(Coverage$Completeness, horiz = TRUE, 
          xlim = c(0,100), las=1, col="#00BFC4",  names.arg = Names_plot$T, cex.names=1.3)
  
  abline(v=100, lty="dotted")
  abline(v=50, lty="dotted")
  
  x <- mean(Coverage$Completeness)
  print(paste("Mean trait coverage across all traits:", x, "%"))
  
  return(Coverage)
  
}
  
  # # # # 

traits.completeness <- readRDS("../../Results/Traits_to_map/traits_completeness_V2.rds")
  
lapply(traits.completeness, FUN = function(x){return(length(unique(x$Best_guess_binomial)))}) %>%  unlist
  
# mean trait completeness across species
lapply(traits.completeness, function(x) {return(mean(x$completeness, na.rm=TRUE))}) %>%  unlist

# species with 0% completeness

# amphibians
x <- traits.completeness$Amphibians %>% 
  dplyr::filter(completeness==0) %>% 
  nrow()
x/nrow(traits.completeness$Amphibians)*100

## example of 0% completeness: Rhinella centralis
Amphibians <- read.csv("../../Data/Trait_data/Amphibians.csv")
Amphibians  %>% filter(is.na(Maturity_d)) %>%  filter(!is.na(Max_longevity_d))
y <- Amphibians  %>% filter(is.na(Body_length_mm)) %>%  filter(!is.na(Body_mass_g))
traits.completeness$Amphibians %>%  filter(Best_guess_binomial %in% y$Best_guess_binomial)

TraitsA <- c("Body_mass_g",
            "Max_longevity_d",
            "Litter_size", 
            "Habitat_breadth_IUCN",
            "Specialisation",
            "Trophic_level",
            "Diel_activity",
            "Maturity_d", 
            "Body_length_mm")
Completeness0 <- function(TraitDF, Traits) {
  # completeness
  TraitDF$Percent <- apply(TraitDF[,Traits], 1, function(y) sum(!is.na(y)))
  TraitDF$Percent <- TraitDF$Percent / length(Traits) * 100 
  
  Completeness_0 <- TraitDF$Best_guess_binomial[TraitDF$Percent==0]
  print(length(Completeness_0))
  return(Completeness_0)
}

Completeness0(Amphibians, TraitsA)

y <- subset(Amphibians, Amphibians$Best_guess_binomial=="Rhinella centralis")


# reptiles
x <- traits.completeness$Reptiles %>% 
  dplyr::filter(completeness==0) %>% 
  nrow()
x/nrow(traits.completeness$Reptiles)*100 
# birds
x <- traits.completeness$Birds %>% 
  dplyr::filter(completeness==0) %>% 
  nrow()
x/nrow(traits.completeness$Birds)*100
# mammals
x <- traits.completeness$Mammals %>% 
  dplyr::filter(completeness==0) %>% 
  nrow()
x/nrow(traits.completeness$Mammals)*100

# NB: mean trait coverage across species = mean trait completeness across traits = number of NAs across species and traits / numbers of species * number of traits (# measurements).

# look at skewness
skewness(traits.completeness$Amphibians$completeness)
skewness(traits.completeness$Reptiles$completeness)
skewness(traits.completeness$Mammals$completeness)
skewness(traits.completeness$Birds$completeness)

# quantiles
quantile(traits.completeness$Amphibians$completeness)
quantile(traits.completeness$Reptiles$completeness)

quantile(traits.completeness$Mammals$completeness)
quantile(traits.completeness$Birds$completeness)

# 80%-100% range completeness
nrow(traits.completeness$Mammals[traits.completeness$Mammals$completeness>=80,])  / nrow(traits.completeness$Mammals) *100
nrow(traits.completeness$Birds[traits.completeness$Birds$completeness>=80,])  / nrow(traits.completeness$Birds) *100

# 0%-40% range completeness
nrow(traits.completeness$Amphibians[traits.completeness$Amphibians$completeness<=50,])  / nrow(traits.completeness$Amphibians) *100
nrow(traits.completeness$Reptiles[traits.completeness$Reptiles$completeness<=50,])  / nrow(traits.completeness$Reptiles) *100




## Coverage for one dataset
Plot_coverage_V2 <- function(traits, FontSize){

  
  GGPoptions <- theme_classic()+ theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=FontSize, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=FontSize), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=FontSize),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=FontSize)) 
  
  Get_cov <- function(df, Class) {
    
    df <- df %>%
      dplyr::select(-completeness, -Family, -Best_guess_binomial, -Genus, -Order)
    
    Coverage <- colnames(df) %>%
      as.data.frame() %>%
      setNames(., "Original")
    
    Coverage$FP <- c("Body size", "Life span", "Litter/clutch size", "Trophic level", "Diel activity", "Habitat breadth", "Use of artificial habitats")
    
    Coverage$Cov <- apply(df, 2,  function(y) sum(!is.na(y))) 
    Coverage$Cov <- Coverage$Cov/nrow(df)*100
    Coverage$Class <- Class
    
    return(Coverage)
  }
  
  MammalCov <- Get_cov(traits[["Mammals"]], "Mammals")
  BirdCov <- Get_cov(traits[["Birds"]], "Birds")
  AmphibianCov <- Get_cov(traits[["Amphibians"]], "Amphibians")
  ReptileCov <- Get_cov(traits[["Reptiles"]], "Reptiles")
  
  Coverage <- rbind(MammalCov, BirdCov, AmphibianCov, ReptileCov)
  
  # order by highest to lowest mean coverage
  Order <- Coverage %>%  group_by(FP) %>% 
    summarise(Mean=mean(Cov)) %>% 
    as.data.frame()
  Order <- Order[order(Order$Mean, decreasing = FALSE),]
  Order <- Order$FP
  
  Plotdata <- Coverage %>%  group_by(FP) %>% mutate(position=rank(Cov))
  
  
  browser()
  # plot for all classes together
  p <- ggplot(Plotdata, aes(FP, Cov, fill=Class, group=position)) + 
    geom_bar(stat="identity", position="dodge", width=0.7) +
    GGPoptions +  
    scale_x_discrete(limits=Order)+
    coord_flip() + ylab("Coverage (%)") + xlab("") +
    scale_fill_viridis_d() + 
    geom_hline(yintercept=50, lty="dashed")
  
  return(p)
}

## Function for KWallis test
KTest <- function(TraitsM, TraitsB, TraitsR, TraitsA, Pairwise){
  
  ResultsM <- TraitsM %>%
    mutate(Class="Mammalia")
  
  ResultsR <- TraitsR %>%
    mutate(Class="Reptiles")
  
  ResultsA <- TraitsA %>%
    mutate(Class="Amphibians")
  
  ResultsB <- TraitsB %>%
    mutate(Class="Birds")
  
  Results <- rbind(ResultsB, ResultsA, ResultsR, ResultsM)
  
  Results$Class <- as.factor(Results$Class)
  
  if(!Pairwise) {
    return(kruskal.test(Results$completeness, Results$Class))
  }
  
  else{
    return(pairwise.wilcox.test(Results$completeness, Results$Class,
                                p.adjust.method = "BH"))
  }
  
}

traits.completeness <- readRDS("../../Results/Traits_to_map/traits_completeness_V2.rds")


# plotting MANUSCRIPT FIGURE

p1 <- Plot_completeness(traits.completeness, 13, 4, FALSE) + xlab("Completeness (%)") + 
  theme(plot.title = element_text(size=13, face="bold")) + scale_fill_viridis_d() + scale_color_viridis_d()

p2 <- Plot_coverage_V2(traits.completeness, 13)

p2$data %>%  group_by(Class) %>%  summarise(Mean=mean(Cov), Median=median(Cov))

p <- ggpubr::ggarrange(p2 + ggtitle("(a)"),p1+ ggtitle("(b)"), common.legend = TRUE, widths = c(0.54, 0.46), legend="bottom")
ggsave(p, filename="../../Results/Coverage_completeness/Figure1_revised.pdf", width = 10, height = 4)

## Coverage for PREDICTS -- ms functional diversity
Regional.stdFD <- read.csv("../../../1.2.Trait_imputations_FD_using_Global_Gaps_data/Results/dbFD_indices_for_analysis/dbFD_region_std.csv")
Studies <- unique(Regional.stdFD$SS)
Predicts <- readRDS("../../../1.2.Trait_imputations_FD_using_Global_Gaps_data/Data/PredictsVertebrates.rds") %>% 
  dplyr::filter(SS %in% Studies)
Species <- unique(Predicts[, c("Class", "Best_guess_binomial")])
Species %>%  group_by(Class) %>% summarise(C=n())

traits.completeness.Predicts <- lapply(traits.completeness, function(x){
  x <- subset(x, Best_guess_binomial %in% Species$Best_guess_binomial)
  return(x)})

p2 <- Plot_coverage_V2(traits.completeness.Predicts, 13)


# Pairwise Wilcoxon test 
KTest(TraitsM=traits.completeness$Mammals, 
      TraitsR=traits.completeness$Reptiles, 
      TraitsB=traits.completeness$Birds, 
      TraitsA=traits.completeness$Amphibians, 
      Pairwise = TRUE)

# Kruskal wallis (not pairwise)
KTest(TraitsM=traits.completeness$Mammals, 
      TraitsR=traits.completeness$Reptiles, 
      TraitsB=traits.completeness$Birds, 
      TraitsA=traits.completeness$Amphibians, 
      Pairwise = FALSE)


# # # # FIGURE 1 WITH UNCORRECTED TAXONOMY

UncorTraits <- readRDS("../../Results/Traits_to_map/traits_completeness_V2_UNCORRECTED_TAXONOMY.rds")

# # # Trait completeness

p1un <- Plot_completeness(UncorTraits, 13, 4, FALSE) + xlab("Completeness (%)") + 
  theme(plot.title = element_text(size=13, face="bold")) + scale_fill_viridis_d() + scale_color_viridis_d()

p2un <- Plot_coverage_V2(UncorTraits, 13)

pun <- ggpubr::ggarrange(p2un + ggtitle("(a)"),p1un+ ggtitle("(b)"), common.legend = TRUE, widths = c(0.54, 0.46), legend="bottom")
ggsave(pun, filename="../../Results/Coverage_completeness/Figure1_revised_UNCORRECTED_TAXONOMY.pdf", width = 10, height = 4)


# Pairwise Wilcoxon test 
KTest(TraitsM=UncorTraits$Mammals, 
      TraitsR=UncorTraits$Reptiles, 
      TraitsB=UncorTraits$Birds, 
      TraitsA=UncorTraits$Amphibians, 
      Pairwise = TRUE)

# Kruskal wallis (not pairwise)
KTest(TraitsM=UncorTraits$Mammals, 
      TraitsR=UncorTraits$Reptiles, 
      TraitsB=UncorTraits$Birds, 
      TraitsA=UncorTraits$Amphibians, 
      Pairwise = FALSE)

# 100-80% range in completeness
y <-  p1$data

y %>% 
  filter(Class=="Mammals") %>% 
  filter(completeness<=50) %>% 
  nrow()

361/nrow(y[y$Class=="Mammals",])*100

 ## Comparing coverage before VS after taxonomic corrections

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

# # # # # # # # # #  # # # # # # # # # # # # # #

Plot.Delta.Cov <- function(TraitData1, TraitData2, Traits_name, Main ){#Order=c("By1", "By2")) {

  
  Names <- c("Body_size", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm","Body_length_mm",
             "Svl_length_mm","Life_span_proxy", "Longevity_d", "Maturity_d", "AFR_d",
             "Litter_size", "Diel_activity", "Trophic_level","Diet_breadth","Primary_diet","Specialisation","Habitat_breadth_IUCN") %>%
    as.data.frame() %>%
    setNames(., "Original")
  
  Names$FP <- c("Body size", "Svl length", "Forearm length", "Head length", "Body length", "Svl length",
                "Life span proxy", "Longevity", "Sexual maturity age", "Age 1st reproduction",
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

pdf(file="../../Results/Coverage_completeness/DeltaCoverage.pdf", width=8, height=7, family="Times", pointsize=15)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(7,2,1,1))
par(mfrow=c(2,2))

Plot.Delta.Cov(AmphibiansU, traits.completeness$Amphibians, Traits_name = c("Body_size",
                                                                        "Life_span_proxy",
                                                                        "Litter_size",
                                                                        "Trophic_level",
                                                                        "Diel_activity",
                                                                        "Habitat_breadth_IUCN",
                                                                        "Specialisation"), Main="Amphibians")


Plot.Delta.Cov(BirdsU, traits.completeness$Birds, Traits_name = c("Body_size",
                                                                            "Life_span_proxy",
                                                                            "Litter_size",
                                                                            "Trophic_level",
                                                                            "Diel_activity",
                                                                            "Habitat_breadth_IUCN",
                                                                            "Specialisation"), Main="Birds")

Plot.Delta.Cov(MammalsU, traits.completeness$Mammals, Traits_name = c("Body_size",
                                                                  "Life_span_proxy",
                                                                  "Litter_size",
                                                                  "Trophic_level",
                                                                  "Diel_activity",
                                                                  "Habitat_breadth_IUCN",
                                                                  "Specialisation"), Main="Mammals")

Plot.Delta.Cov(ReptilesU, traits.completeness$Reptiles, Traits_name = c("Body_size",
                                                                        "Life_span_proxy",
                                                                        "Litter_size",
                                                                        "Trophic_level",
                                                                        "Diel_activity",
                                                                        "Habitat_breadth_IUCN",
                                                                        "Specialisation"), Main="Reptiles")


par(xpd=NA)
legend(x=-110, y=-6, title="Trait coverage",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("#F8766D", "#00BFC4"),  bty="n")
dev.off()


###### OLDER SCRIPTS


# # # # Trait coverage
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(3,2))
# Plot_coverage(traits.completeness[[3]]); title("(a) Mammals", adj=0)
# Plot_coverage(traits.completeness[[2]]); title("(b) Birds", adj=0)
# Plot_coverage(traits.completeness[[1]]); title("(c) Amphibians", adj=0)
# Plot_coverage(traits.completeness[[4]]); title("(d) Reptiles", adj=0)
# mtext(at=50, line=-13, "Trait coverage (%)", cex=0.8)
# mtext(at=-110, line=-13, "Trait coverage (%)", cex=0.8)
# 
# plot.new()
# vps <- baseViewports()
# pushViewport(vps$figure)
# vp1 <-plotViewport(c(0,2.5,0.5,-10)) ## create new vp with margins, you play with this values 
# print(p,vp = vp1)
# 
# box("outer")
# dev.off()

####

# # # # Trait completeness
# # Plot_completeness(traits.completeness, 12, 5, TRUE)
# # Plot_completeness(traits.completeness, 12, 3, FALSE, "Mean")
# p <- Plot_completeness(traits.completeness, 12, 4, FALSE) + xlab("Completeness (%)")
# ggsave(p, filename="../../Results/Coverage_completeness/Completeness.pdf", height=3, width=6)
# ggsave(p, filename="../../Results/Coverage_completeness/Completeness.png", height=3, width=6)
# 
#      
# # # # Trait coverage
# 
# pdf(file="../../Results/Coverage_completeness/Coverage.pdf", width=8, height=5.5, family="Times", pointsize=15)
# #png(filename="../../Results/Coverage_completeness/Coverage.png", width=8, height=5.5, units="in", family="Times", pointsize=15, res=400)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(2,2))
# Plot_coverage(traits.completeness[[3]]); title("(a) Mammals", adj=0)
# Plot_coverage(traits.completeness[[2]]); title("(b) Birds", adj=0)
# Plot_coverage(traits.completeness[[1]]); title("(c) Amphibians", adj=0)
# Plot_coverage(traits.completeness[[4]]); title("(d) Reptiles", adj=0)
# mtext(at=50, line=-10, "Coverage (%)", cex=0.8)
# mtext(at=-155, line=-10, "Coverage (%)", cex=0.8)
# box("outer")
# dev.off()



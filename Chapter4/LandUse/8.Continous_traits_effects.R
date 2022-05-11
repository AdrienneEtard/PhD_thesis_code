## continous traits: effects of land use & use intensity on slopes of relationships between traits and occurrence probability

##############################################################################################################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(StatisticalModels)
library(scales)
library(viridis)
library(dplyr)
library(lme4)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

## loading models
Birds <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Birds_Diet.rds")
Mammals <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Mammals_Diet.rds")
Amphibians <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Amphibians_Diet.rds") ## no primary diet
Reptiles <- readRDS("../Results/GLMER_Models/Class_specific/Use_intensity/Reptiles_Diet.rds")

## function to get slopes by resampling in the models coefficients

Get_slopes_resampling <- function(Model, PD, Trait, Class){
  
  
  ## create new data 
  if(PD){
  Newdata <- expand.grid(LandUseGrouped=levels(Model@frame$LandUseGrouped ), 
                         Use_intensity=levels(Model@frame$Use_intensity), 
                         Specialisation =levels(Model@frame$Specialisation),  
                         Diel_activity =levels(Model@frame$Diel_activity),  
                         Primary_diet =levels(Model@frame$Primary_diet),  
                         log10_Body_mass_g=mean(Model@frame$log10_Body_mass_g),  
                         log10_Litter_size=mean(Model@frame$log10_Litter_size),  
                         log10_Lifespan_proxy=mean(Model@frame$ log10_Lifespan_proxy),  
                         log10_Range_area=mean(Model@frame$log10_Range_area),  
                         sqrt_Habitat_breadth_IUCN=mean(Model@frame$sqrt_Habitat_breadth_IUCN),
                         sqrt_Diet_breadth =mean(Model@frame$ sqrt_Diet_breadth),
                         Occurrence=1)
  } else{
    Newdata <- expand.grid(LandUseGrouped=levels(Model@frame$LandUseGrouped ), 
                           Use_intensity=levels(Model@frame$Use_intensity), 
                           Specialisation =levels(Model@frame$Specialisation),  
                           Diel_activity =levels(Model@frame$Diel_activity),  
                           #Primary_diet =levels(Model@frame$Primary_diet),  
                           log10_Body_mass_g=mean(Model@frame$log10_Body_mass_g),  
                           log10_Litter_size=mean(Model@frame$log10_Litter_size),  
                           log10_Lifespan_proxy=mean(Model@frame$ log10_Lifespan_proxy),  
                           log10_Range_area=mean(Model@frame$log10_Range_area),  
                           sqrt_Habitat_breadth_IUCN=mean(Model@frame$sqrt_Habitat_breadth_IUCN),
                           sqrt_Diet_breadth =mean(Model@frame$ sqrt_Diet_breadth),
                           Occurrence=1)
  }
  

  # subsetting model's coefficients and variance-covariance matrix for the trait
  Fixed_effects <- fixef(Model) 
  Fixed_effects_slopes <- Fixed_effects[grepl(Trait, names(Fixed_effects))]
  
  V_Cov <- vcov(Model)
  # subset relevant rows
  V_Cov_sub <- V_Cov[grepl(Trait, rownames(V_Cov)),] 
  # subset relevant columns
  Col_retain <- which(grepl(Trait, colnames(V_Cov_sub)))
  
  # sanity check
  if(length(Col_retain)!=nrow(V_Cov_sub)){
    print("Problem!")
    stop()
  } else{
    print("Sanity check ok")
  }
    
  V_Cov_sub <- V_Cov_sub[, Col_retain]
  
  # sanity check n2
  if(nrow(V_Cov_sub)!=ncol(V_Cov_sub)){
    print("Problem!")
          stop()
  } else{
    print("Sanity check number 2 ok")
  }

  
  print("Resampling in progress...")
  ## getting slopes and 95% CI by re sampling
  preds <- sapply(X=1:5000, FUN=function(i){
    coefs <- mvrnorm(n=1, mu=Fixed_effects_slopes, Sigma=V_Cov_sub)
    mm <- model.matrix(terms(Model), Newdata)
    ## drop coefs that shouldn't be estimated (the ones that correspond to intercepts)
    #print(setdiff(colnames(mm), names(coefs)))
    to_drop <- setdiff(colnames(mm), names(coefs))
    if(length(to_drop)!=0){
      mm <- as.data.frame(mm)
      mm <- mm[, -which(colnames(mm) %in% to_drop)]
      mm <- as.matrix(mm)
    }
    ## get slope estimates
    y <- mm %*% coefs
    return(y)
  })
  
  preds_Slopes <- data.frame(Median=apply(X=preds, MARGIN=1, FUN=median),
                             Upper=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.975),
                             Lower=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.025))
  
  preds_Slopes <- cbind(preds_Slopes, Newdata)
  
  preds_Slopes$Class <- Class
  preds_Slopes$Trait <- Trait

  return(preds_Slopes)
  
}

Run_on_traits <- function(Model, PD, Traits, Class) {
  
  ResList <- list()
  
  for (t in Traits){
    cat("Running for", t, "\n", "\n")
    ResList[[t]] <- Get_slopes_resampling(Model, PD, Trait = t, Class)
  }
  
  return(data.table::rbindlist(ResList))
  
}

##############################################################################################################################

## CONTINUOUS TRAITS:

Cont_traits <- c("log10_Body_mass_g",
                 "log10_Litter_size",
                 "log10_Lifespan_proxy",
                 "log10_Range_area",
                 "sqrt_Habitat_breadth_IUCN",
                 "sqrt_Diet_breadth")


## SLOPES for birds
SlopesBirds <- Run_on_traits(Model = Birds, PD = TRUE, Traits = Cont_traits, Class="Birds")
SlopesBirds <- SlopesBirds %>% 
  dplyr::select(-Primary_diet)

## SLOPES for mammals
SlopesMammals <- Run_on_traits(Model = Mammals, PD = TRUE, Traits = Cont_traits, Class="Mammals")
SlopesMammals <- SlopesMammals %>% 
  dplyr::select(-Primary_diet)

## SLOPES for amphibians
SlopesAmphibians <- Run_on_traits(Model = Amphibians, PD = FALSE, Traits = Cont_traits, Class="Amphibians")

## SLOPES for reptiles
SlopesReptiles <- Run_on_traits(Model = Reptiles, PD = FALSE, Traits = Cont_traits, Class="Reptiles")

## all slopes
Slopes <- rbind(SlopesAmphibians, SlopesBirds, 
                SlopesMammals, SlopesReptiles)

## setting CI for reference to NA ## minimally used PV
Slopes$Lower[Slopes$LandUseGrouped=="Primary vegetation" & Slopes$Use_intensity=="Minimal use"] <- NA
Slopes$Upper[Slopes$LandUseGrouped=="Primary vegetation" & Slopes$Use_intensity=="Minimal use"] <- NA


## removing certain effects before plotting because of massive error bars

## removing all traits, amphibians, urban -- all null effects
Slopes$Median[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Lifespan_proxy" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Range_area" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Habitat_breadth_IUCN" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Body_mass_g" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

Slopes$Median[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Lower[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA
Slopes$Upper[Slopes$Trait=="log10_Litter_size" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Amphibians"] <- NA

## removing diet breadth urban mammals
Slopes$Median[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Mammals"] <- NA
Slopes$Lower[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Mammals"] <- NA
Slopes$Upper[Slopes$Trait=="sqrt_Diet_breadth" & Slopes$LandUseGrouped=="Urban" & Slopes$Class=="Mammals"] <- NA


############################################# Plotting slopes


Plot_slopes_fun_Allfacet <- function(SlopesDF) {
  
  SlopesDF$H_line <- NA
  
  for(comb in levels(SlopesDF$Comb)){
    
    ref <- SlopesDF$Median[SlopesDF$Comb==comb & 
                             SlopesDF$LandUseGrouped=="Primary vegetation" &
                             SlopesDF$Use_intensity=="Minimal use"]
    SlopesDF$H_line[SlopesDF$Comb==comb] <- ref
  }
  
  # ## assigning lines of reference
  # Vals <- SlopesDF %>% 
  #   dplyr::group_by(Comb) %>% 
  #   dplyr::filter(LandUseGrouped=="Primary vegetation") %>% 
  #   dplyr::filter(Use_intensity=="Minimal use") %>% 
  #   dplyr::distinct(Median) %>% 
  #   as.data.frame()
  # 
  # #Vals <- Vals[order(match(Vals$Comb, levels(SlopesDF$Comb))),]
  # 
  # Vals <- Vals$Median
  # 
  # Mult <- SlopesDF %>%  
  #   group_by(Comb) %>%
  #   summarise(C=n())
  # 
  # Mult <- Mult$C
  # 
  # SlopesDF$H_line <- NA
  # 
  # Seq <- 1
  # for(i in 1:length(Mult)){
  #   a <- Seq
  #   b <- Seq + Mult[i]-1
  #   Seq <- Seq + Mult[i]
  #   SlopesDF$H_line[a:b] <- Vals[i]
  #   
  # }
  
  ggplot(SlopesDF,
         aes(LandUseGrouped, Median, ymin = Lower, ymax = Upper, group=Use_intensity, shape=Use_intensity, col=Use_intensity)) +
    geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    geom_rect(xmin=8.5, xmax=9.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
    ylab("") + xlab("") +
    geom_hline(aes(yintercept=H_line, group=Use_intensity),linetype="dashed") +
    geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
    geom_point(size=2, position=position_dodge(width = 0.7)) +
    GGPoptions +
    ggtitle("") +
    theme(panel.spacing = unit(0, "lines")) +
    theme( strip.text.x = element_text(size = 12, face = "bold"),
           strip.text.y = element_text(size = 12, face = "bold")) +
    scale_colour_manual(values=c("black", "#FF7F0E", "#D55E00"), name="Use intensity")+
    theme(legend.position = "right") +
    scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
    theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
    ylab("Slope estimate (+/- 95% CI)") +
    theme(legend.position = "top") + facet_wrap(~Comb, scale="free_y", ncol=4)+
    scale_x_discrete(labels=c("Primary", "Secondary","Plantation", "Agricultural", "Urban")) +
    theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
          axis.text.x = element_text(size = 10, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")))
  
}

Slopes$Comb <- paste(Slopes$Class, Slopes$Trait)
Slopes$Comb <- as.factor(Slopes$Comb)
Slopes$Comb %>%  levels()
Names <- expand.grid(c("BM", "LP","LCS","RA", "DB", "HB"), c("Amphibians", "Birds", "Mammals", "Reptiles"))
Names <- paste0(Names$Var2, ": ",Names$Var1)
Levels  <- expand.grid(c("Amphibians", "Birds", "Mammals", "Reptiles"), c("BM", "LP","LCS","RA", "DB", "HB"))
Levels <- paste0(Levels$Var1, ": ",Levels$Var2)
Slopes$Comb <- factor(Slopes$Comb, labels=Names)
Slopes$Comb <- factor(Slopes$Comb, levels=Levels)

p <- Plot_slopes_fun_Allfacet(Slopes)
p

ggsave(p, 
       filename="E:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Figures/Slopes_LU.pdf",
       width = 8, height = 10)



















# Plot_slopes_fun_Classfacet <- function(SlopesDF, trait) {
#   
#   SlopesDF <- subset(SlopesDF,SlopesDF$Trait==trait)
#   
#   ## assigning lines of reference
#   Ref_A <- SlopesDF$Median[SlopesDF$LandUse=="Primary vegetation" & 
#                              SlopesDF$Use_intensity=="Minimal use" & SlopesDF$Class %in% c("Amphibians")]
#   # Ref_A2 <- SlopesDF$Median[SlopesDF$LandUse=="Primary vegetation" & 
#   #                            SlopesDF$Use_intensity=="Light use" & SlopesDF$Class %in% c("Amphibians")]
#   # Ref_A3 <- SlopesDF$Median[SlopesDF$LandUse=="Primary vegetation" & 
#   #                            SlopesDF$Use_intensity=="Intense use" & SlopesDF$Class %in% c("Amphibians")]
#   Ref_B <- SlopesDF$Median[SlopesDF$LandUse=="Primary vegetation" & 
#                              SlopesDF$Use_intensity=="Minimal use" & SlopesDF$Class %in% c("Birds")]
#   Ref_M <- SlopesDF$Median[SlopesDF$LandUse=="Primary vegetation" & 
#                              SlopesDF$Use_intensity=="Minimal use" & SlopesDF$Class %in% c("Mammals")]
#   Ref_R <- SlopesDF$Median[SlopesDF$LandUse=="Primary vegetation" & 
#                              SlopesDF$Use_intensity=="Minimal use" & SlopesDF$Class %in% c("Reptiles")]
#   
#   SlopesDF$H_line <- c(rep(Ref_A, 5*3),
#                         rep(Ref_B, 5*3),
#                         rep(Ref_M, 5*3),
#                         rep(Ref_R, 4*3))
# 
#   #browser()
#   
#     ggplot(SlopesDF,
#          aes(LandUseGrouped, Median, ymin = Lower, ymax = Upper, group=Use_intensity, shape=Use_intensity, col=Use_intensity)) +
#     geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
#     geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
#     geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
#     geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
#     geom_rect(xmin=8.5, xmax=9.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
#     ylab("") + xlab("") +
#     geom_hline(aes(yintercept=H_line, group=Use_intensity),linetype="dashed") +
#     geom_errorbar(width=.2, size=0.6, position=position_dodge(width = 0.7), stat="identity") + #aes(linetype=Significance)
#     geom_point(size=2, position=position_dodge(width = 0.7)) +
#     GGPoptions +
#     ggtitle("") +
#     theme(panel.spacing = unit(0, "lines")) +
#     theme( strip.text.x = element_text(size = 12, face = "bold"),
#            strip.text.y = element_text(size = 12, face = "bold")) +
#     scale_colour_manual(values=c("black", "red", "coral"), name="Use intensity")+
#     #scale_colour_viridis_d() +
#     theme(legend.position = "right") +
#     scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
#     theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
#     ylab("Slope estimate (+/- 95% CI)") +
#     theme(legend.position = "bottom") + facet_wrap(~Class, scale="free_y", nrow=6) + 
#     ggtitle(trait)
#   
# }


## plotting for all continuous traits

# Slopes_continuous_traits <- 
# 
# ggarrange(
# Plot_slopes_fun_Classfacet(Slopes, trait="log10_Body_mass_g")+ ggtitle("(a) Body mass") + theme(axis.text.x = element_blank()),
# Plot_slopes_fun_Classfacet(Slopes, trait="log10_Litter_size")+ ggtitle("(b) Litter/clutch size") + theme(axis.text.x = element_blank()) + theme(axis.title.y=element_blank()),
# Plot_slopes_fun_Classfacet(Slopes, trait="log10_Lifespan_proxy")+ ggtitle("(c) Lifespan proxy") + theme(axis.text.x = element_blank()) + theme(axis.title.y=element_blank()),
# Plot_slopes_fun_Classfacet(Slopes, trait="sqrt_Diet_breadth")+ ggtitle("(d) Diet breadth") +
#   scale_x_discrete(labels=c("Primary", 
#                             "Secondary",
#                             "Plantation", 
#                             "Agricultural",
#                             "Urban")), 
# Plot_slopes_fun_Classfacet(Slopes, trait="sqrt_Habitat_breadth_IUCN")+ ggtitle("(e) Habitat breadth") + theme(axis.title.y=element_blank()) +
#   scale_x_discrete(labels=c("Primary", 
#                             "Secondary",
#                             "Plantation", 
#                             "Agricultural",
#                             "Urban")),
# Plot_slopes_fun_Classfacet(Slopes, trait="log10_Range_area")+ ggtitle("(f) Range area") + theme(axis.title.y=element_blank())+
#   scale_x_discrete(labels=c("Primary", 
#                             "Secondary",
#                             "Plantation", 
#                             "Agricultural",
#                             "Urban")),
# nrow=2, ncol=3, common.legend = TRUE, heights = c(0.40, 0.60))
# # 
# 
# ggsave(Slopes_continuous_traits, filename="D:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Figures/Slopes_LU.pdf",
#        width = 8, height = 12)


## different facetting



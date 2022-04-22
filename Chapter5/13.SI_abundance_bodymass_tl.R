## model special issue -- abundance
## checking changes in size-spectrum

library(dplyr)
library(lme4)
library(ggplot2)
library(StatisticalModels)
library(performance)
library(ggpubr)
library(viridis)

## for plotting effects

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

Predict_effects <- function(Newdata, Model, rescale, Cont_or_Cat, seMultiplier, LU_n) {
  
  # which coefs are not estimated?
  coefs <- mvrnorm(n=1, mu=fixef(Model), Sigma=vcov(Model))
  mm <- model.matrix(terms(Model), Newdata)
  print(setdiff(colnames(mm), names(coefs)))
  
  # predictions, for categories
  if(Cont_or_Cat=="categorical"){
    
    preds <- sapply(X=1:1000, FUN=function(i){
      
      #browser()
      mm <- model.matrix(terms(Model), Newdata)
      to_drop <- setdiff(colnames(mm), names(coefs))
      
      if(length(to_drop!=0)){
        mm <- as.data.frame(mm)
        mm <- mm[, -which(colnames(mm) %in% to_drop)]
        mm <- as.matrix(mm)}
      
      coefs <- mvrnorm(n=1, mu=fixef(Model), Sigma=vcov(Model))
      y <- mm %*% coefs
      
      # backtranforming
      y <- exp(y)
      
      # rescaling
      if(rescale){
        
        # inittialisation
        seq <- 1:LU_n
        y[seq] <- y[seq]/y[seq[1]]*100
        
        # for loop to rescale all values
        for(i in 1:(nrow(Newdata)/LU_n-1)){
          seq <- seq + LU_n
          y[seq] <- y[seq]/y[seq[1]]*100
        }
      }
      
      return(y)
    })
    
    preds <- data.frame(Median=apply(X=preds, MARGIN=1, FUN=median),
                        Upper=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.975),
                        Lower=apply(X=preds, MARGIN=1, FUN=quantile, probs=0.025))
    
    preds <- cbind(preds, Newdata)
    return(preds)
  }
  
  if(Cont_or_Cat=="continuous"){
    
    mm <- model.matrix(terms(Model), Newdata)
    
    # predictions
    preds <- mm %*% fixef(Model)
    
    # backtransform
    Newdata$Estimate <- exp(preds)
    
    # getting the errors around the predicted values
    # variance covariance matrix and matrix multiplication
    VarCoVar <- as.matrix(vcov(Model))
    VarCoVar <- base::tcrossprod(VarCoVar, mm)
    
    # matrix multiplication to get estimates and take diagonal matrix
    pvar <- diag(mm %*% VarCoVar)
    
    # estimate error of predictions and backtransform
    Lower <- preds - seMultiplier*sqrt(pvar)
    Upper <- preds + seMultiplier*sqrt(pvar)
    
    Newdata$Lower <- exp(Lower)
    Newdata$Upper <- exp(Upper)

    return(Newdata)
    
  }
}


Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")


Labels=c("PV", "SV", "PF", "PA", "CR", "UR")

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols)


######################################################################################################################

## loading and checking PREDICTS data
PredictsTraits <- "../Results/PredictsTraits.rds" %>% 
  readRDS() 

levels(PredictsTraits$LandUse)
levels(PredictsTraits$Use_intensity)

## abundance data
AbSub <- PredictsTraits[PredictsTraits$Occurrence==1,] 
AbSub <- AbSub[AbSub$Diversity_metric=="abundance",] 
AbSub$logAb <- log(AbSub$Effort_corrected_measurement+1)  
hist(AbSub$logAb)
hist(log(AbSub$Effort_corrected_measurement))

## body mass
AbSub$logBM <- log(AbSub$Body_mass_g)

## Running abundance model -- grouping secondary vegetation
Start <- Sys.time()
print(Start)
Model <- lme4::lmer(logAb ~
                       LandUse +
                       Use_intensity +
                       Trophic_level +
                       logBM +
                       LandUse:Use_intensity +
                       LandUse:Trophic_level +
                       LandUse:logBM +
                       Use_intensity:Trophic_level +
                       Use_intensity:logBM +
                       Trophic_level:logBM +
                       LandUse:Trophic_level:logBM +
                       Use_intensity:Trophic_level:logBM +
                       (1|SS) +
                       (1|SSBS) +
                       (1|Best_guess_binomial),
                     data=AbSub)

End <- Sys.time()
print(Start-End) ## 12 seconds...super quick!


######################################################################################################################

## plotting effects: slopes of the relationship
## getting slopes by resampling in the models coefficients

Newdata <- expand.grid(LandUse=levels(Model@frame$LandUse), 
                       Use_intensity=levels(Model@frame$Use_intensity), 
                       Trophic_level=levels(Model@frame$Trophic_level),  
                       logBM=1,  
                       logAb=1)

# subsetting model's coefficients and variance-covariance matrix
Fixed_effects <- fixef(Model) 
Fixed_effects_slopes <- Fixed_effects[grepl("logBM", names(Fixed_effects))]

V_Cov <- vcov(Model)
# subset relevant rows
V_Cov_sub <- V_Cov[grepl("logBM", rownames(V_Cov)),] 
# subset relevant columns
Col_retain <- which(grepl("logBM", colnames(V_Cov_sub)))
length(Col_retain)==nrow(V_Cov_sub)
V_Cov_sub <- V_Cov_sub[, Col_retain]
nrow(V_Cov_sub)==ncol(V_Cov_sub)

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

## plotting slopes
scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[c(2:5,7)]
cbPalette <- c("#000000", cols)
scales::show_col(cbPalette)

# setting levels for land use
preds_Slopes$LandUse <- factor(preds_Slopes$LandUse, 
                               labels = c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban"))

# for adding horizontal line
preds_Slopes$h_line <- rep(preds_Slopes$Median[preds_Slopes$LandUse=="Primary" &
                                                 preds_Slopes$Use_intensity=="Minimal use"], each=18)

# remove CIs for Primary vegetation, minimal use
preds_Slopes$Lower[preds_Slopes$LandUse=="Primary"&
                     preds_Slopes$Use_intensity=="Minimal use"] <- NA

preds_Slopes$Upper[preds_Slopes$LandUse=="Primary"&
                     preds_Slopes$Use_intensity=="Minimal use"] <- NA

# adding in whether the slope differ within land use intensities for each land use
Signif_UI <- function(preds_Slopes, TL) {
  
  sub <- subset(preds_Slopes, Trophic_level==TL)
  
  # split Sub by land use
  SubSplit <- split(x=sub, f=sub$LandUse)
  
  # apply function to each element
  FunToApply <- function(X) {
    ValueMinimal <- X$Median[X$Use_intensity=="Minimal use"]
    X$Significance <- NA
    X$Significance[2] <- !dplyr::between(ValueMinimal, X$Lower[2], X$Upper[2])      
    X$Significance[3] <- !dplyr::between(ValueMinimal, X$Lower[3], X$Upper[3])      
    return(X)
  }
  
  SubSplitRes <- lapply(SubSplit, FunToApply)
  SubSplitRes <- data.table::rbindlist(SubSplitRes)
  
  return(SubSplitRes)
}

preds_SlopesC <- Signif_UI(preds_Slopes, "Carnivore")
preds_SlopesO <- Signif_UI(preds_Slopes, "Omnivore")
preds_SlopesH <- Signif_UI(preds_Slopes, "Herbivore")

preds_Slopes <- rbind(preds_SlopesC, preds_SlopesO, preds_SlopesH)
preds_Slopes$Significance[is.na(preds_Slopes$Significance)] <- TRUE

preds_Slopes$Significance <- factor(preds_Slopes$Significance, levels=c("TRUE", "FALSE"))

Slopes_plot <-
  ggplot(preds_Slopes,
         aes(LandUse, Median, ymin = Lower, ymax = Upper, group=Use_intensity, shape=Use_intensity, col=Significance)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=8.5, xmax=9.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  ylab("") + xlab("") +
  geom_hline(aes(yintercept=h_line), col="#808080") +
  #geom_hline(yintercept=0, col="darkgrey") +
  geom_hline(yintercept=preds_Slopes$Estimate[preds_Slopes$LandUse=="Primary" & preds_Slopes$Use_intensity=="Minimal use"],linetype="dashed") +
  #geom_hline(yintercept=0,linetype="dotted", col="blue") +
  geom_errorbar(width=.2, size=0.6, position=position_dodge(width = 0.7), stat="identity") + #aes(linetype=Significance)
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=c("black", "#A0A0A0"),
                      name="Within land uses, is slope \nsignificantly different from \nminimal use-intensity level") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  #guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Slope estimate (+/- 95% CI)")
#scale_linetype_manual(values=c(1,8), name="Within land uses, is slope \nsignificantly different from \nminimal use-intensity level")

Slopes_plot


######################################################################################################################

## plotting continuous effects
Newdata <- expand.grid(LandUse=levels(Model@frame$LandUse), 
                       Use_intensity=levels(Model@frame$Use_intensity), 
                       Trophic_level=levels(Model@frame$Trophic_level),  
                       logBM=seq(from=min(Model@frame$logBM),
                                                to=max(Model@frame$logBM),
                                                length.out=100),  
                       logAb=1)


## predictions
Predictions <- Predict_effects(Newdata = Newdata,
                               Model = Model,
                               rescale = FALSE,
                               Cont_or_Cat = "continuous",
                               seMultiplier = 1.96,
                               LU_n = 6)

Predictions$LandUse <- factor(Predictions$LandUse, 
                              labels = c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban"))



# plot for minimal uses for carnivores and omnivores, and for minimal and light uses for herbivores (in Primary, secondary, plantation)
Predictions_MU <- subset(Predictions, Use_intensity %in% c("Minimal use"))
# Predictions_LU <- subset(Predictions, Use_intensity =="Light use" &
#                            LandUse %in% c("Primary", "Secondary", "Plantation", "Pasture") &
#                            Trophic_level=="Carnivore")
Predictions_IU1 <- subset(Predictions, Use_intensity =="Intense use" &
                           LandUse %in% c("Primary") &
                           Trophic_level=="Omnivore")
Predictions_IU2 <- subset(Predictions, Use_intensity =="Intense use" &
                           LandUse %in% c("Primary", "Secondary", "Plantation", "Pasture", "Cropland") &
                           Trophic_level=="Carnivore")


PredsToplot <- rbind(Predictions_MU, Predictions_IU1, Predictions_IU2)

PredsToplot$Significance_from_PV <- NA
PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Omnivore" &
                                   PredsToplot$LandUse %in% c("Primary","Pasture", "Cropland", "Urban")] <- TRUE

PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Herbivore" &
                                   PredsToplot$LandUse %in% c("Urban")] <- TRUE

PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Carnivore" & 
                                   PredsToplot$LandUse %in% c("Cropland") &
                                   PredsToplot$Use_intensity=="Intense use"] <- TRUE

PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Carnivore" & 
                                   PredsToplot$LandUse %in% c("Plantation") &
                                   PredsToplot$Use_intensity=="Minimal use"] <- TRUE

PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Carnivore" & 
                                   PredsToplot$LandUse %in% c("Pasture","Urban")] <- TRUE

PredsToplot$Significance_from_PV[PredsToplot$Trophic_level=="Carnivore" & 
                                   PredsToplot$LandUse %in% c("Secondary") &
                                   PredsToplot$Use_intensity!="Minimal use"] <- TRUE

PredsToplot$Significance_from_PV[is.na(PredsToplot$Significance_from_PV)] <- FALSE
PredsToplot$Significance_from_PV[PredsToplot$LandUse=="Primary"] <- TRUE
PredsToplot$Significance_from_PV <- factor(PredsToplot$Significance_from_PV, levels=c("TRUE", "FALSE"))


p_continuous_effects <- 
  ggplot(Predictions[Predictions$Use_intensity=="Intense use",], 
         aes(logBM, Estimate, col=Use_intensity, fill=Use_intensity, ymin=Lower, ymax=Upper)) + 
  geom_line(size=0.8) + 
  # scale_fill_manual(values=c("black", "black"), name="Use intensity") +
  # scale_colour_manual(values=c("black", "darkgrey"), name="Use intensity") +
  geom_ribbon(alpha=0.05, col=NA) + 
  facet_grid(LandUse~Trophic_level, scales="free")+
  xlab("Body mass (log g)") + ylab("Abundance") +
  GGPoptions +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  #scale_y_continuous(breaks = seq(0, 0.8, by = 0.25)) +
  scale_linetype_manual(values=c(1,2), name="Within land-use intensity, is \nrelationship significantly \ndifferent from primary")
p_continuous_effects


## Arranging the plot

p_manuscript <- 
  ggarrange(  Slopes_plot +
                ggtitle("(a) Slope of the relationships between log-abundance and log-body mass"), 
              p_continuous_effects +
                ggtitle("(b) Predicted effects of log-body mass on log-abundance"),
              
              nrow=2,
              #common.legend = TRUE,
              legend = "right",
              heights = c(0.50, 0.50))



ggsave(p_manuscript, 
       filename = "../Results/Figures/13.Abundance_BM/Figure_predictions_slopes.pdf", height=9, width=11)

ggsave(p_manuscript, 
       filename = "../Results/Figures/13.Abundance_BM/Figure_predictions_slopes.png", height=9, width=11)

######################################################################################################################

## discretinz effects of body mass


logBM <- unique(Model@frame[, c("Best_guess_binomial", "logBM")])
Low <- quantile(logBM$logBM)["25%"]
Median <- quantile(logBM$logBM)["50%"]
High <- quantile(logBM$logBM)["75%"]

Newdata <- expand.grid(LandUse=levels(Model@frame$LandUse), 
                       logBM=c(Low, Median, High),  
                       Use_intensity=levels(Model@frame$Use_intensity), 
                       Trophic_level=levels(Model@frame$Trophic_level),  
                       logAb=1)

Effects_discrete <- Predict_effects(Newdata = Newdata,
                                    Model = Model,
                                    rescale = TRUE, 
                                    Cont_or_Cat = "categorical", 
                                    seMultiplier = 1.96, 
                                    LU_n = 18)


Effects_discrete$Median <- Effects_discrete$Median - 100
Effects_discrete$Lower <- Effects_discrete$Lower - 100
Effects_discrete$Upper <- Effects_discrete$Upper - 100

Effects_discrete$logBM <- factor(Effects_discrete$logBM, 
                                    labels=c("Low (25% quantile)",
                                             "Median",
                                             "High (75% quantile)"))

Effects_discrete$Signif <- TRUE
Effects_discrete$Signif[Effects_discrete$Trophic_level=="Carnivore" & Effects_discrete$LandUse=="Cropland"] <- FALSE


pAb <- 
  ggplot(Effects_discrete,
         aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=logBM, group=logBM)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="#E0E0E0", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed", col="darkgrey") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  GGPoptions +
  ggtitle("") +
  facet_grid(Use_intensity~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_shape_manual(values=c(16,3,4), name="Body mass") +
  theme(legend.position = "right") +
  scale_colour_manual(values=c("#A0A0A0", "black"), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Abundance\n(% difference from primary)") +
  scale_x_discrete(labels=c("Primary", "Secondary", "Plantation", "Pasture",
                            "Cropland", "Urban"))
  

pAb

ggsave(pAb,
       filename = "../Results/Figures/13.Abundance_BM/Figure_Abundance_difference.pdf", height=6, width=8)

ggsave(pAb,
       filename = "../Results/Figures/13.Abundance_BM/Figure_Abundance_difference.png", height=6, width=8)


## Arranging plots



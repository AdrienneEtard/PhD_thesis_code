## Investigating effects of land use and use intensity on total, assemblage-level tRMR, within trophic levels and diet


##################################################################################################################

library(StatisticalModels)
library(scales)
library(viridis)
library(performance)
library(dplyr)
library(ggpubr)

## for plotting

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")
Labels=c("PV", "SV", "PF", "PA", "CR", "UR")

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[c(2:5,7)]
cbPalette <- c("#000000", cols)
scales::show_col(cbPalette)

Predict_effects <- function(Newdata, Model, rescale, Cont_or_Cat, seMultiplier, LU_n) {
  
  # which coefs are not estimated?
  coefs <- mvrnorm(n=1, mu=fixef(Model), Sigma=vcov(Model))
  mm <- model.matrix(terms(Model), Newdata)
  print(setdiff(names(coefs), colnames(mm)))
  
  # predictions, for categories
  if(Cont_or_Cat=="categorical"){
    
    preds <- sapply(X=1:1000, FUN=function(i){
      
      coefs <- mvrnorm(n=1, mu=fixef(Model), Sigma=vcov(Model))
      mm <- model.matrix(terms(Model), Newdata)
      
      # # drop coefs that couldn't be estimated
      print(setdiff(colnames(mm), names(coefs)))
      to_drop <- setdiff(colnames(mm), names(coefs))
      
      if(length(to_drop)!=0){
        mm <- as.data.frame(mm)
        mm <- mm[, -which(colnames(mm) %in% to_drop)]
        mm <- as.matrix(mm)
      }
      
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
    Newdata$Estimate <- logit2prob(preds)
    
    # getting the errors around the predicted values
    # variance covariance matrix and matrix multiplication
    VarCoVar <- as.matrix(vcov(Model))
    VarCoVar <- base::tcrossprod(VarCoVar, mm)
    
    # matrix multiplication to get estimates and take diagonal matrix
    pvar <- diag(mm %*% VarCoVar)
    
    # estimate error of predictions and backtransform
    Lower <- preds - seMultiplier*sqrt(pvar)
    Upper <- preds + seMultiplier*sqrt(pvar)
    
    Newdata$Lower <- logit2prob(Lower)
    Newdata$Upper <- logit2prob(Upper)
    
    return(Newdata)
    
  }
}


###########################################################################################
###########################################################################################
###########################################################################################


## fitting model to explain total assemblage-level RMR (hypothesis 1)


###########################################################################################
###########################################################################################
###########################################################################################


## load PREDICTS data
PredictsTraits <- "../Results/PredictsTraits.rds" %>% 
  readRDS() 

levels(PredictsTraits$LandUse)
levels(PredictsTraits$Predominant_land_use)
levels(PredictsTraits$Use_intensity)
length(unique(PredictsTraits$SS))

## addind in Diet data
DietData <- readRDS("D:/3.Explanatory_traits/Results/Model_data_allvertebrates.rds")
DietData <- unique(DietData[, c("Best_guess_binomial", "sqrt_Diet_breadth", "Primary_diet")])
DietData$Primary_diet %>%  table()
PredictsTraits <- left_join(PredictsTraits, DietData)

## Occurrence data in PREDICTS
PredictsTraits$Occurrence <- ifelse(PredictsTraits$Measurement > 0, 1, 0)
PredictsTraits$Class <- droplevels(PredictsTraits$Class)

## filter for occurring species, for which Abundance was measured, and remove NA land uses and use intensities
# with grouped secondary vegetation
Predicts_SV <- PredictsTraits %>% 
  filter(Occurrence!=0) %>% 
  filter(Diversity_metric=="abundance") %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity))

# # with all stages for secondary vegetation
PredictsTraits <- PredictsTraits %>%
  filter(Occurrence!=0) %>%
  filter(Diversity_metric=="abundance") %>%
  filter(!is.na(Predominant_land_use), !is.na(Use_intensity))

any(is.na(PredictsTraits$Predominant_land_use))
any(is.na(Predicts_SV$LandUse))

any(is.na(PredictsTraits$Use_intensity))
any(is.na(Predicts_SV$Use_intensity))

## checking sample sizes in sites and studies (bigger when grouped SV)
length(unique(PredictsTraits$SSBS))
length(unique(PredictsTraits$SS))

length(unique(Predicts_SV$SSBS))
length(unique(Predicts_SV$SS))

Predicts_SV$Best_guess_binomial %>%  unique %>%  length()
Predicts_SV$Primary_diet %>%  unique()

## checking sample sizes by diet
Predicts_SV_diet <- unique(Predicts_SV[, c("Class", "Best_guess_binomial", "Primary_diet")]) 
table(Predicts_SV_diet$Class, Predicts_SV_diet$Primary_diet)

## grouping diet
Predicts_SV$Primary_diet <- as.character(Predicts_SV$Primary_diet)
Predicts_SV$Primary_diet[Predicts_SV$Primary_diet=="SE"] <- "PL|SE"
Predicts_SV$Primary_diet[Predicts_SV$Primary_diet=="PL"] <- "PL|SE"
Predicts_SV$Primary_diet[Predicts_SV$Primary_diet=="NE"] <- "FR|NE"
Predicts_SV$Primary_diet[Predicts_SV$Primary_diet=="FR"] <- "FR|NE"

###########################################################################################

## get total RMR (assemblage level) -- by trophic level

## total assemblage level RMR (mass-dependent) as sum(RMR*Abundance)
tRMR_SV_TL <- Predicts_SV %>% 
  group_by(SS, SSB, SSBS, LandUse, Use_intensity, Trophic_level, Annual_mean_temperature, Latitude) %>%
  summarise(tRMR=sum(exp(log_BMR)*(Measurement)))

## checking distribution & log-transform
hist(tRMR_SV_TL$tRMR)
hist(log(tRMR_SV_TL$tRMR))
tRMR_SV_TL$log_tRMR <- log(tRMR_SV_TL$tRMR)

## check sample sizes
tRMR_SV_TL$SS %>%  unique %>%  length()
tRMR_SV_TL$SSBS %>%  unique %>% length()
tRMR_SV_TL %>%
  group_by(LandUse, Use_intensity) %>%
  summarise(c=n()) %>% as.data.frame() %>%
  arrange(desc(c))

## get total RMR (assemblage level) -- by diet groups
tRMR_SV_PD <- Predicts_SV %>% 
  group_by(SS, SSB, SSBS, LandUse, Use_intensity, Primary_diet, Annual_mean_temperature, Latitude) %>%
  summarise(tRMR=sum(exp(log_BMR)*(Measurement)))

tRMR_SV_PD$log_tRMR <- log(tRMR_SV_PD$tRMR)


###########################################################################################
## model selection: effects of land use, use intensity, trophic level or diet on tRMR


########################################################################################### with Trophic level
## with secondary vegetation grouped
ModelSelect <- GLMERSelect(modelData=tRMR_SV_TL,
                           responseVar = "log_tRMR",
                           fitFamily = "gaussian",
                           fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
                           fitInteractions = TRUE,
                           randomStruct = "(1|SS)",
                           verbose = TRUE)

## checking wheher Annual mean temperautre is retained
ModelSelect <- GLMERSelect(modelData=tRMR_SV_TL,
                           responseVar = "log_tRMR",
                           fitFamily = "gaussian",
                           fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
                           fixedTerms = list(Annual_mean_temperature=1),
                           fitInteractions = TRUE,
                           randomStruct = "(1|SS)",
                           verbose = TRUE)

# main effect of annual mean temperature -- not significant 
ModelSelect$model %>%  summary()

## checking latitude -- is it retained?
ModelSelect <- GLMERSelect(modelData=tRMR_SV_TL,
                           responseVar = "log_tRMR",
                           fitFamily = "gaussian",
                           fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
                           fixedTerms = list(Latitude=1),
                           fitInteractions = TRUE,
                           randomStruct = "(1|SS)",
                           verbose = TRUE)

# main effect of latitude -- not significant 
ModelSelect$model %>%  summary()

## checking whether model with 3-way interactions fits better than full model
Full_model_TL <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
                        LandUse:Use_intensity+
                        LandUse:Trophic_level+
                        Use_intensity:Trophic_level+
                        (1|SS),
                      data=tRMR_SV_TL)

Model_3_way_TL <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
                      LandUse:Use_intensity+
                      LandUse:Trophic_level+
                      Use_intensity:Trophic_level+
                      LandUse:Use_intensity:Trophic_level+
                      (1|SS),
                    data=tRMR_SV_TL)

anova(Full_model_TL, Model_3_way_TL)

## for trophic levels, model with 3-way interactions fits better

## testing formally if adding annual mean temperature would improve fit
m3 <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+Annual_mean_temperature+
             LandUse:Use_intensity+
             LandUse:Trophic_level+
             Use_intensity:Trophic_level+
             LandUse:Use_intensity:Trophic_level+
             (1|SS),
           data=tRMR_SV_TL)

m3 <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+Annual_mean_temperature+
             LandUse:Use_intensity+
             LandUse:Trophic_level+
             Use_intensity:Trophic_level+
             LandUse:Use_intensity:Trophic_level+
             (1|SS),
           data=m3@frame)

m2 <-  lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
              LandUse:Use_intensity+
              LandUse:Trophic_level+
              Use_intensity:Trophic_level+
              LandUse:Use_intensity:Trophic_level+
              (1|SS),
            data=m3@frame)

anova(m2, m3) ## adding annual mean temperature doesn't improve fit here

########################################################################################### with diet

## with secondary vegetation grouped
ModelSelect <- GLMERSelect(modelData=tRMR_SV_PD,
                           responseVar = "log_tRMR",
                           fitFamily = "gaussian",
                           fixedFactors = c("LandUse","Use_intensity", "Primary_diet"),
                           fitInteractions = TRUE,
                           randomStruct = "(1|SS)",
                           verbose = TRUE)
## only two interactions are retained here


## checking wheher Annual mean temperautre is retained / singificant as main effect
ModelSelect <- GLMERSelect(modelData=tRMR_SV_PD,
                           responseVar = "log_tRMR",
                           fitFamily = "gaussian",
                           fixedFactors = c("LandUse","Use_intensity", "Primary_diet"),
                           fixedTerms = list(Annual_mean_temperature=1),
                           fitInteractions = TRUE,
                           randomStruct = "(1|SS)",
                           verbose = TRUE)

"log_tRMR~LandUse+Primary_diet+poly(Annual_mean_temperature,1)+Use_intensity+
LandUse:Primary_diet+
LandUse:poly(Annual_mean_temperature,1)+
Use_intensity:Primary_diet+
Use_intensity:poly(Annual_mean_temperature,1)+
Primary_diet:poly(Annual_mean_temperature,1)+(1|SS)"

# main effect of annual mean temperature -- significant effects
ModelSelect$model %>%  summary()

## checking latitude -- is it retained / significant
ModelSelect <- GLMERSelect(modelData=tRMR_SV_PD,
                           responseVar = "log_tRMR",
                           fitFamily = "gaussian",
                           fixedFactors = c("LandUse","Use_intensity", "Primary_diet"),
                           fixedTerms = list(Latitude=1),
                           fitInteractions = TRUE,
                           randomStruct = "(1|SS)",
                           verbose = TRUE)

# main effect of latitude -- not significant 
ModelSelect$model %>%  summary()

## checking whether model with 3-way interactions fits better than full model
Full_model_PD <- lmer(log_tRMR~LandUse+Use_intensity+Primary_diet+
                        LandUse:Use_intensity+
                        LandUse:Primary_diet+
                        Use_intensity:Primary_diet+
                        (1|SS),
                      data=tRMR_SV_PD)


Model_3_way_PD <- lmer(log_tRMR~LandUse+Use_intensity+Primary_diet+
                         LandUse:Use_intensity+
                         LandUse:Primary_diet+
                         Use_intensity:Primary_diet+
                         LandUse:Use_intensity:Primary_diet+
                         (1|SS),
                       data=tRMR_SV_PD)

anova(Full_model_TL, Model_3_way_TL)

## for trophic levels, model with 3-way interactions fits better

## checking effects of temperature
summary(lmer(log_tRMR~LandUse+Use_intensity+Primary_diet+Annual_mean_temperature+
               LandUse:Use_intensity+
               LandUse:Primary_diet+
               Use_intensity:Primary_diet+
               LandUse:Use_intensity:Primary_diet+
               (1|SS),
             data=tRMR_SV_PD))

## no significant effect of annual mean temperature here (main effect)

## checking formally: does adding annual mean temperature to the full model improve fit?
m2 <- lmer(log_tRMR~LandUse+Use_intensity+Primary_diet+Annual_mean_temperature+
             LandUse:Use_intensity+
             LandUse:Primary_diet+
             Use_intensity:Primary_diet+
             LandUse:Use_intensity:Primary_diet+
             (1|SS),
           data=tRMR_SV_PD)

m2 <- lmer(log_tRMR~LandUse+Use_intensity+Primary_diet+Annual_mean_temperature+
             LandUse:Use_intensity+
             LandUse:Primary_diet+
             Use_intensity:Primary_diet+
             LandUse:Use_intensity:Primary_diet+
             (1|SS),
           data=m2@frame)

m1 <- lmer(log_tRMR~LandUse+Use_intensity+Primary_diet+
             LandUse:Use_intensity+
             LandUse:Primary_diet+
             Use_intensity:Primary_diet+
             LandUse:Use_intensity:Primary_diet+
             (1|SS),
           data=m2@frame)

anova(m1,m2) ## adding annual mean temperature doesn't improve fit here


######################################################## PLOTTING EFFECTS -- trophic level

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])

## for trophic level

NewdataTL <- expand.grid(LandUse=levels(Model_3_way_TL@frame$LandUse), 
                       Use_intensity=levels(Model_3_way_TL@frame$Use_intensity), 
                       Trophic_level=levels(Model_3_way_TL@frame$Trophic_level),  
                       log_tRMR=1)

Effects_TL <-  Predict_effects(NewdataTL, 
                               Model_3_way_TL, 
                               TRUE,
                               "categorical", 
                               seMultiplier = NULL,
                               18)

# no error bars for ref level
Effects_TL$Upper[Effects_TL$LandUse=="Primary vegetation" & Effects_TL$Use_intensity=="Minimal use"] <- NA
Effects_TL$Lower[Effects_TL$LandUse=="Primary vegetation"& Effects_TL$Use_intensity=="Minimal use"] <- NA

# recentrering predictions
Effects_TL$Median <- Effects_TL$Median-100
Effects_TL$Lower <- Effects_TL$Lower-100
Effects_TL$Upper <- Effects_TL$Upper-100

## plotting
pTL <- ggplot(Effects_TL,
              aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=Limits, labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level log-total RMR \n(% difference from primary)")

ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.pdf", width=9, height=4)
ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.png", width=9, height=4)

## model diagnostics
check_model(Model_3_way_TL, check = c("qq", "normality"))
r2(Model_3_way_TL)

## anova on model terms
Anova <- anova(Model_3_way_TL)
Anova_ss <- Anova$"Sum Sq"
Anova_ss <- cbind(Anova,PctExp=Anova_ss/sum(Anova_ss)*100)
print(Anova_ss[order(Anova_ss$PctExp),])


######################################################## PLOTTING EFFECTS -- trophic level

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])

## for trophic level

NewdataTL <- expand.grid(LandUse=levels(Model_3_way_TL@frame$LandUse), 
                         Use_intensity=levels(Model_3_way_TL@frame$Use_intensity), 
                         Trophic_level=levels(Model_3_way_TL@frame$Trophic_level),  
                         log_tRMR=1)

Effects_TL <-  Predict_effects(NewdataTL, 
                               Model_3_way_TL, 
                               TRUE,
                               "categorical", 
                               seMultiplier = NULL,
                               18)

# no error bars for ref level
Effects_TL$Upper[Effects_TL$LandUse=="Primary vegetation" & Effects_TL$Use_intensity=="Minimal use"] <- NA
Effects_TL$Lower[Effects_TL$LandUse=="Primary vegetation"& Effects_TL$Use_intensity=="Minimal use"] <- NA

# recentrering predictions
Effects_TL$Median <- Effects_TL$Median-100
Effects_TL$Lower <- Effects_TL$Lower-100
Effects_TL$Upper <- Effects_TL$Upper-100

## plotting
pTL <- ggplot(Effects_TL,
              aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=Limits, labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level log-total RMR \n(% difference from primary)")

ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.pdf", width=9, height=4)
ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.png", width=9, height=4)

## model diagnostics
check_model(Model_3_way_TL, check = c("qq", "normality"))
r2(Model_3_way_TL)

## anova on model terms
Anova <- anova(Model_3_way_TL)
Anova_ss <- Anova$"Sum Sq"
Anova_ss <- cbind(Anova,PctExp=Anova_ss/sum(Anova_ss)*100)
print(Anova_ss[order(Anova_ss$PctExp),])


######################################################## PLOTTING EFFECTS -- trophic level

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])

## for trophic level

NewdataTL <- expand.grid(LandUse=levels(Model_3_way_TL@frame$LandUse), 
                         Use_intensity=levels(Model_3_way_TL@frame$Use_intensity), 
                         Trophic_level=levels(Model_3_way_TL@frame$Trophic_level),  
                         log_tRMR=1)

Effects_TL <-  Predict_effects(NewdataTL, 
                               Model_3_way_TL, 
                               TRUE,
                               "categorical", 
                               seMultiplier = NULL,
                               18)

# no error bars for ref level
Effects_TL$Upper[Effects_TL$LandUse=="Primary vegetation" & Effects_TL$Use_intensity=="Minimal use"] <- NA
Effects_TL$Lower[Effects_TL$LandUse=="Primary vegetation"& Effects_TL$Use_intensity=="Minimal use"] <- NA

# recentrering predictions
Effects_TL$Median <- Effects_TL$Median-100
Effects_TL$Lower <- Effects_TL$Lower-100
Effects_TL$Upper <- Effects_TL$Upper-100

## plotting
pTL <- ggplot(Effects_TL,
              aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=Limits, labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level log-total RMR \n(% difference from primary)")

ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.pdf", width=9, height=4)
ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.png", width=9, height=4)

## model diagnostics
check_model(Model_3_way_TL, check = c("qq", "normality"))
r2(Model_3_way_TL)

## anova on model terms
Anova <- anova(Model_3_way_TL)
Anova_ss <- Anova$"Sum Sq"
Anova_ss <- cbind(Anova,PctExp=Anova_ss/sum(Anova_ss)*100)
print(Anova_ss[order(Anova_ss$PctExp),])

######################################################## PLOTTING EFFECTS -- trophic level

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])

## for trophic level

NewdataTL <- expand.grid(LandUse=levels(Model_3_way_TL@frame$LandUse), 
                         Use_intensity=levels(Model_3_way_TL@frame$Use_intensity), 
                         Trophic_level=levels(Model_3_way_TL@frame$Trophic_level),  
                         log_tRMR=1)

Effects_TL <-  Predict_effects(NewdataTL, 
                               Model_3_way_TL, 
                               TRUE,
                               "categorical", 
                               seMultiplier = NULL,
                               18)

# no error bars for ref level
Effects_TL$Upper[Effects_TL$LandUse=="Primary vegetation" & Effects_TL$Use_intensity=="Minimal use"] <- NA
Effects_TL$Lower[Effects_TL$LandUse=="Primary vegetation"& Effects_TL$Use_intensity=="Minimal use"] <- NA

# recentrering predictions
Effects_TL$Median <- Effects_TL$Median-100
Effects_TL$Lower <- Effects_TL$Lower-100
Effects_TL$Upper <- Effects_TL$Upper-100

## plotting
pTL <- ggplot(Effects_TL,
              aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=Limits, labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level log-total RMR \n(% difference from primary)")

ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.pdf", width=9, height=4)
ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.png", width=9, height=4)

## model diagnostics
check_model(Model_3_way_TL, check = c("qq", "normality"))
r2(Model_3_way_TL)

## anova on model terms
Anova <- anova(Model_3_way_TL)
Anova_ss <- Anova$"Sum Sq"
Anova_ss <- cbind(Anova,PctExp=Anova_ss/sum(Anova_ss)*100)
print(Anova_ss[order(Anova_ss$PctExp),])


######################################################## PLOTTING EFFECTS -- diet

NewdataPD <- expand.grid(LandUse=levels(Model_3_way_PD@frame$LandUse), 
                         Use_intensity=levels(Model_3_way_PD@frame$Use_intensity), 
                         Primary_diet=levels(Model_3_way_PD@frame$Primary_diet),  
                         log_tRMR=1)

Effects_PD <-  Predict_effects(NewdataPD, 
                               Model_3_way_PD, 
                               TRUE,
                               "categorical", 
                               seMultiplier = NULL,
                               18)

# no error bars for ref level
Effects_PD$Upper[Effects_PD$LandUse=="Primary vegetation" & Effects_PD$Use_intensity=="Minimal use"] <- NA
Effects_PD$Lower[Effects_PD$LandUse=="Primary vegetation"& Effects_PD$Use_intensity=="Minimal use"] <- NA

# recentrering predictions
Effects_PD$Median <- Effects_PD$Median-100
Effects_PD$Lower <- Effects_PD$Lower-100
Effects_PD$Upper <- Effects_PD$Upper-100

## plotting
pPD <- ggplot(Effects_PD,
              aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=Limits, labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(~Primary_diet, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level log-total RMR \n(% difference from primary)")

ggsave(pPD, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV_PD.pdf", width=12, height=4)
ggsave(pPD, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV_PD.png", width=12, height=4)

## model diagnostics
check_model(Model_3_way_TL, check = c("qq", "normality"))
r2(Model_3_way_TL)

## anova on model terms
Anova <- anova(Model_3_way_TL)
Anova_ss <- Anova$"Sum Sq"
Anova_ss <- cbind(Anova,PctExp=Anova_ss/sum(Anova_ss)*100)
print(Anova_ss[order(Anova_ss$PctExp),])







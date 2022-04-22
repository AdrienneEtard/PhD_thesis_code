## Investigating effects of land use and use intensity on total, assemblage-level tRMR, within trophic levels


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


## Occurrence data in PREDICTS
PredictsTraits$Occurrence <- ifelse(PredictsTraits$Measurement > 0, 1, 0)

## filter for occurring species, for which Abundance was measured, and remove NA land uses and use intensities

# for grouped secondary vegetation
Predicts_SV <- PredictsTraits %>% 
  filter(Occurrence!=0) %>% 
  filter(Diversity_metric=="abundance") %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity))

# with all stages for secondary vegetation
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


###########################################################################################

## get total RMR (assemblage level) 

## total assemblage level RMR (mass-dependent) as sum(RMR*Abundance)
tRMR <- PredictsTraits %>% 
  group_by(SS, SSB, SSBS, Predominant_land_use, Use_intensity, Trophic_level) %>%
  summarise(tRMR=sum(exp(log_BMR)*(Measurement)))

tRMR_SV <- Predicts_SV %>% 
  group_by(SS, SSB, SSBS, LandUse, Use_intensity, Trophic_level, Annual_mean_temperature) %>%
  summarise(tRMR=sum(exp(log_BMR)*(Measurement)))

Predicts_SV$Best_guess_binomial %>%  unique %>%  length()

## checking distribution
hist(tRMR$tRMR)
hist(log(tRMR$tRMR))

hist(tRMR_SV$tRMR)
hist(log(tRMR_SV$tRMR))

## log-transformation of tRMR
tRMR$log_tRMR <- log(tRMR$tRMR)
tRMR_SV$log_tRMR <- log(tRMR_SV$tRMR)

# ## sample sizes
# tRMR$SS %>%  unique %>%  length()
# tRMR$SSBS %>%  unique %>% length()
# tRMR %>% 
#   group_by(Predominant_land_use, Use_intensity) %>% 
#   summarise(c=n()) %>% as.data.frame() %>% 
#   arrange(desc(c))


tRMR_SV$SS %>%  unique %>%  length()
tRMR_SV$SSBS %>%  unique %>% length()

###########################################################################################
## model selection: effects of land use, use intensity, trophic level on tRMR


###########################################################################################
## with secondary vegetation grouped
ModelSelect <- GLMERSelect(modelData=tRMR_SV,
                            responseVar = "log_tRMR",
                            fitFamily = "gaussian",
                            fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
                            fitInteractions = TRUE,
                            randomStruct = "(1|SS)",
                            verbose = TRUE)

# ## checking wheher Annual mean temperautre is retained
# ModelSelect <- GLMERSelect(modelData=tRMR_SV,
#                            responseVar = "log_tRMR",
#                            fitFamily = "gaussian",
#                            fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
#                            fixedTerms = list(Annual_mean_temperature=1),
#                            fitInteractions = TRUE,
#                            randomStruct = "(1|SS)",
#                            verbose = TRUE)
# 
# # temperature retianed, but main effect of annual mean temperature -- not significant 
# ModelSelect$model %>%  summary()
# 
# ## checking latitude -- is it retained?
# ModelSelect <- GLMERSelect(modelData=tRMR_SV,
#                            responseVar = "log_tRMR",
#                            fitFamily = "gaussian",
#                            fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
#                            fixedTerms = list(Latitude=1),
#                            fitInteractions = TRUE,
#                            randomStruct = "(1|SS)",
#                            verbose = TRUE)
# # main effect of latitude retained but not significant 
# ModelSelect$model %>%  summary()


## best-fitting model includes all 3 interaction terms
"log_tRMR~LandUse+Use_intensity+Trophic_level+
LandUse:Use_intensity+
LandUse:Trophic_level+
Use_intensity:Trophic_level+
(1|SS)"

Best_model_tRMR_SV <- ModelSelect$model
Best_model_tRMR_SV %>%  summary()
r2(Best_model_tRMR_SV)
AIC(Best_model_tRMR_SV)

Anova <- anova(Best_model_tRMR_SV)
Anova_ss <- Anova$"Sum Sq"
Anova_ss <- cbind(Anova,PctExp=Anova_ss/sum(Anova_ss)*100)
print(Anova_ss[order(Anova_ss$PctExp),])

## does including the 3-way interaction improve model?
Model_Best <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
                      LandUse:Use_intensity+
                      LandUse:Trophic_level+
                      Use_intensity:Trophic_level+
                      (1|SS),
                    data=tRMR_SV)

Model_3_way <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
                      LandUse:Use_intensity+
                      LandUse:Trophic_level+
                      Use_intensity:Trophic_level+
                      LandUse:Use_intensity:Trophic_level+
                      (1|SS),
                    data=tRMR_SV)

saveRDS(Model_3_way, "../Results/6.tRMR/Model3way.rds")


## test which model fits better
AnovaTest <- anova(Model_Best, Model_3_way) ## model with 3-way interaction fits better
AnovaTest$`Pr(>Chisq)`

## model diagnostics
check_model(Model_3_way, check = c("qq", "normality"))
rm(ModelSelect)
r2(Model_3_way)

## plotting effects for model with SV grouped (model with 3-way interactions)

Model_3_way <- readRDS( "../Results/6.tRMR/Model3way.rds")

Newdata <- expand.grid(LandUse=levels(Model_3_way@frame$LandUse), 
                       Use_intensity=levels(Model_3_way@frame$Use_intensity), 
                       Trophic_level=levels(Model_3_way@frame$Trophic_level),  
                       log_tRMR=1)

Effects_SV <-  Predict_effects(Newdata, 
                               Model_3_way, 
                               TRUE,
                               "categorical", 
                               seMultiplier = NULL,
                               18)

# no error bars for ref level
Effects_SV$Upper[Effects_SV$LandUse=="Primary vegetation" & Effects_SV$Use_intensity=="Minimal use"] <- NA
Effects_SV$Lower[Effects_SV$LandUse=="Primary vegetation"& Effects_SV$Use_intensity=="Minimal use"] <- NA

# recentrering predictions
Effects_SV$Median <- Effects_SV$Median -100
Effects_SV$Lower <- Effects_SV$Lower -100
Effects_SV$Upper <- Effects_SV$Upper -100

## plotting
scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])



pSV <- ggplot(Effects_SV,
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
  facet_wrap(~Trophic_level, nrow=1) +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level total RMR \n(% difference from primary)")

ggsave(pSV, filename="c:/Users/adrie/OneDrive/Desktop/Thesis/figures/Chapter5/Figure3.pdf", width=9, height=4)


ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.pdf", width=9, height=4)
ggsave(pSV, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR_SV.png", width=9, height=4)
saveRDS(pSV, "../Results/Figures/6.Model_total_BMR/Plot_data_SV.rds")

#############################################################################################

Effects_SV$Y_lines <- rep(Effects_SV$Median[Effects_SV$LandUse=="Primary vegetation"], each=6)

ggplot(Effects_SV,
       aes(LandUse, Median, ymin = Lower, ymax = Upper)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +
  geom_hline(aes(yintercept = Y_lines), linetype="dashed")+
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=Limits, labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(Use_intensity~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level total RMR \n(% difference from primary)")




###########################################################################################

## checking changes in median RMR in occurring species 

PredictsTraitsAll <- "../Results/PredictsTraits.rds" %>% 
  readRDS() %>% 
  filter(Occurrence!=0) %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity))

PredictsTraitsAll$SS %>% unique %>% length()
PredictsTraitsAll$SSBS %>% unique %>% length()

length(unique(PredictsTraitsAll$Best_guess_binomial[PredictsTraitsAll$Thermoregulation=="Endotherms"]))
length(unique(PredictsTraitsAll$Best_guess_binomial[PredictsTraitsAll$Thermoregulation=="Ectotherms"]))
length(unique(PredictsTraitsAll$Best_guess_binomial))

tRMR_Median <- PredictsTraitsAll %>% 
  group_by(SS, SSB, SSBS, LandUse, Use_intensity, Trophic_level) %>%
  summarise(RMR_Median=median(exp(log_BMR)*(Occurrence)))

hist(tRMR_Median$RMR_Median)
hist(log(tRMR_Median$RMR_Median))
tRMR_Median$log_RMR_Median <- log(tRMR_Median$RMR_Median)

tRMR_Median$SS %>%  unique %>%  length()
tRMR_Median$SSBS %>%  unique %>% length()
tRMR %>%
  group_by(Predominant_land_use, Use_intensity) %>%
  summarise(c=n()) %>% as.data.frame() %>%
  arrange(desc(c))

ModelSelect <- GLMERSelect(modelData=tRMR_Median,
                           responseVar = "log_RMR_Median",
                           fitFamily = "gaussian",
                           fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
                           fitInteractions = TRUE,
                           randomStruct = "(1|SS)",
                           verbose = TRUE)

##  all 3 interactions retained
Model_full_medianRMR <- lmer(log_RMR_Median~
                                LandUse+Use_intensity+Trophic_level+
                                LandUse:Use_intensity+
                                LandUse:Trophic_level+
                                Use_intensity:Trophic_level+
                                (1|SS),
                              data=tRMR_Median)

Model_3_way_medianRMR <- lmer(log_RMR_Median~
                                LandUse+Use_intensity+Trophic_level+
                                LandUse:Use_intensity+
                                LandUse:Trophic_level+
                                Use_intensity:Trophic_level+
                                LandUse:Use_intensity:Trophic_level+
                                (1|SS),
                    data=tRMR_Median)

anova(Model_full_medianRMR, Model_3_way_medianRMR)

## plotting effects for model on median RMR

Newdata <- expand.grid(LandUse=levels(Model_3_way_medianRMR@frame$LandUse), 
                       Use_intensity=levels(Model_3_way_medianRMR@frame$Use_intensity), 
                       Trophic_level=levels(Model_3_way_medianRMR@frame$Trophic_level),  
                       log_RMR_Median=1)

Effects_SV <-  Predict_effects(Newdata, 
                               Model_3_way_medianRMR, 
                               TRUE,
                               "categorical", 
                               seMultiplier = NULL,
                               18)

# no error bars for ref level
Effects_SV$Upper[Effects_SV$LandUse=="Primary vegetation" & Effects_SV$Use_intensity=="Minimal use"] <- NA
Effects_SV$Lower[Effects_SV$LandUse=="Primary vegetation"& Effects_SV$Use_intensity=="Minimal use"] <- NA

# recentrering predictions
Effects_SV$Median <- Effects_SV$Median -100
Effects_SV$Lower <- Effects_SV$Lower -100
Effects_SV$Upper <- Effects_SV$Upper -100

## plotting
scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])

pMedian <- ggplot(Effects_SV,
              aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=LandUse)) +
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
  ylab("Assemblage-level median RMR \n(% difference from primary)")


##################################################################### Arranging figure for manuscript

Total <- pSV$data
Total$Which <- "tRMR"
Median <- pMedian$data
Median$Which <- "mRMR"
colnames(Total)[7] <- "RMR"
colnames(Median)[7] <- "RMR"

ToPlot <- rbind(Total, Median) 
ToPlot$Which <- factor(ToPlot$Which, levels = c("tRMR", "mRMR"))

PlotMedianTotal <- 
ggplot(ToPlot,
       aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=LandUse)) +
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
  facet_grid(Which~Trophic_level, scales="free") +
  theme(panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(1, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level RMR (% difference from primary)")


ggsave(PlotMedianTotal,
       filename="../Results/Figures/6.Model_total_BMR/Median_total_BMR_SV.pdf", width=9, height=5)

ggsave(PlotMedianTotal,
       filename="../Results/Figures/6.Model_total_BMR/Median_total_BMR_SV.png", width=9, height=5)

saveRDS(PlotMedianTotal, "../Results/Figures/6.Model_total_BMR/Plot_data_Total_Median.rds")





#ggarrange(pSV, pMedian, common.legend = TRUE)

# ## checking changes in abundance among species
# 
# tAbundance <- Predicts_SV %>% 
#   group_by(SS, SSB, SSBS, LandUse, Use_intensity, Trophic_level) %>%
#   summarise(tAbundance=sum(Measurement))
# 
# hist(tAbundance$tAbundance)
# hist(log(tAbundance$tAbundance))
# tAbundance$log_abundance <- log(tAbundance$tAbundance)
# 
# Model_abundance <- lmer(log_abundance~
#                           LandUse+Use_intensity+Trophic_level+
#                           LandUse:Use_intensity+
#                           LandUse:Trophic_level+
#                           Use_intensity:Trophic_level+
#                           LandUse:Use_intensity:Trophic_level+
#                           (1|SS),
#                         data=tAbundance)
# 
# ## plotting effects for model on abundance
# 
# Newdata <- expand.grid(LandUse=levels(Model_abundance@frame$LandUse), 
#                        Use_intensity=levels(Model_abundance@frame$Use_intensity), 
#                        Trophic_level=levels(Model_abundance@frame$Trophic_level),  
#                        log_abundance=1)
# 
# Effects_ab <-  Predict_effects(Newdata, 
#                                Model_3_way_medianRMR, 
#                                TRUE,
#                                "categorical", 
#                                seMultiplier = NULL,
#                                18)
# 
# # no error bars for ref level
# Effects_ab$Upper[Effects_ab$LandUse=="Primary vegetation" & Effects_ab$Use_intensity=="Minimal use"] <- NA
# Effects_ab$Lower[Effects_ab$LandUse=="Primary vegetation"& Effects_ab$Use_intensity=="Minimal use"] <- NA
# 
# # recentrering predictions
# Effects_ab$Median <- Effects_ab$Median -100
# Effects_ab$Lower <- Effects_ab$Lower -100
# Effects_ab$Upper <- Effects_ab$Upper -100
# 
# ## plotting
# scales::show_col(viridis(option = "plasma", n=8))
# cols <- viridis(option = "plasma", n=8)[1:7]
# cbPalette <- c("#000000", cols[c(2:5,7)])
# 
# pAbundance <- ggplot(Effects_ab,
#                   aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=LandUse)) +
#   geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
#   ylab("") + xlab("") +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
#   geom_point(size=2, position=position_dodge(width = 0.7)) +
#   scale_x_discrete(limits=Limits, labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
#   GGPoptions +
#   ggtitle("") +
#   facet_grid(~Trophic_level, scales="free") +
#   theme(panel.spacing = unit(0, "lines")) +
#   theme( strip.text.x = element_text(size = 12, face = "bold"),
#          strip.text.y = element_text(size = 12, face = "bold")) +
#   scale_colour_manual(values=cbPalette, name="Use intensity") +
#   theme(legend.position = "right") +
#   scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
#   guides(colour=FALSE) +
#   theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
#   ylab("Assemblage-level log-total abundance \n(% difference from primary)")
# 
# 
# 
# 
# ###########################################################################################
# 
# # ## with secondary vegetation not grouped (3 stages)
# # ModelSelect <- GLMERSelect(modelData=tRMR,
# #                            responseVar = "log_tRMR",
# #                            fitFamily = "gaussian",
# #                            fixedFactors = c("Predominant_land_use","Use_intensity", "Trophic_level"),
# #                            fitInteractions = TRUE,
# #                            randomStruct = "(1|SS)",
# #                            verbose = TRUE)
# # 
# # ## best-fitting model includes all 3 interaction terms
# # "log_tRMR~Predominant_land_use+Use_intensity+Trophic_level+
# # Predominant_land_use:Use_intensity+
# # Predominant_land_use:Trophic_level+
# # Use_intensity:Trophic_level+
# # (1|SS)"
# # 
# # Best_model_tRMR <- ModelSelect$model
# # Best_model_tRMR %>%  summary()
# # r2(Best_model_tRMR)
# # AIC(Best_model_tRMR)
# # 
# # Anova <- anova(Best_model_tRMR)
# # Anova_ss <- Anova$"Sum Sq"
# # Anova_ss <- cbind(Anova,PctExp=Anova_ss/sum(Anova_ss)*100)
# # print(Anova_ss[order(Anova_ss$PctExp),])
# # 
# # ## model diagnostics
# # # check_model(Best_model_tRMR, check = c("qq", "normality"))
# # 
# # 
# # ###########################################################################################
# # ###########################################################################################
# # ###########################################################################################
# # 
# # 
# # # ## plotting effects for model with SV in 3 stages
# # 
# # Newdata <- expand.grid(Predominant_land_use=levels(Best_model_tRMR@frame$Predominant_land_use), 
# #                        Use_intensity=levels(Best_model_tRMR@frame$Use_intensity), 
# #                        Trophic_level=levels(Best_model_tRMR@frame$Trophic_level),  
# #                        log_tRMR=1)
# # 
# # Effects <-  Predict_effects(Newdata, 
# #                                Best_model_tRMR, 
# #                                TRUE,
# #                                "categorical", 
# #                                seMultiplier = NULL, 
# #                                24)
# # 
# # # no error bars for ref level
# # Effects$Upper[Effects$LandUse=="Primary vegetation" & Effects$Use_intensity=="Minimal use"] <- NA
# # Effects$Lower[Effects$LandUse=="Primary vegetation"& Effects$Use_intensity=="Minimal use"] <- NA
# # 
# # # recentrering predictions
# # Effects$Median <- Effects$Median -100
# # Effects$Lower <- Effects$Lower -100
# # Effects$Upper <- Effects$Upper -100
# # 
# # ## plotting
# # scales::show_col(viridis(option = "plasma", n=8))
# # cols <- viridis(option = "plasma", n=8)[1:7]
# # cbPalette <- c("#000000", cols)
# # scales::show_col(cbPalette)
# # 
# # p_all <- ggplot(Effects,
# #               aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
# #   geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
# #   geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
# #   geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
# #   geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
# #   geom_rect(xmin=8.5, xmax=9.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
# #   ylab("") + xlab("") +
# #   geom_hline(yintercept = 0, linetype="dashed") +
# #   geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
# #   geom_point(size=2, position=position_dodge(width = 0.7)) +
# #   scale_x_discrete(
# #                    labels=c("Primary",
# #                             "Secondary (mature)", 
# #                             "Secondary (intermediate)",
# #                             "Secondary (young)",
# #                             "Plantation", 
# #                             "Pasture",
# #                             "Cropland", 
# #                             "Urban")) +
# #   GGPoptions +
# #   ggtitle("") +
# #   facet_grid(~Trophic_level, scales="free") +
# #   theme(panel.spacing = unit(0, "lines")) +
# #   theme( strip.text.x = element_text(size = 12, face = "bold"),
# #          strip.text.y = element_text(size = 12, face = "bold")) +
# #   scale_colour_manual(values=cbPalette, name="Use intensity") +
# #   theme(legend.position = "right") +
# #   scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
# #   guides(colour=FALSE) +
# #   theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
# #   ylab("Assemblage-level log-total RMR \n(% difference from primary)")
# # 
# # 
# # ggsave(p_all, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR.pdf", width=9, height=4)
# # ggsave(p_all, filename="../Results/Figures/6.Model_total_BMR/Model_total_BMR.png", width=9, height=4)
# # 
# 
# 

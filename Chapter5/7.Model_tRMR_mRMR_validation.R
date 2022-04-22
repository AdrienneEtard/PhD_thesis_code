## Script for fitting models explaining tRMR and mRMR by landuse and landuse intesnity
## validation using complete RMR data

library(StatisticalModels)
library(scales)
library(viridis)
library(performance)
library(dplyr)

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

Predict_Effects <- function(Newdata, Model, rescale, Cont_or_Cat, seMultiplier, LU_n) {

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
    #Newdata$Estimate <- logit2prob(preds)
    
    # getting the errors around the predicted values
    # variance covariance matrix and matrix multiplication
    VarCoVar <- as.matrix(vcov(Model))
    VarCoVar <- base::tcrossprod(VarCoVar, mm)
    
    # matrix multiplication to get estimates and take diagonal matrix
    pvar <- diag(mm %*% VarCoVar)
    
    # estimate error of predictions and backtransform
    Lower <- preds - seMultiplier*sqrt(pvar)
    Upper <- preds + seMultiplier*sqrt(pvar)
    
    # Newdata$Lower <- logit2prob(Lower)
    # Newdata$Upper <- logit2prob(Upper)
    
    return(Newdata)
    
  }
}


###########################################################################################
###########################################################################################
###########################################################################################


## fitting model to explain total assemblage-level RMR: MODEL VALIDATION ON COMPLETE DATA

BMR_data <- read.csv("../Results/2.BMR_data_to_impute.csv") %>% 
  dplyr::select(-Thermoregulation, -Body_mass_g, -Imputed, -Order,-Family, -Class)
colnames(BMR_data)[1] <- "Best_guess_binomial"

PredictsTraits_complete <- "../Results/PredictsTraits.rds" %>% 
  readRDS()

## join
PredictsTraits_complete <- left_join(PredictsTraits_complete, BMR_data)

## filter for non-0 abundance
Predicts_complete <- PredictsTraits_complete %>% 
  mutate(Occurrence=ifelse(Measurement > 0, 1, 0)) %>% 
  filter(Occurrence!=0) %>% 
  filter(Diversity_metric=="abundance") %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity)) %>% 
  dplyr::select(-Residual_BMR_log_log) %>% 
  filter(!is.na(BMR_ml_O2_per_h))

Predicts_complete_m <- PredictsTraits_complete %>% 
  mutate(Occurrence=ifelse(Measurement > 0, 1, 0)) %>% 
  filter(Occurrence!=0) %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity)) %>% 
  dplyr::select(-Residual_BMR_log_log) %>% 
  filter(!is.na(BMR_ml_O2_per_h))

# levels(PredictsTraits_complete$LandUse)
# levels(PredictsTraits_complete$Use_intensity)
# levels(Predicts_complete$LandUse)
# levels(Predicts_complete$Use_intensity)

## assess how many species

# for tRMR
Predicts_complete$Best_guess_binomial %>%
  unique %>%  length

Predicts_complete %>%  
  group_by(Class, Best_guess_binomial) %>%
  summarise(C=n()) %>% 
  group_by(Class) %>% 
  summarise(C=n())

# for mRMR
Predicts_complete_m$Best_guess_binomial %>%
  unique %>%  length

Predicts_complete_m %>%  
  group_by(Class, Best_guess_binomial) %>%
  summarise(C=n()) %>% 
  group_by(Class) %>% 
  summarise(C=n())


## calculating total BMR
tRMR_complete <- Predicts_complete %>% 
  group_by(SS, SSB, SSBS, LandUse, Use_intensity, Trophic_level) %>%
  summarise(tRMR=sum(BMR_ml_O2_per_h*(Measurement)))

hist(tRMR_complete$tRMR)
hist(log(tRMR_complete$tRMR))
tRMR_complete$log_tRMR <- log(tRMR_complete$tRMR)

## calculating median BMR
mRMR_complete <- Predicts_complete_m %>% 
  group_by(SS, SSB, SSBS, LandUse, Use_intensity, Trophic_level) %>%
  summarise(mRMR=median(BMR_ml_O2_per_h*Occurrence))

hist(mRMR_complete$mRMR)
hist(log(mRMR_complete$mRMR))
mRMR_complete$log_mRMR <- log(mRMR_complete$mRMR)


## model selection for tRMR and mRMR
Model_complete_Select <- GLMERSelect(modelData=tRMR_complete,
                                         responseVar = "log_tRMR",
                                         fitFamily = "gaussian",
                                         fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
                                         fitInteractions = TRUE,
                                         randomStruct = "(1|SS)",
                                         verbose = TRUE)

Model_complete_Select2 <- GLMERSelect(modelData=mRMR_complete,
                                     responseVar = "log_mRMR",
                                     fitFamily = "gaussian",
                                     fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
                                     fitInteractions = TRUE,
                                     randomStruct = "(1|SS)",
                                     verbose = TRUE)

## best fitting model includes all three interactions for tRMR - like the model on the imputed data
## and two interactions for mRMR

## does including the 3-way interaction improve models? -- for tRMR
Model_Best_tRMR <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
                     LandUse:Use_intensity+
                     LandUse:Trophic_level+
                     Use_intensity:Trophic_level+
                     (1|SS),
                   data=tRMR_complete)

Model_3_way_tRMR <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
                      LandUse:Use_intensity+
                      LandUse:Trophic_level+
                      Use_intensity:Trophic_level+
                      LandUse:Use_intensity:Trophic_level+
                      (1|SS),
                    data=tRMR_complete)

## test which model fits better
anova(Model_Best_tRMR, Model_3_way_tRMR) ## model with 3-way interaction fits better (but 2 coefficients cannot be estimated)


## for tRMR, compare models coefficients between models on complete and imputed data
ModelImputed <- readRDS("../Results/6.tRMR/Model3way.rds")
ModelComplete <- Model_3_way_tRMR

Summary_Imputed <- as.data.frame(summary(ModelImputed)$coefficients)
Summary_Complete <- as.data.frame(summary(ModelComplete)$coefficients)
Summary_Imputed$Which <- "Collected and imputed"
Summary_Complete$Which <- "Collected"
Summary_Imputed$Coef <- rownames(Summary_Imputed)
Summary_Complete$Coef <- rownames(Summary_Complete)

## filter out coefficients that couldn't be estimated with complete model in summary imputed
Match <- intersect(Summary_Imputed$Coef, Summary_Complete$Coef)
Summary_Imputed <- Summary_Imputed %>% 
  filter(Coef %in% Match)
Summary_Complete$Coef==Summary_Imputed$Coef


## plot coefficients

pdf(file = "../Results/Figures/7.Model_total_BMR_validation/SI_Coefficients_Complete_against_Imputed.pdf", width=5, height=7, pointsize = 13)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,1), oma=c(2,1,2,1))
par(mfrow=c(2,1))

plot(Summary_Imputed$Estimate~Summary_Complete$Estimate, pch=19, 
     ylab="Main model", 
     xlab="Model on empirical RMR values only")
abline(a=0, b=1, lty="dashed", col="red")
title("(a) Coefficient estimates", adj = 0, line = 1)
plot(Summary_Imputed$`Std. Error`~Summary_Complete$`Std. Error`, pch=19,
     ylab="Main model", 
     xlab="Model on empirical RMR values only")
abline(a=0, b=1, lty="dashed", col="red")
title("(b) Standard error estimates", adj = 0, line = 1)
dev.off()

cor(Summary_Imputed$Estimate, Summary_Complete$Estimate)
cor(Summary_Imputed$`Std. Error`, Summary_Complete$`Std. Error`)



########################################################################################################################################
## does including the 3-way interaction improve models? -- for mRMR
Model_all_mRMR <- lmer(log_mRMR~LandUse+Use_intensity+Trophic_level+
                          LandUse:Use_intensity+
                          LandUse:Trophic_level+
                          Use_intensity:Trophic_level+
                          (1|SS),
                        data=mRMR_complete)

Model_3_way_mRMR <- lmer(log_mRMR~LandUse+Use_intensity+Trophic_level+
                           LandUse:Use_intensity+
                           LandUse:Trophic_level+
                           Use_intensity:Trophic_level+
                           LandUse:Use_intensity:Trophic_level+
                           (1|SS),
                         data=mRMR_complete)

## test which model fits better
anova(Model_all_mRMR, Model_3_way_mRMR) ## model with 3-way interaction fits better (but 2 coefficients cannot be estimated)

## in both cases, adding the three-way interaction improves model fit.
## models with 3-way interaction fit better (but 2 coefficients cannot be estimated)

r2(Model_3_way_tRMR)
r2(Model_3_way_mRMR)

# af <- anova(Model_3_way)
# afss <- af$"Sum Sq"
# afss <- cbind(af,PctExp=afss/sum(afss)*100)
# print(afss[order(afss$PctExp),])

## model diagnostics
check_model(Model_3_way_mRMR, check = c("qq", "normality"))
check_model(Model_3_way_tRMR, check = c("qq", "normality"))


###########################################################################################
###########################################################################################
###########################################################################################

## plotting Effects_tRMR

############################# for the best fitting model -- tRMR
Newdata <- expand.grid(LandUse=levels(Model_3_way_tRMR@frame$LandUse), 
                       Use_intensity=levels(Model_3_way_tRMR@frame$Use_intensity), 
                       Trophic_level=levels(Model_3_way_tRMR@frame$Trophic_level),  
                       log_tRMR=1)

Effects_tRMR <-  Predict_Effects(Newdata, Model_3_way_tRMR, TRUE,"categorical", seMultiplier = NULL, 18)
Effects_tRMR$Median <- Effects_tRMR$Median -100
Effects_tRMR$Lower <- Effects_tRMR$Lower -100
Effects_tRMR$Upper <- Effects_tRMR$Upper -100

## These two Effects_tRMR can't be estimated
# "LandUseCropland:Use_intensityIntense use:Trophic_levelOmnivore" 
# "LandUseUrban:Use_intensityIntense use:Trophic_levelOmnivore"

Effects_tRMR$Median[Effects_tRMR$LandUse=="Cropland" & Effects_tRMR$Use_intensity=="Intense use" & Effects_tRMR$Trophic_level=="Omnivore"] <- NA
Effects_tRMR$Median[Effects_tRMR$LandUse=="Urban" & Effects_tRMR$Use_intensity=="Intense use" & Effects_tRMR$Trophic_level=="Omnivore"] <- NA

Effects_tRMR$Lower[Effects_tRMR$LandUse=="Cropland" & Effects_tRMR$Use_intensity=="Intense use" & Effects_tRMR$Trophic_level=="Omnivore"] <- NA
Effects_tRMR$Lower[Effects_tRMR$LandUse=="Urban" & Effects_tRMR$Use_intensity=="Intense use" & Effects_tRMR$Trophic_level=="Omnivore"] <- NA

Effects_tRMR$Upper[Effects_tRMR$LandUse=="Cropland" & Effects_tRMR$Use_intensity=="Intense use" & Effects_tRMR$Trophic_level=="Omnivore"] <- NA
Effects_tRMR$Upper[Effects_tRMR$LandUse=="Urban" & Effects_tRMR$Use_intensity=="Intense use" & Effects_tRMR$Trophic_level=="Omnivore"] <- NA

## plotting

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])

p_complete_tRMR <- ggplot(Effects_tRMR,
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
  ylab("Assemblage-level log-tRMR \n(% difference from primary)")

#p_complete <- p_complete_tRMR$data
#saveRDS(p_complete, "../Results/Figures/7.Model_total_BMR_validation/Plot_data.rds")
# ggsave(p_complete, filename="../Results/Figures/7.Model_total_BMR_validation/Model_complete_RMR_SV.pdf", width=9, height=4)
# ggsave(p_complete, filename="../Results/Figures/7.Model_total_BMR_validation/Model_complete_RMR_SV.png", width=9, height=4)

############################# for the best fitting model -- mRMR
Newdata <- expand.grid(LandUse=levels(Model_3_way_mRMR@frame$LandUse), 
                       Use_intensity=levels(Model_3_way_mRMR@frame$Use_intensity), 
                       Trophic_level=levels(Model_3_way_mRMR@frame$Trophic_level),  
                       log_mRMR=1)

Effects_mRMR <-  Predict_Effects(Newdata, Model_3_way_mRMR, TRUE,"categorical", seMultiplier = NULL, 18)
Effects_mRMR$Median <- Effects_mRMR$Median -100
Effects_mRMR$Lower <- Effects_mRMR$Lower -100
Effects_mRMR$Upper <- Effects_mRMR$Upper -100

## These two Effects_tRMR can't be estimated
# "LandUseCropland:Use_intensityIntense use:Trophic_levelOmnivore" 
# "LandUseUrban:Use_intensityIntense use:Trophic_levelOmnivore"

Effects_mRMR$Median[Effects_mRMR$LandUse=="Cropland" & Effects_mRMR$Use_intensity=="Intense use" & Effects_mRMR$Trophic_level=="Omnivore"] <- NA
Effects_mRMR$Median[Effects_mRMR$LandUse=="Urban" & Effects_mRMR$Use_intensity=="Intense use" & Effects_mRMR$Trophic_level=="Omnivore"] <- NA

Effects_mRMR$Lower[Effects_mRMR$LandUse=="Cropland" & Effects_mRMR$Use_intensity=="Intense use" & Effects_mRMR$Trophic_level=="Omnivore"] <- NA
Effects_mRMR$Lower[Effects_mRMR$LandUse=="Urban" & Effects_mRMR$Use_intensity=="Intense use" & Effects_mRMR$Trophic_level=="Omnivore"] <- NA

Effects_mRMR$Upper[Effects_mRMR$LandUse=="Cropland" & Effects_mRMR$Use_intensity=="Intense use" & Effects_mRMR$Trophic_level=="Omnivore"] <- NA
Effects_mRMR$Upper[Effects_mRMR$LandUse=="Urban" & Effects_mRMR$Use_intensity=="Intense use" & Effects_mRMR$Trophic_level=="Omnivore"] <- NA

## plotting

scales::show_col(viridis(option = "plasma", n=8))
cols <- viridis(option = "plasma", n=8)[1:7]
cbPalette <- c("#000000", cols[c(2:5,7)])

p_complete_mRMR <- ggplot(Effects_mRMR,
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
  ylab("Assemblage-level log-mRMR \n(% difference from primary)")


####### Assembling plots together
Total <- p_complete_tRMR$data
Total$Which <- "tRMR"
Median <- p_complete_mRMR$data
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
       filename="../Results/Figures/7.Model_total_BMR_validation/Median_total_BMR_SV.pdf", width=9, height=5)

ggsave(PlotMedianTotal,
       filename="../Results/Figures/7.Model_total_BMR_validation/Median_total_BMR_SV.png", width=9, height=5)

saveRDS(PlotMedianTotal$data, "../Results/Figures/7.Model_total_BMR_validation/Plot_data_Total_Median.rds")


######################################################################################################################
## plotting complete against imputed

data_complete <- readRDS("../Results/Figures/7.Model_total_BMR_validation/Plot_data.rds")
data_imputed <- readRDS("../Results/Figures/6.Model_total_BMR/Plot_data_SV.rds")$data

data_complete$which <- "Empirical RMR values only"
data_imputed$which <- "Empirical and imputed RMR values"
to_plot <- rbind(data_complete, data_imputed)

collected_vs_imputed_tRMR <- 
  ggplot(to_plot,
       aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=which)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=c(),
                   labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")) +
  GGPoptions +
  ggtitle("") +
  facet_grid(Use_intensity~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 8), name="RMR data used to fit the model") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("Assemblage-level total RMR \n(% difference from primary)")

ggsave(collected_vs_imputed_tRMR, filename="../Results/Figures/7.Model_total_BMR_validation/Complete_VS_imputed_total.pdf",
       width=9, height=9)


# ggsave(p_collected_vs_imputed, filename="../Results/Figures/7.Model_total_BMR_validation/Complete_VS_imputed_median_total.pdf",
#        width=7, height=10)

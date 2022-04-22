## checking wheher Annual mean temperautre is retained

# function for predicting effects
Predict_effects <- function(Newdata, Model, rescale, Cont_or_Cat, seMultiplier, LU_n) {
  
  # browser()
  
  # function for backtransforming
  logit2prob <- function(logit){
    return(exp(logit)/(1+exp(logit)))
  }
  
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
      #y <- logit2prob(y)
      
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
    Newdata$Estimate <- preds
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
    
    Newdata$Lower <- Lower
    Newdata$Upper <- Upper
    
    # Newdata$Lower <- logit2prob(Lower)
    # Newdata$Upper <- logit2prob(Upper)
    
    return(Newdata)
    
  }
}

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

#####################################################
PredictsTraits <- "../Results/PredictsTraits.rds" %>% 
  readRDS() 

## Occurrence data in PREDICTS
PredictsTraits$Occurrence <- ifelse(PredictsTraits$Measurement > 0, 1, 0)

## filter for occurring species, for which Abundance was measured, and remove NA land uses and use intensities
# for grouped secondary vegetation
Predicts_SV <- PredictsTraits %>% 
  filter(Occurrence!=0) %>% 
  filter(Diversity_metric=="abundance") %>% 
  filter(!is.na(LandUse), !is.na(Use_intensity))

## get total RMR (assemblage level) 
## total assemblage level RMR (mass-dependent) as sum(RMR*Abundance)
tRMR_SV <- Predicts_SV %>% 
  group_by(SS, SSB, SSBS, LandUse, Use_intensity, Trophic_level, Annual_mean_temperature) %>%
  summarise(tRMR=sum(exp(log_BMR)*(Measurement)))

## log-transformation of tRMR
tRMR$log_tRMR <- log(tRMR$tRMR)
tRMR_SV$log_tRMR <- log(tRMR_SV$tRMR)


# ModelSelect <- GLMERSelect(modelData=tRMR_SV,
#                            responseVar = "log_tRMR",
#                            fitFamily = "gaussian",
#                            fixedFactors = c("LandUse","Use_intensity", "Trophic_level"),
#                            fixedTerms = list(Annual_mean_temperature=1),
#                            fitInteractions = TRUE,
#                            randomStruct = "(1|SS)",
#                            verbose = TRUE)
# 
# ModelSelect$final.call
# 
# 
# Best_model <-  lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+Annual_mean_temperature+
#                       LandUse:Use_intensity+
#                       LandUse:Trophic_level+
#                       LandUse:Annual_mean_temperature+
#                       Use_intensity:Trophic_level+
#                       Use_intensity:Annual_mean_temperature+
#                       Trophic_level:Annual_mean_temperature+
#                       (1|SS),
#                     data=tRMR_SV)


Model_3_way_temp <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+Annual_mean_temperature+
                      LandUse:Use_intensity+
                      LandUse:Trophic_level+
                      Use_intensity:Trophic_level+
                      LandUse:Use_intensity:Trophic_level+
                      (1|SS),
                    data=tRMR_SV)

Data <- Model_3_way_temp@frame

Model_3_way_temp <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+Annual_mean_temperature+
                           LandUse:Use_intensity+
                           LandUse:Trophic_level+
                           Use_intensity:Trophic_level+
                           LandUse:Use_intensity:Trophic_level+
                           (1|SS),
                         data=Data)

Model_3_way <- lmer(log_tRMR~LandUse+Use_intensity+Trophic_level+
                           LandUse:Use_intensity+
                           LandUse:Trophic_level+
                           Use_intensity:Trophic_level+
                           LandUse:Use_intensity:Trophic_level+
                           (1|SS),
                         data=Data)
X <- summary(Model_3_way_temp)$coefficients %>%  as.data.frame()

anova(Model_3_way, Model_3_way_temp)

Model_3_way_temp@frame$SS %>%  unique %>%  length

## plotting continuous effects
Newdata <- expand.grid(LandUse=levels(ModelSelect$data$LandUse), 
                       Use_intensity=levels(ModelSelect$data$Use_intensity), 
                       Trophic_level=levels(ModelSelect$data$Trophic_level),  
                       Annual_mean_temperature=seq(from=min(ModelSelect$data$Annual_mean_temperature),
                                                to=max(ModelSelect$data$Annual_mean_temperature),
                                                length.out=100),  
                       log_tRMR=1)


## predictions
Predictions <- Predict_effects(Newdata = Newdata,
                               Model = Best_model,
                               rescale = FALSE,
                               Cont_or_Cat = "continuous",
                               seMultiplier = 1.96,
                               LU_n = 6)

ggplot(Predictions,
  aes(Annual_mean_temperature, Estimate, ymin = Lower, ymax = Upper, col=Use_intensity, group=Use_intensity)) +
  ylab("") + xlab("") +
  geom_line() +
  GGPoptions +
  ggtitle("") +
  facet_grid(LandUse~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("tRMR")

ggplot(Predictions,
       aes(Annual_mean_temperature, Estimate, ymin = Lower, ymax = Upper, col=LandUse, group=LandUse)) +
  ylab("") + xlab("") +
  geom_line() +
  GGPoptions +
  ggtitle("") +
  facet_grid(Use_intensity~Trophic_level, scales="free") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=cbPalette, name="Use intensity") +
  theme(legend.position = "right") +
  scale_shape_manual(values=c(19, 17, 8), name="Land use") +
  #guides(colour=FALSE) +
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1)) +
  ylab("tRMR") + xlab("Annual mean temperature")


## discretising effects of Annual mean temperature (low, median, high)

Newdata2 <- expand.grid(LandUse=levels(ModelSelect$data$LandUse), 
                       Use_intensity=levels(ModelSelect$data$Use_intensity), 
                       Trophic_level=levels(ModelSelect$data$Trophic_level),  
                       Annual_mean_temperature=c(-5,20,28),  
                       log_tRMR=1)


## predictions
Predictions2 <- Predict_effects(Newdata = Newdata2,
                               Model = Best_model,
                               rescale = TRUE,
                               Cont_or_Cat = "categorical",
                               seMultiplier = 1.96,
                               LU_n = 6)

Predictions2$Annual_mean_temperature <- as.factor(Predictions2$Annual_mean_temperature)


ggplot(Predictions2, 
       aes(LandUse, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=LandUse)) +
         geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
         geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
         geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
         geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
         ylab("") + xlab("") +
         geom_hline(yintercept = 100, col="black", linetype="dashed") +
         geom_errorbar(width=.2, size=1, position=position_dodge(width = 0.6), stat="identity") +
         geom_point(size=2, position=position_dodge(width = 0.6)) +
         scale_x_discrete(limits=Limits, labels=Labels) +
         GGPoptions +
         scale_colour_viridis_d(option="C") +
         ggtitle("") +
         facet_grid(Annual_mean_temperature~Trophic_level, scales="free") +
         theme(panel.spacing = unit(0, "lines")) +
         theme( strip.text.x = element_text(size = 12, face = "bold"),
                strip.text.y = element_text(size = 12, face = "bold")) +
         scale_colour_manual(values=cbPalette) +
         scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
         guides(colour=FALSE)
       






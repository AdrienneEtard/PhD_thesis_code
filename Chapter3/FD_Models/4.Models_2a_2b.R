################### Models within classes effects (Models 2a, 2b)
library(StatisticalModels)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(viridis)
library(scales)

source("Functions.R")

Classes <- read.csv("../../Results/dbFD_indices_for_analysis/dbFD_classes_std.csv") %>%
  Clean()

# replace FRic and FDis values for SR=1 (not possible to estimate FRic or FDis in that case)
Classes$FRic[Classes$SR==1] <- NA
Classes$FDis[Classes$SR==1] <- NA
Classes$median_sim_FRic[Classes$SR==1] <- NA
Classes$median_sim_FDis[Classes$SR==1] <- NA

# check levels
Classes$Predominant_land_use %>%  unique()
Classes$Predominant_land_use_2 <- as.character(Classes$Predominant_land_use)
Classes$Predominant_land_use_2[Classes$Predominant_land_use_2 %in% c("Mature secondary vegetation",
                                                                     "Intermediate secondary vegetation",
                                                                     #"Secondary vegetation (indeterminate age)",
                                                                     "Young secondary vegetation")]  <- "Secondary vegetation"
Classes$Predominant_land_use_2[Classes$Predominant_land_use_2 %in% c("Cropland",
                                                                     "Pasture")]  <- "Agricultural"
Classes$Predominant_land_use_2 %>%  unique()
Classes$Predominant_land_use_2<- factor(Classes$Predominant_land_use_2, levels=c("Primary vegetation",
                                                                                 "Secondary vegetation",
                                                                                 "Plantation forest",
                                                                                 "Agricultural",
                                                                                 "Urban"))

Classes$Biome
Classes$Use_intensity

# transform values
Classes$asin_sqrt_FRic <- asin(sqrt(Classes$FRic))
Classes$asin_sqrt_FDis <- asin(sqrt(Classes$FDis))


## model 2a, selection

model2a <- GLMERSelect(modelData = Classes,
                       responseVar = "asin_sqrt_FRic",
                       fitFamily = "gaussian",
                       fixedFactors = c("Predominant_land_use_2", "Use_intensity", "Biome", "Class"),
                       fitInteractions=TRUE,
                       randomStruct = "(1|SS)+(1|SSB)",
                       verbose = TRUE)

model2b <- GLMERSelect(modelData = Classes,
                       responseVar = "asin_sqrt_FDis",
                       fitFamily = "gaussian",
                       fixedFactors = c("Predominant_land_use_2", "Use_intensity", "Biome", "Class"),
                       fitInteractions=TRUE,
                       randomStruct = "(1|SS)+(1|SSB)",
                       verbose = TRUE)


model2a <- model2a$model
model2a@call
model2b <- model2b$model
model2b@call

## plotting predictions (figure 3)

Newdata <- expand.grid(Predominant_land_use_2=levels(Classes$Predominant_land_use_2),
                       Use_intensity=levels(Classes$Use_intensity),
                       Biome=levels(Classes$Biome),
                       Class=levels(Classes$Class),
                       asin_sqrt_FDis=0,
                       asin_sqrt_FRic=0)


Predict_effects <- function(Newdata, Model, rescale) {

  preds <- sapply(X=1:1000, FUN=function(i){

    coefs <- mvrnorm(n = 1, mu = fixef(object=Model), Sigma = vcov(object=Model))
    mm <- model.matrix(terms(Model), Newdata)

    # drop coefs that couldn't be estimated
    print(setdiff(colnames(mm), names(coefs)))
    to_drop <- print(setdiff(colnames(mm), names(coefs)))
    mm <- as.data.frame(mm)
    mm <- mm[, -which(colnames(mm) %in% to_drop)]
    mm <- as.matrix(mm)

    y <- mm%*%coefs

    # backtranforming
    y <- sin(y)^2

    # rescaling
    if(rescale){
      # initialisation
      seq <- 1:5
      y[seq] <- y[seq]/y[seq[1]]*100
      # for loop to rescale all values
      for(i in 1:(nrow(Newdata)/5-1)){
        seq <- seq + 5
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

# predictions for FRic
preds_fric <- Predict_effects(Newdata, model2a, TRUE)
preds_fric$Median[preds_fric$Predominant_land_use_2=="Urban" & preds_fric$Class=="Reptiles"] <- NA
preds_fric$Lower[preds_fric$Predominant_land_use_2=="Urban" & preds_fric$Class=="Reptiles"] <- NA
preds_fric$Upper[preds_fric$Predominant_land_use_2=="Urban" & preds_fric$Class=="Reptiles"] <- NA

preds_fric$Median <- preds_fric$Median-100
preds_fric$Lower <- preds_fric$Lower-100
preds_fric$Upper <- preds_fric$Upper-100

# predictions for FDis
preds_fdis <- Predict_effects(Newdata, model2b, TRUE)
preds_fdis$Median[preds_fdis$Predominant_land_use_2=="Urban" & preds_fdis$Class=="Reptiles"] <- NA
preds_fdis$Lower[preds_fdis$Predominant_land_use_2=="Urban" & preds_fdis$Class=="Reptiles"] <- NA
preds_fdis$Upper[preds_fdis$Predominant_land_use_2=="Urban" & preds_fdis$Class=="Reptiles"] <- NA

preds_fdis$Median <- preds_fdis$Median-100
preds_fdis$Lower <- preds_fdis$Lower-100
preds_fdis$Upper <- preds_fdis$Upper-100

## plotting effects

Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Agricultural",
            "Urban")

Labels=c("PV", "SV", "PF", "AGR", "UR")

cols <- scales::show_col(viridis(option = "plasma", n=8))
cbPalette <- c("#000000", viridis(option = "plasma", n=8)[c(2,4,5,6)])
show_col(cbPalette)

p1 <- ggplot(preds_fric, aes(Predominant_land_use_2, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use_2)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  GGPoptions +
  facet_grid(Biome~Class)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")


p2 <- ggplot(preds_fdis[preds_fdis$Biome=="Temperate",], aes(Predominant_land_use_2, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use_2)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=cbPalette) +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  GGPoptions +
  facet_grid(~Class)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(b) FDis (global)") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

figure3 <- ggarrange(p1, p2, common.legend = TRUE, nrow=2, heights = c(3/5,2/5))
ggsave(figure3, filename="../Revisions/PDF_figures/Figure3.pdf",
       height=9, width=10)

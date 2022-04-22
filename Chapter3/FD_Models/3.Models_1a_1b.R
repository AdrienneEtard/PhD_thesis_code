library(StatisticalModels)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(viridis)
library(scales)

source("Functions.R")

## Regional_std is the dataset containing FRic and FDis calculated across terrestrial vertebrates

################### Models across vertebrates classes (Models 1a and 1b)

## loading main results dataset across terrestrial vertebrates
Regional_std <- read.csv("../../Results/dbFD_indices_for_analysis/dbFD_region_std_max.csv") %>%
  Clean()

# replace FRic and FDis values for SR=1 (not possible to estimate FRic or FDis in that case)
Regional_std$FRic[Regional_std$SR==1] <- NA
Regional_std$FDis[Regional_std$SR==1] <- NA
Regional_std$median_sim_FRic[Regional_std$SR==1] <- NA
Regional_std$median_sim_FDis[Regional_std$SR==1] <- NA

# check levels
Regional_std$Predominant_land_use
Regional_std$Biome
Regional_std$Use_intensity

# tranform values
Regional_std$asin_sqrt_FRic <- asin(sqrt(Regional_std$FRic))
Regional_std$asin_sqrt_FDis <- asin(sqrt(Regional_std$FDis))

## sample sizes
Regional_std <- subset(Regional_std,!is.na(Use_intensity))

Regional_std  %>%
  filter(!is.na(FRic)) %>%
  group_by(Predominant_land_use, Use_intensity, Biome) %>%
  summarise(count=n())

Regional_std  %>%
  filter(!is.na(FDis)) %>%
  group_by(Predominant_land_use, Use_intensity, Biome) %>%
  summarise(count=n())

## model 1a selection

model1a <- GLMERSelect(modelData = Regional_std,
                       responseVar = "asin_sqrt_FRic",
                       fitFamily = "gaussian",
                       fixedFactors = c("Predominant_land_use", "Use_intensity", "Biome"),
                       fitInteractions=TRUE,
                       randomStruct = "(1|SS)+(1|SSB)",
                       verbose = TRUE)

model1a <- lme4::lmer(asin_sqrt_FRic~Predominant_land_use+Use_intensity+Biome+
                        Predominant_land_use:Use_intensity+
                        Predominant_land_use:Biome+
                        (1|SS/SSB), data=Regional_std)

## model 1b selection

model1b <- GLMERSelect(modelData = Regional_std,
                       responseVar = "asin_sqrt_FDis",
                       fitFamily = "gaussian",
                       fixedFactors = c("Predominant_land_use", "Use_intensity", "Biome"),
                       fitInteractions=TRUE,
                       randomStruct = "(1|SS)+(1|SSB)",
                       verbose = TRUE)

model1b <- lme4::lmer(asin_sqrt_FDis~Predominant_land_use+Use_intensity+Biome+
                        Predominant_land_use:Use_intensity+
                        (1|SS/SSB), data=Regional_std)

## plotting effects for model 1a and 1b (Figure 2)

## predictions

Newdata <- expand.grid(Predominant_land_use=levels(Regional_std$Predominant_land_use),
                       Use_intensity=levels(Regional_std$Use_intensity),
                       Biome=levels(Regional_std$Biome),
                       asin_sqrt_FDis=0,
                       asin_sqrt_FRic=0)

levels(Newdata$Predominant_land_use)

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
      y[1:8] <- y[1:8]/y[1]*100
      y[9:16] <- y[9:16]/y[9]*100
      y[17:24] <- y[17:24]/y[17]*100
      y[25:32] <- y[25:32]/y[25]*100
      y[33:40] <- y[33:40]/y[33]*100
      y[41:48] <- y[41:48]/y[41]*100
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
preds_fric <- Predict_effects(Newdata, model1a, TRUE)
preds_fric$Median[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA
preds_fric$Lower[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA
preds_fric$Upper[preds_fric$Predominant_land_use=="Mature secondary vegetation" & preds_fric$Use_intensity=="Intense use"] <- NA

preds_fric$Median <- preds_fric$Median-100
preds_fric$Lower <- preds_fric$Lower-100
preds_fric$Upper <- preds_fric$Upper-100

# predictions for FDis
preds_fdis <- Predict_effects(Newdata, model1b, TRUE)
preds_fdis$Median[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA
preds_fdis$Lower[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA
preds_fdis$Upper[preds_fdis$Predominant_land_use=="Mature secondary vegetation" & preds_fdis$Use_intensity=="Intense use"] <- NA

preds_fdis$Median <- preds_fdis$Median-100
preds_fdis$Lower <- preds_fdis$Lower-100
preds_fdis$Upper <- preds_fdis$Upper-100


## plotting

Limits <- c("Primary vegetation",
            "Mature secondary vegetation",
            "Intermediate secondary vegetation",
            "Young secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")

Labels=c("PV", "MSV", "ISV", "YSV", "PF", "PA", "CR", "UR")

cols <- scales::show_col(viridis(option = "plasma", n=8))
cbPalette <- c("#000000", viridis(option = "plasma", n=8)[1:7])
show_col(cbPalette)

p1 <- ggplot(preds_fric, aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
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
  facet_wrap(~Biome)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

p2 <- ggplot(preds_fdis[preds_fdis$Biome=="Temperate",], aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Predominant_land_use)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.6), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.6)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") +
  scale_colour_manual(values=cbPalette) +
  GGPoptions + ggtitle("(b) FDis (Temperate; similar effects for tropical)")+
  guides(colour=FALSE, shape=FALSE)+
  ylab("Relative effect (%)")
  # facet_wrap(~Biome, labeller = labeller(Biome =  c("Temperate" = "Temperate; similar effects for tropical"))) +
  # theme( strip.text.x = element_text(size = 12, face = "bold"),
  #        strip.text.y = element_text(size = 12, face = "bold"))

p2 <- ggarrange(p2, ggplot() + theme_void(), common.legend = TRUE)

## plotting figure 2
fig2 <- ggarrange(p1, p2, common.legend=TRUE, nrow=2)

## saving
ggsave(fig2, filename="../Revisions/PDF_figures/Figure2.pdf",
       height=6, width=10)

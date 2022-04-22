##################################### model 3

library(coin)
library(rstatix)
library(ggpubr)
library(StatisticalModels)
library(ggplot2)
library(dplyr)
library(viridis)
library(scales)

Predicts_site_info <- read.csv("../../Results/Predicts_site_info_ER.csv")
source("Functions.R")

# loading simulation results and testing whether FDis at each site differs from 0-expectations -- function to run the tests
Wilcox.test.site <- function(Data, Studies) {

  pathReg <- "../../Results/dbFD_indices/Sim_Region/Sim_results/Standardised/"

  Results.list <- list()

  for (j in 1:length(Studies)){

    print(paste("Study", j))

    #SimEcoreg <- readRDS(paste0(pathEcoreg, Studies[j], ".rds"))
    SimReg <- readRDS(paste0(pathReg, Studies[j], ".rds"))

    # filter for sites present in results
    Data$SS <- as.character(Data$SS)
    Sites <- as.character(unique(Data$SSBS[Data$SS==Studies[j]]))
    SimReg <- SimReg[Sites]

    Results <- as.data.frame(Sites)

    for (i in 1:length(Sites)) {

      # null expectation distribution for both regional and ecoregional simulations
      Reg.simFRic <- SimReg[[Sites[i]]]$FRic
      Reg.simFDis <- SimReg[[Sites[i]]]$FDis

      # empirical values - regional simulations - first elements of the simulation results
      Reg.empFRic <- Reg.simFRic[1]
      Reg.empFDis <- Reg.simFDis[1]

      #print(Reg.empFRic)

      if(is.null(Reg.empFRic)) {next}

      # wilcoxon signed rank test to test whether null distribution is significantly higher/lower/differs from observed value
      if(!is.na(Reg.empFRic) && length(unique(Reg.simFRic))>1) {

        St.FRic.Greater <- wilcox.test(Reg.simFRic, mu = Reg.empFRic, alternative = "greater", conf.int = TRUE)
        Results$pval.FRic.Greater[Results$Sites==Sites[[i]]] <- St.FRic.Greater$p.value
        Results$pval.FRic.Lower[Results$Sites==Sites[[i]]] <- St.FRic.Greater$conf.int[1]
        Results$pval.FRic.Upper[Results$Sites==Sites[[i]]] <- St.FRic.Greater$conf.int[2]
        Results$FRic[Results$Sites==Sites[[i]]] <- Reg.empFRic
        }

      else{

        Results$pval.FRic.Greater[Results$Sites==Sites[[i]]] <- NA
        Results$pval.FRic.Lower[Results$Sites==Sites[[i]]] <- NA
        Results$pval.FRic.Upper[Results$Sites==Sites[[i]]] <- NA
        Results$FRic[Results$Sites==Sites[[i]]] <- Reg.empFRic

        }

      if(!is.na(Reg.empFDis) && length(unique(Reg.simFDis))>1) {
        St.FDis.Greater <- wilcox.test(Reg.simFDis, mu = Reg.empFDis, alternative = "greater", conf.int = TRUE)

        Results$pval.FDis.Greater[Results$Sites==Sites[[i]]] <- St.FDis.Greater$p.value
        Results$pval.FDis.Lower[Results$Sites==Sites[[i]]] <- St.FDis.Greater$conf.int[1]
        Results$pval.FDis.Upper[Results$Sites==Sites[[i]]] <- St.FDis.Greater$conf.int[2]
        Results$FDis[Results$Sites==Sites[[i]]] <- Reg.empFDis

      } else{
        Results$pval.FDis.Greater[Results$Sites==Sites[[i]]] <- NA
        Results$pval.FDis.Lower[Results$Sites==Sites[[i]]] <- NA
        Results$pval.FDis.Upper[Results$Sites==Sites[[i]]] <- NA
        Results$FDis[Results$Sites==Sites[[i]]] <- Reg.empFDis

      }

    }

    Results.list[[j]] <- Results
  }

  #browser()
  Results.final <- data.table::rbindlist(Results.list)

  return(Results.final)

}

## loading main results dataset across terrestrial vertebrates
Regional_std <- read.csv("../../Results/dbFD_indices_for_analysis/dbFD_region_std_max.csv") %>%
  Clean()

# replace FRic and FDis values for SR=1 (not possible to estimate FRic or FDis in that case)
Regional_std$FRic[Regional_std$SR==1] <- NA
Regional_std$FDis[Regional_std$SR==1] <- NA
Regional_std$median_sim_FRic[Regional_std$SR==1] <- NA
Regional_std$median_sim_FDis[Regional_std$SR==1] <- NA

# test at site-level
Studies <- unique(Regional_std$SS)
ResultsT <- Wilcox.test.site(Regional_std, Studies)

# dataset containing Wilcoxon signed rank tests results (across vertebrates)
colnames(ResultsT)

# merge with site information
colnames(ResultsT)[1] <- "SSBS"
ResultsT <- merge(ResultsT, Predicts_site_info, by="SSBS")
ResultsT <- Clean(ResultsT)

colnames(ResultsT)
ResultsT <- merge(ResultsT, Regional_std[, c("Biome", "SSBS")], by="SSBS")

# interpreting Wilcoxon tests outputs.
# we focus on the test "greater" because we want to check whether the median of the 0-distribution is
# significantly greater than the actual FDis value.
# 1 if null distribution is significantly higher than value, 0 otherwise
ResultsT$OutcomeFDis <- NA
ResultsT$OutcomeFDis[ResultsT$pval.FDis.Greater<0.05] <- 1
ResultsT$OutcomeFDis[ResultsT$pval.FDis.Greater>=0.05] <- 0

# this is equivalent (using the lower confidence interval)
#ResultsT$OutcomeFDis2 <- ifelse(ResultsT$FDis<ResultsT$pval.FDis.Lower, 1,0)

ResultsT %>%  colnames()

# checking levels
ResultsT$Predominant_land_use
ResultsT$Use_intensity
ResultsT$Biome

ResultsT <- subset(ResultsT, !is.na(OutcomeFDis))

# model selection (binomial model)
model3 <- GLMERSelect(modelData = ResultsT,
                      responseVar = "OutcomeFDis",
                      fitFamily = "binomial",
                      fixedFactors = c("Predominant_land_use", "Use_intensity", "Biome"),
                      fitInteractions=TRUE,
                      randomStruct = "(1|SS)+(1|SSB)",
                      verbose = TRUE)

model3 <- glmer(formula = OutcomeFDis ~ Predominant_land_use + Use_intensity + Biome +
                  Predominant_land_use:Use_intensity +
                  Predominant_land_use:Biome + (1|SS/SSB),
                data = ResultsT, family = "binomial",
                control = glmerControl(optCtrl = list(maxfun = 75000)))

###################################################################################

# plotting predictions
Newdata <- expand.grid(Predominant_land_use=levels(ResultsT$Predominant_land_use),
                       Use_intensity=levels(ResultsT$Use_intensity),
                       Biome=levels(Regional_std$Biome),
                       OutcomeFDis=0)


Predict_effects <- function(Newdata, Model, rescale) {

  logit2prob <- function(logit){
    prob <- exp(logit) / (1 + exp(logit))
    return(prob)
  }

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
    y <- logit2prob(y)

    # rescaling
    if(rescale){
      # initialisation
      seq <- 1:8
      y[seq] <- y[seq]-y[seq[1]]
      # for loop to rescale all values
      for(i in 1:(nrow(Newdata)/8-1)){
        seq <- seq + 8
        y[seq] <- y[seq]-y[seq[1]]
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

preds_probs <- Predict_effects(Newdata, model3, TRUE)

preds_probs$Median[preds_probs$Predominant_land_use=="Mature secondary vegetation" & preds_probs$Use_intensity=="Intense use"] <- NA
preds_probs$Lower[preds_probs$Predominant_land_use=="Mature secondary vegetation" & preds_probs$Use_intensity=="Intense use"] <- NA
preds_probs$Upper[preds_probs$Predominant_land_use=="Mature secondary vegetation" & preds_probs$Use_intensity=="Intense use"] <- NA

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
cbPalette <- c("#000000", "#CC0000", "#0000CC")

preds_probs$Expected <- ifelse(preds_probs$Lower>0, "Higher relative probability of occurrence \nof functional under-dispersion\n",
                               ifelse(preds_probs$Upper<0, "Lower relative probability of occurrence \nof functional under-dispersion\n",
                                      "Expected given local species richness\n"))

preds_probs <- subset(preds_probs, !is.na(Expected))

preds_probs$Biome <- factor(preds_probs$Biome, labels=c("(a) Temperate", "(b) Tropical"))

figure4 <- ggplot(preds_probs, aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity, col=Expected)) +
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
  GGPoptions +
  scale_colour_manual(values=cbPalette, name="") + 
  facet_wrap(~Biome)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ylab("Relative effect in p(FDis<FDisnull)")

ggsave(figure4, filename="../Revisions/PDF_figures/Figure4.pdf",
       height=4, width=12)

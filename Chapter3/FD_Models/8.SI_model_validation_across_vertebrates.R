## Model validations (for the models across all terrestrial vertebratres)

library(StatisticalModels)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(viridis)
library(scales)
source("Functions.R")

##############################################################################################################
## Across all vertebrates: robustness to variation in imputed values

Vertebrates_8 <- readRDS("../../Results/dbFD_indices/8_sets.rds")
Vertebrates_8 <- lapply(Vertebrates_8, Clean)
Vertebrates_8 <- lapply(Vertebrates_8, function(x){
  x$FRic[x$SR==1] <- NA
  x$FDis[x$SR==1] <- NA
  # tranform values
  x$asin_sqrt_FRic <- asin(sqrt(x$FRic))
  x$asin_sqrt_FDis <- asin(sqrt(x$FDis))
  # filter out NA use intensity
  x <- x %>%
    filter(!is.na(Use_intensity))
  return(x)
})

# check levels
Vertebrates_8[[1]]$Predominant_land_use
Vertebrates_8[[1]]$Use_intensity

# add tropical versus temperate divide
Verts8 <- read.csv("../../Results/dbFD_indices_for_analysis/dbFD_region_std_max.csv") %>%
  Clean() %>%
  filter(!is.na(Use_intensity)) %>%
  dplyr::select(SSBS, Biome)

Vertebrates_8 <- lapply(Vertebrates_8, function(x){
  x <- left_join(x, Verts8, by="SSBS")
  return(x)
})

## running model 1a and 1b and predictions for each set of metrics

Models_1 <- function(set, N) {

  model1a <- lme4::lmer(asin_sqrt_FRic~Predominant_land_use+Use_intensity+Biome+
                          Predominant_land_use:Use_intensity+
                          Predominant_land_use:Biome+
                          (1|SS/SSB), data=set)

  model1b <- lme4::lmer(asin_sqrt_FDis~Predominant_land_use+Use_intensity+Biome+
                          Predominant_land_use:Use_intensity+
                          (1|SS/SSB), data=set)

  Newdata <- expand.grid(Predominant_land_use=levels(set$Predominant_land_use),
                         Use_intensity=levels(set$Use_intensity),
                         Biome=levels(set$Biome),
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

  preds_fric$Metric <- "FRic"
  preds_fdis$Metric <- "FDis"

  preds_fric$set <- N
  preds_fdis$set <- N

  return(rbind(preds_fric, preds_fdis))
}

Results_8_models <- list()

for(i in 1:8) {
  Results_8_models[[i]] <- Models_1(Vertebrates_8[[i]],i)
}

Results_8_models <- data.table::rbindlist(Results_8_models)


## plotting (Supplementary Information, Figure S20)
Limits <- c("Primary vegetation",
            "Mature secondary vegetation",
            "Intermediate secondary vegetation",
            "Young secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")

Labels=c("PV", "MSV", "ISV", "YSV", "PF", "PA", "CR", "UR")

pFRic <- ggplot(Results_8_models[Results_8_models$Metric=="FRic",],
                aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, col=Use_intensity, group=interaction(set,Use_intensity))) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.8), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.8)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=c("black", "red", "orange")) +
  GGPoptions +
  facet_wrap(~Biome, nrow=2)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

pFDis <- ggplot(Results_8_models[Results_8_models$Metric=="FDis" & Results_8_models$Biome=="Temperate",],
                aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, col=Use_intensity, group=interaction(set,Use_intensity))) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.8), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.8)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=c("black", "red", "orange")) +
  GGPoptions +
  #facet_wrap(~Biome, nrow=2)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(b) FDis (Temperate; similar effects for tropical)") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

p_fig_SI_20 <- ggarrange(pFRic, pFDis, nrow = 2, heights = c(0.65, 0.35))
ggsave(p_fig_SI_20, filename="../Revisions/PDF_figures/Figure_SI_20.pdf",
       height=7, width=9)

rm(pFRic, pFDis)

##############################################################################################################
## Across all vertebrates: robustness to resampling in primary vegetation sites (using the 8th set)

Verts_8 <- Vertebrates_8[[8]]

Resample <- function(results, set_n) {

  ## resample primary vegetation sites (within each biome & use intensity) -- sample size fixed at 50
  resample_function <- function(res, Useintensity, Region){
    res <- res %>%
      dplyr::filter(Predominant_land_use=="Primary vegetation", Biome==Region, Use_intensity==Useintensity) %>%
      dplyr::sample_n(size=50, replace=FALSE)
    return(res)
  }

  Resample_tropical_minimal <- resample_function(results, "Minimal use", "Tropical")
  Resample_tropical_light <- resample_function(results, "Light use", "Tropical")
  Resample_tropical_intense <- resample_function(results, "Intense use", "Tropical")

  Resample_temperate_minimal <- resample_function(results, "Minimal use", "Temperate")
  Resample_temperate_light <- resample_function(results, "Light use", "Temperate")
  Resample_temperate_intense <- resample_function(results, "Intense use", "Temperate")

  results <- results %>%
    dplyr::filter(Predominant_land_use!="Primary vegetation")

  results <- rbind(results,
                   Resample_temperate_minimal, Resample_temperate_light, Resample_temperate_intense,
                   Resample_tropical_minimal, Resample_tropical_light, Resample_tropical_intense)

  results$set <- set_n

  ## run function model
  results_models <- Models_1(results, set_n)
  return(results_models)
  }

## running the resample function 20 times
Resample_results <- list()
for(i in 1:20){
  Resample_results[[i]] <- Resample(Verts_8, i)
  print(i)
}

Resample_results <- data.table::rbindlist(Resample_results)

## plotting figure for supporting information #S21
pFRic <- ggplot(Resample_results[Resample_results$Metric=="FRic",],
                aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, col=Use_intensity, group=interaction(set,Use_intensity))) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 1), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 1)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=c("black", "red", "orange")) +
  GGPoptions +
  facet_wrap(~Biome, nrow=2)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(a) FRic") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

pFDis <- ggplot(Resample_results[Resample_results$Metric=="FDis" & Resample_results$Biome=="Temperate",],
                aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, col=Use_intensity, group=interaction(set,Use_intensity))) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") +
  geom_hline(yintercept = 0, col="black", linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 1), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 1)) +
  scale_x_discrete(limits=Limits, labels=Labels) + xlab("") +
  scale_colour_manual(values=c("black", "red", "orange")) +
  GGPoptions +
  #facet_wrap(~Biome, nrow=2)+
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  ggtitle("(b) FDis (Temperate; similar effects for tropical)") +
  guides(colour=FALSE) +
  ylab("Relative effect (%)")

p_fig_SI_21 <- ggarrange(pFRic, pFDis, nrow = 2, heights = c(0.65, 0.35))
ggsave(p_fig_SI_21, filename="../Revisions/PDF_figures/Figure_SI_21.pdf",
       height=7, width=9)

rm(pFRic, pFDis)

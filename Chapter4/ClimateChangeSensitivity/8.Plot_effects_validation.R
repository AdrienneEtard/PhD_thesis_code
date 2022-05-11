## plotting PGLS predictions (sensitivity ~ traits) ## for models fitted on complete trait data only

library(caper)
library(StatisticalModels)
library(ggplot2)
library(viridis)
library(ggpubr)
library(scales)
library(ggthemes)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


## loading model outputs
Mammals <- readRDS("../Results/10.PGLS_complete_models_results/Mammals.rds")
Reptiles <- readRDS("../Results/10.PGLS_complete_models_results/Reptiles.rds")
Birds <- readRDS( "../Results/10.PGLS_complete_models_results/Birds.rds")
Amphibians <- readRDS("../Results/10.PGLS_complete_models_results/Amphibians.rds")

library(stargazer)
summary(Amphibians)
summary(Mammals)
summary(Reptiles)
summary(Birds)

stargazer(summary(Mammals)$coefficients, summary = FALSE)
stargazer(summary(Birds)$coefficients, summary = FALSE)
stargazer(summary(Amphibians)$coefficients, summary = FALSE)
stargazer(summary(Reptiles)$coefficients, summary = FALSE)


## sample sizes
Mammals$data$data %>%  nrow()
Birds$data$data %>%  nrow()
Amphibians$data$data %>%  nrow()
Reptiles$data$data %>%  nrow()

# ## verifying qqplots (diagnostics)
# pdf("G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/PGLS_diag_mammals.pdf",
#     width=7, height=6, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1.3, oma=c(0.5,0.5,1,0.5), mar=c(4,4,4,4))
# par(mfrow=c(2,2)); plot(Mammals)
# mtext("Mammals: diagnostic plots for the PGLS",
#       side = 3,
#       line = -.1,
#       outer = TRUE,
#       cex = 1.2,
#       adj=0)
# dev.off()
# 
# pdf("G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/PGLS_diag_reptiles.pdf",
#     width=7, height=6, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1.3, oma=c(0.5,0.5,1,0.5), mar=c(4,4,4,4))
# par(mfrow=c(2,2)); plot(Reptiles)
# mtext("Reptiles: diagnostic plots for the PGLS",
#       side = 3,
#       line = -.1,
#       outer = TRUE,
#       cex = 1.2,
#       adj=0)
# dev.off()
# 
# pdf("G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/PGLS_diag_birds.pdf",
#     width=7, height=6, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1.3, oma=c(0.5,0.5,1,0.5), mar=c(4,4,4,4))
# par(mfrow=c(2,2)); plot(Birds)
# mtext("Birds: diagnostic plots for the PGLS",
#       side = 3,
#       line = -.1,
#       outer = TRUE,
#       cex = 1.2,
#       adj=0)
# dev.off()
# 
# pdf("G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/PGLS_diag_amphibians.pdf",
#     width=7, height=6, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1.3, oma=c(0.5,0.5,1,0.5), mar=c(4,4,4,4))
# par(mfrow=c(2,2)); plot(Amphibians)
# mtext("Amphibians: diagnostic plots for the PGLS",
#       side = 3,
#       line = -.1,
#       outer = TRUE,
#       cex = 1.2,
#       adj=0)
# dev.off()




############################################################################################################################################################
## plotting predictions -- for continuous traits -- functions
## plotting prediction intervals around the fitted lines through resampling coefficients


Predict_continuous <- function(Trait, Model) {
  
  colnames(Model$data$data)[colnames(Model$data$data)==Trait] <- "Var"
  Newdata <- data.frame(Var=seq(min(Model$data$data$Var), max(Model$data$data$Var), length.out=100)) 
  
  if(Trait=="log10_Litter_size"){
    Newdata <- data.frame(Var=seq(min(Model$data$data$Var), 2, length.out=100)) 
  }
  
  if(Trait=="log10_Body_mass"){
    Newdata <- data.frame(Var=seq(min(Model$data$data$Var), 3, length.out=100)) 
  }
  
  
  Summary <- as.data.frame(summary(Model)$coefficients)
  Summary$Upper <- Summary$Estimate + 1.96*Summary$`Std. Error`
  Summary$Lower <- Summary$Estimate - 1.96*Summary$`Std. Error`
  
  Resample_predict <- function(Summary, Trait, Newdata){
    
    # significance
    pValsCoefs <- Summary$`Pr(>|t|)`[grepl(Trait, rownames(Summary))]
    
    # resample coefficients -- intercept
    Intercept <- runif(min=Summary$Lower[1], max=Summary$Upper[1], n=1)
    
    # resample coefficients -- slopes
    Sub <- subset(Summary, grepl(Trait, rownames(Summary)))
    Coefs <- apply(Sub, 1, FUN = function(x){runif(min=x[6], max=x[5], n=1)})
    
    # predict
    if(pValsCoefs[1]<=0.05){
      Newdata$preds <- Intercept+Coefs[1]*Newdata$Var
    }
    
    return(Newdata$preds)
  }
  
  ## prediction intervals and estimates from resampling
  
  preds <- sapply(1:1000, FUN=function(i){Resample_predict(Summary, Trait, Newdata)}) %>% 
    as.data.frame()
  
  preds <- data.frame(#Median=apply(X=preds, 1, FUN=median),
    Upper=apply(X=preds, 1, FUN=quantile, probs=0.975),
    Lower=apply(X=preds, 1, FUN=quantile, probs=0.025))
  
  preds$var <- Newdata$Var

  # significance
  pValsCoefs <- Summary$`Pr(>|t|)`[grepl(Trait, rownames(Summary))]
  
  Coefs <- Summary$Estimate[grepl(Trait, rownames(Summary))]
  Intercept <- Summary$Estimate[1]
  
  # predict
  if(pValsCoefs[1]<=0.05){
    preds$Median <- Intercept+Coefs[1]*preds$var
  }

  ##
  return(preds)
}

Plot_continuous <- function(Trait, Classes, Models) {
  results <- list()
  for (i in 1:length(Classes)) {
    Model <- Models[[i]]
    res <- Predict_continuous(Trait, Model)
    res$Class <- Classes[i]
    results[[i]] <- res
  }
  Predictions <- data.table::rbindlist(results)
  p <- ggplot(Predictions, aes(var, Median, col=Class, ymin=Lower, ymax=Upper, group=Class, fill=Class)) +
    geom_line() +
    geom_ribbon(aes(ymin = Lower, ymax = Upper),
                alpha=0.1,
                linetype="dashed",col=NA) +
    ylab("Climate-change sensitivity (log10)") +
    xlab(Trait)+ scale_colour_colorblind() + scale_fill_colorblind() + GGPoptions
  print(p)
  Predictions$Trait <- Trait
  return(Predictions)
}

summary(Amphibians)
summary(Reptiles)
summary(Mammals)
summary(Birds)

## plotting
to_plot <- rbind(Plot_continuous("log10_Body_mass_g", c("Mammals"), list(Mammals)),
                 #Plot_continuous("log10_Lifespan_proxy", c("Birds", "Reptiles"), list(Birds, Reptiles)),
                 Plot_continuous("log10_Litter_size", c("Birds", "Mammals"), list(Birds, Mammals)),
                 Plot_continuous("sqrt_Habitat_breadth_IUCN", c("Birds", "Mammals", "Reptiles"), list(Birds, Mammals, Reptiles)),
                 Plot_continuous("log10_Range_area", c("Amphibians","Birds", "Mammals"), list(Amphibians, Birds, Mammals)))
                 #Plot_continuous("sqrt_Diet_breadth", c("Mammals"), list(Mammals)))

# to_plot$Trait <- factor(to_plot$Trait, levels=c("log10_Body_mass_g", "log10_Lifespan_proxy", "log10_Litter_size",
#                                                 "log10_Range_area", "sqrt_Habitat_breadth_IUCN", "sqrt_Diet_breadth"), 
#                         labels=c("Body mass (g, log10)", "Lifespan (days, log10)", "Litter/clutch size (log10)",
#                                  "Range area (km2, log10)", "Habitat breadth (sqrt)", "Diet breath (sqrt)"))


to_plot$Trait <- factor(to_plot$Trait, levels=c( "log10_Range_area", "sqrt_Habitat_breadth_IUCN",
                                                 "log10_Body_mass_g", "log10_Litter_size"),
                        labels=c( "(a) Range area (km2, log10)", "(b) Habitat breadth (square-root)", 
                                  "(c) Body mass (g, log10)", "(d) Litter/clutch size (log10)"))

p <- ggplot(to_plot, aes(var, Median, col=Class, ymin=Lower, ymax=Upper, group=Class, fill=Class)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lower, ymax = Upper),
              alpha=0.1,
              linetype="dashed",col=NA) +
  ylab("Climate-change sensitivity (log10)") +
  xlab("Trait value") + GGPoptions + 
  facet_wrap(~Trait, scales = "free", nrow =2) +
  theme(strip.text = element_text(face="bold")) + theme(legend.position = "right") +
  viridis::scale_color_viridis(discrete = TRUE,option = "B",  end = 0.7) +
  viridis::scale_fill_viridis(discrete = TRUE,option = "B",  end = 0.7) +
  theme(axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
  theme(legend.position = "top", legend.title=element_blank())


print(p)
ggsave(p, filename="E:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Figures/PGLS_predictions_continuous_validation.pdf",
       width = 8, height = 6)


#ggsave(p, filename="../Results/6.Sensitivity_models_predictions/Continuous_traits.pdf", width=9.5, height=4)

# pmass <- Plot_continuous("log10_Body_mass_g", c("Amphibians", "Birds", "Mammals", "Reptiles"), list(Amphibians, Birds, Mammals, Reptiles)) + xlab("Body mass (log10, g)")  
# plifespan <- Plot_continuous("log10_Lifespan_proxy", c("Birds", "Reptiles"), list(Birds, Reptiles)) + xlab("Lifespan proxy (log10, days)")  
# plcs <- Plot_continuous("log10_Litter_size", c("Amphibians","Birds", "Mammals"), list(Amphibians, Birds, Mammals))+ xlab("Litter/clutch size (log10)")
# pHB <- Plot_continuous("sqrt_Habitat_breadth_IUCN", c("Amphibians","Birds", "Mammals", "Reptiles"), list(Amphibians, Birds, Mammals, Reptiles)) + xlab("Hanitat breadth (sqrt)")
# pRS <- Plot_continuous("log10_Range_area", c("Amphibians","Birds", "Mammals", "Reptiles"), list(Amphibians, Birds, Mammals, Reptiles)) + xlab("Range area (log10, km_sq)")
# pDB <- Plot_continuous("sqrt_Diet_breadth", c("Mammals"), list(MammalsDiet)) + xlab("Diet breadth (sqrt)")

# arranging plots
# ggarrange(pmass, plifespan, plcs, pHB, pRS, pDB)


############################################################################################################################################################
## plotting effect sizes  -- for categorical traits

GetEffectSizes <- function(Model, Trait) {
  colnames(Model$data$data)[colnames(Model$data$data)==Trait] <- "Var"
  Summary <- as.data.frame(summary(Model)$coefficients)
  Sub <- subset(Summary, grepl(Trait, rownames(Summary)))
  Sub <- rbind(Summary[1,], Sub)
  Sub$Estimate[2:nrow(Sub)] <- Sub$Estimate[1]+Sub$Estimate[2:nrow(Sub)] 
  Sub$Upper <- Sub$Estimate + Sub$`Std. Error`*1.96
  Sub$Lower <- Sub$Estimate - Sub$`Std. Error`*1.96
  Sub$Lower[1] <- NA
  Sub$Upper[1] <- NA
  Sub$hline <- Sub$Estimate[1]
  return(Sub)
}

## Diet
# D1 <- GetEffectSizes(Amphibians, "Primary_diet") %>%
#   mutate(Class="Amphibians")
D2 <- GetEffectSizes(Birds, "Primary_diet") %>%
  mutate(Class="Birds")
D3 <- GetEffectSizes(Mammals, "Primary_diet") %>%
  mutate(Class="Mammals")
# D4 <- GetEffectSizes(Reptiles, "Primary_diet") %>%
#   mutate(Class="Reptiles")
Diet <- rbind(D2, D3)

D_NA_amph <- D2
D_NA_amph$Class <- "Amphibians"
D_NA_amph$Estimate <- NA
D_NA_amph$Lower <- NA
D_NA_amph$Upper <- NA
D_NA_amph$hline <- NA

D_NA_rep <- D2
D_NA_rep$Class <- "Reptiles"
D_NA_rep$Estimate <- NA
D_NA_rep$Lower <- NA
D_NA_rep$Upper <- NA
D_NA_rep$hline <- NA

Diet <- rbind(Diet, D_NA_amph, D_NA_rep)

Diet$Trait <- c("FR|NE","IN", "OM", "PL|SE", "VE", "FR|NE", "IN", "OM", "PL|SE", "VE")
Diet$facet <- "Diet"
#Diet$hline <- ifelse(grepl("Intercept", rownames(Diet)), Diet$Estimate, NA) 
Diet$Trait <- factor(Diet$Trait,
                     levels=c( "FR|NE","PL|SE","IN", "VE","OM"), 
                     labels = c("Fruit/nectar", "Plant/seed", "Invertebrate", "Vertebrate","Omnivore"))

Diet$facet <- "(c) Diet"
pDiet <- ggplot(Diet, aes(Trait, Estimate, ymin=Lower, ymax=Upper, group=interaction(Trait, Class), col=Class, shape=Class)) +
  geom_point(position=position_dodge(width=.6), size=2.5) +
  geom_errorbar(width=.3, size=0.5, position=position_dodge(width=.6)) +
  scale_colour_colorblind() + scale_fill_colorblind() + GGPoptions +
  facet_wrap(~facet, scales = "free") + xlab("") +
  theme(strip.text = element_text(face="bold")) + 
  theme(panel.spacing.x = unit(0, "lines")) +
  #facet_wrap(~Class, scales="free", ncol=1) +
  geom_hline(aes(yintercept=hline, col=Class), lty="dashed") +
  coord_flip()  + #ggtitle("Diel activity") +
  scale_y_continuous(n.breaks=4)+ theme(axis.title.x = element_text(vjust=-0.8)) +
  viridis::scale_color_viridis(discrete = TRUE,option = "B",  end = 0.7)  +
  theme(axis.text.y = element_text(size = 11, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
  theme(legend.title = element_blank())

print(pDiet)


## Specialisation
S1 <- GetEffectSizes(Amphibians, "Specialisation") %>%
  mutate(Class="Amphibians")
S2 <- GetEffectSizes(Birds,"Specialisation") %>%
  mutate(Class="Birds")
S3 <- GetEffectSizes(Mammals,"Specialisation") %>%
  mutate(Class="Mammals")
S4 <- GetEffectSizes(Reptiles,"Specialisation") %>%
  mutate(Class="Reptiles")
Spe <- rbind(S1, S2, S3, S4)
Spe$Trait <- rep(c("Artificial habitats user", "Natural habitat specialist"), 4)
Spe$facet <- "Degree of specialisation"
#Spe$hline <- ifelse(grepl("Intercept", rownames(Spe)), Spe$Estimate, NA) 
Spe$Trait <- factor(Spe$Trait, labels = c("Artificial\nhab. user", "Nat. hab.\nspecialist"))

Spe$facet <- "(a) Artificial habitat use"
pSpe <- ggplot(Spe, aes(Trait, Estimate, ymin=Lower, ymax=Upper, group=interaction(Trait, Class), col=Class, shape=Class)) +
  geom_point(position=position_dodge(width=.6), size=2.5) +
  geom_errorbar(width=.3, size=0.5, position=position_dodge(width=.6)) +
  scale_colour_colorblind() + scale_fill_colorblind() + GGPoptions +
  facet_wrap(~facet, scales = "free") + xlab("") +
  theme(strip.text = element_text(face="bold")) + 
  theme(panel.spacing.x = unit(0, "lines")) +
  #facet_wrap(~Class, scales="free", ncol=1) +
  geom_hline(aes(yintercept=hline, col=Class), lty="dashed") +
  coord_flip()  + #ggtitle("Diel activity") +
  scale_y_continuous(n.breaks=4)+ theme(axis.title.x = element_text(vjust=-0.8)) +
  viridis::scale_color_viridis(discrete = TRUE,option = "B",  end = 0.7)  +
  theme(axis.text.y = element_text(size = 11, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"))) +
  theme(legend.title = element_blank())


print(pSpe)

## Diel activity
DA1 <- GetEffectSizes(Amphibians,"Diel_activity") %>%
  mutate(Class="Amphibians")
DA2 <- GetEffectSizes(Birds,"Diel_activity") %>%
  mutate(Class="Birds")
DA3 <- GetEffectSizes(Mammals,"Diel_activity") %>%
  mutate(Class="Mammals")
DA4 <- GetEffectSizes(Reptiles,"Diel_activity") %>%
  mutate(Class="Reptiles")
DA <- rbind(DA1, DA2, DA3, DA4)
DA$Trait <- rep(c("Nocturnal", "Non-\nnocturnal"), 4)
DA$facet <- "Diel activity"
#DA$hline <- ifelse(grepl("Intercept", rownames(DA)), DA$Estimate, NA) 
DA$facet <- "(b) Diel activity"
pDA <- ggplot(DA, aes(Trait, Estimate, ymin=Lower, ymax=Upper, group=interaction(Trait, Class), col=Class, shape=Class)) +
  geom_point(position=position_dodge(width=.6), size=2.5) +
  geom_errorbar(width=.3, size=0.5, position=position_dodge(width=.6)) +
  scale_colour_colorblind() + scale_fill_colorblind() + GGPoptions +
  facet_wrap(~facet, scales = "free") + xlab("") +
  theme(strip.text = element_text(face="bold")) + 
  theme(panel.spacing.x = unit(0, "lines")) +
  #facet_wrap(~Class, scales="free", ncol=1) +
  geom_hline(aes(yintercept=hline, col=Class), lty="dashed") +
  coord_flip()  + #ggtitle("Diel activity") +
  scale_y_continuous(n.breaks=4)+ theme(axis.title.x = element_text(vjust=-0.8)) +
  viridis::scale_color_viridis(discrete = TRUE,option = "B",  end = 0.7)  +
  theme(axis.text.y = element_text(size = 11, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
        axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")))

print(pDA)

## arrange plots
PCat <- ggarrange(
  pSpe + ylab("Climate-change sensitivity (log10)"),
  pDA + ylab("Climate-change sensitivity (log10)"), 
  pDiet + ylab("Climate-change sensitivity (log10)"), common.legend = TRUE, legend="top")

ggsave(PCat, filename="E:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Figures/PGLS_predictions_categories_Validation.pdf",
       width = 8, height = 6)

# Plot <- 
# ggarrange( p + guides(colour=FALSE, fill=FALSE), PCat, 
#          common.legend = TRUE,  widths = c(0.60, 0.40), ncol=2)
# Plot
# 
# ggsave(Plot, filename="E:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Figures/PGLS_predictions2.pdf",
#        width = 11, height = 6)

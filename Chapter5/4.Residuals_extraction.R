library(StatisticalModels)
library(lmerTest)
library(dplyr)


## load data 
BMR_data <- read.csv("../Results/3.BMR_data_imputed.csv")
BMR_data$log_Body_mass <- log(BMR_data$Body_mass_g)
random_str <- "1|Class/Order/Family"

Predicts <- readRDS("D:/4.BMR/Data/PredictsVertebrates.rds")

## are there species in PREDICTS in those that don't have the family? (No.)
NoFam <- BMR_data$Species[is.na(BMR_data$Family)]
length(NoFam)
intersect(NoFam, unique(Predicts$Best_guess_binomial))

## remove data that don't have family in RMR data
BMR_data <- BMR_data %>% 
  filter(!is.na(Family))

## model for extraction of residuals
model1 <- lmer(log_BMR~log(Body_mass_g) + (1|Class/Order/Family), data=BMR_data)
BMR_data$Residual_BMR_log_log <- residuals(model1)
write.csv(BMR_data, "../Results/4.BMR_data_residuals.csv", row.names = FALSE)


GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

ggplot(BMR_data, aes(log(Body_mass_g), Residual_BMR_log_log)) +
  geom_point() + GGPoptions + 
  geom_hline(yintercept = 0, lty="dashed", col="darkgrey", size=1.3) + xlim(c(-1,25)) +
  ylab("Residual RMR (mL O2/hour, log)") + xlab("Body mass (log, g)") +
  scale_color_viridis_d(end=.7) + theme(legend.position = "bottom") +
geom_segment(x = 14, y = 0, xend = 14, yend = 2.2, col="red", size=1.3,
             arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(x = 14, y = 0, xend = 14, yend = -2.2, col="blue", size=1.3,
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(x=15, y=2, label="Positive deviations:
            \nRMR higher than expected\nfrom body mass\nand taxonomy", col="red", hjust=0, size=3.5)+
  geom_text(x=15, y=-2, label="Negative deviations:
            \nRMR lower than expected\nfrom body mass\nand taxonomy", col="blue", hjust=0, size=3.5)




#################################################################################################################################################
# # mod.1<- lmer(log_BMR~log(Body_mass_g)+(1|Class/Order/Family),data=BMR_data)
# # mod.2<- lmer(log_BMR~log(Body_mass_g)+ EV_1 +(1|Class/Order/Family),data=BMR_data)
# # mod.3<- lmer(log_BMR~log(Body_mass_g)+ EV_1 + EV_2  +(1|Class/Order/Family),data=BMR_data)
# # mod.4<- lmer(log_BMR~log(Body_mass_g)+ EV_1 + EV_2 + EV_3 +(1|Class/Order/Family),data=BMR_data)
# # mod.5<- lmer(log_BMR~log(Body_mass_g)+ EV_1 + EV_2 + EV_3 + EV_4 +(1|Class/Order/Family),data=BMR_data)
# # mod.6<- lmer(log_BMR~log(Body_mass_g)+ EV_1 + EV_2 + EV_3 + EV_4+ EV_5 +(1|Class/Order/Family),data=BMR_data)
# mod.6bis<- lmer(log_BMR~log(Body_mass_g)+ EV_1 + EV_2 + EV_3 + EV_4+ EV_5+ Thermoregulation +(1|Class/Order/Family),data=BMR_data)
# mod.6bis_alt<- lmer(log_BMR~log(Body_mass_g)+ EV_1 + EV_2 + EV_3 + EV_4+ EV_5 +(1|Class/Order/Family),data=BMR_data)
# mod_alt<- lmer(log_BMR~log(Body_mass_g)+ Thermoregulation +(1|Class/Order/Family),data=BMR_data)
# AIC(mod.6bis)
# AIC(mod.6bis_alt)
# AIC(mod_alt)
# mod.best<- lmer(log_BMR~log(Body_mass_g)+ EV_4+ EV_5+ Thermoregulation +(1|Class/Order/Family),data=BMR_data)
# 
# step_res <- step(mod.6bis)
# step_res <- step(mod.6bis_alt)
# 
# drop1(mod.6bis)
# final <- get_model(step_res)
# anova(final)
# 
# anova(mod.1, mod.2)
# anova(mod.2, mod.3)
# anova(mod.3, mod.4)
# anova(mod.4, mod.5)
# anova(mod.5, mod.6)
# 
# AIC(mod.6bis)
# AIC(mod.best)
# 
# hist(residuals(mod.6bis))
# hist(residuals(mod.best))
# hist(residuals(mod.7))
# 
# ## include the 1st four EV
# residuals_BMR <- residuals(mod.5)
# hist(residuals_BMR)
# plot(residuals_BMR~log(Full_dat$Body_mass_g), pch=19)

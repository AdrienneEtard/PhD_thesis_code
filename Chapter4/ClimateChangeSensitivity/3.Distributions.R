## range sizes versus sensitivity
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggthemes)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=11, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


setwd("F:/PhD/PhD_R_projects/5.Climatic_niche_space/Code/")
Index <- read.csv("../../Range_maps_work/Results/7_GenerateSpeciesIndices/SpeciesIndices.csv")

RS <- read.csv("../Data/Range_sizes.csv")
RS$Species <- paste0("sp", RS$Index)
RS <- RS %>% 
  filter(Range_area_sq_km!=0)

## match with sensitivity dataframe
Sensitivity_50 <- read.csv("../Results/2.CENFA_summary_dataframes/Vertebrates_50km.csv")
Sensitivity_10 <- read.csv("../Results/2.CENFA_summary_dataframes/Vertebrates_10km.csv")
Sensitivity_5 <- read.csv("../Results/2.CENFA_summary_dataframes/Vertebrates_5km.csv")

Sensitivity_50 <- left_join(Sensitivity_50, RS, by="Species")
Sensitivity_10 <- left_join(Sensitivity_10, RS, by="Species")
Sensitivity_5 <- left_join(Sensitivity_5, RS, by="Species")

## plot sensitivity against range size

# ## 50 km
# Sensitivity_50$Below_2500sq <- ifelse(Sensitivity_50$Range_area_sq_km<=(50*50), TRUE, FALSE)
# p50 <- ggplot(Sensitivity_50, aes(log10(Range_area_sq_km), log10(sensitivity), col=Below_2500sq)) + GGPoptions +
#   geom_point() + 
#   facet_wrap(~Class, ncol=4, scales="free_x")+ 
#   #ggtitle(expression(paste("CENFA sensitivity estimation; resolution:", " 50 km"^{2}))) +
#   theme(panel.margin = unit(0, "lines")) + 
#   xlab(expression(paste("Geographical range area (", "km"^{2}, ",  log"[10], ")"))) + 
#   ylab(expression(paste("Climate-change sensitivity ", " (log"[10], ")"))) +
#   guides(col=guide_legend(expression(paste("Geographical range area: <", " 2500 km"^{2})))) +
#   scale_color_colorblind() +
#   theme(legend.position = "top")
#   
# ggsave(p50, filename="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/CENFA_50km.pdf",
#        width = 8, height=3)
# 
# 
# ## 10 km
# Sensitivity_10$Below_100sq <- ifelse(Sensitivity_10$Range_area_sq_km<=(10*10), TRUE, FALSE)
# p10 <- ggplot(Sensitivity_10, aes(log10(Range_area_sq_km), log10(sensitivity), col=Below_100sq)) + GGPoptions +
#   geom_point() + 
#   facet_wrap(~Class, ncol=4, scales="free_x")+ 
#   #ggtitle(expression(paste("CENFA sensitivity estimation; resolution:", " 10 km"^{2}))) +
#   theme(panel.margin = unit(0, "lines")) + 
#   xlab(expression(paste("Geographical range area (", "km"^{2}, ",  log"[10], ")"))) + 
#   ylab(expression(paste("Climate-change sensitivity ", " (log"[10], ")"))) +
#   guides(col=guide_legend(expression(paste("Geographical range area: <", " 100 km"^{2})))) +
#   scale_color_colorblind() +
#   theme(legend.position = "top")
# 
# ggsave(p10, filename="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/CENFA_10km.pdf",
#        width = 8, height=3)
# 
# 
# 
# ## 5 km
# Sensitivity_5$RA_threshold <- ifelse(Sensitivity_5$Range_area_sq_km>=(10*10), "Above 100 km2", 
#                                      ifelse(Sensitivity_5$Range_area_sq_km<=(5*5), "Below 25 km2", "Between 25 and 100 km2"))
# 
# p5 <- ggplot(Sensitivity_5, aes(log10(Range_area_sq_km), log10(sensitivity), col=RA_threshold)) + GGPoptions +
#   geom_point() + 
#   facet_wrap(~Class, ncol=4, scales="free_x")+ 
#   #ggtitle(expression(paste("CENFA sensitivity estimation; resolution:", " 5 km"^{2}))) +
#   theme(panel.margin = unit(0, "lines")) + 
#   xlab(expression(paste("Geographical range area (", "km"^{2}, ",  log"[10], ")"))) + 
#   ylab(expression(paste("Climate-change sensitivity ", " (log"[10], ")"))) +
#   guides(col=guide_legend("Geographical range area:")) + scale_color_colorblind()+
#   theme(legend.position = "top") +
#   scale_color_manual(labels = c(expression(paste(">", " 100 km"^{2})),
#                                 expression(paste("<", " 25 km"^{2})),
#                                 expression(paste(" 25 km"^{2}, "> & <"," 100 km"^{2}))),
#                       values=c("black", "#E69F00", "#56B4E9"))
#                
# p5
# ggsave(p5, filename="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/CENFA_5km.pdf",
#        width = 8, height=3)


## plotting all together
colnames(Sensitivity_5)
colnames(Sensitivity_50)
colnames(Sensitivity_5)

Sensitivity_50$RA_threshold <- ifelse(Sensitivity_50$Range_area_sq_km>=(50*50), "Above 2500 km2",  "Below 2500 km2")
Sensitivity_50$Resolution <- "50 km2"
Sensitivity_50$Vline <- log10(2500)
Sensitivity_50$Vline2 <-NA
Sensitivity_10$RA_threshold <- ifelse(Sensitivity_10$Range_area_sq_km>=(10*10), "Above 100 km2",  "Below 100 km2")
Sensitivity_10$Resolution <- "10 km2"
Sensitivity_10$Vline <- log10(100)
Sensitivity_10$Vline2 <- NA
Sensitivity_5$Resolution <- "5 km2"
Sensitivity_5$Vline <- log10(25)
Sensitivity_5$Vline2 <- log10(100)

Sensitivity_allres <- rbind(Sensitivity_50, Sensitivity_10, Sensitivity_5)
Sensitivity_allres$Resolution <- factor(Sensitivity_allres$Resolution, 
                                        levels=c("50 km2", "10 km2", "5 km2"),
                                        labels = c(expression(paste("Resolution:", " 50 km"^{2})),
                                                   expression(paste("Resolution:", " 10 km"^{2})),
                                                   expression(paste("Resolution:", " 5 km"^{2}))))
pAllres<- 
ggplot(Sensitivity_allres, aes(log10(Range_area_sq_km), log10(sensitivity))) + GGPoptions +
  geom_point() + 
  facet_grid(Resolution~Class, scales="free_x", labeller = label_parsed)+ 
  theme(panel.margin = unit(0, "lines")) + 
  geom_vline(aes(xintercept=Vline), linetype="dashed", col="red") +
  geom_vline(aes(xintercept=Vline2), linetype="dashed", col="blue") +
  xlab(expression(paste("Geographical range area (", "km"^{2}, ",  log"[10], ")"))) + 
  ylab(expression(paste("Climate-change sensitivity ", " (log"[10], ")"))) +
  guides(col=guide_legend("Geographical range area:")) + 
  theme(legend.position = "top") +
  theme(strip.text = element_text(face="bold", size = 12)) 

ggsave(pAllres, filename="G:/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Habitat_vars_predict_vertebrate_sensitivity_better_than_lht/Supporting Information/Figures/CENFA_allres.pdf",
       width = 8, height=7)


###################################################################################################################################
## match with trait data 
Traits <- readRDS("../Data/Imputed_traits_transformed.rds")[[8]]
Sensitivity_5 <- left_join(Sensitivity_5, Index[, c("Index", "Binomial")])
colnames(Sensitivity_5)[8] <- "Best_guess_binomial"

Sensitivity_5_traits <- Sensitivity_5 %>% 
  dplyr::filter(Best_guess_binomial %in% Traits$Best_guess_binomial)
Sensitivity_5_traits <- left_join(Sensitivity_5_traits, Traits, by="Best_guess_binomial")

## plotting species that have range size below 100 square kilometers
Sensitivity_5_traits$Below_100sq <- ifelse(Sensitivity_5_traits$Range_area_sq_km<=(5*5*4), TRUE, FALSE)
ggplot(Sensitivity_5_traits, aes(log(Range_area_sq_km), log(sensitivity), col=Below_100sq)) +
  geom_point() + 
  facet_wrap(~Class.x)+ 
  #geom_vline(xintercept=log(5*5*3), lty="dashed") +
  geom_vline(xintercept=log(5*5*4), lty="dashed")

## getting rid of some species whose range size fall below 5*5 = 25 square kilometers
# Sensitivity_5_traits$Threshold <- 4.5-log(Sensitivity_5_traits$Range_area_sq_km)*0.5
# Sensitivity_5_traits$Below <- ifelse(log(Sensitivity_5_traits$sensitivity)<Sensitivity_5_traits$Threshold, TRUE, FALSE)

Sensitivity_5_traits %>%  group_by(Class.x, Below_100sq) %>%  summarise(C=n())

## filtering out species that have range size below 100 km square
Sensitivity_5_traits_filter <- Sensitivity_5_traits %>% 
  filter(Below_100sq==FALSE)

## sensitivity against body mass
ggplot(Sensitivity_5_traits_filter, aes(log10_Body_mass_g, log10(sensitivity))) +
  geom_point() + 
  facet_wrap(~Class.x)

## sensitivity against litter/clutch size
ggplot(Sensitivity_5_traits_filter, aes(log10_Litter_size, log10(sensitivity))) +
  geom_point() + 
  facet_wrap(~Class.x)

## sensitivity against lifespan
ggplot(Sensitivity_5_traits_filter, aes(log10_Lifespan_proxy, log10(sensitivity))) +
  geom_point() + 
  facet_wrap(~Class.x)

## sensitivity against habitat breadth
ggplot(Sensitivity_5_traits_filter, aes(sqrt_Habitat_breadth_IUCN, log10(sensitivity))) +
  geom_point() + 
  facet_wrap(~Class.x)

# ## sensitivity against diet breadth
ggplot(Sensitivity_5_traits_filter, aes(sqrt_Diet_breadth, log10(sensitivity))) +
  geom_point() +
  facet_wrap(~Class.x)

## sensitivity against trophic levels
ggplot(Sensitivity_5_traits_filter, aes(Trophic_level, log10(sensitivity))) +
  geom_boxplot() + 
  facet_wrap(~Class.x)

## sensitivity against diel activity patterns
ggplot(Sensitivity_5_traits_filter, aes(Diel_activity, log10(sensitivity))) +
  geom_boxplot() + 
  facet_wrap(~Class.x)

## sensitivity against specialisation
ggplot(Sensitivity_5_traits_filter, aes(Specialisation, log10(sensitivity))) +
  geom_boxplot() + 
  facet_wrap(~Class.x)

Sensitivity_5_traits %>%  group_by(Class.x) %>%  summarise(COunt=n())
Sensitivity_5_traits_filter %>%  group_by(Class.x) %>%  summarise(COunt=n())

Traits %>%  group_by(Class) %>%  summarise(COunt=n())

write.csv(Sensitivity_5_traits, "../Results/Traits_CENFA_RS_all.csv", row.names = FALSE)
write.csv(Sensitivity_5_traits_filter, "../Results/Traits_CENFA_RS_filtered.csv", row.names = FALSE)

## distributions
Sensitivity_5_traits_filter <- read.csv( "../Results/Traits_CENFA_RS_filtered.csv")

ggplot(Sensitivity_5_traits_filter, aes(log10(sensitivity))) +
  #geom_boxplot() + 
  geom_histogram(fill="darkblue") +
  facet_wrap(~Class.x, scales="free_y") + GGPoptions + xlab("Sensitivity (log10)")+
  theme(panel.spacing.y = unit(0, "lines")) +
  theme(strip.text = element_text(face="bold"))

Sensitivity_5_traits_filter %>%  group_by(Class.x) %>%  summarise(Count=n())

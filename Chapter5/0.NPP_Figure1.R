## code to make figure 1 (map, distributions and NPP)

library(dplyr)
library(lme4)
library(ggplot2)
library(StatisticalModels)
library(performance)
# library(ggmap)
# library(rnaturalearth)
# library(rnaturalearthdata)
library(ggpubr)

Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")

Labels=c("PV", "SV", "PF", "PA", "CR", "UR")

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

PredictsTraits <-  "../Results/PredictsTraits.rds" %>%  readRDS()

#############################################
## sample sizes and mapping sites
length(unique(PredictsTraits$SS))
length(unique(PredictsTraits$SSBS))
length(unique(PredictsTraits$SS[!is.na(PredictsTraits$Use_intensity) & !is.na(PredictsTraits$LandUse)]))

length(unique(PredictsTraits$Best_guess_binomial))
length(unique(PredictsTraits$Best_guess_binomial[PredictsTraits$Thermoregulation=="Endotherms"]))
length(unique(PredictsTraits$Best_guess_binomial[PredictsTraits$Thermoregulation=="Ectotherms"]))

Predicts_to_map <- unique(PredictsTraits[, c("SS", "SSBS","Longitude", "Latitude", "Predominant_land_use", "Use_intensity")])

Predicts_to_map$Use_intensity <- as.character(Predicts_to_map$Use_intensity)
Predicts_to_map$Use_intensity %>%  unique()
Predicts_to_map$Use_intensity <- ifelse(Predicts_to_map$Use_intensity=="Cannot decide", NA, Predicts_to_map$Use_intensity)
Predicts_to_map$Use_intensity <- factor(Predicts_to_map$Use_intensity, levels=c("Minimal use", "Light use", "Intense use"))
Predicts_to_map$Use_intensity %>%  levels()

Predicts_to_map$Predominant_land_use <- as.character(Predicts_to_map$Predominant_land_use)
Predicts_to_map$Predominant_land_use %>%  unique()
Predicts_to_map$Predominant_land_use[grepl("econdary vegetation", Predicts_to_map$Predominant_land_use)] <- "Secondary vegetation"
Predicts_to_map$Predominant_land_use <- factor(Predicts_to_map$Predominant_land_use,
                                               levels=c("Primary vegetation",
                                                        "Secondary vegetation",
                                                        "Plantation forest",
                                                        "Pasture",
                                                        "Cropland",
                                                        "Urban"))

Model_res <- readRDS("../Results/8.Occurrence_models/Model_occurrence_3_way_LU_UI_TL_Res.rds")

Predicts_to_map$Predominant_land_use %>% levels()
Predicts_to_map <- subset(Predicts_to_map, SSBS %in% unique(Model_res@frame$SSBS))
length(unique(Predicts_to_map$SSBS))

# mapping sites
world <- ne_countries(scale = "medium", returnclass = "sf")
GGPoptions2 <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family = "serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=8), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=8),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=15))

Predicts_to_map_2 <- Predicts_to_map %>% 
  filter(!is.na(Use_intensity)) %>% 
  group_by(SS, Predominant_land_use) %>%  summarise(MeanLong=mean(Longitude), MeanLat=mean(Latitude), Nsites=n())

map <- ggplot2::ggplot(data = world) + xlab("") + ylab("") +
  geom_sf(fill="grey", colour="grey") +
  geom_point(data=Predicts_to_map[!is.na(Predicts_to_map$Longitude) & 
                                    !is.na(Predicts_to_map$Predominant_land_use) &
                                    !is.na(Predicts_to_map$Use_intensity),],
             aes(y=Latitude, x=Longitude, size=Nsites), size=1.5, alpha=.5) +
  GGPoptions2 + theme(panel.grid.major = element_line(colour = 'transparent')) +
  theme(legend.title = element_blank()) + theme(legend.position = c(0.5, -0.125)) + ggtitle("(a)") 

# ggplot2::ggplot(data = world) + xlab("") + ylab("") +
#   geom_sf(fill="grey", colour="grey") +
#   geom_point(data=Predicts_to_map_2[!is.na(Predicts_to_map_2$MeanLong) & 
#                                       !is.na(Predicts_to_map_2$MeanLat),],
#              aes(y=MeanLat, x=MeanLong, size=Nsites), alpha=.5) +
#   theme_classic() + GGPoptions2 + theme(panel.grid.major = element_line(colour = 'transparent')) +
#   theme(legend.title = element_blank()) + theme(legend.position = c(0.5, -0.125)) + ggtitle("(a)") + 
#   facet_wrap(~Predominant_land_use)

#############################################
# mapping land uses distributions

Limits <- c("Primary vegetation",
            "Secondary vegetation",
            "Plantation forest",
            "Pasture",
            "Cropland",
            "Urban")

Labels=paste(Limits, c("(PV)", "(SV)", "(PF)", "(PA)", "(CR)", "(UR)"))
Labels=c("Primary", "Secondary", "Plantation", "Pasture", "Cropland", "Urban")

Counts <-Predicts_to_map[!is.na(Predicts_to_map$Predominant_land_use) &!is.na(Predicts_to_map$Use_intensity),] %>% 
  group_by(Predominant_land_use, Use_intensity) %>% 
  summarise(Count=n())

LUs <- 
ggplot(Counts, aes(Predominant_land_use, Count ,group=Use_intensity, fill=Use_intensity)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  GGPoptions + ylab("Number of sites") + xlab("") + 
  scale_x_discrete(limits=Limits, labels=Labels)+ 
  scale_fill_manual(values=c("#000000", "#A0A0A0", "#C6C6C6"), name="Use intensity") + 
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1))  + ggtitle("(b)")  +
  geom_text(aes(label = Count, x = Predominant_land_use, y = Count+5),
            position = position_dodge(width = 1), size=3, angle=90, hjust=-0.2) + ylim(c(0,1550))


# LUs <- ggplot(Predicts_to_map[!is.na(Predicts_to_map$Predominant_land_use) &!is.na(Predicts_to_map$Use_intensity),],
#               aes(Predominant_land_use, group=Use_intensity, fill=Use_intensity)) +
#   geom_bar(position = position_dodge()) +
#   GGPoptions + 
#   xlab("") + ylab("Number of sites")  + 
#   scale_x_discrete(limits=Limits, labels=Labels)+ 
#   scale_fill_manual(values=c("#000000", "orange", "#DB5C68FF"), name="Use intensity") + 
#   theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1))  + ggtitle("(b)") 
# #theme(legend.position = c(0.85, 0.85))


#############################################
# NPP plot

NPP <- "../Results/Predictions_Figure1_Meghan_rescale2.csv" %>% read.csv()
NPP$Use_intensity <- factor(NPP$Use_intensity, levels=c("Minimal use", "Light use", "Intense use"))
NPP$Ref <- rep(NPP$Median[NPP$Predominant_land_use=="Primary vegetation"], each=6)

p <- ggplot(NPP,
            aes(Predominant_land_use, Median, ymin = Lower, ymax = Upper, shape=Use_intensity)) +
  geom_rect(xmin=0, xmax=1.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=2.5, xmax=3.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=4.5, xmax=5.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  geom_rect(xmin=6.5, xmax=7.5,ymin=-Inf,ymax=Inf, fill="lightgrey", col=NA) +
  ylab("") + xlab("") +  
  #geom_hline(aes(yintercept = Ref,col=Use_intensity), linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_errorbar(width=.2, size=0.5, position=position_dodge(width = 0.7), stat="identity") +
  geom_point(size=2, position=position_dodge(width = 0.7)) +
  scale_x_discrete(limits=Limits, labels=Labels) +
  GGPoptions +
  #scale_colour_viridis_d(option="C") +
  ggtitle("") +
  theme(panel.spacing = unit(0, "lines")) +
  theme( strip.text.x = element_text(size = 12, face = "bold"),
         strip.text.y = element_text(size = 12, face = "bold")) +
  scale_colour_manual(values=c("#000000", "orange", "#DB5C68FF"), name="Use intensity") +
  theme(legend.position = "top") +
  scale_shape_manual(values=c(19, 17, 8), name="Use intensity") + 
  theme(axis.text.x = element_text(angle = 55, vjust = 1.05, hjust=1))+
  #guides(colour=FALSE) +
  ylab("NPP (% difference from primary)") + ggtitle("(c)") + theme(legend.position = "right")

p

bc <- ggarrange(LUs + theme(legend.position = "bottom"),
          p  + theme(legend.position = "bottom"), nrow = 1)

Map <- ggarrange(map, ggplot() + theme_void())
# LUplot <- ggarrange(LUs, ggplot() + theme_void(), widths = c(0.75,0.25))
Figure1 <- ggarrange(Map,
                     bc,
                     nrow = 2, legend = "right", heights =  c(0.40,0.60))
Figure1

# Figure1 <- ggarrange(map,
#                      LUs,
#                      p,
#                      nrow = 3, legend = "right", widths = c(0.75,1,1))
# Figure1
ggsave(Figure1, 
       filename="../Results/Figures/Figure1.pdf",
       width=8, height = 13)

ggsave(Figure1,
       filename="c:/Users/adrie/OneDrive/Desktop/Thesis/figures/Chapter5/Figure1.pdf",
       width=9, height=8.5)




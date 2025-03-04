## code to prepare trait data for further scripts, to calculate functional diversity indices
## check variance inflation factors
## transform and z-score

## For imputed datasets, select traits to keep with VIF stepwise selection procedure: traits should not be too correlated.
## Before this operation traits should be transformed and z-scored.

## VIF stepwise selection process: there are several procedures that exist - ex:
# 1. fit a linear model where all candidate traits are explanatory variable. The dependent variable is a dummy continuous variable
# generated by simulating a normally distributed variable (sampling the normal distribution)
# 2. Calculate VIF for all explanatory variable
# 3. Remove all variables for which VIF > threshold
# 4. Fit model again with the new set of variables, and repeat the procedure until all VIF score for all remaining varaibles are < threshold.

## Prepare trait datasets by class and also overall (either standardised and z-scored across or within classes)

library(dplyr)
library(pedometrics)
library(car)
require(pls)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(mctest)
library(PerformanceAnalytics)
library(Rnalytica)
library(stargazer)
source("Functions.R")
`%nin%` <- Negate(`%in%`)

## Load imputed datasets (8 sets)
Imputed <- readRDS("../../Results/Imputed_sets/List_of_8_sets.rds")

## Add class and change column name (lifespan proxy with which we work)
for(i in 1:8){
Imputed[[i]][["M"]]$Imputed.Dataset$Class <- "Mammals"
colnames(Imputed[[i]][["M"]]$Imputed.Dataset)[6] <- "Lifespan_proxy"

Imputed[[i]][["B"]]$Imputed.Dataset$Class <- "Birds"
colnames(Imputed[[i]][["B"]]$Imputed.Dataset)[6] <- "Lifespan_proxy"

Imputed[[i]][["R"]]$Imputed.Dataset$Class <- "Reptiles"
colnames(Imputed[[i]][["R"]]$Imputed.Dataset)[8] <- "Lifespan_proxy"

Imputed[[i]][["A"]]$Imputed.Dataset$Class <- "Amphibians"
colnames(Imputed[[i]][["A"]]$Imputed.Dataset)[7] <- "Lifespan_proxy"
}

## Bind into a single dataset (put all classes together) for each of the 8 imputed sets and transform and z-score the traits
Imputed_t <- Bind_transform_traits(Imputed)

## By Class: transforming and z-scoring within class
Imputed_t_class <- Bind_transform_traits_BYCLASS(Imputed)

## Saving transformed and zscored datasets
saveRDS(Imputed_t, "../../Results/Imputed_zscored_standardised/All_imputed_sets.rds")
saveRDS(Imputed_t_class, "../../Results/Imputed_zscored_standardised/All_imputed_sets_byClass.rds")

## 8th dataset for studying correlation (number 8 has been randomly chosen)
i <- 8
Selected <- readRDS("../../Results/Imputed_zscored_standardised/All_imputed_sets.rds")[[8]]

## Checking correlations among traits

## 1. First approach: correlations amongst the continuous variables - correlation matrices
Conti <- c("log10_Body_mass_g",
           "log10_Lifespan_proxy",
           "log10_Litter_size",
           "sqrt_Habitat_breadth_IUCN")
Contitraits <- Selected[, Conti]

Conti_cor <- cor(Contitraits, method = "pearson")
Conti_cor %>% det() ## if determinant close to 0 = high degree of multicollinearity
Conti_cor[upper.tri(Conti_cor, diag = FALSE)] <- NA
colnames(Conti_cor) <- c("Body mass", "Longevity", "Litter/clutch size", "Habitat breadth")
rownames(Conti_cor) <- c("Body mass", "Longevity", "Litter/clutch size", "Habitat breadth")
stargazer(Conti_cor, summary = FALSE, digits = 2)
Conti_cor[lower.tri(Conti_cor, diag = FALSE)] <- NA
colnames(Conti_cor) <- c("Body mass", "Longevity", "Litter/clutch size", "Habitat breadth")
rownames(Conti_cor) <- c("Body mass", "Longevity", "Litter/clutch size", "Habitat breadth")


# ## 2. Run stepwise variable selection based on VIF scores on all 8 imputed trait datasets

# # With all variables
# All_traits <- VIF_Selection(Imputed_t, ByClass = FALSE)
# All_traits_M <- VIF_Selection(Imputed_t, ByClass=TRUE, class="Mammalia")
# All_traits_B <- VIF_Selection(Imputed_t, ByClass=TRUE, class="Aves")
# All_traits_R <- VIF_Selection(Imputed_t, ByClass=TRUE, class="Reptilia")
# All_traits_A <- VIF_Selection(Imputed_t, ByClass=TRUE, class="Amphibia")
# #
# # # Just with continuous variables
# # ContVIF <- VIF_Selection_Continuous(Imputed_t, ByClass = FALSE)
# # M_ContVIF <- VIF_Selection_Continuous(Imputed_t, ByClass=TRUE, class="Mammalia")
# # A_ContVIF <- VIF_Selection_Continuous(Imputed_t, ByClass=TRUE, class="Aves")
# # R_ContVIF <- VIF_Selection_Continuous(Imputed_t, ByClass=TRUE, class="Reptilia")
# # A_ContVIF <- VIF_Selection_Continuous(Imputed_t, ByClass=TRUE, class="Amphibia")

# 2. Run stepwise variable selection based on VIF scores on the selected imputed trait dataset
# dummy <-  rnorm(nrow(Selected))
#
# model1 <- lm(dummy ~
#                Selected$log10_Body_mass_g +
#                Selected$log10_Lifespan_proxy +
#                Selected$log10_Litter_size +
#                Selected$sqrt_Habitat_breadth_IUCN +
#                Selected$Specialisation +
#                Selected$Diel_activity +
#                Selected$Trophic_level)
#
# pedometrics::stepVIF(model1, threshold = 5, verbose = TRUE)

VIFVal <- stepwise.vif(dataset=Selected, metrics = c(Conti, "Specialisation", "Diel_activity", "Trophic_level"), vif.threshold = 5, verbose = TRUE)$VIF
colnames(VIFVal) <- "GVIF"
rownames(VIFVal) <- c("Diel activity", "Trophic level", "Specialisation",  "Body mass", "Habitat breadth","Litter/clutch size", "Lifespan")
stargazer(VIFVal, summary = FALSE)

## Results: all lower than the threshold value of 5.

## 4. Fourth approach: Factor Analysis of Mixed Data

## Another approach : Factor analysis of mixed type data with one randomly selected dataset
X <- Selected
X <- X[, c("log10_Body_mass_g",
           "log10_Lifespan_proxy",
           "log10_Litter_size",
           "sqrt_Habitat_breadth_IUCN",
           "Specialisation",
           "Diel_activity",
           "Trophic_level")]
colnames(X) <- c("BM", "L", "LCS", "HB", "Sp", "DA", "TL")

gc()
memory.limit(size=50000)

# run the FAMD

S1 <- Sys.time()
Y <- FAMD (X, ncp = 5, sup.var = NULL, ind.sup = NULL, graph = TRUE)
S2 <- Sys.time()
S1-S2

saveRDS(Y, "../../Results/FAMD_selectedtraits.rds")

# eigenvalues and percent variation in each dimension, and plot contributions of each dimensions in explained variance
Y <- readRDS("../../Results/FAMD_selectedtraits.rds")

eig.val <- get_eigenvalue(Y)
fviz_eig(Y, addlabels = TRUE, ylim = c(0, 50))

# contribution of each variable in explained variance for each dimension
var <- get_famd_var(Y)
var$contrib

# contribution of each variables to the dimensions
par(mfrow=c(4,4))
fviz_contrib(Y, "var", axes = 1)
fviz_contrib(Y, "var", axes = 2)
fviz_contrib(Y, "var", axes = 3)
fviz_contrib(Y, "var", axes = 4)

# visualisation
par(mfrow=c(2,1))
fviz_famd_var(Y, "quanti.var", repel = TRUE,
              col.var = "black")
fviz_famd_var(Y, "quali.var", repel = TRUE,
              col.var = "black")
fviz_famd_var(Y, repel = TRUE,
              col.var = "black")

fviz_famd_ind(Y, geom = "point",  alpha.ind = 0.2)

Ycoord <- as.data.frame(Y$ind$coord)
Ycoord$Class <- Selected$Class
plot(Ycoord$Dim.1, Ycoord$Dim.2, pch=19, col=as.factor(Ycoord$Class))

GGPoptions2 <- theme(
  panel.border = element_rect(colour = "white", fill=NA),
  text = element_text(size=13, family="sans"),
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12),
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=20))

Figpoints <- ggplot(Ycoord, aes(Dim.1, Dim.2, col=Class)) + geom_point(alpha=0.4) + geom_vline(xintercept = 0, lty="dashed") +
   geom_hline(yintercept = 0, lty="dashed")  + GGPoptions2 + theme_bw() + xlab("Dim. 1 (24.1%)") + ylab("Dim. 2 (19.3%)")+
  theme(panel.border = element_rect(colour = "white", fill=NA)) + theme(legend.position = "right", legend.title = element_blank()) + theme(text = element_text(size=20, family="serif"))

ggsave(Figpoints, filename="../../Results/Traits_FAMD.pdf", width=6, height = 4)

### FAMD on Gower distance matrix subsetted for PREDICTS species
Gowerdist <- readRDS("../../Results/Gower_distances/Gower_distance_PREDICTS_8_sets.rds")
Gowerdist <- Gowerdist[[8]]

Predicts <- readRDS("../../Data/PredictsVertebrates.rds")
SP <- unique(Predicts$Best_guess_binomial)
SpC <- unique(Predicts[, c("Best_guess_binomial", "Class")])
Gowr <-  usedist::dist_subset(Gowerdist, SP)
if(!ade4::is.euclid(Gowr)){
  Gowr <- ade4::cailliez(Gowr)
}

x.pco <- ade4::dudi.pco(Gowr, scannf = FALSE, full = TRUE)
Axes <- x.pco$li[,1:2]
Axes$Class <- SpC$Class

plot(Axes$A1, Axes$A2, col=Axes$Class)

Sesub <- subset(Selected, Selected$Best_guess_binomial %in% SP)
Sesub <- cbind(Sesub, Axes)
plot(Sesub$A1, Sesub$A2, col=as.factor(Sesub$Class))

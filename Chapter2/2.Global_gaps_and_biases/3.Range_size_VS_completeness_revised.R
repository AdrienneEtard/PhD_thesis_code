library(ggplot2)
library(dplyr)
library(stargazer)

GGPoptions <- theme_classic()+ theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=16, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=16), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=16),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=16)) 


# # # # # # # # # # # # # # # # # # # # # # # # #


traits.completeness <- readRDS("../Results/Traits_to_map/traits_completeness_V2.rds")

RangeSizes <- read.csv("../Data/Distribution_maps/RangeSizes.csv") 

RangeSizes %>%
  group_by(Class) %>%
  summarise(Count=n())

RangeSizes <- RangeSizes %>% 
  select(-Class, -Index)
colnames(RangeSizes)[1] <- "Best_guess_binomial"

## Add range sizes to completeness datasets

# number of species in trait datasets
lapply(traits.completeness, nrow) 

# number of species for which distribution maps are processed (based on counts before cutting maps here)
lapply(traits.completeness, function(x){return(length(intersect(x$Best_guess_binomial, RangeSizes$Best_guess_binomial)))})

traits.completeness <- lapply(traits.completeness, function(x)
  {return(merge(x, RangeSizes, by="Best_guess_binomial"))})
lapply(traits.completeness, nrow) 

## not completeness, but number of sampled traits
N_sampled <- function(traits) {
  
  Sub <- traits %>% 
    select(RangeSize_sqkm_BeforeCuts, RangeSize_sqkm_AfterCuts)
  traits <- traits %>%
    dplyr::select(-completeness, 
                  -Order,
                  -Family,
                  -Genus,
                  -Best_guess_binomial,
                  -RangeSize_sqkm_BeforeCuts, 
                  -RangeSize_sqkm_AfterCuts)
  
  traits$N_sampled <- apply(
    traits, 1, function(x) 
    {return (length(x[!is.na(x)]))}
  ) 
  return(cbind(traits, Sub))
}

traits.completeness <- lapply(traits.completeness, N_sampled)
traits.completeness$Amphibians$Class <- "Amphibians"
traits.completeness$Birds$Class <- "Birds"
traits.completeness$Mammals$Class <- "Mammals"
traits.completeness$Reptiles$Class <- "Reptiles"
Data <- rbind(traits.completeness$Amphibians, traits.completeness$Birds, traits.completeness$Mammals, traits.completeness$Reptiles)
Data$Class <- factor(Data$Class, levels = c("Mammals", "Birds", "Amphibians", "Reptiles"))


## Fit Poisson model with class as interacting factor

# https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/ for plotting CI


FitPoissonModel <- function(Data, Before_or_After) {
  
  if (Before_or_After=="Before"){
    
    Data <- Data %>% 
      filter(!is.na(RangeSize_sqkm_BeforeCuts)) %>% 
      filter(RangeSize_sqkm_BeforeCuts!=0)
      
    Data$logRS <- log(Data$RangeSize_sqkm_BeforeCuts)
    Model <- glm(as.numeric(N_sampled) ~ logRS + Class + logRS:Class, family="poisson", data=Data)
    
  }
  
  if (Before_or_After=="After"){
    
    Data <- Data %>% 
      filter(!is.na(RangeSize_sqkm_AfterCuts)) %>% 
      filter(RangeSize_sqkm_AfterCuts!=0)
    
    Data$logRS <- log(Data$RangeSize_sqkm_AfterCuts)
    Model <- glm(as.numeric(N_sampled) ~ logRS + Class + logRS:Class, family="poisson", data=Data)
    
  }
  
  
  return(list(Model=Model, ModelData=Data))
  
}

PlotModel <- function(Model, Facet_Wrap) {
  
  ## grad the inverse link function
  ilink <- family(Model$Model)$linkinv
  
  ## some data to predict at: 100 values over the range of logRS
  ndata1 <- with(Model$ModelData, tibble(logRS = seq(min(logRS), max(logRS),length = 100), Class="Amphibians"))
  ndata2 <- with(Model$ModelData, tibble(logRS = seq(min(logRS), max(logRS),length = 100), Class="Reptiles"))
  ndata3 <- with(Model$ModelData, tibble(logRS = seq(min(logRS), max(logRS),length = 100), Class="Birds"))
  ndata4 <- with(Model$ModelData, tibble(logRS = seq(min(logRS), max(logRS),length = 100), Class="Mammals"))
  ndata <- rbind(ndata1, ndata2, ndata3, ndata4)
  
  ## add the fitted values by predicting from the model for the new data
  ndata <- tibble::add_column(ndata, fit = predict(Model$Model, newdata = ndata, type = 'response'))
  
  ## add fit and se.fit on the **link** scale
  ndata <- bind_cols(ndata, setNames(as_tibble(predict(Model$Model, ndata, se.fit = TRUE)[1:2]),c('fit_link','se_link')))
  
  ## create the interval and backtransform
  ndata <- mutate(ndata,
                  fit_resp  = ilink(fit_link),
                  right_upr = ilink(fit_link + (2 * se_link)),
                  right_lwr = ilink(fit_link - (2 * se_link)))

  p <- ggplot(ndata, aes(x = logRS, y = fit,  group=Class, col=Class, fill=Class)) +
    geom_line(size = 1) +
    scale_fill_viridis_d() + scale_color_viridis_d() + 
    geom_ribbon(data = ndata, aes(ymin = right_lwr, ymax = right_upr), alpha = 0.2, colour=NA) +
    GGPoptions + xlab(expression("Range size (log km"^2*")")) + ylab("Numbers of sampled traits")
   # + geom_point(data = Model$ModelData,
   #             aes(x = logRS, y = N_sampled, group=Class, col=Class),
   #             alpha=.1, size=0.1, position=position_jitter(h=.1))
   # 
  if(Facet_Wrap) {
    p <- p + facet_wrap(~Class)
  }
  
  return(p)
}

## Poisson model with RS before cutting by altitudinal limits
MBefore <- FitPoissonModel(Data, "Before")
summary(MBefore$Model)
PlotModel(MBefore, TRUE)
pB <- PlotModel(MBefore, FALSE) + theme(legend.title = element_blank()) 



## Poisson model after cutting by altitudinal limits
MAfter <- FitPoissonModel(Data, "After")
summary(MAfter$Model)
stargazer(summary(MAfter$Model)$coefficients, digits=NA)
PlotModel(MAfter, TRUE)
pA <- PlotModel(MAfter, FALSE) + theme(legend.title = element_blank()) 

## sample sizes
MBefore$Model$data %>% 
  group_by(Class) %>% 
  summarise(Count=n())


MAfter$Model$data %>% 
  group_by(Class) %>% 
  summarise(Count=n())


# save plots
ggsave(pA, filename="../Results/Poisson_model/Plot_AfterCuttingRS.pdf", height=4, width =6)
ggsave(pB, filename="../Results/Poisson_model/Plot_BeforeCuttingRS.pdf", height=4, width =6)
ggsave(pB, filename="../Results/Poisson_model/Plot_BeforeCuttingRS.png", height=4, width =6)

# estimate of overdispersion residual deviance/residual df (oversdispersion if this is significantly greater than 1)
MBefore$Model$deviance / MBefore$Model$df.residual
with(MBefore$Model, cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))

MAfter$Model$deviance / MAfter$Model$df.residual
with(MAfter$Model, cbind(res.deviance = deviance, df = df.residual, p = pchisq(deviance, df.residual, lower.tail=FALSE)))

summary(MAfter$Model)


## plotting differences in range sizes
RangeSizes <- read.csv("../Data/Distribution_maps/RangeSizes.csv") 

RangeSizesPlot <- RangeSizes %>%
  filter(!is.na(RangeSize_sqkm_BeforeCuts)) %>% 
  filter(RangeSize_sqkm_BeforeCuts!=0) %>% 
  filter(RangeSize_sqkm_AfterCuts!=0) 

RangeSizesPlot$Delta <- RangeSizesPlot$RangeSize_sqkm_BeforeCuts - RangeSizesPlot$RangeSize_sqkm_AfterCuts
 
RangeSizesPlot %>%
  group_by(Class) %>%
  summarise(Count=n())

hist(log(RangeSizesPlot$Delta), breaks=50)

p <- ggplot(RangeSizesPlot, aes(log(RangeSize_sqkm_BeforeCuts), log(RangeSize_sqkm_AfterCuts))) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, col="blue") +
  facet_wrap(~Class, scales = "free") + GGPoptions +
  xlab("Range size before cutting (log square km)") + ylab("Range size after cutting (log square km)")

ggsave(p, filename = "../Results/Poisson_model/RangeSizes_before_aftercuts.pdf", width=6, height = 5)

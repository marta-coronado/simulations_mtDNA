#### LIBRARIES AND FUNCTIONS ####

library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(ggrepel)
library(dplyr)
library(reshape2)

theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#FFDF4B","#49B269")), ...)
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}


#### DATA ####

setwd("/run/user/1001/gvfs/sftp:host=158.109.38.201,user=mcoronado/data/home/mcoronado/Simulations_mtDNA_Pipeline/mtDNA_general/Simulations/")

PhiPack_result <- read.table(file="PhiPack_global_result.txt",header=T,sep=" ")

PhiPack_filter<-PhiPack_result %>% group_by(scenarioType,simulation,recRate,rho) %>%
  summarise(percentatgeNSS=sum(pValNSS<0.05)/n(), 
            percentatgeMaxChi=sum(pValMaxChi<0.05)/n(),
            percentatgePHI=sum(pValPHI<0.05)/n())
PhiPack_filter<-melt(PhiPack_filter,id.vars = c("scenarioType","simulation","recRate","rho"),measure.vars = c("percentatgeNSS","percentatgeMaxChi","percentatgePHI"))

PhiPack_filter$variable <- as.character(PhiPack_filter$variable)
PhiPack_filter$variable[PhiPack_filter$variable == "percentatgeNSS"] <- "NSS"
PhiPack_filter$variable[PhiPack_filter$variable == "percentatgeMaxChi"] <- "Max X²"
PhiPack_filter$variable[PhiPack_filter$variable == "percentatgePHI"] <- "PHI (Normal)"
PhiPack_filter$variable <- as.factor(PhiPack_filter$variable)

PhiPack_filter$scenarioType <- as.character(PhiPack_filter$scenarioType)
PhiPack_filter$scenarioType[PhiPack_filter$scenarioType == "1"] <- "Scenario 1"
PhiPack_filter$scenarioType[PhiPack_filter$scenarioType == "2"] <- "Scenario 2"
PhiPack_filter$scenarioType[PhiPack_filter$scenarioType == "3"] <- "Scenario 3"
PhiPack_filter$scenarioType <- as.factor(PhiPack_filter$scenarioType)
PhiPack_filter$value<-PhiPack_filter$value*100

ggplot(PhiPack_filter, aes(x=rho,y=value,colour=variable)) + 
  geom_point() +
  geom_line() + 
  facet_grid(cols = vars(scenarioType)) +
  labs(x = "Population recombination parameter (ρ)", y = "Power (% p-value<0.05)", colour="Recombination tests") + 
  scale_x_continuous(breaks = round(seq(min(PhiPack_filter$rho), max(PhiPack_filter$rho), by = 2),1)) +
  theme_Publication() + scale_colour_Publication()

ggsave("PhiPack_general_3scenarios.png",width = 7.50, height = 5)


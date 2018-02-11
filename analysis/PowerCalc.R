rm(list=ls())
graphics.off()
pardefault <- par()
options(warn=1)



# define path for the dir of current file
rstudioapi::getActiveDocumentContext
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
CurrentPath<-(dirname(rstudioapi::getActiveDocumentContext()$path))


#load libraries
library(pacman)
pacman::p_load(pwr,rstudioapi, stats, dtplyr,plyr,arm,data.table,lme4,lmerTest,ggplot2,corrplot,scatterplot3d,car,directlabels,ggfortify,piecewiseSEM,compare,dplyr,
               weights, tidyr, cowplot, ggpubr,  pROC,stringi,pracma,ggcorrplot,reshape,ggrepel,reshape2,xlsx,akima,sjPlot,RColorBrewer,sandwich,foreign,RCurl,
               eba,stargazer,texreg,car,pscl,broom,extrafont,fontcm,tikzDevice,sos) 


#load data
filefunc<-list.files(CurrentPath,pattern = "ofc", recursive = TRUE)

read.table("ofc_fslmeants_results.txt",header = False)



# pwr.t.test(n = , d = , sig.level = , power = , type = c("two.sample", "one.sample", "paired"))


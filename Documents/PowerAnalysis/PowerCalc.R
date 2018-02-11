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
#these are the results files of fslmeant with mask for ofc or temporal clustered, that Rotem got in her analysis. obtain by these vommant lines:
#for ofc: fslmeants -i cope17.gfeat/cope1.feat/filtered_func_data.nii.gz -o ofc_fslmeants_results.txt -m Masks/left_ofc.nii.gz 
#for temporal: fslmeants -i cope17.gfeat/cope1.feat/filtered_func_data.nii.gz -o temporal_fslmeants_results.txt -m Masks/left_temporal_occipital_after.nii.gz 

AllDataOfc<-read.table("ofc_fslmeants_results.txt")
AllDataTemporal<-read.table("temporal_fslmeants_results.txt")

#take only the experiment group - 36 first subjects
Nsub=36
OfcGroupData<-as.data.frame(AllDataOfc[1:Nsub,])
TemporalGroupData<-as.data.frame(AllDataTemporal[1:Nsub,])

#calculate effect size by cohen d
MeanOfc<-mean(OfcGroupData[,1])
SdOfc<-sd(OfcGroupData[,1])
EsOfc<-MeanOfc/SdOfc

MeanTemporal<-mean(TemporalGroupData[,1])
SdTemporal<-sd(TemporalGroupData[,1])
EsTemporal<-MeanTemporal/SdTemporal

#calculate power
OfcResults<-pwr.t.test(n = , d = EsOfc, sig.level = 0.001, power = .8, type = "one.sample")
TemporalResults<-pwr.t.test(n = , d = EsTemporal, sig.level = 0.001, power = .8, type = "one.sample")


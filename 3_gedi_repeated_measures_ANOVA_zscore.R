# Repeat-measures ANOVA, or mixed-effects model

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 08/18/2024


library(tidyverse)
library(rstatix)
library(performance)
library(nlme)
library(sp)
library(dplyr)
library(broom.mixed)


#############################
# Parameters

# minimum time after fire limit (days)
#minAfter <- 0 # 0 years
#minAfter <- 365 # 1 year
minAfter <- 730 # 2 years

# maximum time after fire limit (days)
#maxAfter <- 365 # 1 year
#maxAfter <- 730 # 2 years
maxAfter <- 1500 # no limit

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

# distance between pairs
distance <- 45

# allowable absolute difference in elevation (m)
delta_elev <- 10

# allowable growth in tAGBD Mg/ha/year
biomass <- 99999 # not used

# date for output file
version <- "240822"

#inDir <- "E:/active/project/calfire_gedi/rq2_paper"
inDir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_wildfire_change_in_footprint_metrics"

# GEDI metrics
metrics <- read.csv(paste0(inDir,"/data/gedi_structure_variables_reduced_pai_240809.csv"))
metrics <- metrics %>% filter(selected == "Y")

inputFile <- paste0(inDir,"/data/master_gedi_fire_difference_",distance,"m_240818.csv")

outputFile <- paste0(inDir,"/data/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")

# read in data
df <- read.csv(inputFile) 

# calculate biomass growth limit based on time difference
df$abgd_limit <- (df$timeDiff/365) * biomass

# filter to time before and after fire limits, absolute change in elevation limit, & biomass growth limit
dfFiltered <- df %>% filter(timeBefore <= maxBefore & timeAfter >= minAfter & timeAfter <= maxAfter & 
                              abs(delta_GEDI_elev) <= delta_elev &
                              delta_strct_tAGBD <= abgd_limit)

# prepare data
postVars <- c(paste("post",metrics$metric,sep="_"),"mtbs_classP","post_x_teale","post_y_teale")
d1 <- dfFiltered %>% select(!!postVars)
colnames(d1)[1:dim(metrics)[1]] <- metrics$metric
ncol <- dim(d1)[2]
d1<- cbind(d1[(ncol-2):ncol], stack(d1[1:(ncol-3)]))
d1$time <- "post-fire"
d1 <- tibble::rowid_to_column(d1, "index") %>% rename(x_teale=post_x_teale,y_teale=post_y_teale)

preVars <- c(paste("pre",metrics$metric,sep="_"),"mtbs_classP","pre_x_teale","pre_y_teale")
d2 <- dfFiltered %>% select(!!preVars)
colnames(d2)[1:dim(metrics)[1]] <- metrics$metric
ncol <- dim(d2)[2]
d2<- cbind(d2[(ncol-2):ncol], stack(d2[1:(ncol-3)]))
d2$time <- "pre-fire"
d2 <- tibble::rowid_to_column(d2, "index") %>% rename(x_teale=pre_x_teale,y_teale=pre_y_teale)

d <- rbind(d1,d2)
d <- d %>% rename(metric = ind) %>% mutate(treatment = ifelse(mtbs_classP == 0,"Control","Fire")) 

d$treatment <- factor(d$treatment , levels=c("Control","Fire"))
d$time <- factor(d$time , levels=c("pre-fire","post-fire"))

# Function to rescale using z score and filter outliers with abs(z) > 4
zscore<-function(df,z=4){
  tmp <- as.matrix(apply(as.matrix(df),2,scale))
  if (ncol(tmp)==1) {tmp[which(abs(tmp)>z)] <- NA}
  if (ncol(tmp)>1) {tmp[which(abs(tmp)>z,arr.ind=T)] <- NA}
  return(tmp)
}

output <- data.frame()
# loop through each metric and create graphs, perform ANOVA
for (strctVar in metrics$metric){
  
  # strctVar = metrics$metric[89] 
  print(strctVar)
  
  # gather data for analysis
  data <- as_tibble(d %>% filter(metric == !!strctVar) %>% mutate_at(vars(index,mtbs_classP,time,treatment), as.factor))
  
  dataFiltered <- data %>% mutate(zscore = zscore(data$values)) %>% na.omit() %>% distinct(x_teale,y_teale, .keep_all=T)

  # filter for paired samples only
  n <- dataFiltered %>% count(index)
  dataFiltered <- left_join(dataFiltered,n,by="index") %>% filter(n == 2)
  
  # sample data to match control
  set.seed(123)
  dataFilteredPre <- dataFiltered %>% filter(time == "pre-fire")
  n <- min(table(dataFilteredPre$treatment))
  dataSampledPre <- dataFilteredPre %>% sample_n_by(treatment, size = n)
  dataSampledPost <- dataFiltered %>% filter(time == "post-fire" & index %in% dataSampledPre$index)
  dataSampled <- rbind(dataSampledPre,dataSampledPost)
  
  tryCatch({
    aspatial <- nlme::lme(zscore ~ treatment * time, data = dataSampled, random = ~ 1|index, control = lmeControl(opt = "optim"))
    spatial <- update(aspatial, correlation = corLin(form=~ x_teale + y_teale,nugget=F))
    
    aspatialResults <- tidy(aspatial)
    aspatialInterceptEst <- aspatialResults$estimate[1]
    aspatialTreatmentEst <- aspatialResults$estimate[2]
    aspatialTimeEst <- aspatialResults$estimate[3]
    aspatialTreatmentTimeEst <- aspatialResults$estimate[4]
    aspatialInterceptPval <- aspatialResults$p.value[1]
    aspatialTreatmentPval <- aspatialResults$p.value[2]
    aspatialTimePval <- aspatialResults$p.value[3]
    aspatialTreatmentTimePval <- aspatialResults$p.value[4]
    aspatialAIC <- AIC(aspatial)
    
    spatialResults <- tidy(spatial)
    spatialInterceptEst <- spatialResults$estimate[1]
    spatialTreatmentEst <- spatialResults$estimate[2]
    spatialTimeEst <- spatialResults$estimate[3]
    spatialTreatmentTimeEst <- spatialResults$estimate[4]
    spatialInterceptPval <- spatialResults$p.value[1]
    spatialTreatmentPval <- spatialResults$p.value[2]
    spatialTimePval <- spatialResults$p.value[3]
    spatialTreatmentTimePval <- spatialResults$p.value[4]
    spatialAIC <- AIC(spatial)
    
    anovaModelPval <- anova(aspatial,spatial)$`p-value`[2]
    
    results <- data.frame(strctVar,n,anovaModelPval, 
                 aspatialInterceptEst,aspatialTreatmentEst,aspatialTimeEst,aspatialTreatmentTimeEst,
                 aspatialInterceptPval,aspatialTreatmentPval,aspatialTimePval,aspatialTreatmentTimePval,aspatialAIC,
                 spatialInterceptEst,spatialTreatmentEst,spatialTimeEst,spatialTreatmentTimeEst,
                 spatialInterceptPval,spatialTreatmentPval,spatialTimePval,spatialTreatmentTimePval,spatialAIC
    )
    
    colnames(results) <- c("metric","n","ModelTestPVal",
                          "aspatialInterceptEst","aspatialTreatmentEst","aspatialTimeEst","aspatialTreatmentTimeEst",
                          "aspatialInterceptPval","aspatialTreatmentPval","aspatialTimePval","aspatialTreatmentTimePval","aspatialAIC",
                          "spatialInterceptEst","spatialTreatmentEst","spatialTimeEst","spatialTreatmentTimeEst",
                          "spatialInterceptPval","spatialTreatmentPval","spatialTimePval","spatialTreatmentTimePval","spatialAIC"
    )
    
    # bind to output table
    output <- rbind(output, results)
  
  }, error = function(e) {
    # Handle error by skipping and printing a message
    cat("Error occurred: Skipping ANOVA")
  })
}

# write out repeated ANOVA table
write.csv(output, outputFile, row.names=F)


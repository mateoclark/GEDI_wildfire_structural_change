# Multivariate spatial regression models

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 08/09/2024

library(spatialreg)
library(sp)
library(spdep)
library(dplyr)

# working directory
indir <- "/GEDI_wildfire_structural_change"
#indir <- "E:/active/project/calfire_gedi/rq2_paper"

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit
#maxAfter <- 730 # 2 years

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

# Number of nearest neighbors
k <- 50  

# Distance of pairs
distance <- 45

# allowable absolute difference in elevation (m)
delta_elev <- 10

# allowable growth in tAGBD Mg/ha/year
biomass <- 99999 # not used

# predictor to use
predictor <- "mtbs_dnbrOW" # post-fire offseted dNBR and weighted for 45-m circle

version <- "240822"

# input .csv file with GEDI metric differences
inputFile <- paste0(indir,"/master_gedi_fire_difference_",distance,"m_240822.csv")

# output .csv file with GEDI metric difference spatial and OLS regression results
outputFile <- paste0(indir,"/master_gedi_fire_difference_spatialreg_interactions_",predictor,"_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")

# GEDI metrics
metrics <- read.csv(paste0(indir,"/gedi_structure_variables_reduced_pai_240927.csv"))
metrics <- metrics %>% filter(selected == "Y")
metrics$metricDelta <- paste0("delta_",c(metrics$metric))
metricsDelta <- metrics$metricDelta
metrics$metricPrefire <- paste0("pre_",c(metrics$metric))
metricsPrefire <- metrics$metricPrefire

# read data
data <- read.csv(inputFile)

# calculate biomass growth limit based on time difference
data$abgd_limit <- (data$timeDiff/365) * biomass

# filter data to include only data in fires (not control)
dataFiltered <- data %>% filter(!mtbs_classP == 0 & 
                                timeBefore <= maxBefore & timeAfter >= minAfter & timeAfter <= maxAfter & 
                                abs(delta_GEDI_elev) <= delta_elev &
                                delta_strct_tAGBD <= abgd_limit
                                ) %>%
  select(all_of(c("GEDI_shot","post_x_teale","post_y_teale","topoSlope","typeForest","beamType",metricsDelta,metricsPrefire,"mtbs_dnbrOW","mtbs_rdnbrW","mtbs_CBIW","timeAfter"))) %>%
  rename(x = post_x_teale,
         y = post_y_teale)

# Function to rescale using z score and set outliers to NA
outlierRemoval <- function(df, z = 4) {
  tmp <- as.matrix(apply(as.matrix(df), 2, scale))
  if (ncol(tmp) == 1) {
    tmp[which(abs(tmp) > z)] <- NA
  } else {
    tmp[which(abs(tmp) > z, arr.ind = TRUE)] <- NA
  }
  return(tmp)
}

# process each metric
output <- data.frame()
for (i in 1:length(metricsDelta)){
  
  response <- metricsDelta[i]
  print(response)
  
  d <- data.frame(response = dataFiltered[[response]],
                  predictor = dataFiltered[[predictor]],
                  classForest = dataFiltered$typeForest,
                  beamType = dataFiltered$beamType,
                  x = dataFiltered$x,
                  y = dataFiltered$y,
                  time = dataFiltered$timeAfter,
                  topoSlope = dataFiltered$topoSlope,
                  prefireMetric = dataFiltered[,gsub("delta","pre",response)]
  )
  d <- d %>% filter(!is.infinite(response))
  d$response <- outlierRemoval(d$response) # z-score transformation and outlier removal
  d$predictor <- outlierRemoval(d$predictor) # z-score transformation and outlier removal
  d$time <- outlierRemoval(d$time) # z-score transformation and outlier removal
  d$topoSlope <- outlierRemoval(d$topoSlope) # z-score transformation and outlier removal
  d$prefireMetric <- outlierRemoval(d$prefireMetric) # z-score transformation and outlier removal
  d <- d %>% na.omit() 
  
  if ((dim(d)[1] > 1000) & (sum(d$response) != 0)){
    
    n <- dim(d)[1]
    
    coordinates(d) <- c("x","y")
    spatial_neighbors <- knearneigh(d,k = k)
    
    # Build spatial weights matrix
    # Weights are based on inverse distance
    W <- knn2nb(spatial_neighbors)
    W <- nb2listw(W,style = "W") # W is row standardized (sums over all links to n),better for large data
    
    # Create a spatial error model LU sparse matrix
    se_model <- errorsarlm(
      formula = response ~ predictor,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # beam type interaction model 
    se_model_beam <- errorsarlm(
      formula = response ~ predictor * beamType,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # forest class interaction model
    se_model_forest <- errorsarlm(
      formula = response ~ predictor * classForest,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # time difference interaction model
    se_model_time <- errorsarlm(
      formula = response ~ predictor * time,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # time difference no interaction model
    se_model_timeNoX <- errorsarlm(
      formula = response ~ predictor + time,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # slope with interaction model
    se_model_topoSlope <- errorsarlm(
      formula = response ~ predictor * topoSlope,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # slope no interaction model
    se_model_topoSlopeNoX <- errorsarlm(
      formula = response ~ predictor + topoSlope,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # slope alone, no dNBR model
    se_model_topoSlopeNoBurn <- errorsarlm(
      formula = response ~ topoSlope,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # prefireMetric with interaction model
    se_model_prefireMetric <- errorsarlm(
      formula = response ~ predictor * prefireMetric,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # prefireMetric no interaction model
    se_model_prefireMetricNoX <- errorsarlm(
      formula = response ~ predictor + prefireMetric,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # prefireMetric alone, no dNBR model
    se_model_prefireMetricNoBurn <- errorsarlm(
      formula = response ~ prefireMetric,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # prefireMetric beam type interaction model 
    se_model_prefireMetricNoBurn_beam <- errorsarlm(
      formula = response ~ prefireMetric * beamType,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    # prefireMetric forest class interaction model
    se_model_prefireMetricNoBurn_forest <- errorsarlm(
      formula = response ~ prefireMetric * classForest,
      data = d,
      listw = W,
      method="LU" # sparse matrix
    )
    
    anova_beam <- anova(se_model,se_model_beam)
    anova_forest <- anova(se_model,se_model_forest)
    anova_time <- anova(se_model,se_model_time)
    anova_timeNoX <- anova(se_model,se_model_timeNoX)
    anova_timeNoXwX <- anova(se_model_timeNoX,se_model_time)
    anova_topoSlope <- anova(se_model,se_model_topoSlope)
    anova_topoSlopeNoX <- anova(se_model,se_model_topoSlopeNoX)
    anova_topoSlopeNoXwX <- anova(se_model_topoSlopeNoX,se_model_topoSlope)
    anova_prefireMetric <- anova(se_model,se_model_prefireMetric)
    anova_prefireMetricNoX <- anova(se_model,se_model_prefireMetricNoX)
    anova_prefireMetricNoXwX <- anova(se_model_prefireMetricNoX,se_model_prefireMetric)
    anova_prefireMetricNoBurn_forest <- anova(se_model_prefireMetricNoBurn,se_model_prefireMetricNoBurn_forest)
    anova_prefireMetricNoBurn_beam <- anova(se_model_prefireMetricNoBurn,se_model_prefireMetricNoBurn_beam)
    
    # OLS model
    ols_model <- lm(response ~ predictor,
                    data = d)
    ols_time_model <- lm(response ~ predictor * time,
                         data = d)
    
    ## Organize regression results
    
    # OLS
    ols_AIC <- AIC(ols_model)
    ols_summary <- summary(ols_model)
    ols_int <- ols_summary$coefficients[1,1]
    ols_slope <- ols_summary$coefficients[2,1]
    ols_pvalslope <- ols_summary$coefficients[1,4]
    ols_pvalint <- ols_summary$coefficients[2,4]
    ols_r2 <- ols_summary$adj.r.squared
    
    olstime_AIC <- AIC(ols_time_model)
    olstime_summary <- summary(ols_time_model)
    olstime_int <- olstime_summary$coefficients[1,1]
    olstime_slope <- olstime_summary$coefficients[2,1]
    olstime_time <- olstime_summary$coefficients[3,1]
    olstime_interaction <- olstime_summary$coefficients[4,1]
    olstime_pvalint <- olstime_summary$coefficients[1,4]
    olstime_pvalslope <- olstime_summary$coefficients[2,4]
    olstime_pvaltime <- olstime_summary$coefficients[3,4]
    olstime_pvalinteraction <- olstime_summary$coefficients[4,4]
    olstime_r2 <- olstime_summary$adj.r.squared
    
    anova_ols_time <- anova(ols_model,ols_time_model)
    ols_anova_pvalue <- anova_ols_time$`Pr(>F)`[2]
    
    # Spatial error 
    model_summary <- summary(se_model,Nagelkerke = TRUE)
    spatial_int <-  model_summary$Coef[1,1]
    spatial_slope <-  model_summary$Coef[2,1]
    spatial_pvalint <- model_summary$Coef[1,4]
    spatial_pvalslope <- model_summary$Coef[2,4]
    spatial_psuedo_r2 <- model_summary$NK
    spatial_AIC <- AIC(se_model)
    
    # Spatial error with beam type
    beam_anova_pvalue <- anova_beam$`p-value`[2]
    model_summary <- summary(se_model_beam,Nagelkerke = TRUE)
    
    beam_int_coverage <-  model_summary$Coef[1,1]
    beam_slope_coverage <-  model_summary$Coef[2,1]
    beam_int_full <-  model_summary$Coef[3,1]
    beam_int_mixed <-  model_summary$Coef[4,1]
    beam_slope_full <-  model_summary$Coef[5,1]
    beam_slope_mixed <-  model_summary$Coef[6,1]
    
    beam_pvalint_coverage <-  model_summary$Coef[1,4]
    beam_pvalslope_coverage <-  model_summary$Coef[2,4]
    beam_pvalint_full <-  model_summary$Coef[3,4]
    beam_pvalint_mixed <-  model_summary$Coef[4,4]
    beam_pvalslope_full <-  model_summary$Coef[5,4]
    beam_pvalslope_mixed <-  model_summary$Coef[6,4]
    
    beam_psuedo_r2 <- model_summary$NK
    beam_AIC <- AIC(se_model_beam)
    
    # Spatial error with forest type
    forest_anova_pvalue <- anova_forest$`p-value`[2]
    model_summary <- summary(se_model_forest,Nagelkerke = TRUE)
    
    forest_int_conifer <-  model_summary$Coef[1,1]
    forest_slope_conifer <-  model_summary$Coef[2,1]
    forest_int_hardwood <-  model_summary$Coef[3,1]
    forest_int_mixed <-  model_summary$Coef[4,1]
    forest_slope_hardwood <-  model_summary$Coef[5,1]
    forest_slope_mixed <-  model_summary$Coef[6,1]
    
    forest_pvalint_conifer <-  model_summary$Coef[1,4]
    forest_pvalslope_conifer <-  model_summary$Coef[2,4]
    forest_pvalint_hardwood <-  model_summary$Coef[3,4]
    forest_pvalint_mixed <-  model_summary$Coef[4,4]
    forest_pvalslope_hardwood <-  model_summary$Coef[5,4]
    forest_pvalslope_mixed <-  model_summary$Coef[6,4]
    
    forest_psuedo_r2 <- model_summary$NK
    forest_AIC <- AIC(se_model_forest)
    
    
    # Spatial error with time
    time_anova_pvalue <- anova_time$`p-value`[2]
    model_summary <- summary(se_model_time,Nagelkerke = TRUE)
    
    time_int <-  model_summary$Coef[1,1]
    time_slope <-  model_summary$Coef[2,1]
    time_time <-  model_summary$Coef[3,1]
    time_interaction <-  model_summary$Coef[4,1]
    time_pvalint <-  model_summary$Coef[1,4]
    time_pvalslope <-  model_summary$Coef[2,4]
    time_pvaltime <-  model_summary$Coef[3,4]
    time_pvalinteraction <-  model_summary$Coef[4,4]
    
    time_psuedo_r2 <- model_summary$NK
    time_AIC <- AIC(se_model_time)
    
    # Spatial error with time (no interaction)
    timeNoX_anova_pvalue <- anova_timeNoX$`p-value`[2]
    model_summary <- summary(se_model_timeNoX,Nagelkerke = TRUE)
    
    timeNoX_int <-  model_summary$Coef[1,1]
    timeNoX_slope <-  model_summary$Coef[2,1]
    timeNoX_time <-  model_summary$Coef[3,1]
    timeNoX_pvalint <-  model_summary$Coef[1,4]
    timeNoX_pvalslope <-  model_summary$Coef[2,4]
    timeNoX_pvaltime <-  model_summary$Coef[3,4]
    
    timeNoX_psuedo_r2 <- model_summary$NK
    timeNoX_AIC <- AIC(se_model_timeNoX)
    
    timeNoXwX_anova_pvalue <- anova_timeNoXwX$`p-value`[2]
    
    # Spatial error with topoSlope
    topoSlope_anova_pvalue <- anova_topoSlope$`p-value`[2]
    model_summary <- summary(se_model_topoSlope,Nagelkerke = TRUE)
    
    topoSlope_int <-  model_summary$Coef[1,1]
    topoSlope_slope <-  model_summary$Coef[2,1]
    topoSlope_topoSlope <-  model_summary$Coef[3,1]
    topoSlope_interaction <-  model_summary$Coef[4,1]
    topoSlope_pvalint <-  model_summary$Coef[1,4]
    topoSlope_pvalslope <-  model_summary$Coef[2,4]
    topoSlope_pvaltopoSlope <-  model_summary$Coef[3,4]
    topoSlope_pvalinteraction <-  model_summary$Coef[4,4]
    
    topoSlope_psuedo_r2 <- model_summary$NK
    topoSlope_AIC <- AIC(se_model_topoSlope)
    
    # Spatial error with topoSlope (no interaction)
    topoSlopeNoX_anova_pvalue <- anova_topoSlopeNoX$`p-value`[2]
    model_summary <- summary(se_model_topoSlopeNoX,Nagelkerke = TRUE)
    
    topoSlopeNoX_int <-  model_summary$Coef[1,1]
    topoSlopeNoX_slope <-  model_summary$Coef[2,1]
    topoSlopeNoX_topoSlope <-  model_summary$Coef[3,1]
    topoSlopeNoX_pvalint <-  model_summary$Coef[1,4]
    topoSlopeNoX_pvalslope <-  model_summary$Coef[2,4]
    topoSlopeNoX_pvaltopoSlope <-  model_summary$Coef[3,4]
    
    topoSlopeNoX_psuedo_r2 <- model_summary$NK
    topoSlopeNoX_AIC <- AIC(se_model_topoSlopeNoX)
    
    topoSlopeNoXwX_anova_pvalue <- anova_topoSlopeNoXwX$`p-value`[2]
    
    # Spatial error with topoSlope, no burn index
    model_summary <- summary(se_model_topoSlopeNoBurn,Nagelkerke = TRUE)
    
    topoSlopeNoBurn_int <-  model_summary$Coef[1,1]
    topoSlopeNoBurn_slope <-  model_summary$Coef[2,1]
    topoSlopeNoBurn_pvalint <-  model_summary$Coef[1,4]
    topoSlopeNoBurn_pvalslope <-  model_summary$Coef[2,4]
    
    topoSlopeNoBurn_psuedo_r2 <- model_summary$NK
    topoSlopeNoBurn_AIC <- AIC(se_model_topoSlopeNoBurn)
    
    # Spatial error with prefireMetric
    prefireMetric_anova_pvalue <- anova_prefireMetric$`p-value`[2]
    model_summary <- summary(se_model_prefireMetric,Nagelkerke = TRUE)
    
    prefireMetric_int <-  model_summary$Coef[1,1]
    prefireMetric_slope <-  model_summary$Coef[2,1]
    prefireMetric_prefireMetric <-  model_summary$Coef[3,1]
    prefireMetric_interaction <-  model_summary$Coef[4,1]
    prefireMetric_pvalint <-  model_summary$Coef[1,4]
    prefireMetric_pvalslope <-  model_summary$Coef[2,4]
    prefireMetric_pvalprefireMetric <-  model_summary$Coef[3,4]
    prefireMetric_pvalinteraction <-  model_summary$Coef[4,4]
    
    prefireMetric_psuedo_r2 <- model_summary$NK
    prefireMetric_AIC <- AIC(se_model_prefireMetric)
    
    # Spatial error with prefireMetric (no interaction)
    prefireMetricNoX_anova_pvalue <- anova_prefireMetricNoX$`p-value`[2]
    model_summary <- summary(se_model_prefireMetricNoX,Nagelkerke = TRUE)
    
    prefireMetricNoX_int <-  model_summary$Coef[1,1]
    prefireMetricNoX_slope <-  model_summary$Coef[2,1]
    prefireMetricNoX_prefireMetric <-  model_summary$Coef[3,1]
    prefireMetricNoX_pvalint <-  model_summary$Coef[1,4]
    prefireMetricNoX_pvalslope <-  model_summary$Coef[2,4]
    prefireMetricNoX_pvalprefireMetric <-  model_summary$Coef[3,4]
    
    prefireMetricNoX_psuedo_r2 <- model_summary$NK
    prefireMetricNoX_AIC <- AIC(se_model_prefireMetricNoX)
    
    prefireMetricNoXwX_anova_pvalue <- anova_prefireMetricNoXwX$`p-value`[2]
    
    # Spatial error with prefireMetric, no burn index
    model_summary <- summary(se_model_prefireMetricNoBurn,Nagelkerke = TRUE)
    
    prefireMetricNoBurn_int <-  model_summary$Coef[1,1]
    prefireMetricNoBurn_slope <-  model_summary$Coef[2,1]
    prefireMetricNoBurn_pvalint <-  model_summary$Coef[1,4]
    prefireMetricNoBurn_pvalslope <-  model_summary$Coef[2,4]
    
    prefireMetricNoBurn_psuedo_r2 <- model_summary$NK
    prefireMetricNoBurn_AIC <- AIC(se_model_prefireMetricNoBurn)
    
    # Spatial error with prefireMetric, no burn index & beam type
    prefireMetricNoBurn_beam_anova_pvalue <- anova_prefireMetricNoBurn_beam$`p-value`[2]
    model_summary <- summary(se_model_prefireMetricNoBurn_beam,Nagelkerke = TRUE)
    
    prefireMetricNoBurn_beam_int_coverage <-  model_summary$Coef[1,1]
    prefireMetricNoBurn_beam_slope_coverage <-  model_summary$Coef[2,1]
    prefireMetricNoBurn_beam_int_full <-  model_summary$Coef[3,1]
    prefireMetricNoBurn_beam_int_mixed <-  model_summary$Coef[4,1]
    prefireMetricNoBurn_beam_slope_full <-  model_summary$Coef[5,1]
    prefireMetricNoBurn_beam_slope_mixed <-  model_summary$Coef[6,1]
    
    prefireMetricNoBurn_beam_pvalint_coverage <-  model_summary$Coef[1,4]
    prefireMetricNoBurn_beam_pvalslope_coverage <-  model_summary$Coef[2,4]
    prefireMetricNoBurn_beam_pvalint_full <-  model_summary$Coef[3,4]
    prefireMetricNoBurn_beam_pvalint_mixed <-  model_summary$Coef[4,4]
    prefireMetricNoBurn_beam_pvalslope_full <-  model_summary$Coef[5,4]
    prefireMetricNoBurn_beam_pvalslope_mixed <-  model_summary$Coef[6,4]
    
    prefireMetricNoBurn_beam_psuedo_r2 <- model_summary$NK
    prefireMetricNoBurn_beam_AIC <- AIC(se_model_prefireMetricNoBurn_beam)
    
    # Spatial error with prefireMetric, no burn index & forest type
    prefireMetricNoBurn_forest_anova_pvalue <- anova_prefireMetricNoBurn_forest$`p-value`[2]
    model_summary <- summary(se_model_prefireMetricNoBurn_forest,Nagelkerke = TRUE)
    
    prefireMetricNoBurn_forest_int_conifer <-  model_summary$Coef[1,1]
    prefireMetricNoBurn_forest_slope_conifer <-  model_summary$Coef[2,1]
    prefireMetricNoBurn_forest_int_hardwood <-  model_summary$Coef[3,1]
    prefireMetricNoBurn_forest_int_mixed <-  model_summary$Coef[4,1]
    prefireMetricNoBurn_forest_slope_hardwood <-  model_summary$Coef[5,1]
    prefireMetricNoBurn_forest_slope_mixed <-  model_summary$Coef[6,1]
    
    prefireMetricNoBurn_forest_pvalint_conifer <-  model_summary$Coef[1,4]
    prefireMetricNoBurn_forest_pvalslope_conifer <-  model_summary$Coef[2,4]
    prefireMetricNoBurn_forest_pvalint_hardwood <-  model_summary$Coef[3,4]
    prefireMetricNoBurn_forest_pvalint_mixed <-  model_summary$Coef[4,4]
    prefireMetricNoBurn_forest_pvalslope_hardwood <-  model_summary$Coef[5,4]
    prefireMetricNoBurn_forest_pvalslope_mixed <-  model_summary$Coef[6,4]
    
    prefireMetricNoBurn_forest_psuedo_r2 <- model_summary$NK
    prefireMetricNoBurn_forest_AIC <- AIC(se_model_prefireMetricNoBurn_forest)
    
    
    results <- data.frame(
      response,predictor,n,
      ols_AIC,ols_int,ols_slope,ols_pvalslope,ols_pvalint,ols_r2,
      ols_anova_pvalue,olstime_AIC,olstime_int,olstime_slope,olstime_time,olstime_interaction,olstime_pvalint,olstime_pvalslope,olstime_pvaltime,olstime_pvalinteraction,olstime_r2,
      spatial_AIC,spatial_int,spatial_slope,spatial_pvalint,spatial_pvalslope,spatial_psuedo_r2,
      beam_anova_pvalue,beam_int_coverage,beam_slope_coverage,beam_int_full,beam_int_mixed,beam_slope_full,
      beam_slope_mixed,beam_pvalint_coverage,beam_pvalslope_coverage,beam_pvalint_full,beam_pvalint_mixed,beam_pvalslope_full,
      beam_pvalslope_mixed,beam_psuedo_r2,beam_AIC,
      forest_anova_pvalue,forest_int_conifer,forest_slope_conifer,forest_int_hardwood,forest_int_mixed,forest_slope_hardwood,
      forest_slope_mixed,forest_pvalint_conifer,forest_pvalslope_conifer,forest_pvalint_hardwood,forest_pvalint_mixed,forest_pvalslope_hardwood,
      forest_pvalslope_mixed,forest_psuedo_r2,forest_AIC,
      time_anova_pvalue,time_int,time_slope,time_time,time_interaction,time_pvalint,time_pvalslope,time_pvaltime,time_pvalinteraction,time_psuedo_r2,time_AIC,
      timeNoX_anova_pvalue,timeNoX_int,timeNoX_slope,timeNoX_time,timeNoX_pvalint,timeNoX_pvalslope,timeNoX_pvaltime,timeNoX_psuedo_r2,timeNoX_AIC,
      timeNoXwX_anova_pvalue,
      topoSlope_anova_pvalue,topoSlope_int,topoSlope_slope,topoSlope_topoSlope,topoSlope_interaction,topoSlope_pvalint,topoSlope_pvalslope,topoSlope_pvaltopoSlope,topoSlope_pvalinteraction,topoSlope_psuedo_r2,topoSlope_AIC,
      topoSlopeNoX_anova_pvalue,topoSlopeNoX_int,topoSlopeNoX_slope,topoSlopeNoX_topoSlope,topoSlopeNoX_pvalint,topoSlopeNoX_pvalslope,topoSlopeNoX_pvaltopoSlope,topoSlopeNoX_psuedo_r2,topoSlopeNoX_AIC,
      topoSlopeNoXwX_anova_pvalue,
      topoSlopeNoBurn_int,topoSlopeNoBurn_slope,topoSlopeNoBurn_pvalint,topoSlopeNoBurn_pvalslope,topoSlopeNoBurn_psuedo_r2,topoSlopeNoBurn_AIC,
      prefireMetric_anova_pvalue,prefireMetric_int,prefireMetric_slope,prefireMetric_prefireMetric,prefireMetric_interaction,prefireMetric_pvalint,prefireMetric_pvalslope,prefireMetric_pvalprefireMetric,prefireMetric_pvalinteraction,prefireMetric_psuedo_r2,prefireMetric_AIC,
      prefireMetricNoX_anova_pvalue,prefireMetricNoX_int,prefireMetricNoX_slope,prefireMetricNoX_prefireMetric,prefireMetricNoX_pvalint,prefireMetricNoX_pvalslope,prefireMetricNoX_pvalprefireMetric,prefireMetricNoX_psuedo_r2,prefireMetricNoX_AIC,
      prefireMetricNoXwX_anova_pvalue,
      prefireMetricNoBurn_int,prefireMetricNoBurn_slope,prefireMetricNoBurn_pvalint,prefireMetricNoBurn_pvalslope,prefireMetricNoBurn_psuedo_r2,prefireMetricNoBurn_AIC,
      prefireMetricNoBurn_beam_anova_pvalue,prefireMetricNoBurn_beam_int_coverage,prefireMetricNoBurn_beam_slope_coverage,prefireMetricNoBurn_beam_int_full,prefireMetricNoBurn_beam_int_mixed,prefireMetricNoBurn_beam_slope_full,
      prefireMetricNoBurn_beam_slope_mixed,prefireMetricNoBurn_beam_pvalint_coverage,prefireMetricNoBurn_beam_pvalslope_coverage,prefireMetricNoBurn_beam_pvalint_full,prefireMetricNoBurn_beam_pvalint_mixed,prefireMetricNoBurn_beam_pvalslope_full,
      prefireMetricNoBurn_beam_pvalslope_mixed,prefireMetricNoBurn_beam_psuedo_r2,prefireMetricNoBurn_beam_AIC,
      prefireMetricNoBurn_forest_anova_pvalue,prefireMetricNoBurn_forest_int_conifer,prefireMetricNoBurn_forest_slope_conifer,prefireMetricNoBurn_forest_int_hardwood,prefireMetricNoBurn_forest_int_mixed,prefireMetricNoBurn_forest_slope_hardwood,
      prefireMetricNoBurn_forest_slope_mixed,prefireMetricNoBurn_forest_pvalint_conifer,prefireMetricNoBurn_forest_pvalslope_conifer,prefireMetricNoBurn_forest_pvalint_hardwood,prefireMetricNoBurn_forest_pvalint_mixed,prefireMetricNoBurn_forest_pvalslope_hardwood,
      prefireMetricNoBurn_forest_pvalslope_mixed,prefireMetricNoBurn_forest_psuedo_r2,prefireMetricNoBurn_forest_AIC
    )
    output <- rbind(output,results)
    
    remove(W)
  }
}

# Write out statistics file
write.csv(output,outputFile,row.names = F)


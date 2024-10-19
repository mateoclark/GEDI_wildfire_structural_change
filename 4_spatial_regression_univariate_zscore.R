# Univariate spatial regression models

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 08/09/2024


library(spatialreg)
library(sp)
library(spdep)
library(dplyr)
library(sf)
library(foreach)
library(doParallel)

# working directory
indir <- "G:/Shared drives/CALFIRE_GEDI/Manuscripts/RQ2_wildfire_change_in_footprint_metrics"
#indir <- "E:/active/project/calfire_gedi/rq2_paper"

# Detect the number of available cores
#numCores <- detectCores() - 1
numCores <- 7

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit
#maxAfter <- 730 # 2 years

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

# distance between pairs
distance <- 45

# allowable absolute difference in elevation (m)
delta_elev <- 10

# allowable growth in tAGBD Mg/ha/year
biomass <- 99999 # not used

# version date
version <- "241012"

# Number of nearest neighbors
k <- 50  

# predictors to use (continuous)
predictors <- c("mtbs_dnbrOW", "topoSlope","timeAfter","prefireMetric","vpd","windSpeed","et")

# input .csv file with GEDI metric differences
inputFile <- paste0(indir,"/data/master_gedi_fire_difference_",distance,"m_240822.csv")

# GEDI metrics
metrics <- read.csv(paste0(indir,"/data/gedi_structure_variables_reduced_pai_240927.csv"))
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
dataFiltered <- data %>% filter(timeBefore <= maxBefore & timeAfter >= minAfter & timeAfter <= maxAfter  &
                                                     abs(delta_GEDI_elev) <= delta_elev &
                                                    delta_strct_tAGBD <= abgd_limit &
                                                    !mtbs_classP == 0) %>% 
                  select(any_of(c("GEDI_shot","post_x_teale","post_y_teale","typeForest","beamType",
                  metricsDelta,metricsPrefire,"mtbs_dnbrOW","mtbs_rdnbrW","mtbs_CBIW","timeAfter","timeDiff",
                  "topoSlope", "vpd","et","windSpeed"))) %>%
                  rename(x = post_x_teale,
                         y = post_y_teale)

# Function to rescale using z score and set outliers to NA
outlierRemoval<-function(df,z=4){
  tmp <- as.matrix(apply(as.matrix(df),2,scale))
  if (ncol(tmp)==1) {tmp[which(abs(tmp)>z)] <- NA}
  if (ncol(tmp)>1) {tmp[which(abs(tmp)>z,arr.ind=T)] <- NA}
  return(tmp)
}

# Register the parallel backend
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Define a function to process each predictor and metric
process_predictor_metric <- function(predictor, metricsDelta, dataFiltered, distance, maxBefore, minAfter, maxAfter, indir) {
  output <- data.frame()
  for (i in 1:length(metricsDelta)) {
    response <- metricsDelta[i]
    if (predictor == "prefireMetric") {
      p <- gsub("delta_strct","pre_strct",response)
      d <- data.frame(response = dataFiltered[[response]],
                      predictor = dataFiltered[[p]],
                      classForest = dataFiltered$typeForest,
                      beamType = dataFiltered$beamType,
                      x = dataFiltered$x,
                      y = dataFiltered$y
      )
    } else {
      d <- data.frame(response = dataFiltered[[response]],
                      predictor = dataFiltered[[predictor]],
                      classForest = dataFiltered$typeForest,
                      beamType = dataFiltered$beamType,
                      x = dataFiltered$x,
                      y = dataFiltered$y
      )
    }
    d <- d %>% filter(!is.infinite(response))
    d$response <- outlierRemoval(d$response) # z-score transformation and outlier removal
    d$predictor <- outlierRemoval(d$predictor) # z-score transformation and outlier removal
    d <- d %>% na.omit() 
    if ((dim(d)[1] > 1000) & (sum(d$response) != 0)) {
      n <- dim(d)[1]
      coordinates(d) <- c("x","y")
      spatial_neighbors <- knearneigh(d,k = k)
      W <- knn2nb(spatial_neighbors)
      W <- nb2listw(W,style = "W")
      se_model <- errorsarlm(formula = response ~ predictor, data = d, listw = W, method="LU")
      se_model_beam <- errorsarlm(formula = response ~ predictor * beamType, data = d, listw = W, method="LU")
      se_model_forest <- errorsarlm(formula = response ~ predictor * classForest, data = d, listw = W, method="LU")
      anova_beam <- anova(se_model,se_model_beam)
      anova_forest <- anova(se_model,se_model_forest)
      ols_model <- lm(response ~ predictor, data = d)
      ols_AIC <- AIC(ols_model)
      ols_summary <- summary(ols_model)
      ols_int <- ols_summary$coefficients[1,1]
      ols_slope <- ols_summary$coefficients[2,1]
      ols_pvalslope <- ols_summary$coefficients[1,4]
      ols_pvalint <- ols_summary$coefficients[2,4]
      ols_r2 <- ols_summary$adj.r.squared
      model_summary <- summary(se_model,Nagelkerke = TRUE)
      spatial_int <-  model_summary$Coef[1,1]
      spatial_slope <-  model_summary$Coef[2,1]
      spatial_stderrint <- model_summary$Coef[1,2]
      spatial_stderrslope <- model_summary$Coef[2,2]
      spatial_pvalint <- model_summary$Coef[1,4]
      spatial_pvalslope <- model_summary$Coef[2,4]
      spatial_psuedo_r2 <- model_summary$NK
      spatial_AIC <- AIC(se_model)
      beam_anova_pvalue <- anova_beam$`p-value`[2]
      model_summary <- summary(se_model_beam,Nagelkerke = TRUE)
      beam_int_coverage <-  model_summary$Coef[1,1]
      beam_slope_coverage <-  model_summary$Coef[2,1]
      beam_int_full <-  model_summary$Coef[3,1]
      beam_int_mixed <-  model_summary$Coef[4,1]
      beam_slope_full <-  model_summary$Coef[5,1]
      beam_slope_mixed <-  model_summary$Coef[6,1]
      beam_stderrint_coverage <-  model_summary$Coef[1,2]
      beam_stderrslope_coverage <-  model_summary$Coef[2,2]
      beam_stderrint_full <-  model_summary$Coef[3,2]
      beam_stderrint_mixed <-  model_summary$Coef[4,2]
      beam_stderrslope_full <-  model_summary$Coef[5,2]
      beam_stderrslope_mixed <-  model_summary$Coef[6,2]
      beam_pvalint_coverage <-  model_summary$Coef[1,4]
      beam_pvalslope_coverage <-  model_summary$Coef[2,4]
      beam_pvalint_full <-  model_summary$Coef[3,4]
      beam_pvalint_mixed <-  model_summary$Coef[4,4]
      beam_pvalslope_full <-  model_summary$Coef[5,4]
      beam_pvalslope_mixed <-  model_summary$Coef[6,4]
      beam_psuedo_r2 <- model_summary$NK
      beam_AIC <- AIC(se_model_beam)
      forest_anova_pvalue <- anova_forest$`p-value`[2]
      model_summary <- summary(se_model_forest,Nagelkerke = TRUE)
      forest_int_conifer <-  model_summary$Coef[1,1]
      forest_slope_conifer <-  model_summary$Coef[2,1]
      forest_int_hardwood <-  model_summary$Coef[3,1]
      forest_int_mixed <-  model_summary$Coef[4,1]
      forest_slope_hardwood <-  model_summary$Coef[5,1]
      forest_slope_mixed <-  model_summary$Coef[6,1]
      forest_stderrint_conifer <-  model_summary$Coef[1,2]
      forest_stderrslope_conifer <-  model_summary$Coef[2,2]
      forest_stderrint_hardwood <-  model_summary$Coef[3,2]
      forest_stderrint_mixed <-  model_summary$Coef[4,2]
      forest_stderrslope_hardwood <-  model_summary$Coef[5,2]
      forest_stderrslope_mixed <-  model_summary$Coef[6,2]
      forest_pvalint_conifer <-  model_summary$Coef[1,4]
      forest_pvalslope_conifer <-  model_summary$Coef[2,4]
      forest_pvalint_hardwood <-  model_summary$Coef[3,4]
      forest_pvalint_mixed <-  model_summary$Coef[4,4]
      forest_pvalslope_hardwood <-  model_summary$Coef[5,4]
      forest_pvalslope_mixed <-  model_summary$Coef[6,4]
      forest_psuedo_r2 <- model_summary$NK
      forest_AIC <- AIC(se_model_forest)
      results <- data.frame(
        response,predictor,n,
        ols_AIC,ols_int,ols_slope,ols_pvalslope,ols_pvalint,ols_r2,
        spatial_AIC,spatial_int,spatial_slope,spatial_stderrint,spatial_stderrslope,spatial_pvalint,spatial_pvalslope,spatial_psuedo_r2,
        beam_anova_pvalue,beam_int_coverage,beam_slope_coverage,beam_int_full,beam_int_mixed,beam_slope_full,beam_slope_mixed,
        beam_stderrint_coverage,beam_stderrslope_coverage,beam_stderrint_full,beam_stderrint_mixed,beam_stderrslope_full,beam_stderrslope_mixed,
        beam_pvalint_coverage,beam_pvalslope_coverage,beam_pvalint_full,beam_pvalint_mixed,beam_pvalslope_full,beam_pvalslope_mixed,
        beam_psuedo_r2,beam_AIC,
        forest_anova_pvalue,forest_int_conifer,forest_slope_conifer,forest_int_hardwood,forest_int_mixed,forest_slope_hardwood,forest_slope_mixed,
        forest_stderrint_conifer,forest_stderrslope_conifer,forest_stderrint_hardwood,forest_stderrint_mixed,forest_stderrslope_hardwood,
        forest_pvalint_conifer,forest_pvalslope_conifer,forest_pvalint_hardwood,forest_pvalint_mixed,forest_pvalslope_hardwood,
        forest_pvalslope_mixed,forest_psuedo_r2,forest_AIC
      )
      output <- rbind(output,results)
      remove(W)
    }
  }
  outputFile <- paste0(indir,"/data/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictor,"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")
  write.csv(output,outputFile,row.names = F)
}

# Run predictors with parallelization
foreach(predictor = predictors, .packages = c("dplyr", "sf", "sp", "spatialreg", "spdep")) %dopar% {
  process_predictor_metric(predictor, metricsDelta, dataFiltered, distance, maxBefore, minAfter, maxAfter, indir)
}

# Stop the cluster
stopCluster(cl)
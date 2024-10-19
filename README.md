# GEDI_wildfire_structural_change
 Code and data for paper "Changes in GEDI-based measures of forest structure after large California wildfires relative to pre-fire conditions"

Matthew L. Clark, Ph.D.  
Department of Geography, Environment, and Planning  
Sonoma State University, California USA  
matthew.clark@sonoma.edu  

**Code**  

*1_prepare_master_data_file.R*  
- Loads in GEDI structure x fire data .csv and outputs a master .csv file with information needed for analysis. Also removes shots that are in treatments or in overlapping fires.

*2_footprint_gedi_metric_difference.R*  
- Pair GEDI footprints and calculate post- minus pre-fire metric differences.

*3_gedi_repeated_measures_ANOVA_zscore.R*  
- Repeat-measures ANOVA, or mixed-effects model

*4_spatial_regression_univariate_zscore.R*  
- Univariate spatial regression models

*5_spatial_regression_multivariate_interactions_zscore_predictors.R*  
- Multivariate spatial regression models

*6_graphs.R*  
- Graphs for manuscript

 
**Data**  

*calfire_gedi_treatments.gdb.zip*  
- Treatment polygons from CALMAPPER and FACTS GIS data  
 
*gedi_structure_variables_reduced_pai_240927.csv*  
- Full GEDI metrics list, including a "selected" column to indicate which metrics to analyze.

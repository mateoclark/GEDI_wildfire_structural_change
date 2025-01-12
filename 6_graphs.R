# Graphs for paper "Changes in GEDI-based measures of forest structure after large California wildfires 1 relative to pre-fire conditions"
# Remote Sensing of Environment

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 01/12/2025


library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
library(forcats)
library(rstatix)
library(scales)


# working directory
indir <- "/GEDI_wildfire_structural_change"

# distance between pairs
distance <- 45

# allowable absolute difference in elevation (m)
delta_elev <- 10

# allowable growth in tAGBD Mg/ha/year
biomass <- 99999 # not used when set very high

# lookup
metricLookup <- data.frame(
  metric = c("strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
             "strct_mPAI_b10","strct_prop_int_btw_0_10m","strct_prop_int_below_10m","strct_RH_10","strct_RH_25",
             "strct_RH_50","strct_RH_75","strct_RH_98"),
  labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI","mPAI0to10","RE0to10m","REBelow10m",
             "RH10","RH25","RH50","RH75","RH98")
)

# GEDI metrics to use
metricList <- read.csv(paste0(indir,"/gedi_structure_variables_reduced_pai_240927.csv"))  %>% filter(selected == "Y")
metrics <- left_join(metricList, metricLookup, by="metric") %>% select(metric, labels)

###############################################################################
### Plot footprint burn severity and time since fire figure
###############################################################################

# version date
version <- "240822"

pal <- RColorBrewer::brewer.pal(8, "Set2")

# read change metrics file
file <- paste0(indir,"/master_gedi_fire_difference_",distance,"m_",version,".csv")
df <- read.csv(file)

# get wildfire burn severity classes
mtbsCodes <- data.frame(
  mtbs_classP = c(0, 1, 2, 3, 4, 5),
  FireSeverityMTBS = c("Control","Unchanged", "Low","Moderate","High","Unchanged") # treat class 5 (growth) as unchanged as small numbers
)
df <- left_join(df,mtbsCodes,by="mtbs_classP")
df$FireSeveritydNBRMillerThode <- ifelse(df$mtbs_dnbrOW >= 367, "High",ifelse(df$mtbs_dnbrOW >= 177, "Moderate",ifelse(df$mtbs_dnbrOW > 41, "Low", "Unchanged"))) 
df$FireSeveritydNBRMillerThode <- ifelse(df$FireSeverityMTBS == "Control", "Control",df$FireSeveritydNBRMillerThode)
                                                                                                                                                               
# maximum time after fire limit (days)
maxAfter <- 1500 # no limit

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

# calculate biomass growth limit based on time difference
df$abgd_limit <- (df$timeDiff/365) * biomass

data <- df %>% filter(timeBefore <= maxBefore & timeAfter >= minAfter & timeAfter <= maxAfter &
                        !is.na(mtbs_dnbrOW) &  
                        abs(delta_GEDI_elev) <= delta_elev &
                        delta_strct_tAGBD <= abgd_limit)
data$FireSeverityOrdered <- factor(data$FireSeveritydNBRMillerThode, levels = c("Control","Unchanged", "Low","Moderate","High"))
data$typeForestOrdered <- factor(data$typeForest, levels = c("conifer","hardwood","mixed"), labels = c("Conifer","Hardwood","Mixed"))
data1 <- data %>% group_by(typeForestOrdered) %>% count(FireSeverityOrdered) %>% filter(!is.na(FireSeverityOrdered))

# count of paired footprints by burn severity class and forest type
plotBurnSeverityForestClass <- ggplot(data1, aes(x=FireSeverityOrdered, y=n, fill = typeForestOrdered)) +
  geom_bar(alpha=0.8, stat="identity") +
  theme_bw() +  
  scale_fill_manual(values=pal[1:3]) +
  labs(x="dNBR Burn Severity Class",y="Count") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 5)) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.key.width = unit(.1, "in")
  ) 
#dev.new();print(plotBurnSeverityForestClass)

# count of paired footprints by time after fire, color by region
data2 <- data
data2$beamTypeOrdered <- factor(data2$beamType, levels = c("coverage","full","mixed"), labels = c("Coverage","Full Power","Mixed"))
data2$regionOrdered <- factor(data2$region, levels = c("CC","NC","S"), labels = c("Central Coast","North Coast","Sierra Nevada"))


plotTimeAfter <- ggplot(data2, aes(timeAfter,fill=typeForestOrdered)) +
  geom_histogram(bins = 50) +
  labs(x="Time Since Fire (days)",y="Count") +
  theme_bw() +
  scale_fill_manual(values=pal[1:3]) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 4)) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10,angle=45,hjust=1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.key.width = unit(.1, "in")
  ) +
  geom_vline(xintercept = c(365,730), linetype = "dashed", color = "black")
#dev.new();print(plotTimeAfter)

# combine graphs for figure
combinedFootprint <- ggarrange(plotBurnSeverityForestClass,plotTimeAfter, 
                          labels = c("A", "B"),
                          #hjust = -1,
                          common.legend = TRUE,
                          legend = "bottom",
                          ncol = 2,
                          nrow = 1,
                          align = "h")
dev.new();print(combinedFootprint)

outFile <- paste0(indir,"/gedi_paired_footprints_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",version,".png")
ggsave(outFile, width = 6, height = 3, units = "in", dpi=600)

###############################################################################
### Tables for results section
###############################################################################

data2$timeSinceFire <- ifelse(data2$timeAfter <= 365,"<= 1 year",ifelse(data2$timeAfter <= 730, "1-2 years", "> 2 years"))
table(data2$timeSinceFire)
table(data2$typeForestOrdered,data2$timeSinceFire)

table(data2$FireSeveritydNBRMillerThode)
table(data2$typeForestOrdered,data2$FireSeveritydNBRMillerThode)

table(data2$region)/sum(table(data2$region)) * 100
table(data2$beamType,data2$FireSeveritydNBRMillerThode)

summary(data2$mtbs_dnbrOW)

mean(data2$distance)
sd(data2$distance)

# filter data to include only data needed
dataFiltered <- data %>% select(all_of(c("topoSlope","vpd","windSpeed","et","typeForest","beamType",paste0("delta_",metrics$metric),paste0("pre_",metrics$metric),"mtbs_dnbrOW","FireSeveritydNBRMillerThode","timeAfter")))

dataFiltered <- dataFiltered %>% 
  rename(
  dNBR = mtbs_dnbrOW,
  H = pre_strct_H,
  nmode = pre_strct_nmode,
  VDR = pre_strct_VDR,
  Cover = pre_strct_cover,
  tAGBD = pre_strct_tAGBD,
  tPAI = pre_strct_tPAI,          
  mPAI0to10 = pre_strct_mPAI_b10,
  REBelow10m = pre_strct_prop_int_below_10m,
  RH25 = pre_strct_RH_25,
  RH50 = pre_strct_RH_50,
  RH75 = pre_strct_RH_75,
  RH98 = pre_strct_RH_98
  )

# Function to calculate min, max, average, and standard deviation for selected columns
summarizeColumns <- function(data, columns) {
  summary <- data %>%
    select(all_of(columns)) %>%
    summarise_all(list(
      min = ~min(., na.rm = TRUE),
      max = ~max(., na.rm = TRUE),
      mean = ~mean(., na.rm = TRUE),
      sd = ~sd(., na.rm = TRUE)
    ))
  
  # Reshape the summary to have min, max, mean, and sd in columns
  summary_long <- summary %>%
    pivot_longer(cols = everything(), 
                 names_to = c("variable", "stat"), 
                 names_sep = "_") %>%
    pivot_wider(names_from = stat, values_from = value)
  
  return(summary_long)
}

dataFiltered <- dataFiltered %>% rename(
  FHD = H,
  aVDR = VDR
)

# statistics for both control and burned area
summaryResult <- summarizeColumns(dataFiltered, c("dNBR","timeAfter","topoSlope","et","vpd","windSpeed",
                                                  "FHD","nmode","aVDR","Cover","tAGBD","tPAI","mPAI0to10",          
                                                  "REBelow10m","RH25","RH50","RH75","RH98"))
# Print the summary
print(summaryResult)

# within burned area statistics
dataFilteredFire <- dataFiltered %>% filter(FireSeveritydNBRMillerThode != "Control")
summaryResultFire <- summarizeColumns(dataFilteredFire, c("dNBR","timeAfter","topoSlope","et","vpd","windSpeed",
                                                          "FHD","nmode","aVDR","Cover","tAGBD","tPAI","mPAI0to10",          
                                                          "REBelow10m","RH25","RH50","RH75","RH98"))

# Print the summary
print(summaryResultFire)

###############################################################################
### Plot NLME SPATIAL REPEAT-MEASURES ANOVA 
###############################################################################

# version date
version <- "240822"

pvalue = 0.01
plot_list = list()

# change in elevation limit
# delta_elev = 10
delta_elev = 10

# y-axis limits
ylimits = c(-0.5,0.5)

# Under 1 year since fire
maxAfter <- 365 
minAfter <- 0
maxBefore <- 730 # 2 years
file <- paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")
df <- read.csv(file) %>% filter(metric %in% metrics$metric)

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df1 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df1) <- c("metric","ind","pval","pval-type","zscore","effect")
df1$significant <- ifelse(df1$pval <= pvalue,"Sig.","Not Sig.")
df1$effectFactor <- factor(df1$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df1$metricFactor <- factor(df1$metric, 
                           levels = c("strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI","strct_mPAI_b10",
                                      "strct_prop_int_below_10m","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98"))
pal <- RColorBrewer::brewer.pal(3, "Set2")

if (length(unique(df1$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot1 <- ggplot(df1, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(
    position = "dodge", stat = "identity",
    pattern_density = 0.1,
    pattern_fill = "black",
    pattern_angle = 45,
    pattern_spacing = 0.02,
  ) +
  theme_bw() +
  labs(title = "A. Within 1 year of fire", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",expression(mPAI["0-10m"]),
                              expression(RE["<10m"]),"RH25","RH50","RH75","RH98")) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.key.width = unit(.1, "in")
  )

#dev.new();print(plot1)
plot_list[[1]] = plot1

# 1-2 years since fire
maxAfter <- 730
minAfter <- 365
maxBefore <- 730 # 2 years
file <- paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")
df <- read.csv(file) %>% filter(metric %in% metrics$metric)

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df2 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df2) <- c("metric","ind","pval","pval-type","zscore","effect")
df2$significant <- ifelse(df2$pval <= pvalue,"Sig.","Not Sig.")
df2$effectFactor <- factor(df2$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df2$metricFactor <- factor(df2$metric, 
                           levels = c("strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI","strct_mPAI_b10",
                                      "strct_prop_int_below_10m","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98"))
if (length(unique(df2$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot2 <- ggplot(df2, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  theme_bw() +
  labs(title = "B. 1-2 years since fire", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",expression(mPAI["<10m"]),
                              expression(RE["<10m"]),"RH25","RH50","RH75","RH98")) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text( size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot2)
plot_list[[2]] = plot2

# Over 2 years since fire
maxAfter <- 1500
minAfter <- 730
maxBefore <- 730 # 2 years
file <- paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")
df <- read.csv(file) %>% filter(metric %in% metrics$metric)

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df3 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df3) <- c("metric","ind","pval","pval-type","zscore","effect")
df3$significant <- ifelse(df3$pval <= pvalue,"Sig.","Not Sig.")
df3$effectFactor <- factor(df3$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df3$metricFactor <- factor(df3$metric, 
                           levels = c("strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI","strct_mPAI_b10",
                                      "strct_prop_int_below_10m","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98"))
if (length(unique(df3$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot3 <- ggplot(df3, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  theme_bw() +
  labs(title = "C. 2-3 years since fire", y = "Normalized Effect") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",expression(mPAI["<10m"]),
                              expression(RE["<10m"]),"RH25","RH50","RH75","RH98")) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text( size = 12),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot3)
plot_list[[3]] = plot3

# All time since fire
maxAfter <- 1500 
minAfter <- 0
maxBefore <- 730 # 2 years
file <- paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")
df <- read.csv(file) %>% filter(metric %in% metrics$metric)

df$effectTreatment <- df$spatialInterceptEst + df$spatialTreatmentEst
df$effectTime <- df$spatialInterceptEst + df$aspatialTimeEst
df$effectInteraction <- df$spatialInterceptEst + df$spatialTreatmentEst + df$aspatialTimeEst + df$spatialTreatmentTimeEst

df4 <- cbind(stack(df[1]), stack(df[18:20]),stack(df[22:24]))
colnames(df4) <- c("metric","ind","pval","pval-type","zscore","effect")
df4$significant <- ifelse(df4$pval <= pvalue,"Sig.","Not Sig.")
df4$effectFactor <- factor(df4$effect, levels = c("effectTreatment", "effectTime","effectInteraction"),
                           labels = c("Burned", "Timing","Interaction"))
df4$metricFactor <- factor(df4$metric, 
                           levels = c("strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI","strct_mPAI_b10",
                                      "strct_prop_int_below_10m","strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98"))
if (length(unique(df4$significant)) == 1) pattern = "none" else pattern = c("stripe", "none")

# with hatch for significant
plot4 <- ggplot(df4, aes(x = metricFactor, y = zscore, fill = effectFactor, pattern = significant)) +
  geom_bar_pattern(position = "dodge", stat = "identity",
                   pattern_density = 0.1,
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_spacing = 0.02) +
  theme_bw() +
  labs(title = "D. All years since fire") +
  scale_pattern_manual(values = pattern) + 
  scale_fill_manual(values=c(pal[2],pal[1],pal[3])) +
  guides(pattern = "none",
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  coord_cartesian(ylim = ylimits) +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",expression(mPAI["<10m"]),
                              expression(RE["<10m"]),"RH25","RH50","RH75","RH98")) +
  theme(
    plot.title = element_text( size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  theme(legend.position = "none")

#dev.new();print(plot4)
plot_list[[4]] = plot4

combined <- ggarrange(plotlist = plot_list,
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 2,
                      align = "hv")

dev.new();print(combined)

outFile <- paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",pvalue,"_",version,".png")
ggsave(outFile, width = 7, height = 6, units = "in", dpi=600)

###############################################################################
### Gather significance stats for ANOVA
###############################################################################
df1$period <- "<1 year"
df2$period <- "1-2 years"
df3$period <- "2-3 years"
df4$period <- "All years"
df_combined <- rbind(df1,df2,df3,df4)
outFile <- paste0(indir,"/nlme_repeated_anova_spatial_test_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",pvalue,"_",version,".csv")
write.csv(df_combined, outFile, row.names=F)

###############################################################################
### DENSITY -- Difference in metrics between pre- and post-fire
###############################################################################

# version date
version <- "240822"

# read change metrics file
file <- paste0(indir,"/master_gedi_fire_difference_",distance,"m_",version,".csv")
df <- read.csv(file)

# get wildfire burn severity classes
mtbsCodes <- data.frame(
  mtbs_classP = c(0, 1, 2, 3, 4),
  FireSeverityMTBS = c("Control","Unchanged", "Low","Moderate","High")
)
df <- left_join(df,mtbsCodes,by="mtbs_classP")
df$FireSeveritydNBRMillerThode <- ifelse(df$mtbs_dnbrOW >= 367, "High",ifelse(df$mtbs_dnbrOW >= 177, "Moderate",ifelse(df$mtbs_dnbrOW > 41, "Low", "Unchanged"))) 
df$FireSeveritydNBRMillerThode <- ifelse(df$FireSeverityMTBS == "Control", "Control",df$FireSeveritydNBRMillerThode)

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

# allowable absolute difference in elevation (m)
delta_elev <- 10

# allowable growth in tAGBD Mg/ha/year
biomass <- 99999

# calculate biomass growth limit based on time difference
df$abgd_limit <- (df$timeDiff/365) * biomass

# filter data
data <- df %>% filter(timeBefore <= maxBefore & timeAfter >= minAfter & timeAfter <= maxAfter &
                      !is.na(mtbs_dnbrOW) &
                      abs(delta_GEDI_elev) <= delta_elev &
                      delta_strct_tAGBD <= abgd_limit
                      )

data$FireSeverityOrdered <- factor(data$FireSeveritydNBRMillerThode, levels = c("Control","Unchanged", "Low","Moderate","High"))
data$typeForestOrdered <- factor(data$typeForest, levels = c("conifer","hardwood","mixed"), labels = c("Conifer","Hardwood","Mixed"))

# GEDI metrics to use
metricList <- read.csv(paste0(indir,"/gedi_structure_variables_reduced_pai_240927.csv"))  %>% filter(selected == "Y")
metrics <- left_join(metricList, metricLookup, by="metric") %>% select(metric, labels)

# prepare data
postVars <- c(paste("post",metrics$metric,sep="_"),"typeForestOrdered","FireSeverityOrdered","GEDI_utmX","GEDI_utmY")
d1 <- data %>% select(!!postVars)
colnames(d1)<- c(metrics$labels,"typeForestOrdered","FireSeverityOrdered","GEDI_utmX","GEDI_utmY")
ncol <- dim(d1)[2]
d1<-  cbind(d1[(ncol-3):ncol], stack(d1[1:(ncol-4)]))
d1$time <- "Post-fire"

preVars <- c(paste("pre",metrics$metric,sep="_"),"typeForestOrdered","FireSeverityOrdered","pre_GEDI_utmX","pre_GEDI_utmY")
d2 <- data %>% select(!!preVars)
colnames(d2)<- c(metrics$labels,"typeForestOrdered","FireSeverityOrdered","GEDI_utmX","GEDI_utmY")
ncol <- dim(d2)[2]
d2<- cbind(d2[(ncol-3):ncol], stack(d2[1:(ncol-4)]))
d2$time <- "Pre-fire"

d <- rbind(d1,d2)
d <- d %>% rename(metric = ind)

d$timeOrdered <- factor(d$time,levels = c("Pre-fire","Post-fire"))

# Function to rescale using z score and filter outliers with abs(z) > 4
zscore<-function(df,z=4){
  tmp <- as.matrix(apply(as.matrix(df),2,scale))
  if (ncol(tmp)==1) {tmp[which(abs(tmp)>z)] <- NA}
  if (ncol(tmp)>1) {tmp[which(abs(tmp)>z,arr.ind=T)] <- NA}
  return(tmp)
}

detailedLabel <- function(m){
  l <- m
  if (m == "mPAI0to10") {
    l = expression(paste(mPAI["0-10m"], " (Mean Plant Area Index below 10-m height, or low-stature fuels)"))
  }else if (m == "REBelow10m") {
    l = expression(paste(RE["<10m"], " (Proportion of energy below 10-m height)"))
  }else if (m == "FHD") {
    l = "FHD (Foliage height diversity)"
  }else if (m == "Cover") {
    l = "Cover (Total canopy cover, percent)"
  }else if (m == "RH98") {
    l = "RH98 (Top-of-canopy height, meters)"
  }else if (m == "RH75") {
    l = "RH75 (Relative height at 75% of returned energy, meters)"
  }else if (m == "RH50") {
    l = "RH50 (Relative height at 50% of returned energy, meters)"
  }else if (m == "RH25") {
    l = "RH25 (Relative height at 25% of returned energy, meters)"
  }else if (m == "aVDR") {
    l = "aVDR (aboveground Vertical Distribution Ratio)"
  }else if (m == "nmode") {
    l = "nmode (Number of detected modes in a waveform)"
  }else if (m == "tPAI") {
    l = "tPAI (Total Plant Area Index)"
  }else if (m == "tAGBD") {
    l = "tAGBD (Total aboveground biomass density, Mg/ha)"
  }else {
    l = m
  }
}

# make plot for select metrics
metricsToPlot <- c("FHD","Cover","RH98","mPAI0to10","REBelow10m")

pal <- RColorBrewer::brewer.pal(5, "Set2")

plot_list <- list()
for (i in 1:length(metricsToPlot)) {
  m <- metricsToPlot[i]
  l <- detailedLabel(m)
  dataSubset = d %>% filter(metric == m)
  dataSubsetFiltered <- dataSubset %>% mutate(zscore = zscore(dataSubset$values)) %>% na.omit() %>% distinct(GEDI_utmX,GEDI_utmY, .keep_all=T)
  if (m == "Cover") dataSubsetFiltered$values = dataSubsetFiltered$values * 100
  max <- summary(dataSubsetFiltered$values)[6]
  if (i == 1){
    p <- ggplot(dataSubsetFiltered , aes(values, fill=timeOrdered)) + 
      facet_wrap(~FireSeverityOrdered,ncol=5) + 
      theme_bw() +
      scale_fill_manual(values=c(pal[3],pal[2])) +
      geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 3)) +
      labs(x=l,y="Density") +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            legend.text = element_text(size = 12),
            legend.title = element_blank(),
            strip.text = element_text(size = 10,margin = margin(t = 2, b = 2))
      )
  } else {
    if (max >= 5){
      p <- ggplot(dataSubsetFiltered , aes(values, fill=timeOrdered)) + 
        facet_wrap(~FireSeverityOrdered,ncol=5) + 
        theme_bw() +
        scale_fill_manual(values=c(pal[3],pal[2])) +
        geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
        scale_x_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01),breaks = pretty_breaks(n = 3)) +
        labs(x=l,y="Density") +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              strip.text = element_blank()
        )
    } else {
      p <- ggplot(dataSubsetFiltered , aes(values, fill=timeOrdered)) + 
        facet_wrap(~FireSeverityOrdered,ncol=5) + 
        theme_bw() +
        scale_fill_manual(values=c(pal[3],pal[2])) +
        geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
        scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 3)) +
        labs(x=l,y="Density") +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              strip.text = element_blank()
        )
    }
  }
 

  #dev.new();print(p)
  plot_list[[i]] = p

}

# combine plots into final figure
combined <- ggarrange(plotlist = plot_list,
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 1, nrow = length(metricsToPlot),
                      heights = c(1.2,1,1,1,1,1),
                      align = "v")

dev.new();print(combined)

outFile <- paste0(indir,"/severity_difference_select_metrics_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",version,".png")
ggsave(outFile, width = 7, height = 6, units = "in", dpi=600)


# make plot for supplemental extra metrics
metricsToPlot <- c("aVDR","nmode","tPAI","tAGBD")

pal <- RColorBrewer::brewer.pal(5, "Set2")

plot_list <- list()
for (i in 1:length(metricsToPlot)) {
  m <- metricsToPlot[i]
  l <- detailedLabel(m)
  dataSubset = d %>% filter(metric == m)
  dataSubsetFiltered <- dataSubset %>% mutate(zscore = zscore(dataSubset$values)) %>% na.omit() %>% distinct(GEDI_utmX,GEDI_utmY, .keep_all=T)
  if (m == "Cover") dataSubsetFiltered$values = dataSubsetFiltered$values * 100
  max <- summary(dataSubsetFiltered$values)[6]
  if (i == 1){
    p <- ggplot(dataSubsetFiltered , aes(values, fill=timeOrdered)) + 
      facet_wrap(~FireSeverityOrdered,ncol=5) + 
      theme_bw() +
      scale_fill_manual(values=c(pal[3],pal[2])) +
      geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 3)) +
      labs(x=l,y="Density") +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            legend.text = element_text(size = 12),
            legend.title = element_blank(),
            strip.text = element_text(size = 10,margin = margin(t = 2, b = 2))
      )
  } else {
    if (max >= 5){
      p <- ggplot(dataSubsetFiltered , aes(values, fill=timeOrdered)) + 
        facet_wrap(~FireSeverityOrdered,ncol=5) + 
        theme_bw() +
        scale_fill_manual(values=c(pal[3],pal[2])) +
        geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
        scale_x_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01),breaks = pretty_breaks(n = 3)) +
        labs(x=l,y="Density") +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              strip.text = element_blank()
        )
    } else {
      p <- ggplot(dataSubsetFiltered , aes(values, fill=timeOrdered)) + 
        facet_wrap(~FireSeverityOrdered,ncol=5) + 
        theme_bw() +
        scale_fill_manual(values=c(pal[3],pal[2])) +
        geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
        scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 3)) +
        labs(x=l,y="Density") +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              strip.text = element_blank()
        )
    }
  }
  
  
  #dev.new();print(p)
  plot_list[[i]] = p
  
}

# combine plots into final figure
combined <- ggarrange(plotlist = plot_list,
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 1, nrow = length(metricsToPlot),
                      heights = c(1.2,1,1,1,1,1),
                      align = "v")

dev.new();print(combined)

outFile <- paste0(indir,"/severity_difference_extra_metrics_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",version,".png")
ggsave(outFile, width = 7, height = 6, units = "in", dpi=600)


# make plot for RH metrics
metricsToPlot <- c("RH25","RH50","RH75","RH98")
#metricsToPlot <- c("RH10","RH25","RH50","RH75","RH98")

pal <- RColorBrewer::brewer.pal(5, "Set2")
plot_list <- list()
for (i in 1:length(metricsToPlot)) {
  m <- metricsToPlot[i]
  l <- detailedLabel(m)
  dataSubset = d %>% filter(metric == m)
  #print(summary(dataSubset$values))
  dataSubsetFiltered <- dataSubset %>% mutate(zscore = zscore(dataSubset$values)) %>% na.omit() %>% distinct(GEDI_utmX,GEDI_utmY, .keep_all=T)
  max <- summary(dataSubsetFiltered$values)[6]
  p <- ggplot(dataSubsetFiltered , aes(values, fill=timeOrdered)) + 
      facet_wrap(~FireSeverityOrdered,ncol=5) + 
      theme_bw() +
      scale_fill_manual(values=c(pal[3],pal[2])) +
      geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
      scale_x_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 5)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01),breaks = pretty_breaks(n = 3)) +
      labs(x=l,y="Density") +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            legend.text = element_text(size = 12),
            legend.title = element_blank(),
            strip.text = element_text(size = 10,margin = margin(t = 2, b = 2))
      )
            
  #dev.new();print(p)
  plot_list[[i]] = p
  
}

# combine plots into final figure
combined <- ggarrange(plotlist = plot_list,
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 1, nrow = length(metricsToPlot),
                      heights = c(1.2,1,1,1,1,1),
                      align = "v")

dev.new();print(combined)

outFile <- paste0(indir,"/severity_difference_rh_metrics_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",version,".png")
ggsave(outFile, width = 7, height = 6, units = "in", dpi=600)

###############################################################################
### Plots of pre-fire metric distributions by forest type 
###############################################################################

metricsToPlot <- unique(d$metric)
pal <- RColorBrewer::brewer.pal(5, "Set2")
plot_list <- list()
for (i in 1:length(metricsToPlot)) {
  m <- metricsToPlot[i]
  l <- m
  if (m == "mPAI0to10") {
    l = expression(mPAI["0-10m"])
  }else if (m == "RE0to10m") {
    l = expression(WI["0-10m"])
  }else if (m == "REBelow10m") {
    l = expression(RE["<10m"])
  }
  
  dataSubset = d %>% filter(metric == m & timeOrdered == "Pre-fire") %>% na.omit()
  max <- summary(dataSubset$values)[6]
  print(m)
  print(max)
  print(table(dataSubset$typeForestOrdered))

  if (i < 3){
    p <- ggplot(dataSubset, aes(values, fill=typeForestOrdered)) + 
      #facet_wrap(~typeForestOrdered,ncol=3) + 
      theme_bw() +
      scale_fill_manual(values=pal) +
      geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
      scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 2)) +
      labs(x=l,y="Density") +
      theme(axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            legend.text = element_text(size = 12),
            legend.title = element_blank(),
            strip.text = element_text(size = 10,margin = margin(t = 2, b = 2))
      )
  } else {
    if (max >= 5){
      p <- ggplot(dataSubset, aes(values, fill=typeForestOrdered)) + 
        #facet_wrap(~typeForestOrdered,ncol=3) + 
        theme_bw() +
        scale_fill_manual(values=pal) +
        geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
        scale_x_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.01),breaks = pretty_breaks(n = 2)) +
        labs(x=l,y="Density") +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              strip.text = element_blank()
        )
    } else {
      p <- ggplot(dataSubset, aes(values, fill=typeForestOrdered)) + 
        #facet_wrap(~typeForestOrdered,ncol=3) + 
        theme_bw() +
        scale_fill_manual(values=pal) +
        geom_density(aes(y=after_stat(density)), color='gray50', alpha=0.50, position = "identity") +
        scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 2)) +
        labs(x=l,y="Density") +
        theme(axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 8),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_blank(),
              legend.text = element_text(size = 12),
              legend.title = element_blank(),
              strip.text = element_blank()
        )
    }
  }
  
  
  #dev.new();print(p)
  plot_list[[i]] = p
  
}

# combine plots into final figure
combined <- ggarrange(plotlist = plot_list,
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 6,
                      align = "v")

dev.new();print(combined)

outFile <- paste0(indir,"/forest_type_pre-fire_metrics_distributions_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",version,".png")
ggsave(outFile, width = 7, height = 7, units = "in", dpi=600)


###########################################################################################################
### Spatial regression of change in pre- and post-fire univariate predictor r2 and slopes
##########################################################################################################

# version date
version <- "240822"

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit
#maxAfter <- 365 # 1 year

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years
#maxBefore <- 365 # 1 year

# elevation change limit
delta_elev <- 10

# biomass threshold
biomass <- 99999

# get regression results table
predictors <- c("mtbs_dnbrOW", "topoSlope","timeAfter","prefireMetric","vpd","windSpeed","et")
dfdNBR <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[1],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dftopoSlope <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[2],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dftimeAfter <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[3],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfpreMetric <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[4],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfvpd <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[5],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfwindSpeed <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[6],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfet <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[7],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
                                           
df <- rbind(dfdNBR,dftimeAfter,dftopoSlope,dfvpd,dfet,dfwindSpeed,dfpreMetric)

df <- df %>% filter(response %in% paste0("delta_",metrics$metric))

df$metricFactor <- factor(df$response, 
                          levels = c("delta_strct_H","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                                     "delta_strct_mPAI_b10","delta_strct_prop_int_below_10m",
                                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98"))

#df2 <- df %>% select(ols_slope,spatial_slope,spatial_stderrslope,spatial_pvalslope,predictor,metricFactor)
df2 <- df %>% select(ols_slope,spatial_slope,spatial_pvalslope,predictor,metricFactor)
df2$significant <- ifelse(df2$spatial_pvalslope <= 0.01,"Sig.","Not Sig.")
df2$modelFactor <- factor(df2$predictor, 
                          levels = c("mtbs_dnbrOW","timeAfter","topoSlope","prefireMetric","et","vpd","windSpeed"),
                          labels = c("Landsat burn severity (dNBR)","Elapsed Time","Topographic Slope","Pre-fire Structure","Evapotranspiration","Vapor Pressure Deficit","Wind Speed"))
pal <- RColorBrewer::brewer.pal(8, "Set2")
pal1 <- c(pal[1:3],pal[8],pal[4:5],pal[7])

df2$metricFactor <- fct_rev(df2$metricFactor)

plotSlope <- ggplot(df2,aes(x=spatial_slope,y=metricFactor,shape = significant,color = modelFactor)) +
  theme_bw() +
  scale_color_manual(values=pal1) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(breaks = seq(-0.8, 0.4, by = 0.2)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2), size = 1.5) +
  labs(x = "Normalized slope") +
  guides(shape = "none") +
  scale_y_discrete(labels = c("RH98","RH75","RH50","RH25",
                              expression(RE["<10m"]),expression(mPAI["0-10m"]),
                              "tPAI","tAGBD","Cover","aVDR","nmode","FHD")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")


#dev.new();print(plotSlope)

# Plot of r2
df3 <- df %>% select(ols_r2,spatial_psuedo_r2,predictor,metricFactor)
df3$modelFactor <- factor(df3$predictor, 
                          levels = c("mtbs_dnbrOW","timeAfter","topoSlope","prefireMetric","et","vpd","windSpeed"),
                          labels = c("Landsat burn severity (dNBR)","Elapsed Time","Topographic Slope","Pre-fire Structure","Evapotranspiration","Vapor Pressure Deficit","Wind Speed"))
plotR2 <- ggplot(df3,aes(x=metricFactor,y=spatial_psuedo_r2,fill = modelFactor)) +
  geom_bar(stat="identity",position=position_dodge()) +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1)) +
  scale_fill_manual(values=pal1) +
  ylab(expression("Pseudo r "^2)) +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                              expression(mPAI["0-10m"]),expression(RE["<10m"]),
                              "RH25","RH50","RH75","RH98")) +
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(.1, "in") 
  )

#dev.new();print(plotR2)


# combine plots into final figure
combined <- ggarrange(plotR2,plotSlope,
                      labels = c("A", "B"),
                      common.legend =TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 1)

dev.new();print(combined)

outFile <- paste0(indir,"/regression_models_univariate_r2_slopes_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".png")
ggsave(outFile, width = 8, height = 5, units = "in", dpi=600)

outFile <- paste0(indir,"/regression_models_univariate_statistics_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv")
write.csv(df,outFile, row.names = F)

###############################################################################
### Plot of spatial regression plot slopes
###############################################################################

# version date
version <- "240822"

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

# get regression results table
df <- read.csv(paste0(indir,"/master_gedi_fire_difference_spatialreg_interactions_mtbs_dnbrOW_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
df <- df %>% filter(response %in% paste0("delta_",metrics$metric))

df$metricFactor <- factor(df$response, 
                          levels = c("delta_strct_H","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                                     "delta_strct_mPAI_b10", "delta_strct_prop_int_below_10m",
                                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98")
                          )
pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:3])

# Plot of slopes for dNBR and dNBR + forest type model
df4 <- df %>% select(spatial_slope,forest_slope_conifer,forest_slope_hardwood,forest_slope_mixed,metricFactor)
df4$forest_slope_hardwood <- df4$forest_slope_hardwood + df4$forest_slope_conifer
df4$forest_slope_mixed <- df4$forest_slope_mixed + df4$forest_slope_conifer

df4_stacked <- stack(df4,select = -metricFactor) %>% rename(slope = values,model = ind)
df4_stacked$response <- rep(df4$metricFactor,4)
df4_stacked$modelFactor <- factor(df4_stacked$model, 
                                  levels = c("spatial_slope","forest_slope_conifer","forest_slope_hardwood","forest_slope_mixed"),
                                  labels = c("Baseline model","Conifer Forest","Hardwood Forest","Mixed Forest"))

df5 <- df %>% select(spatial_pvalslope,forest_pvalslope_conifer,forest_pvalslope_hardwood,forest_pvalslope_mixed,metricFactor)
df5_stacked <- stack(df5,select = -metricFactor) %>% rename(pval = values,model = ind)
df5_stacked$response <- rep(df5$metricFactor,4)
df5_stacked$significant <- ifelse(df5_stacked$pval <= 0.01,"Sig.","Not Sig.")

df_stacked <- cbind(df4_stacked,df5_stacked$significant) %>% rename(significant = "df5_stacked$significant")
df_stacked$response <- fct_rev(df_stacked$response)

plotForestSlope <- ggplot(df_stacked,aes(x=slope,y=response,shape = significant,color = modelFactor)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2), size = 1.5) +
  theme_bw() +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = c("RH98","RH75","RH50","RH25",
                              expression(RE["<10m"]),expression(mPAI["0-10m"]),
                              "tPAI","tAGBD","Cover","aVDR","nmode","FHD")) +
  labs(x = expression("Normalized effect on " * Delta * " Structure"),
       title = "dNBR") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

#dev.new();print(plotForestSlope)

# Plot of slopes for dNBR and dNBR + beam type model
df6 <- df %>% select(spatial_slope,beam_slope_coverage,beam_slope_full,beam_slope_mixed,metricFactor)
df6$beam_slope_full <- df6$beam_slope_full + df6$beam_slope_coverage
df6$beam_slope_mixed <- df6$beam_slope_mixed + df6$beam_slope_coverage

df6_stacked <- stack(df6,select = -metricFactor) %>% rename(slope = values,model = ind)
df6_stacked$response <- rep(df6$metricFactor,4)
df6_stacked$modelFactor <- factor(df6_stacked$model, 
                                  levels = c("spatial_slope","beam_slope_coverage","beam_slope_full","beam_slope_mixed"),
                                  labels = c("Baseline model","Coverage Beam","Full Power Beam","Mixed Beam"))

df7 <- df %>% select(spatial_pvalslope,beam_pvalslope_coverage,beam_pvalslope_full,beam_pvalslope_mixed,metricFactor)
df7_stacked <- stack(df7,select = -metricFactor) %>% rename(pval = values,model = ind)
df7_stacked$response <- rep(df7$metricFactor,4)
df7_stacked$significant <- ifelse(df7_stacked$pval <= 0.01,"Sig.","Not Sig.")

df_stacked <- cbind(df6_stacked,df7_stacked$significant) %>% rename(significant = "df7_stacked$significant")
df_stacked$response <- fct_rev(df_stacked$response)

plotBeamSlope <- ggplot(df_stacked,aes(x=slope,y=response,shape = significant,color = modelFactor)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2), size = 1.5) +
  theme_bw() +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = c("RH98","RH75","RH50","RH25",
                              expression(RE["<10m"]),expression(mPAI["0-10m"]),
                              "tPAI","tAGBD","Cover","aVDR","nmode","FHD")) +
  labs(x = expression("Normalized effect on " * Delta * " Structure"),
       title = "dNBR") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

#dev.new();print(plotBeamSlope)


# Plot of slopes for Prefire Metrics + forest type model
df6 <- df %>% select(prefireMetricNoBurn_slope,prefireMetricNoBurn_forest_slope_conifer,prefireMetricNoBurn_forest_slope_hardwood,prefireMetricNoBurn_forest_slope_mixed,metricFactor)
df6$prefireMetricNoBurn_forest_slope_hardwood <- df6$prefireMetricNoBurn_forest_slope_hardwood + df6$prefireMetricNoBurn_forest_slope_conifer
df6$prefireMetricNoBurn_forest_slope_mixed <- df6$prefireMetricNoBurn_forest_slope_mixed + df6$prefireMetricNoBurn_forest_slope_conifer

df6_stacked <- stack(df6,select = -metricFactor) %>% rename(slope = values,model = ind)
df6_stacked$response <- rep(df6$metricFactor,4)
df6_stacked$modelFactor <- factor(df6_stacked$model, 
                                  levels = c("prefireMetricNoBurn_slope","prefireMetricNoBurn_forest_slope_conifer","prefireMetricNoBurn_forest_slope_hardwood","prefireMetricNoBurn_forest_slope_mixed"),
                                  labels = c("Baseline model","Conifer Forest","Hardwood Forest","Mixed Forest"))


df7 <- df %>% select(prefireMetricNoBurn_pvalslope,prefireMetricNoBurn_forest_pvalslope_conifer,prefireMetricNoBurn_forest_pvalslope_hardwood,prefireMetricNoBurn_forest_pvalslope_mixed,metricFactor)
df7_stacked <- stack(df7,select = -metricFactor) %>% rename(pval = values,model = ind)
df7_stacked$response <- rep(df7$metricFactor,4)
df7_stacked$significant <- ifelse(df7_stacked$pval <= 0.01,"Sig.","Not Sig.")

df_stacked <- cbind(df6_stacked,df7_stacked$significant) %>% rename(significant = "df7_stacked$significant")
df_stacked$response <- fct_rev(df_stacked$response)

plotPrefireForestSlope <- ggplot(df_stacked,aes(x=slope,y=response,shape = significant,color = modelFactor)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2), size = 1.5) +
  theme_bw() +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = c("RH98","RH75","RH50","RH25",
                              expression(RE["<10m"]),expression(mPAI["0-10m"]),
                              "tPAI","tAGBD","Cover","aVDR","nmode","FHD")) +
  labs(x = expression("Normalized effect on " * Delta * " Structure"),
       title = "Pre-fire Structure") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

#dev.new();print(plotPrefireForestSlope)

# Plot of slopes for Prefire Metrics + beam type model
df8 <- df %>% select(prefireMetricNoBurn_slope,prefireMetricNoBurn_beam_slope_coverage,prefireMetricNoBurn_beam_slope_full,prefireMetricNoBurn_beam_slope_mixed,metricFactor)
df8$prefireMetricNoBurn_beam_slope_full <- df8$prefireMetricNoBurn_beam_slope_full + df8$prefireMetricNoBurn_beam_slope_coverage
df8$prefireMetricNoBurn_beam_slope_mixed <- df8$prefireMetricNoBurn_beam_slope_mixed + df8$prefireMetricNoBurn_beam_slope_coverage

df8_stacked <- stack(df8,select = -metricFactor) %>% rename(slope = values,model = ind)
df8_stacked$response <- rep(df8$metricFactor,4)
df8_stacked$modelFactor <- factor(df8_stacked$model, 
                                  levels = c("prefireMetricNoBurn_slope","prefireMetricNoBurn_beam_slope_coverage","prefireMetricNoBurn_beam_slope_full","prefireMetricNoBurn_beam_slope_mixed"),
                                  labels = c("Baseline model","Coverage Beam","Full Power Beam","Mixed Beam"))


df9 <- df %>% select(prefireMetricNoBurn_pvalslope,prefireMetricNoBurn_beam_pvalslope_coverage,prefireMetricNoBurn_beam_pvalslope_full,prefireMetricNoBurn_beam_pvalslope_mixed,metricFactor)
df9_stacked <- stack(df9,select = -metricFactor) %>% rename(pval = values,model = ind)
df9_stacked$response <- rep(df9$metricFactor,4)
df9_stacked$significant <- ifelse(df9_stacked$pval <= 0.01,"Sig.","Not Sig.")

df_stacked <- cbind(df8_stacked,df9_stacked$significant) %>% rename(significant = "df9_stacked$significant")
df_stacked$response <- fct_rev(df_stacked$response)

plotPrefireBeamSlope <- ggplot(df_stacked,aes(x=slope,y=response,shape = significant,color = modelFactor)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.2), size = 1.5) +
  theme_bw() +
  scale_color_manual(values=pal2) +
  scale_shape_manual(values = c(1, 16)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(labels = c("RH98","RH75","RH50","RH25",
                              expression(RE["<10m"]),expression(mPAI["0-10m"]),
                              "tPAI","tAGBD","Cover","aVDR","nmode","FHD")) +
  labs(x = expression("Normalized effect on " * Delta * " Structure"),
       title = "Pre-fire Structure") +
  guides(shape = "none") +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

#dev.new();print(plotPrefireBeamSlope)

# combine forest and toposlope plots into final figure
combined1 <- ggarrange(plotPrefireForestSlope,plotForestSlope,
                      labels = c("A", "B"),
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 1,
                      align = "h")

dev.new();print(combined1)

outFile <- paste0(indir,"/regression_models_dNBR_prefireMetric_forest_slopes_",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".png")
ggsave(outFile, width = 6, height = 4, units = "in", dpi=600)

# combine beam type and toposlope plots into final figure
combined2 <- ggarrange(plotPrefireBeamSlope, plotBeamSlope,
                      labels = c("A", "B"),
                      common.legend = TRUE,
                      legend = "bottom",
                      ncol = 2, nrow = 1,
                      align = "h")

dev.new();print(combined2)

outFile <- paste0(indir,"/regression_models_dNBR_prefireMetric_beam_slopes_",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".png")
ggsave(outFile, width = 6, height = 4, units = "in", dpi=600)

###############################################################################
### Plot of spatial regression plot slopes for dNBR * Prefire Metrics
###############################################################################

# version date
version <- "240822"
#version <- "241225"

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit
#maxAfter <- 365 # 1 year

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years
#maxBefore <- 365 # 1 year

# get regression results table
df <- read.csv(paste0(indir,"/master_gedi_fire_difference_spatialreg_interactions_mtbs_dnbrOW_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
df <- df %>% filter(response %in% paste0("delta_",metrics$metric))

df$metricFactor <- factor(df$response, 
                          levels = c("delta_strct_H","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                                     "delta_strct_mPAI_b10", "delta_strct_prop_int_below_10m",
                                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98")
                          )
pal <- RColorBrewer::brewer.pal(8, "Set2")
pal2 <- c(pal[8],pal[1:3])

# Plot of slopes for dNBR * Prefire Metrics interaction 
df10 <- df %>% select(prefireMetric_slope,prefireMetric_prefireMetric,prefireMetric_interaction, metricFactor)

df10_stacked <- stack(df10,select = -metricFactor) %>% rename(slope = values,model = ind)
df10_stacked$response <- rep(df10$metricFactor,3)
df10_stacked$modelFactor <- factor(df10_stacked$model, 
                                  levels = c("prefireMetric_slope","prefireMetric_prefireMetric","prefireMetric_interaction"),
                                  labels = c("dNBR","Pre-fire Structure","Interaction"))


df11 <- df %>% select(prefireMetric_pvalslope,prefireMetric_pvalprefireMetric,prefireMetric_pvalinteraction, metricFactor)
df11_stacked <- stack(df11,select = -metricFactor) %>% rename(pval = values,model = ind)
df11_stacked$response <- rep(df11$metricFactor,3)
df11_stacked$significant <- ifelse(df11_stacked$pval <= 0.01,"Sig.","Not Sig.")

df_stacked <- cbind(df10_stacked,df11_stacked$significant) %>% rename(significant = "df11_stacked$significant")

df_stacked$response <- fct_rev(df_stacked$response)

plotdNBRPreFireSlope <- ggplot(df_stacked,aes(x=slope,y=response,shape = significant,color = modelFactor)) +
  geom_point(size = 1.5) +
  theme_bw() +
  scale_color_manual(values=pal) +
  #scale_shape_manual(values = c(1, 16)) +
  scale_shape_manual(values = c(16)) +
  labs(x = "Normalized slope") +
  guides(shape = "none") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1),breaks = pretty_breaks(n = 5)) +
  scale_y_discrete(labels = c("RH98","RH75","RH50","RH25",
                              expression(RE["<10m"]),expression(mPAI["0-10m"]),
                              "tPAI","tAGBD","Cover","aVDR","nmode","FHD")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right"
  ) + 
  geom_vline(xintercept = c(0), linetype = "dashed", color = "black")

dev.new();print(plotdNBRPreFireSlope)

outFile <- paste0(indir,"/regression_models_dNBR_prefireMetric_interaction_zscore_slopes_",distance,"m_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".png")
ggsave(outFile, width = 5, height = 4, units = "in", dpi=600)

###############################################################################
### Plot of sample size for regression analysis
###############################################################################

# version date
version <- "240822"

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit
#maxAfter <- 365 # 1 year

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years
#maxBefore <- 365 # 1 year

# get regression results table

predictors <- c("mtbs_dnbrOW", "topoSlope","timeAfter","prefireMetric","vpd","windSpeed","et")
dfdNBR <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[1],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dftopoSlope <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[2],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dftimeAfter <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[3],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfpreMetric <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[4],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfvpd <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[5],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfwindSpeed <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[6],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))
dfet <- read.csv(paste0(indir,"/univariate_spatialreg_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_",predictors[7],"_zscore_maxBefore",maxBefore,"_minAfter",minAfter,"_maxAfter",maxAfter,"_",version,".csv"))

df <- rbind(dfdNBR,dftimeAfter,dftopoSlope,dfvpd,dfet,dfwindSpeed,dfpreMetric)
df <- df %>% filter(response %in% paste0("delta_",metrics$metric))

df$metricFactor <- factor(df$response, 
                          levels = c("delta_strct_H","delta_strct_nmode","delta_strct_VDR","delta_strct_cover","delta_strct_tAGBD","delta_strct_tPAI",
                                     "delta_strct_mPAI_b10","delta_strct_prop_int_below_10m",
                                     "delta_strct_RH_25","delta_strct_RH_50","delta_strct_RH_75","delta_strct_RH_98"))

df$predictorFactor <- factor(df$predictor, 
                             levels = c("mtbs_dnbrOW","timeAfter","topoSlope","prefireMetric","et","vpd","windSpeed"),
                             labels = c("dNBR","Elapsed Time","Topographic Slope","Pre-fire Structure","Evapotranspiration","Vapor Pressure Deficit","Wind Speed"))
df10 <- df %>% select(n,metricFactor)
df10_stacked <- stack(df10,select = -metricFactor) %>% rename(n = values)
df10_stacked$response <- df10$metricFactor
df10_stacked$predictor <- df$predictorFactor

plotSampleSize <- ggplot(df10_stacked,aes(x=response,y=n)) +
  facet_wrap(~predictor,ncol=2) +
  geom_bar(stat="identity",position=position_dodge(),fill = pal[3]) +
  theme_bw() +
  ylab("Count") +
  coord_cartesian(ylim = c(32000, 34000)) +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                              expression(mPAI["0-10m"]),expression(RE["<10m"]),
                              "RH25","RH50","RH75","RH98")) +
  theme(axis.text.x = element_text(size = 10,angle=90,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width = unit(.1, "in") 
  ) 

dev.new();print(plotSampleSize)
outFile <- paste0(indir,"/regression_models_sample_size_",distance,"m_",version,".png")
ggsave(outFile, width = 6, height = 5, units = "in", dpi=600)

###############################################################################
### Plot of sample size for ANOVA analysis
###############################################################################

# version date
version <- "240822"

df <- read.csv(paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore730_minAfter0_maxAfter365_",version,".csv"))
df11 <- df %>% select(n,metric) %>% mutate(time = "<1 year") %>% filter(metric %in% metrics$metric)
df <- read.csv(paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore730_minAfter365_maxAfter730_",version,".csv"))
df12 <- df %>% select(n,metric) %>% mutate(time = "1-2 years") %>% filter(metric %in% metrics$metric)
df <- read.csv(paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore730_minAfter730_maxAfter1500_",version,".csv"))
df13 <- df %>% select(n,metric) %>% mutate(time = "2-3 years") %>% filter(metric %in% metrics$metric)
df <- read.csv(paste0(indir,"/nlme_repeated_anova_spatial_zscore_abgd",biomass,"Mghayr_elev",delta_elev,"m_dist",distance,"m_maxBefore730_minAfter0_maxAfter1500_",version,".csv"))
df14 <- df %>% select(n,metric) %>% mutate(time = "All years") %>% filter(metric %in% metrics$metric)

df15 <- rbind(df11,df12,df13,df14)

df15$metricFactor <- factor(df11$metric, 
                            levels = c("strct_H","strct_nmode","strct_VDR","strct_cover","strct_tAGBD","strct_tPAI",
                                       "strct_mPAI_b10","strct_prop_int_below_10m",
                                       "strct_RH_25","strct_RH_50","strct_RH_75","strct_RH_98")
                            )

pal <- RColorBrewer::brewer.pal(5, "Set2")
plotSampleSizeAnova <- ggplot(df15,aes(x=metricFactor,y=n, fill = time)) +
  geom_bar(stat="identity",position=position_dodge()) +
  #facet_wrap(~metricFactor,ncol=2) +
  theme_bw() +
  ylab("Count") +
  scale_fill_manual(values=pal) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 5)) +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                              expression(mPAI["0-10m"]),expression(RE["<10m"]),
                              "RH25","RH50","RH75","RH98")) +
  theme(axis.text.x = element_text(size = 10,angle=45,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width = unit(.1, "in") 
  ) 

dev.new();print(plotSampleSizeAnova)
outFile <- paste0(indir,"/anova_sample_size_",distance,"m_",version,".png")
ggsave(outFile, width = 6, height = 3, units = "in", dpi=600)

###############################################################################
### Plot footprints per fire
###############################################################################

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

file <- paste0(indir,"/master_gedi_fire_difference_",distance,"m_240822.csv")
df <- read.csv(file)

# calculate biomass growth limit based on time difference
df$abgd_limit <- (df$timeDiff/365) * biomass

data <- df %>% filter(timeBefore <= maxBefore & timeAfter >= minAfter & timeAfter <= maxAfter &
                        !is.na(mtbs_dnbrOW) &  
                        abs(delta_GEDI_elev) <= delta_elev &
                        delta_strct_tAGBD <= abgd_limit) %>%
                select(fire, region)

fireCount <- data %>% group_by(fire) %>% summarize(count = n())

file <- paste0(indir,"/fire_data_renamed.csv")
fireLookup <- read.csv(file) %>% rename(fire = fireName)

fireCount <- left_join(fireCount,fireLookup, by="fire")

pal <- RColorBrewer::brewer.pal(5, "Set2")
plotSampleSizeFires <- ggplot(fireCount,aes(x=name,y=count, fill = region)) +
  geom_bar(stat="identity",position=position_dodge()) +
  theme_bw() +
  ylab("Count") +
  scale_fill_manual(values=pal) +
  scale_y_continuous(labels = scales::number_format(accuracy = 1),breaks = pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(size = 9,angle=90,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.width = unit(.1, "in") 
  ) 

dev.new();print(plotSampleSizeFires)
outFile <- paste0(indir,"/fire_sample_size_",distance,"m_",version,".png")
ggsave(outFile, width = 6, height = 3, units = "in", dpi=600)

outFile <- paste0(indir,"/fire_sample_size_",distance,"m_",version,".csv")
write.csv(fireCount, outFile, row.names = F)

###################################################################################
## Plot of Change in tAGBD vs Change in RE10m metrics, colored by dNBR
###################################################################################

# version date
version <- "240822"

library(viridis)

# read change metrics file
file <- paste0(indir,"/master_gedi_fire_difference_",distance,"m_",version,".csv")
df <- read.csv(file)

# get wildfire burn severity classes
mtbsCodes <- data.frame(
  mtbs_classP = c(0, 1, 2, 3, 4),
  FireSeverityMTBS = c("Control","Unchanged", "Low","Moderate","High")
)
df <- left_join(df,mtbsCodes,by="mtbs_classP")
df$FireSeveritydNBRMillerThode <- ifelse(df$mtbs_dnbrOW >= 367, "High",ifelse(df$mtbs_dnbrOW >= 177, "Moderate",ifelse(df$mtbs_dnbrOW > 41, "Low", "Unchanged"))) 
df$FireSeveritydNBRMillerThode <- ifelse(df$FireSeverityMTBS == "Control", "Control",df$FireSeveritydNBRMillerThode)

# maximum time after fire limit (days)
maxAfter <- 1500 # no limit

# minimum time after fire limit (days)
minAfter <- 0 # 0 years

# maximum time before fire limit (days)
maxBefore <- 730 # 2 years

# elevation change limit
delta_elev <- 10

# elevation change limit
abgd_limit <- 99999

# get data
data <- df %>% filter(timeBefore <= maxBefore & timeAfter >= minAfter & timeAfter <= maxAfter &
                        !is.na(mtbs_dnbrOW) &  
                        abs(delta_GEDI_elev) <= delta_elev &
                        delta_strct_tAGBD <= abgd_limit)
data$FireSeverityOrdered <- factor(data$FireSeveritydNBRMillerThode, levels = c("Control","Unchanged", "Low","Moderate","High"))
data$typeForestOrdered <- factor(data$typeForest, levels = c("conifer","hardwood","mixed"), labels = c("Conifer","Hardwood","Mixed"))

# GEDI metrics to use
metricList <- read.csv(paste0(indir,"/gedi_structure_variables_reduced_pai_240927.csv"))  %>% filter(selected == "Y")
metrics <- left_join(metricList, metricLookup, by="metric") %>% select(metric, labels)

# select only data from fires
d16 <- data %>% filter(FireSeverityOrdered != "Control")

vars <- c(paste0("pre_",metrics$metric),paste0("delta_",metrics$metric))
d17 <- d16 %>% select(all_of(c("mtbs_dnbrOW",vars))) %>%
  mutate(across(c("mtbs_dnbrOW", vars), ~ zscore(.), .names = "{.col}_zscore")) %>%
  na.omit() 

### Plot Intensity<10m vs. Biomass
x <- d17$pre_strct_tAGBD
y <- d17$pre_strct_prop_int_below_10m
c <- d17$delta_strct_prop_int_below_10m

d18 <- data.frame(x=x,y=y,c=c) %>%
  filter(is.finite(x), is.finite(y), is.finite(c))

# Fit a linear model
model <- lm(y ~ x, data = d18)

# Get the R-squared value
r_squared <- summary(model)$r.squared

plot9 <- ggplot(d18,aes(x=x,y=y)) +
  stat_summary_hex(aes(z = c), fun = mean, bins = 50) + 
  theme_bw() +
  xlab("") +
  ylab("") +
  coord_cartesian(xlim = c(0, 800), ylim = c(0,100)) +
  labs(
    x = "Pre-fire tAGBD",
    y =  expression("Pre-fire RE"["<10m"]),
    fill = c) +
  scale_fill_distiller(palette="RdYlGn", direction= 1, guide = guide_colorbar(
    title = expression("Decrease    " * Delta * " RE"["<10m"] * "    Increase"),
    title.position = "top", 
    title.hjust = 0.5,
    label = TRUE,
    barwidth = 10,
    barheight = 1 
  )) +
  geom_smooth(method = "lm",formula = "y~x", color = "black", size = 0.5) + 
  annotate("text", x = 500, y = 90,
           label = paste("r? =", round(r_squared, 2)),
           hjust = 0, vjust = 0, size = 4, color = "black") 


x <- d17$delta_strct_tAGBD
y <- d17$delta_strct_prop_int_below_10m
c <- d17$mtbs_dnbrOW

d18 <- data.frame(x=x,y=y,c=c)

# Fit a linear model
model <- lm(y ~ x, data = d18)

# Get the R-squared value
r_squared <- summary(model)$r.squared

plot10 <- ggplot(d18,aes(x=x,y=y)) +
  stat_summary_hex(aes(z = c), fun = mean, bins = 50) + 
  theme_bw() +
  xlab("") +
  ylab("") +
  coord_cartesian(xlim = c(-600,600), ylim = c(-100,100)) +
  labs(
    x = expression(Delta * " tAGBD"),
    y = expression(Delta * " RE"["<10m"]),
    fill = x) +
  scale_fill_distiller(palette = "RdYlGn", guide = guide_colorbar(
    title = "Low        Avg dNBR        High",  
    title.position = "top",  
    label = TRUE, 
    title.hjust = 0.5,
    barwidth = 10,  
    barheight = 1  
  )) +
  geom_smooth(method = "lm",formula = "y~x", color = "black", size = 0.5) + 
  annotate("text", x = 150, y = 80,
           label = paste("r? =", round(r_squared, 2)),
           hjust = 0, vjust = 0, size = 4, color = "black") 


# Arrange plots with a different common legend
arrange <- ggarrange(plot9, plot10,
                       labels = c("A", "B"),
                       common.legend = F,
                       legend = "bottom",
                       ncol = 2, nrow = 1,
                       align = "h")

dev.new();print(arrange)

outFile <- paste0(indir,"/intensity_graphs_",distance,"m_",version,".png")
ggsave(outFile, width = 7, height = 6, units = "in", dpi=600)

###################################################################################
## Histogram of distances between pairs
###################################################################################

meanDistance <- mean(data$distance)

histogramDistance <- ggplot(data, aes(x = distance)) +
  geom_histogram(bins = 20, color = "black", fill = "gray") +
  theme_bw() +
  xlab("Distance (m)") +
  ylab("Count") +
  geom_vline(aes(xintercept = meanDistance), color = "black", linetype = "dashed", size = 1) 

dev.new();print(histogramDistance)
outFile <- paste0(indir,"/pair_distance_histogram_",distance,"m_",version,".png")
ggsave(outFile, width = 6, height = 6, units = "in", dpi=600)

############ 3D Plot ################
library(plotly)

x <- d16$pre_strct_prop_int_below_10m
y <- d16$delta_strct_prop_int_below_10m
z <- d16$mtbs_dnbrOW
c <- d16$delta_strct_cover *100
d18 <- data.frame(x=x,y=y,c=c)

plot11 <- plot_ly(
  data = d18,
  x = ~x,
  y = ~y,
  z = ~z,
  color = ~c,
  colors = colorRamp(c("blue", "red")),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Pre-Fire Intensity"),
      yaxis = list(title = "Change in Intensity"),
      zaxis = list(title = "dNBR")
    ),
    title = "3D Plot of Pre-Fire Intensity, Change in Intensity, and dNBR"
  )

plot11

###################################################################################
## Correlation matrices
###################################################################################

library(reshape2)  # for the melt function

# Function to rescale using z score and filter outliers with abs(z) > 4
zscore<-function(df,z=4){
  tmp <- as.matrix(apply(as.matrix(df),2,scale))
  if (ncol(tmp)==1) {tmp[which(abs(tmp)>z)] <- NA}
  if (ncol(tmp)>1) {tmp[which(abs(tmp)>z,arr.ind=T)] <- NA}
  return(tmp)
}

##### correlation of change metrics #####
vars <- paste0("delta_",metrics$metric)

d19 <- d16 %>% select(all_of(c("distance","timeAfter","mtbs_dnbrOW",vars))) %>% mutate(across(c("distance","timeAfter","mtbs_dnbrOW",vars), zscore)) %>% na.omit() 

cn <- gsub("delta_", "", vars) 
cn <- recode(cn, !!!setNames(metricLookup$labels, metricLookup$metric))
colnames(d19) <- c("Pair distance","Time Since Fire","dNBR",cn)

d19 <- d19 %>% select(c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                          "mPAI0to10","REBelow10m",
                          "RH25","RH50","RH75","RH98","dNBR","Pair distance","Time Since Fire"))
                            
# Calculate the correlation matrix
corr_matrix <- cor(d19)

# Convert the matrix to a tidy format
melted_corr_matrix <- melt(corr_matrix)

# Plot the correlation matrix using ggplot2
plot12 <- ggplot(data = melted_corr_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size = 10)) +
  coord_fixed() +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                              expression(mPAI["0-10m"]),expression(RE["<10m"]),
                              "RH25","RH50","RH75","RH98","dNBR","Pair distance","Time Since Fire")) +
  scale_y_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                              expression(mPAI["0-10m"]),expression(RE["<10m"]),
                              "RH25","RH50","RH75","RH98","dNBR","Pair distance","Time Since Fire")) +
  labs(title = "Correlation of Change in Structural Metrics",
       x = "",
       y = "") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 2)
  
dev.new();print(plot12)
outFile <- paste0(indir,"/correlation_graph_change_",distance,"m_",version,".png")
ggsave(outFile, width = 6, height = 6, units = "in", dpi=600,bg = "transparent")

##### correlation pre-fire & post-fire metrics #####

vars <- paste0("pre_",metrics$metric)
d23 <- d16 %>% select(all_of(c(vars))) %>% mutate(across(c(vars), zscore)) %>% na.omit() 
vars <- paste0("post_",metrics$metric)
d24 <- d16 %>% select(all_of(c(vars))) %>% mutate(across(c(vars), zscore)) %>% na.omit() 

cn <- gsub("post_", "", vars) 
cn <- recode(cn, !!!setNames(metricLookup$labels, metricLookup$metric))
colnames(d23) <- cn
colnames(d24) <- cn

d25 <- rbind(d23,d24)

d25 <- d25 %>% select(c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                        "mPAI0to10","REBelow10m",
                        "RH25","RH50","RH75","RH98"))

labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",expression(mPAI["<10m"]),
           expression(RE["<10m"]),"RH25","RH50","RH75","RH98")

# Calculate the correlation matrix
corr_matrix <- cor(d25)

# Convert the matrix to a tidy format
melted_corr_matrix <- melt(corr_matrix)

# Plot the correlation matrix using ggplot2
plot13 <- ggplot(data = melted_corr_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size = 10)) +
  coord_fixed() +
  scale_x_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                              expression(mPAI["0-10m"]),expression(RE["<10m"]),
                              "RH25","RH50","RH75","RH98")) +
  scale_y_discrete(labels = c("FHD","nmode","aVDR","Cover","tAGBD","tPAI",
                              expression(mPAI["0-10m"]),expression(RE["<10m"]),
                              "RH25","RH50","RH75","RH98")) +
  labs(title = "Correlation of all Structural Metrics (Pre- and Post-fire)",
       x = "",
       y = "") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 2)

dev.new();print(plot13)
outFile <- paste0(indir,"/correlation_graph_pre-post_",distance,"m_",version,".png")
ggsave(outFile, width = 6, height = 6, units = "in", dpi=600,bg = "transparent")



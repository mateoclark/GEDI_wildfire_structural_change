# Pair GEDI footprints and calculate post- minus pre-fire metric differences

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 08/18/2024

library(dplyr)
library(sf)
library(lubridate)

# input directory
inDir <- "E:/active/project/calfire_gedi/rq2_paper"

# search distance to find matched points
searchDistance <- 45 # meters

# months for analysis
analysisMonths <- c(5:10) # May through October

# minimum number of shots to include fire in analysis
numShots <- 50

# input master .csv file with GEDI and fire information
inputFile <- paste0(inDir,"/data/master_gedi_fire_data_240809.csv")

# output .csv file with GEDI metric differences by fire
outputFile <- paste0(inDir,"/data/master_gedi_fire_difference_45m_240818.csv")

# input .csv file with fire information
fireFile <- paste0(inDir,"/data/fire_data.csv")

# fire weather auxiliary data
auxFile <- paste0(inDir,"/data/fire_csvs/2024_3_28/4_t_w.csv")

# load input data
df <- read.csv(inputFile)
df <- df %>% filter(Month %in% analysisMonths)
df <- tibble::rowid_to_column(df, "index")

# create a mean PAI below 10m metric
df <- df %>% rowwise() %>% mutate(strct_mPAI_b10 = mean(c(strct_PAI_inc_00_05,strct_PAI_inc_05_10)))

# get GEDI structure metrics
metrics <- colnames(df)[grep(colnames(df),pattern ="strct")]
metrics <- c("GEDI_elev",metrics)

# read fire information
fireInfo <- read.csv(fireFile)

# Convert date
fireInfo$alarmDate <- as.Date(fireInfo$alarmDate, format = "%m/%d/%Y")
fireInfo$containmentDate <- as.Date(fireInfo$containmentDate, format = "%m/%d/%Y")

# select pre- and post-fire footprints by fire
fireNames <- unique(df$fire)
finaldt <- data.frame()
for (fireName in fireNames){
  
  print(fireName)
  
  # select data for the fire
  dfSelect <- df %>% filter(fire == fireName)
                
  # set time before or after fire
  f <- fireInfo[fireInfo$fireName == fireName,]
  
  if (nrow(f) > 0 ){
    dfSelect$timeBefore <- as.numeric(abs(difftime(dfSelect$Date,f$alarmDate,units = "days")))
    dfSelect$timeAfter <- as.numeric(abs(difftime(dfSelect$Date,f$containmentDate,units = "days")))
    
    # select out pre- and post-fire footprints and create spatial layers
    prefire <- dfSelect %>% filter(fire == fireName & GEDI_fireTiming == "pre")
    prefire$timeAfter <- NULL
    postfire <- dfSelect %>% filter(fire == fireName & GEDI_fireTiming == "post")
    postfire$timeBefore <- NULL
    
    zone <- prefire$GEDI_utmZ[1]
    if (zone == 10) crs = 32610 else crs = 32611
    
    prefireLayer <- st_as_sf(prefire, coords = c("GEDI_utmX", "GEDI_utmY"), crs = crs, remove = FALSE)
    postfireLayer <- st_as_sf(postfire, coords = c("GEDI_utmX", "GEDI_utmY"), crs = crs)
    
    # create a buffer search area around post-fire footprints
    postfireLayerBuffered <- st_buffer(postfireLayer, searchDistance)
    
    # find pre-fire footprints that overlap post-fire buffers and also have same CALVEG class
    prefireLayerPostBuffer <- st_join(prefireLayer,postfireLayerBuffered) %>% filter(typeForest.x == typeForest.y)
    
    if (dim( prefireLayerPostBuffer)[1] > 0){
      # summarize mean values of GEDI metrics for all pre-fire data within a post-fire buffer
      metrics.x <- paste(metrics,"x",sep=".")
      data <- data.frame(prefireLayerPostBuffer %>% select(all_of(c("index.y","GEDI_shot.x","GEDI_utmX","GEDI_utmY","mtbs_dnbrOW.x","GEDI_beam.x","timeBefore","Date.x",metrics.x)))) %>% select (-geometry)
      data$timeDiff <- as.numeric(abs(difftime(prefireLayerPostBuffer$Date.x,prefireLayerPostBuffer$Date.y,units = "days")))
      
      # when multiple pre-fire footprints exist, select the minimum time difference pair; if more than 1, chose one at random
      dataMinTimeDiff <- data %>% group_by(index.y) %>% filter(timeDiff == min(timeDiff)) %>% slice_sample(n = 1)
      dataSummary <- dataMinTimeDiff %>% group_by(index.y) %>% slice_sample(n = 1)
      
      premetrics <- paste("pre",metrics,sep="_")
      colnames(dataSummary) <- c("index","pre_GEDI_shot","pre_GEDI_utmX", "pre_GEDI_utmY","pre_mtbs_dnbrOW","pre_Beam","timeBefore","pre_Date",premetrics,"timeDiff")
      
      # join post- and pre-fire data
      dt <- left_join(postfire,dataSummary, by="index") %>% filter(timeBefore != "NA")
      
      # calculate difference between pre- and post-fire metrics
      post <- dt %>% select(all_of(metrics))
      pre <- dt %>% select(all_of(premetrics))
      
      # calculate change
      delta <- post - pre
      colnames(delta) <- paste("delta",metrics,sep="_")

      # create master data frame with all variables
      dt_master <- cbind(dt,delta)
      
      # replace column names for post-fire metrics
      colNames <- colnames(dt_master)
      postmetrics <- paste("post",metrics,sep="_")
      colnames(dt_master) <- replace(colNames, which(colNames %in% metrics), postmetrics)
      dt_master <- dt_master %>% rename("post_Beam" = "GEDI_beam")
      
      # bind to final data table
      finaldt <- rbind(finaldt,dt_master)
      
      # remove data
      remove(prefire)
      remove(postfire)
      remove(prefireLayer)
      remove(postfireLayer)
      remove(postfireLayerBuffered)
      remove(prefireLayerPostBuffer)
      remove(dt)
      remove(pre)
      remove(post)
      remove(delta)
      gc()
    }
  }
}

# filter for selected fires, remove wetlands and remove pairs with pre-fire RH98 height below 3 m (eg shrublands, grasslands)
dt <- finaldt %>% filter(typeForest != "wetland" & pre_strct_RH_98 > 3)

# find fires that have total number of shots above numShots
selectedFires <- dt %>% count(fire) %>% filter(n >= numShots)
dt <- dt %>% filter(fire %in% selectedFires$fire)

# function to select beam type (if both full and coverage, then labeled Mixed)
beam <- function(x){
  beamTypes = unique(x)
  if (length(beamTypes) > 1) beam = "mixed"
  else beam = beamTypes
  return(beam)
}

# add a field called beamType that indicates if pre- and post-fire beams are both coverage, full or mixed
dt <- dt %>% rowwise() %>% mutate(beamType = beam(c(pre_Beam,post_Beam)))

# add distance
postxy <- data.frame(x=dt$GEDI_utmX,y=dt$GEDI_utmY)
prexy <- data.frame(x=dt$pre_GEDI_utmX,y=dt$pre_GEDI_utmY)
euclidean_distance <- function(p1, p2) {
  sqrt(sum((p1 - p2)^2))
}
dt$distance <- apply(cbind(postxy, prexy), 1, function(row) {
  euclidean_distance(row[1:2], row[3:4])
})

# add post-fire footprint auxiliary predictors based on GEDIshot
aux <- read.csv(auxFile)
auxSelect <- aux %>% select(c("GEDI_shot","topo_slope", "w_vpd_d_max","w_ET_b4_mean","w_vs_d_max")) 
auxMean<- auxSelect %>% group_by(GEDI_shot) %>% summarize(
  topoSlope = mean(topo_slope, na.rm = TRUE),
  vpd = mean(w_vpd_d_max, na.rm = TRUE),
  windSpeed = mean(w_vs_d_max, na.rm = TRUE),
  et = mean(w_ET_b4_mean, na.rm = TRUE)
)
dt <- left_join(dt, auxMean, by = "GEDI_shot")

# Function to create a CRS string for UTM based on the zone
crsUTM <- function(zone) {
  paste0("+proj=utm +zone=", zone, " +datum=WGS84 +units=m +no_defs")
}

# Iterate through each UTM zone and transform to California Albers
# Need a uniform projection for distance matrix
dtSF <- do.call(rbind, lapply(split(dt, dt$GEDI_utmZ), function(x) {
  st_as_sf(x, coords = c("GEDI_utmX","GEDI_utmY"), crs = crsUTM(x$GEDI_utmZ[1]), remove = F) %>%
    st_transform(crs = 3310)  # California Albers (Teale Albers) NAD83 EPSG:3310
}))

# Add converted coordinates and convert back to data frame
dt1 <- as.data.frame(dtSF) %>% 
  mutate(post_x_teale = st_coordinates(geometry)[,1],
         post_y_teale = st_coordinates(geometry)[,2]) %>%
  select(-geometry)

# Iterate through each UTM zone and transform to California Albers
# Need a uniform projection for distance matrix
dtSF <- do.call(rbind, lapply(split(dt1, dt1$GEDI_utmZ), function(x) {
  st_as_sf(x, coords = c("pre_GEDI_utmX","pre_GEDI_utmY"), crs = crsUTM(x$GEDI_utmZ[1]), remove = F) %>%
    st_transform(crs = 3310)  # California Albers (Teale Albers) NAD83 EPSG:3310
}))

# Add converted coordinates and convert back to data frame
dtFinal <- as.data.frame(dtSF) %>% 
  mutate(pre_x_teale = st_coordinates(geometry)[,1],
         pre_y_teale = st_coordinates(geometry)[,2]) %>%
  select(-geometry)

# write out .csv file of difference statistics
write.csv(dtFinal, outputFile, row.names = FALSE)

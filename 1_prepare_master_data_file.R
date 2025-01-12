# loads in GEDI structure x fire data .csv and outputs a master .csv file with information needed for analysis.
# removes shots that are in treatments or in overlapping fires

# Matthew L. Clark, Ph.D.
# Department of Geography, Environment, and Planning
# Sonoma State University, California USA
# matthew.clark@sonoma.edu
# Version 08/09/2024

library(tidyverse)
library(lubridate)
library(sf)

# working directory
indir <- "/GEDI_wildfire_structural_change"

# input file name with GEDI structure (zipped in GitHub)
fileName <- paste0(indir,"/2_strct_sev_2024_8_9.csv")

excludeFires <- c("BLUEJAY")

# file geodatabase with vector layers (zipped in GitHub)
fgdb <- paste0(indir,"/calfire_gedi_treatments.gdb")

# load fire GEDI data
df1 <- read.csv(fileName)

# remove shots that were not within the fire perimeter
df2 <- filter(df1, typeForest != "NA" & !fire %in% excludeFires)

# get month and year, accounting for any discrepancies in date format in csvs
df2$Date <-as.Date(df2$GEDI_date)
df2$Month <- month(as.POSIXlt(df2$Date, format="%Y/%m/%d",tz="UTC"))
df2$Year <- year(as.POSIXlt(df2$Date, format="%Y/%m/%d",tz="UTC"))
df3 <- df2 %>% select(-GEDI_date) %>% filter(!is.na(Date))

# find shots that are in treatments

## UTM zone 10
treatments <- st_read(fgdb,layer = "facts_calmapper_utmz10")
treatmentsPoly <- st_cast(treatments, "POLYGON")
treatmentsPoly$Treatment <- "Y"
dfUTMz10 <- df3 %>% filter(GEDI_utmZ == 10)
shotsUTMz10 <- st_as_sf(dfUTMz10, coords = c("GEDI_utmX", "GEDI_utmY"),remove=F)
shotsUTMz10 <- st_set_crs(shotsUTMz10, 26910)
shotsUTMz10Treatments <- st_join(shotsUTMz10, treatmentsPoly,left = TRUE)
shotsUTMz10Treatments <- data.frame(shotsUTMz10Treatments) %>% select(-c(Shape_Length,Shape_Area))

## UTM zone 11
treatments <- st_read(fgdb,layer = "facts_calmapper_utmz11")
treatmentsPoly <- st_cast(treatments, "POLYGON")
dfUTMz11 <- df3 %>% filter(GEDI_utmZ == 11)
treatmentsPoly$Treatment <- "Y"
shotsUTMz11 <- st_as_sf(dfUTMz11, coords = c("GEDI_utmX", "GEDI_utmY"),remove=F)
shotsUTMz11 <- st_set_crs(shotsUTMz11, 26911)
shotsUTMz11Treatments <- st_join(shotsUTMz11, treatmentsPoly,left = TRUE)
shotsUTMz11Treatments <- data.frame(shotsUTMz11Treatments) %>% select(-c(Shape_Length,Shape_Area))

df4 <- rbind(shotsUTMz10Treatments,shotsUTMz11Treatments)

# remove duplicate shots due to overlapping fires or overlapping treatments
df5 <- distinct(df4, GEDI_shot, .keep_all= TRUE)

# set null start and end dates to pre-GEDI or post-GEDI dates, respectively
df5$StartDate <- ifelse(is.na(df5$StartDate),"2018-01-01",df5$StartDate) # before GEDI
df5$EndDate <- ifelse(is.na(df5$EndDate),"2030-01-01",df5$EndDate) # after current GEDI 
df5$Treatment <- ifelse(is.na(df5$Treatment),"N",df5$Treatment) # set shots outside of treatment polygons to "N"

# remove shots that are in treatments based on timing of shot and pre- vs post-fire setting
df6 <- df5 %>% filter(
  (Treatment == "Y" & GEDI_fireTiming == "pre" & Date > EndDate)| # in treatment polygon, keep pre-fire shots after treatment ends
  (Treatment == "Y" & GEDI_fireTiming == "post" & Date < StartDate)| # in treatment polygon, keep post-fire shots before treatment starts
   Treatment == "N" # keep shots outside of treatment polygons
) %>% select(-geometry)

df7 <- df6 %>% select(!c("StartDate","EndDate"))

# write out GEDI fire master .csv file
outfile <- paste0(indir,"/master_gedi_fire_data_240809.csv")
write.csv(df7, outfile, row.names = FALSE)


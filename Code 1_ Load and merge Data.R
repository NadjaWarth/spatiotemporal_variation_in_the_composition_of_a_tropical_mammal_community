#CODE 1 

### Optional: install devtools and camtraptor ----
# Analysis Camera trapping data 
# Instructions on the ctdp package from: https://wec.wur.nl/r/ctdp/
# Camtraptor: https://github.com/inbo/camtraptor/blob/main/R/read_wi.R

#install.packages('devtools')
#library(devtools)
#install.packages("sass", type = "source", repos = "https://cran.r-project.org/")
#install.packages("openssl", repos = "https://cran.r-project.org/")
#install.packages("writexl")  
#devtools::install_github("frictionlessdata/frictionless-r")
#devtools::install_github("inbo/camtraptor")
#devtools::install_gitlab(repo = "camtrap/ctdp", host = "git.wur.nl")

### Load Packages ----
library(vegan)
library(tidyr)
library(dplyr)    # For data manipulation
library(gridExtra)
library(camtraptor)
library(ctdp)
library(readxl)

### Variables and paths used for this script ----
# Sources
agouti_data_path <- "C:/NadjaMA/data/original/26.10_bcnm-mammal-monitoring-team-20231026142012.zip"
wildID_data_path <- "C:/NadjaMA/data/original/WildID_data_08.01.24/"
deployment_history_path <-"C:/NadjaMA/data/original/Camera_deployment_history.xlsx"

# Outputs
agouti_data_output_path <- "C:/NadjaMA/data/processed/merged_Agouti_31.10.rds"
wildID_data_output_path <- "C:/NadjaMA/data/processed/merged_WildID_01.08.rds"
taxonomy_table_output_path <- "C:/NadjaMA/data/processed/Species_list.csv"
camera_data_output <- "C:/NadjaMA/data/processed/data_cam_dd_23.11.rds"

# 1) Load and prepare Agouti data  ----------------------------------------------------------------------
  
# Load the Agouti Data and remove empty deployments: 
Agouti<- read_ctdp(path=agouti_data_path ,tz = "Etc/GMT+5")
Agouti$deployments <- Agouti$deployments[!is.na(Agouti$deployments$locationID), ]
    
# Export species taxonomy table (for WildID data)
taxonomy <- Agouti$taxonomy %>%
  select(scientificName ,vernacularNames.en,class, order)
write_excel_csv(taxonomy, taxonomy_table_output_path)
  
# Merge relevant data from Agouti (only use data from 2018 and later, because WildID is used for the years before 2018) 
# ASSUMPTION: Deployment interval = (time interval from setup) - (pick up or last image taken by a camera before it stopped functioning)
# 
# This pipe merges the subtibbles of the ctdp file and only keeps the columns needed for the analysis.
# From the deployment interval it extracts the Start and End Date of a deployment.
# From the sequence interval it extracts the Date on which the interval was recorded (needed later).
# Then it filters the data to only contain deployments starting from 01.01.2018.
# 
# returns a dataframe with the following columns:
# "locationName","deploymentID","StartDate","EndDate","count","scientificName","vernacularNames.en","class","longitude","latitude","SequDate"

merged_Agouti <- 
  #Create a big data table with all the columns by merging the "subtibbles"
  merge_tibbles(Agouti, dropMedia = TRUE) %>% 
  
  #Select relevant Columns
  select(locationName, deploymentID, deployment_interval,
         sequence_interval, count, scientificName, vernacularNames.en, class,longitude, latitude) %>%
  
  #Make a new column for the Setup date (key to join with deployment history data) and End date of each deployment 
  separate(deployment_interval, into = c("StartDate", "EndDate"), sep = "--", convert = TRUE) %>%  
  mutate(StartDate = as.Date(StartDate), EndDate = as.Date(EndDate)) %>% 
  
  #Regenerate deployment Interval (if needed): 
  #mutate(deployment_interval = paste(StartDate, EndDate, sep = "--")) %>%
  
  #Create a column with a Sequence data instead of interval (needed later to filter out sequences)
  mutate(SequDate = as.POSIXct(str_split(sequence_interval, " -- ", simplify = TRUE)[, 1])) %>% 
  select(-c("sequence_interval")) %>% 
  
  #Only use data from 2018-now
  filter(EndDate>"2018-01-01") %>% 
  
  #Drop all rows that contain NA, could e.g. be observations that were categorized as "unkown"
  drop_na()

# Save output
saveRDS(merged_Agouti,agouti_data_output_path)

# 2) Load and prepare WildID data ---------------
# 2.1) Wild ID export might contain >1 image file. To use the read.wi function the image files in the WildID data folder first need to be merged.
image_files <- list.files(wildID_data_path, pattern = "^images.*\\.csv$", full.names = TRUE)
if (length(image_files) > 1) {
  # Loop through the files, read each file, and store it in the list
  image_data_list <- list()
  for (file in image_files) {
    print(file)
    image_data_list[[length(image_data_list) + 1]] <- read_csv(file)
  }
  # Combine all into one dataframe
  images_combined <- bind_rows(image_data_list)
  
  # Write the combined data frame to a new CSV file in the same folder
  write_excel_csv(images_combined, paste0(wildID_data_path, "images.csv"))
}
    
# 2.2) Load WildID data and taxonomy df
data_WildID <- read_wi(wildID_data_path)
taxonomy <- read.csv(taxonomy_table_output_path)
    
#2.3) Combine relevant data from WILD, prepare it to match the format of merged_Agouti above
merged_WildID <- 
  #### Join deployment and observations data in one table 
  data_WildID$data$deployments %>% 
  left_join(data_WildID$data$observations, by = "deploymentID") %>% 
  
  #### General Adjustments to fit the format of Agouti table 
  # Select relevant columns
  select(locationName, deploymentID, start, end,
         timestamp, count, scientificName, longitude, latitude) %>% 
  
  # Rename columnnames in "merged Agouti" and make sure the format is correct 
  mutate(StartDate = as.Date(start), 
         EndDate = as.Date(end), 
         SequDate = as.Date(timestamp),
         locationName = sub("CT-", "", locationName)) %>%
  
  #Add columns from Agouti Ctdp taxnomoy table that are missing in the WildID output
  left_join(taxonomy%>% select("scientificName", "vernacularNames.en", "class"), by ="scientificName") %>%
  
  #Only use data from 2011-2017
  filter(EndDate<"2018-01-01") %>% 

  
  #### Create image sequences: 
  #a) Convert timestamp from character to date format: 
  mutate(timestamp = as.POSIXct(timestamp, format = "%Y/%m/%d %H:%M:%S")) %>% 

  #b) Devide into sequences: 
  group_by(deploymentID) %>% 
  # Sort by time 
  arrange(timestamp) %>%
  # estimate difference between timestamps
  mutate(time_diff = timestamp -lag(timestamp)) %>% 
  # First value for time_diff is NA.Rows with timestamps <120s appart or NA get the value 1, 
  # others the value 0. Each 1 marks the start of a new sequence 
  mutate(Sequ_N = ifelse(time_diff > 120 | is.na(time_diff) , 1, 0)) %>% 
  # Use cumsum function to create increasing numbers: 
  mutate(Sequ_N = cumsum (Sequ_N)) %>% 
  # Create Sequence ID through combination of Deployment ID and Sequ N 
  mutate(SequenceID = paste(deploymentID,"_",Sequ_N)) %>% 
  ungroup() %>% 

  #c) Only keep rows with max count per species per sequence
  group_by(SequenceID, scientificName) %>% 
  # If there are multiple rows, only keep row with max count value (normally all counts should be the same within one sequence but just to be sure)
  slice(if (n() > 1) which.max(count) else 1) %>% 
  ungroup() %>% 

  #### Drop unneccesarry parts: 
  #Drop old columns
  select(-c("start", "end", "time_diff", "Sequ_N", "SequenceID")) %>% 
  #Drop all rows that contain NA, could e.g. be observations that were categorized as "unkown"
  drop_na()

#2.4.) Save Output: 
saveRDS(merged_WildID, wildID_data_output_path)
    
# 3) Combine Agouti and Wild ID data ------------------
#Load needed data 
merged_WildID<- readRDS(wildID_data_output_path)
merged_Agouti <- readRDS(agouti_data_output_path)

#Combine data from WILD ID and Agouti 
Camera_data <- bind_rows(merged_Agouti, merged_WildID) %>%
  # Filter locationName for BCI-1-..locations
  filter(grepl("BCI-1-", locationName)) %>%
  #Remove timestamp variable: 
  select(-c("timestamp"))
    
# 4) Load and prepare Detection Distance --------------------------
    
#Load deployment history data that includes detection distance 
depl_hist_raw <- read_excel(deployment_history_path)

#Adjustments in basic deployment sheet: 

depl_hist<- depl_hist_raw %>%
  #a) rename column locationID to  LocationName to match with the camera trap data
  rename(locationName = LocationID) %>%  
  #b) fix individual deviations from standard LocationID format "BCI-X-XX"
  mutate(locationName = ifelse(locationName == "BCI-3-11.2", "BCI-3-11", locationName)) %>% 
  #c) select relevant columns only (disregard comments etc.)
  select(locationName, SetupDate, Detection_Distance) %>% 
  #d) Only keep rows where Detection_Distance is not NA & filter out rows where Detection_Distance contains non-numeric characters.
  filter(!is.na(Detection_Distance) & !str_detect(Detection_Distance, "[^0-9.]")) %>%   
  #e) Converting detection distance to numeric class and SetupDate to date class, and rename to "StartDate" to match Ctdp data 
  mutate(Detection_Distance = as.numeric(Detection_Distance), 
         StartDate = as.Date(SetupDate)) %>% 
  select(-c("SetupDate")) %>% 
  #f) one value was at 1050 --> probably a missed comma
  mutate(Detection_Distance = ifelse(Detection_Distance == 1050, 10.50, Detection_Distance))  

#Estimate average detection distance per location: 
#Needed to fill gaps of missing detection distance in certain years 

Mean_dd_data <- depl_hist %>% 
  group_by(locationName) %>% 
  summarise(Mean_dd = mean(Detection_Distance))

# 5) Include Detection Distance in Camera data sheet #####

data_cam_dd <- 
    #Add detection distance from deployment history data 
    left_join(Camera_data, depl_hist, by = c("locationName", "StartDate")) %>% 
    left_join(Mean_dd_data, by = "locationName") %>%
    mutate(Detection_Distance = ifelse(is.na(Detection_Distance), Mean_dd, Detection_Distance)) %>%
    select(-Mean_dd)  # Drop the mean_distance column

saveRDS(data_cam_dd, camera_data_output)
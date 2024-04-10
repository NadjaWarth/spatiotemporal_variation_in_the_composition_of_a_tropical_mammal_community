# -----------------------------------------------------------------------------------------
# Script Name: Code 2_Capture_rate_function.R
# Author: Nadja Warth

# Description:
# This script is pivotal for filtering the camera trapping data to align with the specific 
# research period from December 16 to April 21. It meticulously crops deployments to ensure
# that both the start and end dates fall within this designated timeframe. The end date is 
# determined by the last day the camera was operational, either due to retrieval or the last 
# captured image in the event of a malfunction. Any observations recorded outside the research 
# period are systematically excluded. Additionally, the script filters out any deployments 
# lasting less than 14 days due to their insufficiency for the study's needs.
#
# Following the data filtering, the script then embarks on calculating the effort and capture 
# rates, which are essential metrics for the study.
#
# Usage:
# - Install and load all required packages
# - Update the paths to input files and desired output location as needed.
# - Adjust the parameters of the capture rate function to match the research period (years and timeframes)
#
# NOTE: Ensure the dataset prepared by "Code 1_Load and Merge Data" is correctly formatted and accessible.
# -----------------------------------------------------------------------------------------

### Load Packages ----
library(vegan)
library(tidyr)
library(dplyr)

### Variables used for this script ----

# Sources
data_cam_path <- "C:/NadjaMA/processed/data_cam_dd_23.11.rds"

# Outputs
capture_rates_output_path <- "C:/NadjaMA/processed/Capture_Rate_Data_new_dates.rds"

# Load camera data ----
data_cam_dd<- readRDS(data_cam_path)

# 1) define time constraint and estimate effort function -----------

# estimate_effort_within_timeframe function: 

#' Function to estimate the effort within specific timeframes in a year.
#' If a timeframe covers two years, important is the start of the timeframe.
#' (e.g timeframe: ('12-01', '01-15'), years: (2015, 2018, ...) -> timeframes to analyze: (2015-12-01, 2016-01-15), (2018-12-01, 2019-01-15) ...)
#' 
#' If the Start / End Date of a deployment lie outside the timeframe boundaries, the start / end will be adjusted to be the first / last day of the period. 
#' Sequences will be filtered as well (based on the startDate of the sequence) and the Timespan of a deployment will be estimated based on the adjusted dates. 
#' Deployments that lie entierly outside of the given period will be disregarded
#'
#' Effort = detection distance (m) * timespan (days)
#'
#' @param data Input data: Dataframe with the following collumns expected: StartDate, EndDate, SequDate, Detection_Distance.
#' @param years list of years to analyze within the data. If no years are specified, all occurring years in the data are taken
#' @param timeframes list of timeframes to analyze. A timeframe consists of a (start,end)-tupel in the format ('%m-%d', '%m-%d'). Timeframes should not overlap.
#' 
#' @return list of results for each year
#'
#' @examples
#' \dontrun{
#' # define test_data
#' test_data <- data.frame(
#'   deploymentID = c("06ee9610-93c3-4de5-ade7-0a507b94ad71","06ee9610-93c3-4de5-ade7-0a507b94ad71","044c34ac-7082-4240-9457-ac865aa55fcc", "044c34ac-7082-4240-9457-ac865aa55fcc", "0009d0dc-6337-47d8-8bb4-10da7bb786b4"),
#'   StartDate = as.Date(c("2015-01-15", "2015-03-05", "2016-03-10", "2019-08-20", "2020-02-05"), format = "%Y-%m-%d"),
#'   EndDate = as.Date(c("2015-02-28", "2015-05-07","2016-11-30", "2019-12-31", "2020-12-31"), format = "%Y-%m-%d"),
#'   SequDate = c(as.POSIXct(strptime("2015-01-17 01:30:00", "%Y-%m-%d %H:%M:%S")), as.POSIXct(strptime("2015-03-05 06:30:00", "%Y-%m-%d %H:%M:%S")), as.POSIXct(strptime("2016-04-12 01:30:00", "%Y-%m-%d %H:%M:%S")),as.POSIXct(strptime("2019-08-20 01:30:00", "%Y-%m-%d %H:%M:%S")),as.POSIXct(strptime("2020-02-27 01:30:00", "%Y-%m-%d %H:%M:%S"))),
#'   Detection_Distance = c(10, 5, 15, 10, 20)
#'   )
#' # define timeframes
#' timeframes <- list(c('02-01', '04-30'), c('05-01', '08-31'), c('09-01', '12-31'))
#' 
#' # (optional) define years
#' years <- c(2016,2017,2018)
#' 
#' results <- estimate_effort_within_timeframe(data, years, timeframes)
#' }
#' 
estimate_effort_within_timeframe <- function(data, timeframes, years = NULL) {
  # TODO: Add 0's in capture rate instead of null -> jedes Deployment sollte capture Rate für ALLE Tiere haben
  # TODO: Checks (optional)
  # Check data
  #check if data is in the right format (e.g if StartDate and EndDate exists)
  
  # Check timeframes
  # check if list is not empty, (so at least one timeframe is given) -> error or default c("01-01","12-31")?
  # for all timeframes check:
  # is in right format? (eg. %m-%y) -> error
  # valid dates? (e.g c("40-12", "03-30")) is wrong -> error
  # overlapping timeframes? -> error
  
  # Check years
  # If no years are specified, all take all occurring years in the data
  if (is.null(years)) {
    # Extract the year component from the startDate and endDate columns
    start_years <- as.integer(format(data$StartDate, "%Y"))
    end_years <- as.integer(format(data$EndDate, "%Y"))
    # Combine the years and remove duplicates
    years <- unique(c(start_years, end_years))
    # Sort the years in ascending order
    years <- sort(years)
  }
  
  results <- list()  # Create an empty list to store results for each year
  for (year in years) {
    results_per_timeframe <- list()
    for (timeframe in timeframes) {
      timeframe_start <- as.Date(paste(year, timeframe[1], sep = "-"))
      timeframe_end <- as.Date(paste(year, timeframe[2], sep = "-"))
      
      if (timeframe_start > timeframe_end) {
        timeframe_end <- as.Date(paste(year + 1, timeframe[1], sep = "-"))#Warum +1?? nicht -1?????
      }
      
      truncate_deployments_to_timeframe <- data %>%
        #snap deployment start and end dates to the dates of the timeframe
        mutate(
          StartDate = ifelse(StartDate < timeframe_start & EndDate > timeframe_start, timeframe_start, StartDate),
          EndDate = ifelse(StartDate < timeframe_end & EndDate > timeframe_end, timeframe_end, EndDate)
        ) %>%
        mutate(
          StartDate = as.Date(StartDate),
          EndDate = as.Date(EndDate), 
          SequDate = as.Date(SequDate)
        ) %>% 
        #filter out all deployments outside the timeframe
        filter(StartDate >= timeframe_start & StartDate <= timeframe_end & 
                 EndDate >= timeframe_start & EndDate <= timeframe_end) %>%
        #keep only sequences inside the timeframe
        filter(SequDate >= StartDate & SequDate <= EndDate)
      # Calculate effort per deployment based on first/last sequence and deployment duration
      
      effort_per_deployment <- mutate(truncate_deployments_to_timeframe,
                                      time_span_days = as.numeric(difftime(EndDate, StartDate, units = "days")),
                                      effort = time_span_days * Detection_Distance
      ) %>%
        # Remove all NA rows (false detections)----- still necessary?
        drop_na()
      #filter unneccessary columns?
      
      results_per_timeframe[[paste(timeframe[1], "-", timeframe[2])]] <- effort_per_deployment
    }
    results[[paste0(year)]] <- results_per_timeframe
  }
  return(results)
}



# 2) define function to calculate the capture rate -----------------------------
  
#' Function: capture_rate_fun
#'
#' 
#' Unweighted: Number of observations per species in a given location divided by the number of observations 
#' capture rte / effort (overall, does not consider different lengths of different deployments)
#' Weighted Capture rate: capture rate per species per location based on effort (= time_span_days * Detection_Distance) of each deployment 
#' -->capturerate per day per m 
#'
#' @param data Input data: Dataframe with the following collumns expected: StartDate, EndDate, SequDate, Detection_Distance.
#' @param years list of years to analyze within the data. If no years are specified, all occurring years in the data are taken
#' @param timeframes list of timeframes to analyze. A timeframe consists of a (start,end)-tupel in the format ('%m-%d', '%m-%d'). Timeframes should not overlap.
#' 
#'
#' @return data frame containing the calculated weighted and unweighted capture rates 

capture_rate_fun <- function(data, timeframes, years) {
  results <- list()  # Create an empty list to store results for each year
  taxonomy <- as.data.frame(read.csv("C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Species_list.csv"))
  
  
  #load cropped data that contains estimated effort 
  effort_data <- estimate_effort_within_timeframe(data,timeframes, years) 
  
  for (year in years) {
    capture_rate <- effort_data[[as.character(year)]][['12-16 - 04-21']] %>% 
      select(-c("SequDate"))
  
    deployments <- capture_rate %>% 
      select(-c("scientificName",  "vernacularNames.en", "count", "class")) %>% 
      distinct()
    
    combinations_of_depl_taxonomy <- deployments %>% 
      crossing(taxonomy) 
    
    data <- capture_rate %>% 
      select(c("deploymentID", "scientificName", "count")) %>% 
      
      #Count the total amount of observations per Species in each deployment 
      group_by(deploymentID, scientificName) %>%
      mutate(count_per_spec_depl = sum(count)) %>%
      ungroup() %>% 
      #Remove redundant columns 
      select(-c("count")) %>%
      #Remove identical Rows 
      distinct()
      
      
      #complete function ensures that all combinations of deploymentID and scientificName are present in the result, filling in missing combinations with a count of 0.
      #complete(deploymentID, scientificName,  fill = list(count_per_spec_depl = 0)) 
      
    
    capture_rate_data <- left_join(combinations_of_depl_taxonomy,data, by =c("deploymentID", "scientificName")) %>% 
      mutate(count_per_spec_depl = ifelse(is.na(count_per_spec_depl), 0, count_per_spec_depl)) %>% 
      
      # Calculate weighted capture rate per species per location based on effort 
      #(effort1 * cr_depl1 + effort2 * cr_dpl2...)/total_effort 
      # = (effort1 * (count1 / effort1) + effort2 * (count2 / effort2)...)/total_effort
      # = (count1+count2)/total effort 
      
      #Calculate summed effort and timespan per location 
      group_by(locationName) %>% 
      mutate(sum_effort = sum(unique(effort)), 
             sum_time = sum(unique(time_span_days))) %>% 
      ungroup() %>% 
      
      #Calculate summmed count per spec per location 
      group_by(locationName, scientificName) %>% 
      mutate(count_per_spec=sum(count_per_spec_depl)) %>%
      #Calculate cr_weighted base on the count of a species at a given location and the overall effort at this location
      mutate(cr_weighted = count_per_spec/sum_effort) %>% 
      ungroup() %>% 
      
      # remove all redundant columns (Detection_Distance","time_span_days","effort" are removed because they are different for each deployment, but we now ignore deployments and look at all deployments in a period)
      select(c("locationName","longitude", "latitude","scientificName","vernacularNames.en","class",
               "count_per_spec","sum_effort","sum_time","cr_weighted")) %>%
      #drop all rows that have an identical combination of location Name and Scientific species name 
      #different rows exist because of different detection distances etc. in different rows, which leads to multiple rows of one location and a given species 
      distinct(locationName, scientificName, .keep_all= TRUE) %>% 
      drop_na() %>% 
      
      #Add a year column 
      mutate(Year = year) 
    
    results[[as.character(year)]] <- capture_rate_data # Store the results in the list with the year as a key
  }
  #Combine all subtables in one big table 
  capture_data <- bind_rows(results)%>% 
    mutate_all(~replace_na(., 0)) 
  
  return(capture_data)
}

# 3) Execute capture rate function and Save data---- 
years_to_analyze<- c(2011:2022) # Replace with the years you want to analyze
timeframes <- list(c('12-16', '04-21'))
capture_data<- capture_rate_fun(data = data_cam_dd,timeframes = timeframes,   years =years_to_analyze)

#Save data
saveRDS(capture_data, capture_rates_output_path)
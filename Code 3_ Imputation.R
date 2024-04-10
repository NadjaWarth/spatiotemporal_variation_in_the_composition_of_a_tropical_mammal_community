#CODE 3

### Load packages----
library(vegan)
library(tidyr)
library(dplyr)
library(stringr) 
library(mice)
library(ggplot2)

### Variables used for this script ----

# Sources
capture_rate_data_path <- "C:/NadjaMA/processed/Capture_Rate_Data_new_dates.rds"

# Outputs
cart_imputed_data_output_path <- "C:/NadjaMA/processed/Imputed_data_cart.rds"
pmm_imputed_data_output_path <- "C:/NadjaMA/processed/Imputed_data_pmm.rds"
rf_imputed_data_output_path <- "C:/NadjaMA/processed/Imputed_data_rf.rds"

# 1) Load and prepare Camera data ----
Cam_data<- readRDS(capture_rate_data_path)

# Adjust camera data: 
#Species of interest: 
species_list <- c( "Central American agouti",
                   "collared peccary",
                   "white-nosed coati",
                   "Tomes spiny rat",
                   "common opossum",
                   "lowland paca")

#a) filter for desired species
#b)filter out deployments with time > 13 days, 
#c)Filter out years: 2019,2021 (no / to little data available) --> with new data format this means 2018 and 2020
Cam_data<- Cam_data %>% 
  filter(vernacularNames.en %in% species_list, 
         sum_time >13, 
         Year != 2018) %>%
  mutate(vernacularNames.en = str_to_title(vernacularNames.en))


#Bring data in the right format: 
Cam_dat_ST<- Cam_data %>% 
  arrange(Year) %>% 
  mutate(locationName_Year = paste(locationName, Year, sep = "_")) %>%
  select(locationName_Year, everything()) %>%
  select(c(locationName_Year,vernacularNames.en,cr_weighted)) %>% 
  pivot_wider(names_from = vernacularNames.en, values_from = cr_weighted, values_fill = 0) 

#Add missing deployments / Year-location combinations:
#Generate complete list of all Year-location combinations: 
location_Names <- data.frame(locationName = unique(Cam_data$locationName))
Years <- data.frame(Year = c(2011:2017, 2019, 2021))
all_combinations <- crossing(location_Names, Years) %>%
  mutate(locationName_Year = paste(locationName, Year, sep = "_"))

#Combine complete list with Cam data and generate NA's where data is missing: 
ST_data3<- merge(Cam_dat_ST, all_combinations , by = "locationName_Year", all ="TRUE") %>% 
  rename(Coati = `White-Nosed Coati`)
colnames(ST_data3) <- gsub(" ", "_", colnames(ST_data3))

# 2) Create imputed dataset performing Imputations using PMM, RF, and CART -------

# Define the methods to compare, excluding 'lasso.norm'
methods <- c("pmm", "rf", "cart")

# Initialize a list to store imputed datasets
imputed_datasets <- list()

# Perform imputation for each method
for (method in methods) {
  temp_data <- mice(ST_data3, method = method, m = 1, maxit = 5, print = FALSE)
  imputed_datasets[[method]] <- complete(temp_data)
}

#Check out warning message
temp_data <- mice(ST_data3, method = methods[1], m = 1, maxit = 5, print = FALSE)
summary(temp_data)

### Compare Summary Statistics
# 3) Define functions to print and compare summary statistics ----

#Functions for comparing the summary statistics of the original data with the imputed data for species_name
compare_summary_statistics <- function(data, imputed_datasets, species_name) {
  # Original data summary
  summary_original <- summary(data[[species_name]])
  # Imputed data summary
  summary_imputed <- sapply(imputed_datasets, function(x) summary(x[[species_name]]))
  
  print(paste0("Original Summary for ", species_name))
  print(summary_original)
  print(paste0("Imputed Summary for ", species_name))
  print(summary_imputed)
}

### Visual Comparison of Distributions
#This function creates visualizations to compare the distributions of the original and imputed data.
visual_comparision_of_Distributions <- function(data, imputed_datasets, species_name) {
  return (ggplot() +
    geom_density(data = data, aes_string(x = species_name), fill = "blue", alpha = 0.5, color = "blue") +
    geom_density(data = imputed_datasets[["pmm"]], aes_string(x = species_name), alpha = 0.5, color = "red") +
    geom_density(data = imputed_datasets[["rf"]], aes_string(x = species_name), alpha = 0.5, color = "green") +
    geom_density(data = imputed_datasets[["cart"]], aes_string(x = species_name), alpha = 0.5, color = "orange") +
    labs(title = paste0("Distribution of '",species_name, "' Capture Rates"), x = "Capture Rate", y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = c("Original" = "blue", "PMM" = "red", "RF" = "green", "CART" = "yellow")) +
    guides(fill = guide_legend(title = "Method")))
}
# 4) Print and compare summary statistics for all relevant species ----
#1 Paca
compare_summary_statistics(ST_data3, imputed_datasets, "Lowland_Paca")
visual_comparision_of_Distributions(ST_data3, imputed_datasets, "Lowland_Paca")

#2 "Central_American_Agouti"
compare_summary_statistics(ST_data3, imputed_datasets, "Central_American_Agouti")
visual_comparision_of_Distributions(ST_data3, imputed_datasets, "Central_American_Agouti")

#3 Common_Opossum
compare_summary_statistics(ST_data3, imputed_datasets, "Common_Opossum")
visual_comparision_of_Distributions(ST_data3, imputed_datasets, "Common_Opossum")
    
#4 Coati
compare_summary_statistics(ST_data3, imputed_datasets, "Coati")
visual_comparision_of_Distributions(ST_data3, imputed_datasets, "Coati")

#5 Collared_Peccary
compare_summary_statistics(ST_data3, imputed_datasets, "Collared_Peccary")
visual_comparision_of_Distributions(ST_data3, imputed_datasets, "Collared_Peccary")

#6 Tomes_Spiny_Rat
compare_summary_statistics(ST_data3, imputed_datasets, "Tomes_Spiny_Rat")
visual_comparision_of_Distributions(ST_data3, imputed_datasets, "Tomes_Spiny_Rat")

# 5) Save imputed data for all methods --------
#Cart:
imputed_data_cart <- mice(ST_data3, method = "cart", m = 1, maxit = 5)
completed_data <- complete(imputed_data_cart)
saveRDS(completed_data, cart_imputed_data_output_path)


#PMM: 
imputed_data_pmm <- mice(ST_data3, method = "pmm", m = 1, maxit = 5)
completed_data_pmm <- complete(imputed_data_pmm)
saveRDS(completed_data_pmm, pmm_imputed_data_output_path)


#RF
imputed_data_rf <- mice(ST_data3, method = "rf", m = 1, maxit = 5)
completed_data_rf <- complete(imputed_data_rf)
saveRDS(completed_data_rf, rf_imputed_data_output_path)
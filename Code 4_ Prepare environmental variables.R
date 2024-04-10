#CODE 4
#prepare Environmental Data 

### Load packages ----
library(readxl)
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)

### Variables needed for this script ----
# Sources
environmental_variables_bcinm_source <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/Input for code/Environmental_variables_BCINM.xlsx"
dry_mass_data_source <-"C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/Input for code/Plot50ha.csv"
four_spp_source <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/Input for code/Plot50haFourSpp.csv"
rain_raw <- read_excel("C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/Input for code/Monthly summaries_BCI_horizontal.xlsx", 
                      sheet = "Rain_Seasonal", 
                      skip = 1, col_names = TRUE)
# Outputs
spatial_data_output_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Spatial_data.rds"
temporal_data_output_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Temporal_data.rds"

#1) Load Spatial data ----------
#1.1) Load spatial variables
Spat_dat_raw <-read_excel(environmental_variables_bcinm_source)

#1.2) Create clean dataset with relevant columns 

#Select all numeric variables
numeric_columns <- Spat_dat_raw%>%
  select_if(is.numeric) %>%
  colnames()

Spat_dat_processed<- Spat_dat_raw %>% 
  #Fill all NA in numeric variables with 0
  mutate(across(all_of(numeric_columns), ~ replace_na(., 0))) %>% 
  #Select all relevant variables 
  select(-c(`Slope_Class (?)`, Descibtion_Soil,`NM_Border_Dist[m]`,`Path_Dist[m]_first version`, `Path_Dist[m]`,
            BCNM,`Forest Age`, forest_age_estimation )) %>% 
  rename(Forest_Border_Dist_m = `Forest_Border_Dist[m]`,
         Water_distance = `Water distance`, 
         AST_trees = Astrocaryum,
         ATTA_trees = Attalea, 
         DIP_trees = Dipteryx, 
         GUS_trees = Gustavia, 
         Forest_age = `forest_ age`) %>% 
  # filter locationName for BCI-1-..locations
  filter(grepl("BCI-1-", locationName)) %>% 
  drop_na()

#Convert Soil and Region columns to binary columns
binary_columns_soil <- model.matrix(~ Soil - 1, data = Spat_dat_processed)
binary_columns_region <- model.matrix(~ Region - 1, data = Spat_dat_processed)

Spat_dat <- cbind(Spat_dat_processed, binary_columns_soil, binary_columns_region) %>% 
        mutate(East = rowSums(select(.,c("RegionFrijoles", "RegionBuena Vista", "RegionBohio"))), 
       West =  rowSums(select(.,c("RegionGigante", "RegionPena Blanca"))), 
       BCI = RegionBCI) %>% 
       select(-c("array","Region", "RegionBCI", "RegionGigante", "RegionPena Blanca","RegionFrijoles", "RegionBuena Vista", "RegionBohio" )) %>% 
  mutate(Region = case_when(
    West == 1 ~ "West",
    East == 1 ~ "East",
    BCI == 1 ~ "BCI",
    TRUE ~ NA_character_
  ))
unique(Spat_dat$Soil)    

#Save the data: 
saveRDS(Spat_dat,file = spatial_data_output_path)

#2) Load Temporal data---------
# 2.1.) Load temporal data on fruit counts and rainfall ----

# Data including total seed plus fruit dry mass summed over all species from 2013 to the present
# Explanation on data: 
#   Census – an integer to identify individual census 
#   Value – dry mass of seeds plus fruit summed over all species for each value of census
#   Mean_date – the mean date weighted by number of traps censuses each day for the census
dry_mass <- read_csv(dry_mass_data_source)

# Data on 4 species
# Explanation on data:
#   sp - species abbreviation: 
#   ASTS = Astrocaryum 
#   ATTB = Attalea 
#   DIPP = Dipteryx 
#   GUSS = Gustavia
#   fecha - census date 
#   part - Integer to identify reproductive structures (relevant here: 1)
FourSpp_count <- read_csv(four_spp_source)

# Rainfall per month
colnames(rain_raw) <- gsub("\\.\\.\\..*", "", colnames(rain_raw))
Rain <- rain_raw %>% select(1:4)

# 2.2) Filter and combine Temporal data ----
#   a) Fruit data ----

#' Data_fun_30_30
#'
#' This function crops the data from dry_mass and four_spp_count 
#' to the period from the End of last years period (30.03.year-1) to the end of this years period (30.03.year)
#' Then estimates the average mean value of overall fruit dry mass
#' and the overall count of seeds for the 4 tree species
#'
#' @param dry_mass 
#' @param FourSpp_count 
#' @param years
#'
#' @return dataframe with the average mean value of overall fruit dry mass and the overall count of seeds for the 4 tree species
Data_fun_30_30 <- function(dry_mass, FourSpp_count, years) {
  results <- list()  # Create an empty list to store results for each year
  
  for (year in years) {
    Start_this_year <- as.Date(paste0(year, "-12-16"))
    End_this_year <- as.Date(paste0(year+1, "-04-21"))
    End_last_year  <- as.Date(paste0(year, "-04-21"))

    
    mean_dry_mass <- dry_mass %>% 
        # a) Mean value from 30.03-30.03
        mutate(mean_date = as.Date(mean_date, format = "%m/%d/%Y")) %>% #Somehow mean_dat ends up in the format Y-m-d, which I want, but I dont know how it gets there on it's own 
          filter(mean_date >= End_last_year, mean_date <= End_this_year) %>%
          summarise(mean_dry_mass = mean(VALUE, na.rm = TRUE))%>%
        mutate(Year = year)
      
    sum_spp_count <- FourSpp_count %>% 
        mutate(mean_date = as.Date(fecha, format = "%m/%d/%Y")) %>% #Somehow mean_dat ends up in the format Y-m-d, which I want, but I dont know how it gets there on it's own 
        filter(mean_date >= End_last_year, mean_date <= End_this_year) %>% 
        filter(part == 1) %>% 
        group_by(sp) %>% 
        summarise(sum_quantity = sum(quantity, na.rm = TRUE)) %>% 
        ungroup() %>% 
        pivot_wider(names_from = sp, values_from = sum_quantity, values_fill = 0) %>% 
        mutate(Year = year)
    
    fruitdata <- merge(mean_dry_mass, sum_spp_count)
    results[[as.character(year)]] <- fruitdata
  }
  combined_table <- bind_rows(results) %>% 
    mutate_all(~replace_na(., 0))
  
  return(combined_table)
}
  

years <- c(2010:2022) # one year earlier than needed 
Fruits_30_30 <- Data_fun_30_30(dry_mass, FourSpp_count, years) %>% 
rename( AST_fruit ="ASTS",  
        ATTA_fruit="ATTB", 
        DIP_fruit="DIPP",
        GUS_fruit="GUSS") %>% 
select(-c(mean_dry_mass)) %>% 
mutate(Total_count = rowSums(select(., AST_fruit, DIP_fruit, GUS_fruit, ATTA_fruit)), 
       Fruit_change = Total_count-lag(Total_count)) %>% 
filter(Year %in% c(2010:2022))

summary(Fruits_30_30)

#   b) Rain data ----
years_to_analyze<- c(2010:2022) # Replace with the years you want to analyze
test_timeframes <- list(c('12-16', '04-21'))
Rain_data <- Rain %>% 
  mutate(Year = as.numeric(year),
         #Rain total = rain of last wet season + current dry season 
         Rain_total = total,
         #Rain change = rain in current dry season minus rain in past season 
         Rain_change = dry - lag(dry)) %>% 
  filter(year %in% years_to_analyze) %>% 
  select(c(Year, Rain_total, Rain_change))

#   c) combine Rain and Fruit data and Save the result ----
Temp_dat <- merge(Rain_data, Fruits_30_30) %>% 
  mutate(Year = as.character(Year))

saveRDS(Temp_dat,file = temporal_data_output_path)

#3) Overview Graphs ------------
# 3.1) Temporal Data---------
#   a) Rain per year (Jan-Mar) ----
Rain_per_year <- ggplot(Temp_dat, aes(x = Year, y = Rain)) +
  geom_bar(stat = "identity", 
           fill = "lightgrey", 
           color = "black") +
  labs(title = "Precipitation per research period",
       x = "Year",
       y = "Precipitation [mm]") +
  theme_minimal()+
  theme(axis.text.x = element_text( hjust = 0.5, vjust = 5), 
    axis.text.y = element_text( hjust = 0.5, vjust = 0.5),
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))+  # Adjust the text size) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5)  # Adjust size and center the title
  )

print(Rain_per_year)

#   b) Fruit count per year ----

# Exclude the "Rain" variable
temp_dat_subset <- select(Temp_dat,c("Year", "Dipterix", "Gustavio", "Astrocaryum" ,"Attalea"))

# Reshape the data to long format
long_data <- gather(temp_dat_subset, key = Variable, value = Value, -Year)

# Use the ggplot function to create a histogram
Fruit_count_per_year <- ggplot(long_data, 
                               aes(x = Year, y = Value, 
                                   fill = factor(Variable, levels = c("Attalea", "Astrocaryum", "Gustavio", "Dipterix")))) +
  geom_bar(stat = "identity") +
  labs(title = "Seed / Fruit count per species per Year",
       x = "Year",
       y = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text( hjust = 0.5), 
    text = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  guides(fill = guide_legend(title = "Tree species", title.theme = element_text(size = 16)))+
  scale_fill_manual(values = c("Dipterix" = "lightgreen", "Gustavio" = "lightblue","Astrocaryum" = "orange", "Attalea" = "yellow"))

print(Fruit_count_per_year)

#   c) Mean Fruit dry mass per year ----
Fruit_dm_per_year <- ggplot(Temp_dat, aes(x = Year, y = mean_dry_mass)) +
  geom_bar(stat = "identity", 
           fill = "lightgrey", 
           color = "black") +
  labs(title = "Mean fruit dry mass per year",
       x = "Year",
       y = "Mean dry mass") +
  theme_minimal()+
  theme(axis.text.x = element_text( hjust = 0.5, vjust = 5), 
        axis.text.y = element_text( hjust = 0.5, vjust = 0.5),
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+  # Adjust the text size) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5)  # Adjust size and center the title
  )

print(Fruit_dm_per_year)

# 3.2) Spatial Data -----------------
      
#   a) All Fruit trees ----

# Create a data frame with the selected columns for plotting
plot_data <- Spat_dat %>%
  select(locationName, Astrocaryum, Attalea, Dipterix, Gustavio)

# Reshape the data for plotting (from wide to long format)
plot_data_long <- plot_data %>%
  pivot_longer(cols = c(Astrocaryum, Attalea, Dipterix, Gustavio), 
               names_to = "Variable", values_to = "Value")

# Create the plot
ggplot(plot_data_long, aes(x = locationName, y = Value, fill = Variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of trees per location",
       x = "Location",
       y = "Tree count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14))+  # Adjust the text size) +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5)  # Adjust size and center the title
  )+
  guides(fill = guide_legend(title = "Tree species",   # Customize the legend title
                             title.theme = element_text(size = 16)))+  # Adjust the legend title size
  scale_fill_manual(values = c("Astrocaryum" = "orange", "Attalea" = "yellow", "Dipterix" = "lightgreen", "Gustavio" = "lightblue"))

#   b) Other numeric variables ----
# For example: Path_Dist_m, Forest_Border_Dist_m, Water_distance, Slope_degree, Forest_Area_ha
colnames(Spat_dat)
ggplot(Spat_dat, aes(x = locationName, y = Water_distance)) +
  geom_bar(stat = "identity") +
  labs(title = "Distance to Water",
       x = "Location",
       y = "Water distance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
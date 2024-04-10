#CODE 7
#Maps & Other general graphs 

### Load Packages ----
library(ggplot2)
library(vegan)
library(readxl)
library(ggmap)
library(leaflet)
library(tidyverse)
library(sf)
library(writexl)
library(scatterpie)

### Variables used in this script ----
# Sources
#From Code 2
capture_rates_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Capture_Rate_Data_new_dates.rds"
#From Code 4
spatial_data_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Spatial_data.rds"
temporal_data_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Temporal_data.rds"

island_shp <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/intermediate/croppe map/cropped_map.shp"
Rain_raw<- read_excel("C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/Input for code/Monthly summaries_BCI_horizontal.xlsx", 
                      sheet = "Rain_Seasonal", 
                      skip = 1, col_names = TRUE)
# Outputs
heatmap_effort <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/heatmap_effort_incl2018.png"
capture_rates_boxplot_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Capture_rates_Boxplot_new.png"
fruit_plot_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Fruit plot_minus22.png"
rain_plot_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Rain_plot_minus22.png"
summary_statistics_xlsx <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/summary_stats.xlsx"

#1) Load the data ---------
Spat_dat <- readRDS(spatial_data_path)
Temp_dat_all<- readRDS(temporal_data_path)
Cam_data_raw<- readRDS(capture_rates_path)

#Filter for species of interest in Cam_data_raw
species_list <- c( "Central American agouti",
                   "collared peccary",
                   "white-nosed coati",
                   "Tomes spiny rat",
                   "common opossum",
                   "lowland paca")

#Filter out deployments with time > 13 days, 
#Filter out years: 2019,2021 (because no or to little data available in these years)
Cam_data<- Cam_data_raw %>% 
  filter(vernacularNames.en %in% species_list, 
         sum_time >13,
         Year != 2018) %>% 
  mutate(vernacularNames.en = str_to_title(vernacularNames.en))

Cam_data_all <- Cam_data_raw %>% 
  filter(sum_time >13,
         Year != 2018)

#2) Generate General location Map -----------------------------------

locations <- Cam_data %>% 
  select(c("locationName", "longitude", "latitude")) %>% 
  distinct(locationName, .keep_all = TRUE)

#Load Island shape 
island_shape <- st_read(island_shp)

ggplot() +
  geom_sf(data = island_shape, fill = "lightgreen", color = "white", alpha = 0.5) +
  geom_point(data = locations, aes(x = longitude, y = latitude), color= "darkred", size = 2) +
  geom_text(data = locations, aes(x = longitude, y = latitude, label = locationName), 
            vjust = +1.5, size = 2.5) +
  theme_minimal()

#3) Generate Heatmap effort & Overview over missing deployments  ----------------------------

### Generate Heatmaps ###
heatmap_time_plot <- ggplot(Cam_data_all, aes(x = as.character(Year), y = locationName, fill = sum_time)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "green", limits = c(0, max(Cam_data_all$sum_time))) +
  labs(title = "Deployment time per Location per Year", x = "Year", y = "Location")

heatmap_effort_plot <- ggplot(Cam_data_all, aes(x = as.character(Year), y = locationName, fill = sum_effort)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "green", limits = c(0,max(Cam_data_all$sum_effort))) +
  labs(title = "Effort per Location per Year", x = "Year", y = "Location")

ggsave(heatmap_effort, plot = heatmap_effort_plot, width = 6, height = 6, units = "in")

# Heatmap for specific species (swap species name in plot)
agouti<- Cam_data %>% 
  filter(scientificName=="Dasyprocta punctata")
coati <- agouti<- Cam_data %>% 
  filter(scientificName=="Nasua narica")
pecari <- Cam_data %>% 
  filter(scientificName=="Pecari tajacu")

heatmap_effort_plot <- ggplot(pecari, aes(x = as.character(Year), y = locationName, fill = cr_weighted)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "red", midpoint = median(pecari$cr_weighted),
                       limits = c(min(agouti$cr_weighted), max(pecari$cr_weighted))) +
  labs(title = "Deployment time per Location per Year", x = "Year", y = "Location")

### Getting a Overwiev over missing deployments ###

# Generate all combinations of locationName and Year present in the dataset
all_combinations <- Cam_data_all %>%
  distinct(locationName, Year) %>%
  expand(locationName, Year)

# Left join the complete set of combinations with the original dataset to find missing entries
missing_combinations <- all_combinations %>%
  left_join(Cam_data_all, by = c("locationName", "Year")) %>%
  filter(is.na(sum_time)) %>%
  group_by(locationName, Year) %>%
  summarize(missing_entries = n(), .groups = 'drop')

# Identify combinations with missing sum_time values
missing_combinations <- all_combinations %>%
  anti_join(Cam_data_all, by = c("locationName", "Year")) %>%
  arrange(locationName, Year)

# How many deployments were missing per location?
missing_per_location <- missing_combinations %>%
  group_by(locationName) %>%
  summarize(missing_years_count = n(), .groups = 'drop') %>% 
  arrange(desc(missing_years_count))

# How many deployments were missing per year?
missing_per_year <- missing_combinations %>%
  group_by(Year) %>%
  summarize(missing_years_count = n(), .groups = 'drop') %>% 
  arrange(desc(missing_years_count))

sum(missing_per_year$missing_years_count)

#4) Identify important values on camera traps-------------------

#Count number of observations: (observed animals)
sum(Cam_data$count_per_spec)

Cam_data %>%
  group_by(vernacularNames.en) %>%
  summarize(total_count = sum(count_per_spec))

#How long were Cameras going? 
Time_data <- Cam_data %>% 
  select(c(locationName, sum_time, sum_effort, Year)) %>% 
  filter(!duplicated(.))
sum(Time_data$sum_time)
range(Time_data$sum_time)
mean(Time_data$sum_time)
summary(Time_data)

#How many relevant deployments?
nrow(Time_data)
colnames(Cam_data)

#5) Generate Plots for Variation in species composition------------------------------
#   a) Ordination plot: Variation in species composition-----------
#   ### Prepare data ----
# Remove the years 2010, 2018 and 2020 (lacking data quality)
Temp_dat <- Temp_dat_all %>% 
  filter(!(Year %in% c("2010","2018", "2020")))

# DCA - compare to Code_5 - Ordination Analysis 5)
#Prepare Camera dat
Cam_dat_ST<- Cam_data %>% 
  arrange(Year) %>% 
  mutate(locationName_Year = paste(locationName, Year, sep = "_")) %>%
  select(locationName_Year, everything()) %>%
  select(c(locationName_Year,vernacularNames.en,cr_weighted)) %>% 
  pivot_wider(names_from = vernacularNames.en, values_from = cr_weighted, values_fill = 0) 

dca_result <- decorana(Cam_dat_ST[,-1])

#Prepare environmental data: Make one big table with all environmental variables (Spat + Temp)
#a) Create 10 replicates of the Spat_dat table
replicates <- replicate(10, Spat_dat, simplify = FALSE)#10 is the number of rows of the Tep_data 

#b) Bind replicates together
replicate_spat_data <- bind_rows(replicates, .id = "Replicate") %>% 
  #Create a new column "Year" with values from Temp_dat
  mutate(Year = rep(Temp_dat$Year, each = nrow(Spat_dat))) %>% 
  # Arrange the data based on the Replicate and Year columns
  arrange(Replicate, Year)

#c) Left join Temp_dat with combined_data using "Year" as the key
location_info<-    left_join(replicate_spat_data, Temp_dat, by = "Year") %>% 
  #Create a new column as identifier for each row 
  mutate(locationName_Year = paste(locationName, Year, sep = "_")) %>%
  select(locationName_Year, everything()) %>% 
  select(-c("Year", "locationName","Replicate"))

# Extract the necessary information from the decorana result
dca_data <- data.frame(locationName_Year = Cam_dat_ST$locationName_Year,
                       DCA1 = scores(dca_result, display = "sites")[, 1],
                       DCA2 = scores(dca_result, display = "sites")[, 2],
                       DCA3 = scores(dca_result, display = "sites")[, 3],
                       DCA4 = scores(dca_result, display = "sites")[, 4])

# Now, merge 'dca_data' with 'location_info' based on the common column 'Location'
merged_data <- merge(dca_data, location_info, by = "locationName_Year")

#   ### Plots ----

# Plot the DCA with color based on the "Region" variable
# Convert "Region" to a factor with distinct levels
merged_data$Region <- as.factor(merged_data$Region)
custom_colors <- c("#C3D1F7", "#a5cc6e", "#e68477")

# Match the levels of "Region" to the custom colors
region_colors <- custom_colors[match(levels(merged_data$Region), unique(merged_data$Region))]

# Plot the DCA with manually assigned colors
plot(merged_data$DCA1, merged_data$DCA2, col = region_colors, pch = 16, cex = 1.5, 
     #main = "DCA Plot", 
     xlim = c(-2, 2.5), ylim = c(-0.7, 2), 
     xlab = "DCA Axis 1", ylab = "DCA Axis 2")  # Add axis labels

# Add species names to the plot
text(dca_result, display = "species", cex = 0.9, col = "blue", pos = 3)

# Add legend to the plot
legend("topright", legend = levels(merged_data$Region), col = custom_colors, pch = 16, title = "Region")

# Add arrows to indicate species contributions
arrows(0, 0, scores(dca_result, display = "species")[, 1], scores(dca_result, display = "species")[, 2], col = "black", length = 0.1)

#6) Generate Capture Rates Boxplots: -----------------------

#Single plot for a species of interest 
species_of_interest <-  "lowland paca"

ggplot(data = Cam_data, aes(x = locationName, y = cr_weighted)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  facet_wrap(~vernacularNames.en, scales = "free_y", ncol = 2, nrow = 3) +  # 3 rows, 2 columns
  scale_y_log10() +
  labs(title = "Boxplots of Capture Rates",
       x = "Location", y = "Log Capture Rate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.title.x = element_text(margin = margin(t = 20)))

# Create a boxplots per species over all locations 
p<- ggplot(data = Cam_data, aes(x = locationName, y = cr_weighted)) +
  geom_boxplot(fill = "#C3D1F7", color = "black") +
  facet_wrap(~vernacularNames.en, scales = "free_y", ncol = 2, nrow = 3) +  # 3 rows, 2 columns
  scale_y_log10() +
  labs(title = "Capture Rates",
       x = "Location", y = "Capture Rate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.title.x = element_text(margin = margin(t = 20)))

# Save the plot as a PNG file
ggsave(capture_rates_boxplot_png, p, width = 9.27, height = 11.69, units = "in", dpi = 600)

#7) Generate Plots for Variation in environment -----------------------
#   a) Temporal (Fruits per year, Rain) -------------------------------------

Temp_dat_all <- Temp_dat_all[, c("Year", "Rain_total", "Rain_change", "DIP_fruit","AST_fruit",  "GUS_fruit", "ATTA_fruit", "Fruit_change")]

### Fruits per year ###
# Reshape the data
Temp_dat_long <- Temp_dat_all %>%
  filter(Year != 2022) %>% 
  pivot_longer(cols = c(DIP_fruit, AST_fruit, GUS_fruit, ATTA_fruit), 
               names_to = "Tree_Species", 
               values_to = "Count") %>%
  mutate(Tree_Species = factor(
    case_when(
      Tree_Species == "DIP_fruit" ~ "Dipteryx",
      Tree_Species == "AST_fruit" ~ "Astrocaryum",
      Tree_Species == "GUS_fruit" ~ "Gustavia",
      Tree_Species == "ATTA_fruit" ~ "Attalea"
    ),
    levels = c("Astrocaryum", "Gustavia", "Attalea","Dipteryx")
  ))

# Create the stacked barplot
Fruit_plot <- ggplot(Temp_dat_long, aes(x = Year, y = Count, fill = Tree_Species)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  labs(x = "Year", y = "Fruit and seed count", fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom",  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.title.x = element_text(size = 16), # Increase x-axis title size
        axis.title.y = element_text(size = 16), # Increase y-axis title size
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis text to be vertical
        axis.text.y = element_text(size = 14),  # Increase y-axis text size
        legend.title = element_text(size = 14), # Increase legend title size
        legend.text = element_text(size = 16))+  # Increase legend text size
scale_fill_manual(values = c("Dipteryx" = "#9AEB88", "Astrocaryum" = "#77E3DE" , 
                             "Gustavia" = "lightsalmon", "Attalea" = "#EEE055"))

ggsave(fruit_plot_png, plot = Fruit_plot, width = 6, height = 6, units = "in")

### Rain ###
#Rainfall per month 
colnames(Rain_raw) <- gsub("\\.\\.\\..*", "", colnames(Rain_raw))
Rain_1020<- Rain_raw %>% select(1:4) %>% 
  filter(year %in% c(2010:2021))

# Reshape the data
Rain_long <- Rain_1020 %>%
  select(year, wet, dry) %>%
  pivot_longer(cols = c(wet, dry), 
               names_to = "Season", 
               values_to = "Rainfall") %>%
  mutate(Season = factor(
    Season,
    levels = c("dry","wet")
  ))

# Create the stacked barplot
Rain_plot <- ggplot(Rain_long, aes(x = year, y = Rainfall, fill = Season)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  labs(x = "Year", y = "Rainfall (mm)",
       fill = "Season") +
  theme_minimal() +
  theme(legend.position = "bottom",  
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
              axis.title.x = element_text(size = 16), # Increase x-axis title size
              axis.title.y = element_text(size = 16), # Increase y-axis title size
              axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis text to be vertical
              axis.text.y = element_text(size = 14),  # Increase y-axis text size
              legend.title = element_text(size = 14), # Increase legend title size
              legend.text = element_text(size = 16))+ 
  scale_fill_manual(values = c("wet" = "#344D7E", "dry" = "lightsalmon"))

ggsave(rain_plot_png, plot = Rain_plot, width = 6, height = 6, units = "in")
#   b) Space ---------------------

### Plot Fruit tree map ###
Spat_dat %>% select(locationName, AST_trees, ATTA_trees, DIP_trees, GUS_trees, Trees_Total)

column_sums_soil <- colSums(Spat_dat[, c("SoilTb", "SoilTbo", "SoilTcm","SoilTcv", "SoilTue")])
colSums(Spat_dat[, c("Trees_Total","SoilTb" )])
table(Spat_dat$Forest_age)

locations <- Cam_data %>% 
  select(c("locationName", "longitude", "latitude")) %>% 
  distinct(locationName, .keep_all = TRUE)

# Merge camera locations with species count data
merged_data <- merge(location, Spat_dat, by = "locationName")

# Create data frame in the required format for mapPies function
pie_data <- merged_data %>%
  select(locationName, longitude, latitude, AST_trees, ATTA_trees, DIP_trees, GUS_trees, Trees_Total) 
  
# Load Island shape
island_shape <- st_read(island_shp)

# Plot
ggplot() +
  geom_sf(data = island_shape, fill = "lightgreen", color = "white", alpha = 0.5) +
  geom_point(data = pie_data , aes(x = longitude, y = latitude, size = Trees_Total), color= "darkred", shape = 21)+
  geom_scatterpie(data = pie_data, 
                  aes(x = longitude, y = latitude, r = 0.0001*Trees_Total), 
                  cols = c("AST_trees", "ATTA_trees", "DIP_trees", "GUS_trees")) +
  theme_minimal()

summary(Spat_dat)

#Summary Statistics
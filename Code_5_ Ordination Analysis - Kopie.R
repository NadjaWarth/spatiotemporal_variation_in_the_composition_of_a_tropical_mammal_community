# -----------------------------------------------------------------------------------------
# Script Name: Code 5 Ordination Analysis
# Author: Nadja Warth
# 
# Description:
# This script is designed to delve into the deterministic processes influencing wildlife 
# populations by addressing the third research question of the accompanying thesis: 
# "Which environmental factors are most important in explaining the spatio-temporal variation?". 
# 
# I performed the analysis on three levels. First, I focused on the spatial variation, 
# by estimating an average capture rate for each species across all years at a given location. 
# Second, I focused on the temporal variation, by estimating an average capture rate for 
# each species across all locations, to generate a single mean capture rate per species per year. 
# Lastly, I performed analysis using the values of all locations and years together.
# 
# For all 3 levels I performed indirect and direct gradient analysis, similar to Blake & Loiselle (2018). 
# First, I performed a detrended correspondence analysis (DCA) to obtain the length of the gradients in 
# standard deviations. If the length of the first gradient is larger than 4, unimodal models are appropriate, 
# if the length is smaller than 3, linear models are better suited (Blake & Loiselle, 2018). 
# Since the length of the gradient was >0.2 in all 3 cases, I chose to use a redundancy analysis as linear 
# model, to evaluate how much of the variance in the species composition can be linked to the selected 
# environmental variables. Before performing the RDA, I estimated the variance inflation factor of the variables 
# to check for multicollinearity and removed variables from the model if necessary.
# Afterwards I performed backward selection to test which variables are most relevant in improving the model fit 
# and used those variables in the RDA. 
# To test the significance of the RDA model and the used variables, I performed ANOVAs. 
# Additionally, I generated correlation matrixes using the Kendall rank correlation coefficient, 
# to obtain insights into bivariate relationships between individual species and environmental features. 
#
# Usage:
# - Install and load all required packages
# - Update the paths to input files and desired output location as needed.
# - Ensure all source data files are available and properly formatted before execution.
#
# NOTE: Ensure the imputed data prepared by "Code 3_Imputation" and the spatial and temporal data 
#       prepared by "Code4_Prepare environmental variable" are correctly formatted and accessible.
# -----------------------------------------------------------------------------------------


### Load packages ----

library(vegan)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(corrplot)

### Variables used for this script ----

# Sources
cam_data_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Imputed_data_rf.rds"
spatial_data_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Spatial_data.rds"
temporal_data_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Temporal_data.rds"

# Outputs
correlation_space_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Correlation_space.png"
correlation_time_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Correlation_time.png"
correlation_spt_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Correlation_space+time.png"
ordination_spt_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/SpT_Ord.png"

#1) Load and prepare data ----
Cam_data<- readRDS(cam_data_path)
Spat_dat_raw <- readRDS(spatial_data_path)
Temp_dat_raw <- readRDS(temporal_data_path)

### Scale the Spatial Data:

# Extract the first column
locationName <- Spat_dat_raw[, 1]

# Scale the remaining columns
scaled_data_S <- select(Spat_dat_raw, -c(locationName, Soil, Region)) %>% 
  scale() %>% 
  as.data.frame()

# Combine the first column with the scaled data
Spat_dat <- cbind(locationName, scaled_data_S)
            
### Scale the Temporal data: 

# Extract the first column
Year<- Temp_dat_raw[, 1]

# Scale the remaining columns
scaled_data_T <- Temp_dat_raw[, -1] %>% 
  scale() %>% 
  as.data.frame()

# Combine the first column with the scaled data
Temp_dat <- cbind(Year, scaled_data_T) %>% 
  #Filter out year 2019 and 2021, because they are missing in the camera trapping data. 
  filter(!(Year %in% c("2010","2018", "2020", "2022")))

#2) Quick basic analysis----- 

# Pivot the data to a long format
long_format_data <- Cam_data %>%
  select(-locationName, -Year) %>%  # Exclude non-species columns
  pivot_longer(
    -locationName_Year,             # Exclude the identifier column
    names_to = "Species",
    values_to = "Capture_Rate"
  )

# Calculate stats for each species
species_stats <- long_format_data %>%
  group_by(Species) %>%
  summarise(
    Max_Capture_Rate = max(Capture_Rate, na.rm = TRUE),
    Min_Capture_Rate = min(Capture_Rate, na.rm = TRUE),
    Mean_Capture_Rate = mean(Capture_Rate, na.rm = TRUE),
    Locations_with_Captures = sum(Capture_Rate > 0, na.rm = TRUE)
  )        

#3) Spatial Variation only: ----
#   a) DCA to test  for the appropriate analysis: ----
#Prepare Camera trap data 
Cam_dat_spat <- Cam_data %>%
  group_by(locationName) %>%
  summarise(
    Lowland_Paca = mean(Lowland_Paca, na.rm = TRUE),
    Central_American_Agouti = mean(Central_American_Agouti, na.rm = TRUE),
    Common_Opossum = mean(Common_Opossum, na.rm = TRUE),
    Coati = mean(Coati, na.rm = TRUE),
    Collared_Peccary = mean(Collared_Peccary, na.rm = TRUE),
    Tomes_Spiny_Rat = mean(Tomes_Spiny_Rat, na.rm = TRUE)
  ) %>% 
  ungroup()


#Detrended Correspondence Analysis (DCA)
dca_result_S <- decorana(Cam_dat_spat[,-1])

# Plot DCA
plot(dca_result_S)

# Extract lengths of the segments for each axis
segment_lengths_S <- dca_result_S$evals

# Determine the appropriate model based on the length criteria (<3 --> linear, >4 --> unimodal)
if (max(segment_lengths_S) < 3) {
  model <- "linear"
} else if (max(segment_lengths_S) > 4) {
  model <- "unimodal"
} else {
  model <- "unknown"
}

# Print the selected model
cat("Selected Model:", model, "\n")
#--> shows that the linear model (RDA) should be used!

#   b) linear analyis: RDA-------------------
        
all(unique(Cam_dat_spat$locationName) %in% unique(Spat_dat$locationName))

Spat_data <- inner_join(Cam_dat_spat, Spat_dat, by = "locationName") %>% 
  select(-c("Trees_Total")) %>% 
  rename(Agouti = "Central_American_Agouti",
         AST = AST_trees, 
         ATTA = ATTA_trees, 
         DIP = DIP_trees, 
         GUS = GUS_trees, 
         Forest_border = Forest_Border_Dist_m, 
         Slope = Slope_degree, 
         Forest_Area = Forest_Area_ha)

colnames(Spat_data)

#     Regular RDA: -----------
#Model with all variables: 
rda_S_all <- rda(Spat_data[, 2:ncol(Cam_dat_spat)] ~ ., data = Spat_data[, (ncol(Cam_dat_spat)+1):ncol(Spat_data)])

#Check for colinearity: 
vif.cca(rda_S_all)
corrplot(cor(Spat_data[, (ncol(Cam_dat_spat)+1):ncol(Spat_data)]), method = "color")

#Looks good (Total tree count removed)
###

#Check forward selection 
rda_S_blank <- rda(Spat_data[, 2:ncol(Cam_dat_spat)] ~ 1, data = Spat_data[, (ncol(Cam_dat_spat)+1):ncol(Spat_data)])
# Run the add1() command to see if any of the variables is likely to significantly affect the model
add1(rda_S_blank, scope = formula(rda_S_all), test = "permutation")

# Backward selection:
backward_spat_model <- step(rda_S_all, direction = "backward", test = "permutation")

#Recommanded model after backward selection: 
rda_S_back <- rda(formula = Spat_data[, 2:ncol(Cam_dat_spat)] ~ GUS_trees + Slope_degree + Forest_Area_ha +
                     West, data = Spat_data[, (ncol(Cam_dat_spat) + 1):ncol(Spat_data)])

#4) Analyze output: 
#See: https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
summary(rda_S_back)

#explanatory power??? Is that true? What about Inertia?
RsquareAdj(rda_S_all)$adj.r.squared

anova(rda_S_all) # model is not significant 
anova(rda_S_back)

#The selected environmental variables significantly explain 15.3% (p = 0.001) of the variation in species abundances.

anova.cca(rda_S_all, by = "terms")
#Atta trees, BCNM, Water distance, Slope, Soil Tbo, South, ATTA_fruit 
# almost: Forest_Border_Dist_m

plot(rda_S_all, type = "n", ylim = c(-1.8, 1.8), xlim = c(-2,2), scaling = 1)
text(rda_S_all, display = "bp", col = "blue")
text(rda_S_all, display = "species", col = "red", pch = "+", cex = 1)
points(rda_S_all, display = "sites", cex = 0.8)

#     Bray Curtis RDA--------

rda_S_all0<- capscale(Spat_data[, 2:ncol(Cam_dat_spat)] ~ ., data = Spat_data[, (ncol(Cam_dat_spat)+1):ncol(Spat_data)],
                  distance = "bray")
#Check for colinearity: 
vif.cca(rda_S_all)
corrplot(cor(Spat_data[, (ncol(Cam_dat_spat)+1):ncol(Spat_data)]), method = "color")


rda_S_all<- capscale(formula = Spat_data[, 2:ncol(Cam_dat_spat)] ~ AST + ATTA + DIP + GUS + 
           Forest_border + Water_distance + Slope +
           Forest_Area + Forest_age + SoilTb + SoilTbo + SoilTcm + SoilTcv 
           + East + BCI, data = Spat_data[, (ncol(Cam_dat_spat) +1):ncol(Spat_data)],
           distance = "bray") 


#Backward selection:
backward_spat_model <- step(rda_S_all, direction = "backward", test = "permutation")
rda_S_back0<- capscale(formula = Spat_data[, 2:ncol(Cam_dat_spat)] ~ AST + ATTA + DIP + GUS + Water_distance + 
                        Slope + Forest_Area + SoilTbo + SoilTcm + West, data
                      = Spat_data[, (ncol(Cam_dat_spat) + 1):ncol(Spat_data)], distance = "bray")
rda_S_back<-capscale(formula = Spat_data[, 2:ncol(Cam_dat_spat)] ~ AST + ATTA + DIP + Slope + Forest_Area +
                       SoilTbo + SoilTcm + East + BCI, data = Spat_data[, (ncol(Cam_dat_spat) + 1):ncol(Spat_data)], distance =
                       "bray")
#all variables:
anova.cca(rda_S_all, permutations = how(nperm=100000)) # model is significant 
anova.cca(rda_S_all, by= "terms", permutations = how(nperm=1000))
anova.cca(rda_S_all, permutations = how(nperm=100000), by = "axis")
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_S_all)

summary(rda_S_all)

#Backward selection model: 
anova.cca(rda_S_back, permutations = how(nperm=100000))
anova.cca(rda_S_back, by= "terms", permutations = how(nperm=100000))
anova.cca(rda_S_back, permutations = how(nperm=100000), by = "axis")
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_S_back)

summary(rda_S_back)


#   c) Correlation and Ordination graph ----
### CORRELATION ###
corrplot(cor(Spat_data[,2:ncol(Cam_dat_spat)], Spat_data[, (ncol(Cam_dat_spat)+1):ncol(Spat_data)]), method = "color")
# Calculate correlation matrix
correlation_matrix <- cor(Spat_data[,2:ncol(Cam_dat_spat)], Spat_data[, (ncol(Cam_dat_spat)+1):ncol(Spat_data)], method = "kendall")

# Plot correlation matrix with black text
corrplot(correlation_matrix, method = "color", tl.col = "black", addCoef.col = "black")
png(correlation_space_png, 
    width = 1000, height = 500)  # Adjust width and height as needed
dev.off()
        
####ORDINATION###

## extract scores - these are coordinates in the RDA space
sc_si <- scores(rda_S_back, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(rda_S_back, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(rda_S_back, display="bp", choices=c(1, 2), scaling=1)

# Filter sc_bp for significant variables
significant_vars <- c("AST", "SoilTbo", "SoilTcm", "Slope")
sc_bp_filtered <- sc_bp[significant_vars,]

# Custom triplot function
plotRDA <- function(rda_obj, sc_si, sc_sp, sc_bp_filtered, 
                    xlims = c(-1, 1), ylims = c(-0.5, 0.6), 
                    addSpeciesText = TRUE, addArrowText = TRUE, 
                    TextSize=0.8) {
  perc <- round(100 * (summary(rda_obj)$cont$importance[2, 1:2]), 2)
  
  plot(rda_obj,
       scaling = 2, # set scaling type 
       type = "none", # this excludes the plotting of any points from the results
       frame = TRUE,
       xlim = xlims, 
       ylim = ylims,
       xlab = paste0("RDA1 (", perc[1], "%)"), 
       ylab = paste0("RDA2 (", perc[2], "%)")
  )
  
  points(sc_si, 
         pch = 21, # set shape (here, circle with a fill color)
         col = "white", # outline color
         bg = "#C3D1F7", # fill color
         cex = 1.2) # size
  
  points(sc_sp, 
         pch = 22, # set shape (here, square with a fill color)
         col = "white",
         bg = "#f2bd33", 
         cex = 1.2)
  
  if (addSpeciesText) {
    text(sc_sp + c(0.03, 0), # adjust text coordinates to avoid overlap with points 
         labels = rownames(sc_sp), 
         col = "grey40", 
         font = 1, # normal font
         cex = TextSize)
  }
  
  arrows(0, 0, 
         sc_bp_filtered[,1], sc_bp_filtered[,2], 
         col = "red", 
         lwd = 1, 
         length = 0.1)
  
  if (addArrowText) {
    text(x = sc_bp_filtered[,1], 
         y = sc_bp_filtered[,2] + 0.02, 
         labels = rownames(sc_bp_filtered), 
         col = "red", 
         cex = TextSize, 
         font = 1)
  }
}


Space_Zoom<- plotRDA(rda_obj = rda_S_back, sc_si = sc_si, sc_sp = sc_sp, 
                     sc_bp_filtered = sc_bp_filtered, 
                     xlims = c(-1,1), ylims = c(-0.4, 0.6),
                     TextSize = 0.8)

All <- plotRDA(rda_obj = rda_S_back, sc_si = sc_si, sc_sp = sc_sp, sc_bp_filtered = sc_bp_filtered, 
               xlims = c(-1, 3), ylims = c(-0.5, 3.2), TextSize = 0.4, 
               addSpeciesText = FALSE, addArrowText = FALSE)

      
#4) Temporal variation only: ----
#   a) DCA to test  for the appropriate analysis: ----

#Prepare Camera trap data 
Cam_dat_temp <- Cam_data %>%
  group_by(Year) %>%
  summarise(
    Lowland_Paca = mean(Lowland_Paca, na.rm = TRUE),
    Central_American_Agouti = mean(Central_American_Agouti, na.rm = TRUE),
    Common_Opossum = mean(Common_Opossum, na.rm = TRUE),
    Coati = mean(Coati, na.rm = TRUE),
    Collared_Peccary = mean(Collared_Peccary, na.rm = TRUE),
    Tomes_Spiny_Rat = mean(Tomes_Spiny_Rat, na.rm = TRUE)
  ) %>% 
  mutate(Year = as.character(Year)) %>% 
  ungroup()

#Detrended Correspondence Analysis (DCA)
dca_result_T <- decorana(Cam_dat_temp[,-1])

# Plot DCA
plot(dca_result_T)

# Extract lengths of the segments for each axis
segment_lengths_T <- dca_result_T$evals

# Determine the appropriate model based on the length criteria (<3 --> linear, >4 --> unimodal)
if (max(segment_lengths_T) < 3) {
  model <- "linear"
} else if (max(segment_lengths_T) > 4) {
  model <- "unimodal"
} else {
  model <- "unknown"
}

# Print the selected model
cat("Selected Model:", model, "\n")
#--> linear model (RDA) should be used!

#   b) linear analyis: RDA -------------------

all(unique(Cam_dat_temp$Year) %in% unique(Temp_dat$Year))

Temp_data <- inner_join(Cam_dat_temp, Temp_dat, by = "Year") %>% 
  select(-c(ATTA_fruit, Total_count))


#       Regular RDA: ------

# Model with all variables:
rda_T_all <- rda(Temp_data[, 2:ncol(Cam_dat_temp)] ~ ., data = Temp_data[, (ncol(Cam_dat_temp)+1):ncol(Temp_data)])

#Check for colinearity: 
vif.cca(rda_T_all)
corrplot(cor(Temp_data[, (ncol(Cam_dat_temp)+1):ncol(Temp_data)]), method = "color")

# Forward selection:
rda_T_blank <- rda(Temp_data[, 2:ncol(Cam_dat_temp)] ~ 1, data = Temp_data[, (ncol(Cam_dat_temp)+1):ncol(Temp_data)])
# Run the add1() command to see if any of the variables is likely to significantly affect the model
add1(rda_T_blank, scope = formula(rda_T_all), test = "permutation")

# Backward selection:
backward_temp_model <- step(rda_T_all, direction = "backward", test = "permutation")

# Recommanded model after backward selection: 
rda_T_back <- rda(formula = Temp_data[, 2:ncol(Cam_dat_temp)] ~ Rain_total + Rain_change + AST_fruit + GUS_fruit +
                     Fruit_change, data = Temp_data[, (ncol(Cam_dat_temp) + 1):ncol(Temp_data)])

#Analyze output: 
#See: https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html

summary(rda_T_back)
RsquareAdj(rda_T_back)$adj.r.squared
anova(rda_T_back) # model is significant 

#The selected environmental variables significantly explain 15.3% (p = 0.001) of the variation in species abundances.
anova.cca(rda_T_back, by = "terms")

plot(rda_T_all, type = "n", ylim = c(-1.8, 1.8), xlim = c(-2,2), scaling = 1)
text(rda_T_all, display = "bp", col = "blue")
text(rda_T_all, display = "species", col = "red", pch = "+", cex = 1)
points(rda_T_all, display = "sites", cex = 0.8)

#       Bray Curtis RDA: --------

rda_T_all0<- capscale(Temp_data[, 2:ncol(Cam_dat_temp)] ~ ., data = Temp_data[, (ncol(Cam_dat_temp)+1):ncol(Temp_data)],
                     distance = "bray")
#Check for colinearity: 
vif.cca(rda_T_all)
corrplot(cor(Temp_data[, (ncol(Cam_dat_temp)+1):ncol(Temp_data)]), method = "color")
#high correlation between fruit change and Dip seeds, rain_change and Gus fruits 

rda_T_all<- capscale(Temp_data[, 2:ncol(Cam_dat_temp)] ~ Rain_total+Rain_change+AST_fruit+GUS_fruit+DIP_fruit, data = Temp_data[, (ncol(Cam_dat_temp)+1):ncol(Temp_data)],
                     distance = "bray")

#Backward selection:
backward_temp_model <- step(rda_T_all, direction = "backward", test = "permutation")
#-> no variables found. 

#all variables:
anova.cca(rda_T_all, permutations = how(nperm=100000)) # model is not significant 
anova.cca(rda_T_all, by= "terms", permutations = how(nperm=100000))
anova.cca(rda_T_all, permutations = how(nperm=100000), by = "axis")
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_T_all)

summary(rda_T_all)

#       Correlation graph ----

#Correlation plot
correlation_matrix_time <- cor(Temp_data[,2:ncol(Cam_dat_temp)], Temp_data[, (ncol(Cam_dat_temp)+1):ncol(Temp_data)], method = "kendall")

# Plot correlation matrix with black text
corrplot(correlation_matrix_time, method = "color", tl.col = "black", addCoef.col = "black")

#Save the plot: 
png(correlation_time_png, 
    width = 500, height = 500)  # Adjust width and height as needed
corrplot(correlation_matrix_time, method = "color", tl.col = "black", addCoef.col = "black")
dev.off()

#5) Spatio Temporal Variation: ----

### Prepare Data

#1) Prepare environmental data: Make one big table with all environmental variables (Spat + Temp)

#Create 10 replicates of the Spat_dat table
replicates <- replicate(nrow(Temp_dat), Spat_dat, simplify = FALSE)#9 is the number of rows of the Temp_data 

#Bind replicates together
replicate_spat_data <- bind_rows(replicates, .id = "Replicate") %>% 
  #Create a new column "Year" with values from Temp_dat
  mutate(Year = rep(Temp_dat$Year, each = nrow(Spat_dat))) %>% 
  # Arrange the data based on the Replicate and Year columns
  arrange(Replicate, Year)

#Left join Temp_dat with combined_data using "Year" as the key
Env_data<-    left_join(replicate_spat_data, Temp_dat, by = "Year") %>% 
  #Create a new column as identifier for each row 
  mutate(locationName_Year = paste(locationName, Year, sep = "_")) %>%
  select(locationName_Year, everything()) %>% 
  select(-c("Year", "locationName","Replicate"))

#2) Camera trap data

Cam_dat_ST<- Cam_data %>% 
  select(-c(locationName, Year))

#3) Combine Cam and Env data: 

all(unique(Cam_dat_ST$locationName_Year) %in% unique(Env_data$locationName_Year))

ST_data_0 <- inner_join(Cam_dat_ST, Env_data, by = "locationName_Year") %>% 
  select(-c(ATTA_fruit))

#   a) DCA to test  for the appropriate analysis: ----

#Detrended Correspondence Analysis (DCA)
dca_result <- decorana(Cam_dat_ST[,-1])

# Plot DCA
plot(dca_result)

# Extract lengths of the segments for each axis
segment_lengths <- dca_result$evals

# Determine the appropriate model based on the length criteria (<3 --> linear, >4 --> unimodal)
if (max(segment_lengths) < 3) {
  model <- "linear"
} else if (max(segment_lengths) > 4) {
  model <- "unimodal"
} else {
  model <- "unknown"
}

# Print the selected model
cat("Selected Model:", model, "\n")
#--> linear model (RDA) should be used!


#   b) linear analysis: RDA ----
#       Regular RDA: --------
#Model with all variables: 
rda_ST_all_0 <- rda(ST_data[, 2:ncol(Cam_dat_ST)] ~ ., data = ST_data[, (ncol(Cam_dat_ST)+1):ncol(ST_data)])

#Check for colinearity: 
vif.cca(rda_ST_all_0)
corrplot(cor(ST_data[, (ncol(Cam_dat_ST)+1):ncol(ST_data)]), method = "color")

#Adjusted data: 
ST_data_select <- ST_data_0 %>% 
  select(-c(Trees_Total, DIP_fruit, Total_count))

#new model: 
rda_ST_all <- rda(ST_data_select[, 2:ncol(Cam_dat_ST)] ~ ., data = ST_data_select[, (ncol(Cam_dat_ST)+1):ncol(ST_data_select)])
vif.cca(rda_ST_all)
corrplot(cor(ST_data_select[, (ncol(Cam_dat_ST)+1):ncol(ST_data_select)]), method = "color")


#Forward selection:
rda_ST_blank <- rda(ST_data[, 2:ncol(Cam_dat_ST)] ~ 1, data = ST_data[, (ncol(Cam_dat_ST)+1):ncol(ST_data)])
# Run the add1() command to see if any of the variables is likely to significantly affect the model
add1(rda_ST_blank, scope = formula(rda_ST_all), test = "permutation")

# Backward selection:
backward_temp_model <- step(rda_ST_all, direction = "backward", test = "permutation")

#Recommanded model after backward selection: 
rda_ST_back <- rda(formula = ST_data_select[, 2:ncol(Cam_dat_ST)] ~ AST_trees + ATTA_trees + GUS_trees + Water_distance +
                  Slope_degree + Forest_Area_ha + SoilTb + SoilTcm + SoilTcv + West + Rain_total + Rain_change + ATTA_fruit +
                  Fruit_change, data = ST_data_select[, (ncol(Cam_dat_ST) + 1):ncol(ST_data_select)])

#Analyze output: 
#See: https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
    
summary(rda_ST_all)

#explanatory power: 
RsquareAdj(rda_ST_back)$adj.r.squared

anova(rda_ST_all) # model is significant 

#The selected environmental variables significantly explain 15.3% (p = 0.001) of the variation in species abundances.

anova(rda_ST_all, by = "terms")
#Atta trees, BCNM, Water distance, Slope, Soil Tbo, South, ATTA_fruit 
# almost: Forest_Border_Dist_m

plot(rda_ST_back, type = "n", ylim = c(-1.8, 1.8), xlim = c(-1,1), scaling = 1)
points(rda_ST_back, display = "sites",col="lightblue", cex = 0.8)

text(rda_ST_back, display = "bp", col = "blue")
text(rda_ST_back, display = "species", col = "red", pch = "+", cex = 1)


#       Bray Curtis RDA--------

rda_ST_all0<- capscale(ST_data_0[, 2:ncol(Cam_dat_ST)] ~ ., data = ST_data_0[, (ncol(Cam_dat_ST)+1):ncol(ST_data_0)],
                      distance = "bray")

#Check for colinearity: 
vif.cca(rda_ST_all) 
corrplot(cor(ST_data_0[, (ncol(Cam_dat_ST)+1):ncol(ST_data_0)]), method = "color")
#-> high correlation between fruit change and Dip seeds. and rain_change and Gus fruits 

ST_data <- ST_data_0 %>% 
  select(-c(Trees_Total, Fruit_change, Total_count, SoilTue, West))

rda_ST_all<- capscale(ST_data[, 2:ncol(Cam_dat_ST)] ~ ., data = ST_data[, (ncol(Cam_dat_ST)+1):ncol(ST_data)],
                      distance = "bray")

#Backward selection:
backward_ST_model <- step(rda_ST_all, direction = "backward", test = "permutation")
rda_ST_back <- capscale(formula = ST_data[, 2:ncol(Cam_dat_ST)] ~ AST_trees + ATTA_trees + DIP_trees + GUS_trees + Forest_Border_Dist_m +
                          Slope_degree + Forest_Area_ha + SoilTbo + SoilTcm + East + BCI, data = ST_data[, (ncol(Cam_dat_ST) + 1):ncol(ST_data)], distance = "bray")

#all variables:
anova.cca(rda_ST_all, permutations = how(nperm=10000)) #significant
anova.cca(rda_ST_all, by= "terms", permutations = how(nperm=1000))
anova.cca(rda_ST_all, permutations = how(nperm=1000), by = "axis")
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_ST_all)

summary(rda_ST_all)

#selected variables: 
anova.cca(rda_ST_back, permutations = how(nperm=1000)) # model is not significant 
anova.cca(rda_ST_back, by= "terms", permutations = how(nperm=1000))
anova.cca(rda_ST_back, permutations = how(nperm=1000), by = "axis")
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(rda_ST_back)

summary(rda_ST_back)

        
#   c) Correlation and Ordination Plots ----
#Correlation plot
correlation_matrix_time <- cor(ST_data_0[,2:ncol(Cam_dat_ST)], ST_data_0[, (ncol(Cam_dat_ST)+1):ncol(ST_data_0)], method = "kendall")

# Plot correlation matrix with black text
corrplot(correlation_matrix_time, method = "color", tl.col = "black", addCoef.col = "black")

#Save the plot: 
png(correlation_spt_png, 
    width = 1200, height = 500)  # Adjust width and height as needed
corrplot(correlation_matrix_time, method = "color", tl.col = "black", addCoef.col = "black")
dev.off()
          
#Ordination plot
          
## extract scores - these are coordinates in the RDA space
sc_si <- scores(rda_ST_all, display="sites", choices=c(1,2), scaling=2)
sc_sp <- scores(rda_ST_all, display="species", choices=c(1,2), scaling=2)
sc_bp <- scores(rda_ST_all, display="bp", choices=c(1, 2), scaling=2)

#Which variables were significant?
anova_result <- anova.cca(rda_ST_all, by = "terms", permutations = how(nperm = 1000))

# Convert the anova_result to a data frame
anova_df <- as.data.frame(anova_result)

# Filter for significant variables, e.g., where p-value < 0.05
significant_vars <- anova_df$"Pr(>F)" < 0.05 

# Extract names of significant variables
significant_var_names <- rownames(anova_df)[significant_vars] 
significant_var_names <- as.vector(na.omit(significant_var_names))

# Custom triplot code!
plotRDA <- function(rda_obj, sc_si, sc_sp, sc_bp, multiplier, significant_var_names, xlims = c(-1, 1), ylims = c(-0.5, 0.6), TextSize=0.8, addSpeciesText = TRUE, addArrowText = TRUE) {
  perc <- round(100 * (summary(rda_obj)$cont$importance[2, 1:2]), 2)
  
  plot(rda_obj,
       scaling = 2, # set scaling type 
       type = "none", # this excludes the plotting of any points from the results
       frame = TRUE,
       xlim = xlims, 
       ylim = ylims,
       xlab = paste0("RDA1 (", perc[1], "%)"), 
       ylab = paste0("RDA2 (", perc[2], "%)")
  )
  
  points(sc_si, 
         pch = 21, # set shape (here, circle with a fill colour)
         col = "white", # outline colour
         bg = "#C3D1F7", # fill colour
         cex = 1.2) # size
  
  points(sc_sp, 
         pch = 22, # set shape (here, square with a fill colour)
         col = "white",
         bg = "#f2bd33", 
         cex = 1.2)
  
  # Filter sc_bp scores for significant variables
  sc_bp_filtered <- sc_bp[significant_var_names, , drop = FALSE] # Ensure it doesn't drop to a vector if only one var is significant
  
  arrows(0, 0, 
        sc_bp_filtered[,1]*multiplier, 
        sc_bp_filtered[,2]*multiplier, 
         col = "red", 
         lwd = 1, 
         length = 0.1)
  
  if (addArrowText) {
    text(x = sc_bp_filtered[,1]*multiplier*1.1, 
         y = sc_bp_filtered[,2]*multiplier*1.1, 
         labels = rownames(sc_bp_filtered), 
         col = "red", 
         cex = TextSize, 
         font = 1)
  }
}


# Save the plot: Open a PNG graphics device
png(filename = ordination_spt_png, 
    width = 8, height = 7, units = "in", res = 300)

plotRDA(rda_obj = rda_ST_all, sc_si = sc_si, sc_sp = sc_sp,sc_bp =sc_bp, 
                     multiplier = 6,significant_var_names = significant_var_names,
                     xlims = c(-3,4.5), ylims = c(-4, 2.5),
                     TextSize = 0.8)

# Close the device
dev.off()
          
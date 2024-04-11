# -----------------------------------------------------------------------------------------
# Script Name: Code 6_BD_functions_applied - Bray Curtis
# Author: Nadja Warth
# 
# Description:
# This script conains the Analysis on how the variation in species composition is distributed in 
# space and time, and which role deterministic and stochastic processes play.
# It is integral to addressing key research questions concerning the spatio-temporal distribution 
# of species composition dissimilarity and the relative influence of deterministic versus stochastic 
# assembly processes on community structure. The script utilizes the Bray-Curtis index to quantify 
# beta-diversity, indicating how dissimilarity in species composition varies across spatial and temporal 
# dimensions. Following methodologies outlined by Khattar et al. (2021) and Legendre (2014).
# 
# Sections:
#
# 1) Defines functions from Legendre and Khattar to calculate pairwise dissimilarities and partition
#    the total dissimilarity (Dtotal) into spatial (Dspace), temporal (Dtime), and spatiotemporal (DSpT)
#    components.
# 2) Converts data into the appropriate format for analysis.
# 3) Calculates and plots average beta-diversity across the different dimensions.
# 4) Tests for significant differences in dissimilarity values among the dimensions and plots these differences.
# 5) Assesses the role of stochasticity and determinism in community assembly: 
#    5.a) Generates multiple random community matrices.
#    5.b) Calculates mean and standard deviation for each row of the null dissimilarities.
#    5.c) Calculates Standardized Effect Sizes (SES) to compare observed and null model dissimilarities.
#    5.d) Conducts T-tests to assess significance.
#    5.e) Produces boxplots of SES values.
#    5.f) Graphs to compare means of observed dissimilarities against those of the null models.
# 6) Tests the stochastic and deterministic effects across different years.
# 7) Tests the stochastic and deterministic effects across different locations.
#
# Usage:
# - Install and load all required packages
# - Update the paths to input files and desired output location as needed.
# - Ensure all source data files are available and properly formatted before execution.
#
#  NOTE: Ensure the imputed data prepared by "Code 3_Imputation" are correctly formatted and accessible.
# -----------------------------------------------------------------------------------------


### Load Packages ----
library(tidyverse)
library(vegan)
library(ggplot2)

library(car)
library(dunn.test)

### Variables needed for this script ----

# Sources
cam_data_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Imputed_data_rf.rds"

# Outputs
average_BD_plot_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/average_BD_plot_01.02.png"
dissimilarity_matrix_path <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/dissimilarity_matrix_BrayC_13.03.rds"
dissimilarity_index_plot_png <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/D_index_plot_13.02.png"
null_dissimilarities_table_raw_normalized <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/data/processed/Null_dissimilarities_table_raw_normalized_BrayC_09.03.rds"
mean_SES_plot_horizontal <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Mean_SES_plot_horizontal_13-02.png"
mean_SES_plot <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Mean_SES_plot.png"
mean_SES_plot_years <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Mean_SES_plot_years.png"
mean_SES_plot_location <- "C:/Users/Nadja/OneDrive - Wageningen University & Research/Project Nadja Warth/figures/Graphs/Results/Mean_SES_plot_loc.png"

#1) Define functions from literature ----
#   a)Beta Div Comp function (Legendre 2014) -----------
beta.div.comp <- function(mat, coef="S", quant=TRUE, save.abc=FALSE)
#
# Description --
#
# Podani-family and Baselga-family decompositions of the Jaccard and Sørensen
# dissimilarity coefficients into replacement and richness/abundance difference
# components, for species presence-absence or abundance data, as described
# in Legendre (2014).
#
# Usage --
#
# beta.div.comp(mat, coef="J", quant=FALSE, save.abc=FALSE)
#
# Arguments --
#
# mat : Data in matrix or data.frame form.
# coef : Family of coefficients to be computed --
# "S" or "Sorensen": Podani family, Sørensen-based indices
# "J" or "Jaccard" : Podani family, Jaccard-based indices
# "BS" : Baselga family, Sørensen-based indices
# "BJ" : Baselga family, Jaccard-based indices
# "N" : Podani & Schmera (2011) relativized nestedness index.
# The quantitative form in Sørensen family is the percentage difference.
# The quantitative form in the Jaccard family is the Ruzicka index.
#
# quant=TRUE : Compute the quantitative form of the indices and D.
# =FALSE: Compute the presence-absence form of the coefficients.
# save.abc=TRUE : Save the matrices of parameters a, b and c used in the
# presence-absence calculations.
#
# Details --
#
# For species presence-absence data, the distance coefficients are
# Jaccard=(b+c)/(a+b+c) and Sørensen=(b+c)/(2*a+b+c) with usual abc notation.
#
# For species abundance data, the distance coefficients are
# the Ruzicka index = (B+C)/(A+B+C) and Odum's percentage difference
# (incorrectly called Bray-Curtis) = (B+C)/(2A+B+C), where
# A = sum of the intersections (or minima) of species abundances at two sites,
# B = sum at site 1 minus A, C = sum at site 2 minus A.
#
# The binary (quant=FALSE) and quantitative (quant=TRUE) forms of the S and
# J indices return the same values when computed for presence-absence data.
#

# Value --
#
# repl : Replacement matrix, class = 'dist'.
# rich : Richness/abundance difference or nestedness matrix, class = 'dist'.
# With options "BJ", "BS" and "N", 'rich' contains nestedness indices.
# With option "N", the 'repl' and 'rich' values do not add up to 'D'.
# D : Dissimilarity matrix, class = 'dist'.
# part : Beta diversity partitioning --
# 1. Total beta div. = sum(D.ij)/(n*(n-1)) (Legendre & De Cáceres 2013)
# 2. Total replacement diversity
# 3. Total richness/abundance difference diversity (or nestedness)
# 4. Total replacement div./Total beta div.
# 5. Total richness/abundance diff. div. (or nestedness)/Total beta div.
# Note : Name of the dissimilarity coefficient.
#
# References --
#
# Baselga, A. (2010) Partitioning the turnover and nestedness components of beta
# diversity. Global Ecology and Biogeography, 19, 134–143.
#
# Baselga, A. (2012) The relationship between species replacement, dissimilarity
# derived from nestedness, and nestedness. Global Ecology and Biogeography, 21,
# 1223–1232.
#
# Baselga, A. (2013) Separating the two components of abundance-based
# dissimilarity: balanced changes in abundance vs. abundance gradients. Methods
# in Ecology and Evolution, 4, 552–557.
#
# Carvalho, J.C., Cardoso, P., Borges, P.A.V., Schmera, D. & Podani, J. (2013)
# Measuring fractions of beta diversity and their relationships to nestedness:
# a theoretical and empirical comparison of novel approaches. Oikos, 122,
# 825–834.
#
# Legendre, P. (2014) Interpreting the replacement and richness difference
# components of beta diversity. Global Ecology and Biogeography, 23, xxx–xxx.
#
# Podani, J., Ricotta, C. & Schmera, D. (2013) A general framework for analyzing
# beta diversity, nestedness and related community-level phenomena based on
# abundance data. Ecological Complexity, 15, 52-61.
#
# Podani, J. & Schmera, D. (2011) A new conceptual and methodological framework
# for exploring and explaining pattern in presence-absence data. Oikos, 120,
# 1625–1638.
#
# License: GPL-2
# Author:: Pierre Legendre
{
  coef <- pmatch(coef, c("S", "J", "BS", "BJ", "N"))
  if(coef==5 & quant) stop("coef='N' and quant=TRUE: combination not programmed")
  mat <- as.matrix(mat)
  n <- nrow(mat)
  if(is.null(rownames(mat))) noms <- paste("Site",1:n,sep="")
  else noms <- rownames(mat)
  #
  if(!quant) { # Binary data provided, or make the data binary
    if(coef==1) form="Podani family, Sorensen"
    if(coef==2) form="Podani family, Jaccard"
    if(coef==3) form="Baselga family, Sorensen"
    if(coef==4) form="Baselga family, Jaccard"
    if(coef==5) form="Podani & Schmera (2011) relativized nestedness"
    mat.b <- ifelse(mat>0, 1, 0)
    a <- mat.b %*% t(mat.b)
    b <- mat.b %*% (1 - t(mat.b))
    c <- (1 - mat.b) %*% t(mat.b)
    min.bc <- pmin(b,c)
    #
    if(coef==1 || coef==2) {
      repl <- 2*min.bc # replacement, turnover, beta-3
      rich <- abs(b-c) # nestedness, richness diff., beta-rich
      #
      # Add the denominators
      if(coef==1) { # Sørensen-based components
        repl <- repl/(2*a+b+c)
        rich <- rich/(2*a+b+c)
        D <- (b+c)/(2*a+b+c)
      } else if(coef==2) { # Jaccard-based components
        repl <- repl/(a+b+c)
        rich <- rich/(a+b+c)
        D <- (b+c)/(a+b+c)
      }
    } else if(coef==3) { # Baselga 2010 components based on Sørensen
      D <- (b+c)/(2*a+b+c) # Sørensen dissimilarity
      repl <- min.bc/(a+min.bc) # replacement, turnover
      rich <- D-repl # nestedness-resultant dissimilarity
    } else if(coef==4) { # Baselga 2012 components based on Jaccard
      D <- (b+c)/(a+b+c) # Jaccard dissimilarity
      repl <- 2*min.bc/(a+2*min.bc) # replacement, turnover
      rich <- D-repl # nestedness-resultant dissimilarity
    } else if(coef==5) { # rich = Podani N = nestdness based on Jaccard
      repl <- 2*min.bc/(a+b+c)
      D <- (b+c)/(a+b+c)
      rich <- matrix(0,n,n)
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          aa = a[i,j]; bb = b[i,j]; cc = c[i,j]
          if(a[i,j] == 0) rich[i,j] <- 0
          else rich[i,j] <- (aa + abs(bb-cc))/(aa+bb+cc)
        }
      }
    }
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    D <- as.dist(D)
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    total.div <- sum(D)/(n*(n-1))
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    if(save.abc) {
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form,
                  a=as.dist(a), b=as.dist(b), c=as.dist(c))
    } else {
      res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
    }
    #
  } else { # Quantitative data
    # Calculations based on individuals.within.species
    if(coef==1) form<-"Podani family, percentage difference"
    if(coef==2) form<-"Podani family, Ruzicka"
    if(coef==3) form<-"Baselga family, percentage difference"
    if(coef==4) form<-"Baselga family, Ruzicka"
    # Baselga (2013) notation:
    # A = W = sum of minima in among-site comparisons
    # B = site.1 sum - W = K.1 - W
    # C = site.2 sum - W = K.2 - W
    K <- vector("numeric", n) # site (row) sums
    W <- matrix(0,n,n)
    repl <- matrix(0,n,n)
    rich <- matrix(0,n,n)
    D <- matrix(0,n,n)
    rownames(repl) <- rownames(rich) <- rownames(D) <- noms
    K <- apply(mat,1,sum) # Row sums
    for(i in 2:n) for(j in 1:(i-1)) W[i,j] <- sum(pmin(mat[i,], mat[j,]))
    #
    # Quantitative extensions of the S and J decompositions
    for(i in 2:n) {
      for(j in 1:(i-1)) {
        repl[i,j] <- 2*(min(K[i],K[j])-W[i,j]) # 2*min(B,C)
        rich[i,j] <- abs(K[i]-K[j]) # abs(B-C)
      }
    }
    #
    # Add the denominators
    if(coef==1) { # Sørensen-based (% difference) components
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j]) # 2min(B,C)/(2A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j]) # abs(B-C)/(2A+B+C)
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j]) # (B+C)/(2A+B+C)
        }
      }
    } else if(coef==2) { # Jaccard-based (Ruzicka) components
      for(i in 2:n) {
        for(j in 1:(i-1)) { # Baselga 2013 notation:
          repl[i,j] <- repl[i,j]/(K[i]+K[j]-W[i,j]) # 2min(B,C)/(A+B+C)
          rich[i,j] <- rich[i,j]/(K[i]+K[j]-W[i,j]) # abs(B-C)/(A+B+C)
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j]<-(K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j]) # (B+C)/(A+B+C)
        }
      }
    }
    #
    # Baselga (2013): quantitative extensions of the Baselga (2010) indices
    if(coef==3) { # Baselga (2013) indices decomposing percentage difference
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <- (min(K[i],K[j])-W[i,j])/min(K[i],K[j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/((K[i]+K[j])*min(K[i],K[j]))
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j])
        }
      }
    }
    if(coef==4) { # Decomposing Ruzicka in the spirit of Baselga 2013
      for(i in 2:n) {
        for(j in 1:(i-1)) {
          repl[i,j] <-
            2*(min(K[i],K[j])-W[i,j])/(2*min(K[i],K[j])-W[i,j])
          rich[i,j] <- abs(K[i]-K[j])*W[i,j]/
            ((K[i]+K[j]-W[i,j])*(2*min(K[i],K[j])-W[i,j]))
          # cat(K[i], K[j], W[i,j],"\n")
          D[i,j] <- (K[i]+K[j]-2*W[i,j])/(K[i]+K[j]-W[i,j])
        }
      }
    }
    #
    repl <- as.dist(repl)
    rich <- as.dist(rich)
    D <- as.dist(D)
    repl.div <- sum(repl)/(n*(n-1))
    rich.div <- sum(rich)/(n*(n-1))
    total.div <- sum(D)/(n*(n-1))
    part <- c(total.div,repl.div,rich.div,repl.div/total.div,rich.div/total.div)
    #
    res <- list(repl=repl, rich=rich, D=D, part=part, Note=form)
  }
  res
}

#   b)BD_partition (Khatter 2021)------------------------

BD_partition<-function(x,Component){
  
  if(Component=="Total"){
    Diss_matrix<-beta.div.comp(x,coef="S",quant=T)$D
    Diss_modify_total <- data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_modify_total) <- c("c1", "c2", "distance")
  }
  
  if(Component=="Turnover"){
    Diss_matrix<-beta.div.comp(x,coef="S",quant=T)$repl
    Diss_modify_total <- data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_modify_total) <- c("c1", "c2", "distance")
  }
  
  if(Component=="Nestedness"){
    Diss_matrix<-beta.div.comp(x,coef="S",quant=T)$rich
    Diss_modify_total <- data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_modify_total) <- c("c1", "c2", "distance")
  }
  # Naming the dimension of each entry in Diss_matrix
  Dim<-vector()
  
  for(i in 1:nrow(Diss_modify_total)){
    s1.letter = strsplit(as.character(Diss_modify_total[i,1]), split = "_")[[1]]
    s2.letter = strsplit(as.character(Diss_modify_total[i,2]), split= "_") [[1]]
    if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]){
      Dim[i]<- "space"
    }else{
      if(s1.letter[1] != s2.letter[1] & s1.letter[2] == s2.letter[2]){
        Dim[i]<- "time"
      }else{
        Dim[i]<- "spt"
      }}
  } 
  
  Diss_modify_total$Dim<-Dim
  
  BDtotal<-sum(Diss_modify_total[,3])/(nrow(x)*(nrow(x)-1))
  
  Total_entries_per_dimension<-table(Dim)# Total entries in D representing each dimension
  
  BDspace<-sum(subset(Diss_modify_total,subset=Dim=="space")[,3])/(nrow(x)*(nrow(x)-1))
  BDtime<-sum(subset(Diss_modify_total,subset=Dim=="time")[,3])/(nrow(x)*(nrow(x)-1))
  BDspt<-sum(subset(Diss_modify_total,subset=Dim=="spt")[,3])/(nrow(x)*(nrow(x)-1))
  
  
  BD_partition<-c(BDtotal,BDspace,BDspt,BDtime)
  names(BD_partition)<-c("Total","Space","SPT","Time")
  
  Avg_Cont_to_BDtotal<-c(BDspace,BDspt,BDtime)/Total_entries_per_dimension
  names(Avg_Cont_to_BDtotal)<-c("Avg.Space","Avg.SPT","Avg.Time")
  
  Result<-list(BD_partition,Avg_Cont_to_BDtotal,Total_entries_per_dimension)
  names(Result)<-c("BDpartition","Avg. Cont to BDtotal","Total entries in D per dimension")
  return(Result)
  
}

# Variation of BD_partition function by Khattar (2021), but returns Diss_modify_total needed for own Analysis, instead of summed up values 
BD_partition_matrix <- function(x,Component){
  
  if(Component=="Total"){
    Diss_matrix<-beta.div.comp(x,coef="S",quant=T)$D
    Diss_modify_total <- data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_modify_total) <- c("c1", "c2", "distance")
  }
  
  if(Component=="Turnover"){
    Diss_matrix<-beta.div.comp(x,coef="S",quant=T)$repl
    Diss_modify_total <- data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_modify_total) <- c("c1", "c2", "distance")
  }
  
  if(Component=="Nestedness"){
    Diss_matrix<-beta.div.comp(x,coef="S",quant=T)$rich
    Diss_modify_total <- data.frame(t(combn(as.character(rownames(x)),2)), as.numeric(Diss_matrix))
    names(Diss_modify_total) <- c("c1", "c2", "distance")
  }
  # Naming the dimension of each entry in Diss_matrix
  Dim<-vector()
  
  for(i in 1:nrow(Diss_modify_total)){
    s1.letter = strsplit(as.character(Diss_modify_total[i,1]), split = "_")[[1]]
    s2.letter = strsplit(as.character(Diss_modify_total[i,2]), split= "_") [[1]]
    if(s1.letter[1] == s2.letter[1] & s1.letter[2] != s2.letter[2]){
      Dim[i]<- "space"
    }else{
      if(s1.letter[1] != s2.letter[1] & s1.letter[2] == s2.letter[2]){
        Dim[i]<- "time"
      }else{
        Dim[i]<- "spt"
      }}
  } 
  
  Diss_modify_total$Dim<-Dim
  
  BDtotal<-sum(Diss_modify_total[,3])/(nrow(x)*(nrow(x)-1))
  
  Total_entries_per_dimension<-table(Dim)# Total entries in D representing each dimension
  
  BDspace<-sum(subset(Diss_modify_total,subset=Dim=="space")[,3])/(nrow(x)*(nrow(x)-1))
  BDtime<-sum(subset(Diss_modify_total,subset=Dim=="time")[,3])/(nrow(x)*(nrow(x)-1))
  BDspt<-sum(subset(Diss_modify_total,subset=Dim=="spt")[,3])/(nrow(x)*(nrow(x)-1))
  
  
  BD_partition<-c(BDtotal,BDspace,BDspt,BDtime)
  names(BD_partition)<-c("Total","Space","SPT","Time")
  
  Avg_Cont_to_BDtotal<-c(BDspace,BDspt,BDtime)/Total_entries_per_dimension
  names(Avg_Cont_to_BDtotal)<-c("Avg.Space","Avg.SPT","Avg.Time")
  
  Result<-list(BD_partition,Avg_Cont_to_BDtotal,Total_entries_per_dimension)
  names(Result)<-c("BDpartition","Avg. Cont to BDtotal","Total entries in D per dimension")
  return(Diss_modify_total)
  
}

#2) Convert data into the right format --------
    Cam_data0 <- readRDS(cam_data_path) 

    Cam_data1 <- Cam_data0%>% 
      mutate(Year_Location = paste(Year, locationName, sep = "_")) %>% 
      select(-c(locationName, Year,locationName_Year)) 
    
    #Turn table into matrix: 
    matrix_data <- as.matrix(Cam_data1[, -ncol(Cam_data1)])
    # Set row names to the first column (Year_Location)
    rownames(matrix_data) <- Cam_data1$Year_Location
    
#3) Calculate and plot Average BD ------------------------------------------

Result <-BD_partition(matrix_data,Component = "Total")
    
avg_cont_values<- Result$`Avg. Cont to BDtotal`
    
# Creating a data frame for ggplot
result_data <- data.frame(
  Metric = names(avg_cont_values),
  Value = as.numeric(avg_cont_values))

# Creating a bar plot using ggplot
average_BD_plot<- ggplot(result_data, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", aes(fill = Metric), color = "black", linewidth = 0.5) +  # Update to linewidth
  labs(
    title = expression("Average Contributions to D"[total]),
    x = "",
    y = expression("Average Contributions to D"[total])
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Avg.Space" = "#EDDE9F", "Avg.SPT" = "#C3D1F7", "Avg.Time" = "#344D7E")) +  # Specify fill colors
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none",  # Remove the legend
    text = element_text(size = 16),  # Adjust font size
    axis.title = element_text(size = 16),  # Adjust axis title font size
    axis.text = element_text(size = 16),   # Adjust axis text font size
    plot.title = element_text(size = 20)   # Adjust plot title font size
  ) +
  scale_x_discrete(labels = c(expression("D"[Space]), expression("D"[SpT]), expression("D"[Time])))

#Save the plot: 
ggsave(average_BD_plot_png, plot = average_BD_plot, width = 6, height = 6, units = "in")

#4) Test and plot differences in dissimilarity values: -------------------------------------
dissimilarity_matrix <- BD_partition_matrix(matrix_data, Component = "Total")
saveRDS(dissimilarity_matrix, dissimilarity_matrix_BrayC_path)

#Divide table
D_space <- filter(dissimilarity_matrix, Dim == "space")
D_time <- filter(dissimilarity_matrix, Dim == "time")
D_spt<- filter(dissimilarity_matrix, Dim == "spt")

#Large differences in sample size: 
nrow(D_space)
nrow(D_time)
nrow(D_spt)

summary(D_space)
summary(D_time)
summary(D_spt)

### Kruskal-Wallis test to test for differences: 

# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(distance ~ Dim, data = dissimilarity_matrix)

# Perform post hoc Dunn's test with Bonferroni correction
posthoc_result <- dunn.test(dissimilarity_matrix$distance, dissimilarity_matrix$Dim, method = "bonferroni")

#Visualization: 
dissimilarity_matrix_adj <- dissimilarity_matrix %>%
  mutate(Dim = case_when(
    Dim == "space" ~ "Space",
    Dim == "time" ~ "Time",
    Dim == "spt" ~ "SpT",
    TRUE ~ Dim  # Default case to handle unexpected values
  ))
D_index_plot <- ggplot(dissimilarity_matrix_adj, aes(x = Dim, y = distance, fill = Dim)) +
                geom_boxplot() +
                labs(
                  title = "Dissimilarity index by dimension",
                  x = "",
                  y = "Dissimilarity index"
                ) +
                
                scale_fill_manual(values = c("#EDDE9F",  "#C3D1F7", "#344D7E")) +  # Specify box colors
                theme_minimal() +
                theme(
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  legend.position = "none",  # Remove the legend
                  text = element_text(size = 16),  # Adjust font size
                  axis.title = element_text(size = 16),  # Adjust axis title font size
                  axis.text = element_text(size = 16),   # Adjust axis text font size
                  plot.title = element_text(size = 20)   # Adjust plot title font size
                )
ggsave(dissimilarity_index_plot_png, plot = D_index_plot, width = 6, height = 6, units = "in")   

#5) Test the stochasity / deterministic effect - dimension differences-------------------
### Define functions ----
#Function to generate random communities:
generate_null_community <- function(original_community_data) {
  null_community_data <- original_community_data
  for (col in 1:ncol(null_community_data)) {
    null_community_data[, col] <- sample(null_community_data[, col])
  }
  return(null_community_data)
}
###
#   a) Generate multiple random communities ----
observed_dissimilarities <- BD_partition_matrix(matrix_data, Component = "Total")
num_permutations <- 1000
null_dissimilarities_base<- select(BD_partition_matrix(matrix_data, Component = "Total"), -distance)

# Create an empty table to store dissimilarity values
null_dissimilarities_table <- null_dissimilarities_base

# Generate multiple random communities and add them to null_dissimilarities_table as columns
for (i in 1:num_permutations) {
  null_community_data <- generate_null_community(matrix_data)
  null_dissimilarity_values <- select(BD_partition_matrix(null_community_data, Component = "Total"), distance)
  null_dissimilarities_table <- cbind(null_dissimilarities_table, null_dissimilarity_values)
} 

#Save the data:
saveRDS(null_dissimilarities_table, null_dissimilarities_table_raw_normalized)

#   b) Calculate mean and standard deviation for each row of the null dissimilarities table ----

#(OPTIONAL) Read the data if already available 
# null_dissimilarities_table<- readRDS(null_dissimilarities_table_raw_normalized)

row_means <- apply(null_dissimilarities_table[, -c(1,2,3)], 1, mean)
row_sds <- apply(null_dissimilarities_table[, -c(1, 2,3)], 1, sd)

null_dissimilarities <- cbind(null_dissimilarities_table, Mean = row_means, SD = row_sds) %>% 
  select(c("c1","c2","Dim","Mean","SD")) 

#   c) Standardized Effect Size Calculation ----
SES_total <- (select(observed_dissimilarities, distance) - select(null_dissimilarities, Mean)) / select(null_dissimilarities, SD)
SES_space <- (select(filter(observed_dissimilarities, Dim =="space"), distance) - select(filter(null_dissimilarities, Dim =="space"),Mean)) / select(filter(null_dissimilarities, Dim =="space"),SD)
SES_time  <- (select(filter(observed_dissimilarities, Dim =="time"), distance) - select(filter(null_dissimilarities, Dim =="time"),Mean)) / select(filter(null_dissimilarities, Dim =="time"),SD)
SES_spt   <- (select(filter(observed_dissimilarities, Dim =="spt"), distance) - select(filter(null_dissimilarities, Dim =="spt"),Mean)) / select(filter(null_dissimilarities, Dim =="spt"),SD)

summary(SES_time)
summary(SES_space)
summary(SES_spt)
    
#   d) T-test: ----
t.test(SES_total$distance, mu = 0)# mean -0.002622288 , p-value = 0.6325
wilcox.test(SES_total$distance, mu = 0)#p-value < 2.2e-16

t.test(SES_space, mu = 0)#mean -0.4670475, p-value < 2.2e-16
wilcox.test(SES_space$distance, mu = 0)#p-value < 2.2e-16

t.test(SES_time$distance, mu = 0)#mean of x -0.04511691 p-value = 0.01292
wilcox.test(SES_time$distance, mu = 0)# p-value = 0.1846

t.test(SES_spt, mu = 0)#mean of x 0.01816843, p-value = 0.001794
wilcox.test(SES_spt$distance, mu = 0)#p-value = < 2.2e-16

#Check for normality: 
hist((SES_total$distance), main="Histogramm der Daten", xlab="Werte")
qqnorm(SES_total$distance)
qqline(SES_total$distance)
set.seed(123)
sample_data <- sample(SES_total$distance, 5000)
shapiro.test(sample_data)

#Test for differences among groups: 
# Add a new column to each SES data frame to indicate the dimension
SES_space$Dimension <- 'space'
SES_time$Dimension <- 'time'
SES_spt$Dimension <- 'spt'

#Combine the data frames into one
combined_SES <- rbind(SES_space, SES_time, SES_spt)
glimpse(combined_SES)

#test for normality: 
shapiro.test(combined_SES$distance[combined_SES$Dimension == "space"])
shapiro.test(combined_SES$distance[combined_SES$Dimension ==  "time"])
shapiro.test(combined_SES$distance[combined_SES$Dimension ==  "spt"])

#Homogenity of variance 
leveneTest(distance ~Dimension, data = combined_SES)

#Assumptions for ANOVA not fulfilled --> Kruskal Wallis test instead 
kruskal.test(distance ~ Dimension, data = combined_SES)
dunn.test(combined_SES$distance, combined_SES$Dimension, method="bonferroni")

#   e) Boxplot of SES ----

Mean_SES_plot_horizontal <- 
    ggplot() +
      #geom_boxplot(aes(x = SES_total$distance, y = "SES_total", fill = "SES_total"), color = "black") +
      geom_boxplot(aes(x = SES_space$distance, y = "SES_space", fill = "SES_space"), color = "black") +
      geom_boxplot(aes(x = SES_time$distance, y = "SES_time", fill = "SES_time"), color = "black") +
      geom_boxplot(aes(x = SES_spt$distance, y = "SES_spt", fill = "SES_spt"), color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  # Add a line at x = 0
      labs(title = "Mean SES by dimension", x = "Mean SES", y = "") +
      scale_fill_manual(values = c("#EDDE9F", "#C3D1F7", "#344D7E", "#D0C7A4")) +  # Specify box colors
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        legend.position = "none",  # Remove the legend
        text = element_text(size = 14),  # Adjust font size
        axis.title = element_text(size = 14),  # Adjust axis title font size
        axis.text = element_text(size = 12),   # Adjust axis text font size
        plot.title = element_text(size = 16)   # Adjust plot title font size
      )+ 
      scale_y_discrete(labels = c(expression("BD"[Space]), expression("BD"[SpT]), expression("BD"[Time]),expression("BD"[Total]) ))

Mean_SES_plot <- ggplot() +
    geom_boxplot(aes(y = SES_space$distance, x = "Space", fill = "SES_space"), color = "black") +
    geom_boxplot(aes(y = SES_time$distance, x = "Time", fill = "SES_time"), color = "black") +
    geom_boxplot(aes(y = SES_spt$distance, x = "SpT", fill = "SES_spt"), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Use geom_hline for horizontal line
    labs(title = "SES by dimension", x = "", y = "SES") +
    scale_fill_manual(values = c("#EDDE9F", "#C3D1F7", "#344D7E", "#D0C7A4")) +  # Specify box colors
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      legend.position = "none",  # Remove the legend
      text = element_text(size = 16),  # Adjust font size
      axis.title = element_text(size = 18),  # Adjust axis title font size
      axis.text = element_text(size = 16),   # Adjust axis text font size
      plot.title = element_text(size = 20)   # Adjust plot title font size
    )

ggsave(mean_SES_plot_horizontal, plot = Mean_SES_plot_horizontal, width = 6, height = 6, units = "in")
ggsave(mean_SES_plot, plot = Mean_SES_plot, width = 6, height = 6, units = "in")

#   f) Graph to compare means of observed and null dissimilarity--------------------

Space<- data.frame(c(select(filter(observed_dissimilarities, Dim =="space"), distance), select(filter(null_dissimilarities, Dim =="space"),Mean))) %>% 
  na.omit
Time <- data.frame(c(select(filter(observed_dissimilarities, Dim =="time"), distance), select(filter(null_dissimilarities, Dim =="time"),Mean)))%>% na.omit
 SpT <- data.frame(c(select(filter(observed_dissimilarities, Dim =="spt"), distance), select(filter(null_dissimilarities, Dim =="spt"),Mean)))%>% na.omit 

median(Space$distance)
mean(Time$Mean)

# Combine the data and reshape it
combined_data <- bind_rows(
  mutate(Space, Group = "Space"),
  mutate(Time, Group = "Time"),
  mutate(SpT, Group = "SpT")
) %>% 
  na.omit()

# Calculate mean values for each group and variable
means <- aggregate(. ~ Group, data = combined_data, mean)

# Reshape data to long format
means_long <- tidyr::gather(means, key = "Variable", value = "Mean_Value", -Group)

# Create a grouped barplot
ggplot(means_long, aes(x = Group, y = Mean_Value, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparison of Mean Values by Group",
       x = "Group",
       y = "Mean Value") +
  scale_fill_manual(values = c("Mean" = "blue", "distance" = "green")) +
  theme_minimal()

ggplot(means_long, aes(x = Group, y = Mean_Value, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = ifelse(Variable == "Mean", "NULL", "OBS")),
            position = position_dodge(width = 0.8),
            vjust = 1.5, size = 3) +
  labs(title = "Comparison of Mean Values by Group",
       x = "Group",
       y = "Mean Value") +
  scale_fill_manual(values = c("Mean" = "blue", "distance" = "green")) +
  theme_minimal()

ggplot(means_long, aes(x = Group, y = Mean_Value, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = ifelse(Variable == "Mean", "OBS", "NULL")),
            position = position_dodge(width = 0.8),
            vjust = 1.5, size = 3) +
  labs(title = "Comparison of Mean Values by Group",
       x = "Group",
       y = "Mean Value") +
  scale_fill_manual(values = c("Mean" = "#EDDE9F", "distance" = "#C3D1F7")) +
  theme_minimal()

#6) Test stochasity / deterministic effect - year differences-------
# NOTE:
# Using variable observed_dissimilarities defined in 5a)
# Using variable null_dissimilarities defined in 5b)
    
# Extract years from community identifiers
observed_dissimilarities_years <- observed_dissimilarities %>%
  mutate(year = as.character(substring(c1, 1, 4))) # Assuming year is the first part of the community ID

null_dissimilarities_years <- null_dissimilarities %>%
  mutate(year = as.character(substring(c1, 1, 4))) # Do the same for the null dataset

ses_by_year <- list()
years<- unique(observed_dissimilarities_years$year)

# Calculate SES for each year
for (year_val in years) {
  observed_subset_years <- filter(observed_dissimilarities_years, .data$year == year_val & Dim == "space")
  null_subset_years <- filter(null_dissimilarities_years, .data$year == year_val & Dim == "space")
  
  # Assuming you have a vector of distances in observed_subset and mean, sd values in null_subset
  ses <- (observed_subset_years$distance - null_subset_years$Mean) / null_subset_years$SD
  ses_by_year[[as.character(year_val)]] <- ses
}

# Convert the list to a data frame
ses_data_years <- do.call(rbind, lapply(names(ses_by_year), function(year) {
  data.frame(year = year, SES = ses_by_year[[year]], Dimension = "SES by Year")
})) %>%
  mutate(year = factor(year)) # Ensure year is treated as a categorical variable

# ses_by_year now contains the SES values for each year
sS<- ses_data %>% 
  group_by(year) %>% 
  summarise(mean(SES))

#test for normality: 
shapiro.test(ses_data$SES[ses_data$year == 2011])
shapiro.test(ses_data$SES[ses_data$year == 2012])
shapiro.test(ses_data$SES[ses_data$year == 2013])
shapiro.test(ses_data$SES[ses_data$year == 2014])
shapiro.test(ses_data$SES[ses_data$year == 2015])
shapiro.test(ses_data$SES[ses_data$year == 2016])
shapiro.test(ses_data$SES[ses_data$year == 2019])
shapiro.test(ses_data$SES[ses_data$year == 2021])

#Assumptions for ANOVA not fulfilled --> Kruskal Wallis test instead
kruskal.test(SES ~ year, data = ses_data)
dunn.test(ses_data$SES, ses_data$year, method="bonferroni")

#Generate Graph: 
SES_plot_years <- ggplot(ses_data_years, aes(x = year, y = SES)) +
  geom_boxplot(color = "black", fill = "#EDDE9F") +  # Set fill to a single color
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add a line at y = 0
  labs(title = "SES by Year", x = "Year", y = "SES") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "none",  # Optionally remove the legend if it's not needed
    text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20)
  )
ggsave(mean_SES_plot_years, plot = SES_plot_years, width = 8, height = 7, units = "in")

#7) Test stochasity / deterministic effect - location differences-------
# NOTE:
# Using variable observed_dissimilarities defined in 5a)
# Using variable null_dissimilarities defined in 5b)

# Extract location from community identifiers
# Assuming the location follows the year, e.g., "2011_BCI-1-01"
observed_dissimilarities_sp <- observed_dissimilarities %>%
mutate(location = sub("^[0-9]+_", "", c1)) # Remove the year and underscore

null_dissimilarities_sp <- null_dissimilarities %>%
mutate(location = sub("^[0-9]+_", "", c1)) # Do the same for the null dataset

# Calculate SES for each location
ses_by_location <- list()

locations <- unique(observed_dissimilarities_sp$location)

for (location_val in locations) {
observed_subset_sp <- filter(observed_dissimilarities_sp, .data$location == location_val & Dim == "time")
null_subset_sp <- filter(null_dissimilarities_sp, .data$location == location_val & Dim == "time")
ses <- (observed_subset_sp$distance - null_subset_sp$Mean)/ null_subset_sp$SD
ses_by_location[[location_val]] <- ses
}
# ses_by_location now contains the SES values for each location

# Convert the list to a data frame
ses_data_sp <- do.call(rbind, lapply(names(ses_by_location), function(location) {
data.frame(location = location, SES = unlist(ses_by_location[[location]]), Dimension = "SES by Location")
})) %>%
mutate(location = factor(location)) # Ensure location is treated as a categorical variable

# Correctly summarise SES by location
ses_summary_sp <- ses_data_sp %>% 
group_by(location) %>% 
summarise(Mean_SES = mean(SES, na.rm = TRUE))

# Test for normality by location (example for a few locations, adjust as necessary)
shapiro.tests <- lapply(levels(ses_data_sp$location), function(loc) {
shapiro.test(ses_data_sp$SES[ses_data_sp$location == loc])
})

# Kruskal Wallis test (as assumptions for ANOVA are not met)
kruskal_result <- kruskal.test(SES ~ location, data = ses_data_sp)
dunn.test(ses_data_sp$SES, ses_data_sp$location, method="bonferroni")

# Generate Graph
ses_plot_by_location <- ggplot(ses_data_sp, aes(x = location, y = SES, fill = location)) +
geom_boxplot(color = "black", fill= "#344D7E") +
geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add a line at y = 0
labs(title = "SES by Location", x = "Location", y = "SES") +
theme_minimal() +
theme(
  panel.border = element_rect(color = "black", fill = NA, size = 0.5),
  legend.position = "none",  # Optionally remove the legend if it's not needed
  text = element_text(size = 18),
  axis.title = element_text(size = 18),
  axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for better readability
  axis.text = element_text(size = 16),
  plot.title = element_text(size = 20)
)

ggsave(mean_SES_plot_location, plot =ses_plot_by_location, width = 14, height = 6, units = "in")

### Hayden Gallo ###
### 8/3/2023     ###
### AGB Increment vs. 95% CI Reconstructed Precip ###
### This analysis resulted from Gam findings ###

### clear global environment ###

rm(list=ls())

# load packages #

library(reshape2)
library(tidyr)
library(abind)
library(gratia)
library(ggplot2)
library(Rmisc)
library(ggforce)

### Pick site 

site = 'NRP'

if (site == 'HF'){
  site.alt = 'HARVARD'
  species <- c('Red_Maple', 'Yellow_Birch', 'Red_Oak', 'Hemlock')
} else if (site == 'Goose'){
  site.alt = 'GOOSE'
  species <- c('American_Beech', 'East_White_Pine', 'White_Oak', 'Chestnut_Oak', 'Red_Oak')
}  else if (site == 'Rooster'){
  site.alt = 'ROOSTER'
  species <- c('Red_Maple', 'American_Beech', 'Red_Spruce', 'East_White_Pine', 'Red_Oak')
} else if (site == 'NRP'){
  site.alt = 'NORTHROUND'
  species <- c('Black_Birch', 'East_White_Pine', 'Red_Oak', 'Hemlock')
}


 
dir = paste0('/save/workflows/hgallo')

### loading tree ring data and transforming to species difference in agb per year

load(paste0('/data/dbfiles/sda.obs.',site.alt,'.Rdata'))
tree_ring_data <- data.frame(obs.list$obs.mean)
rownames(tree_ring_data) <- species
final_year <- 1899+ncol(tree_ring_data)
tree_dates <- seq(as.Date("1900",  format = "%Y"), length.out = ncol(tree_ring_data), by = "year")
tree_ring_data <- data.frame(t(tree_ring_data))
tree_ring_data$year <- 1900:final_year
tree_ring_data$date <- tree_dates
for (i in 1:length(species)){
  temp_diff <- diff(tree_ring_data[,i])
  temp_diff <- append(temp_diff,0,0)
  tree_ring_data[,i] <- temp_diff
}


tree_ring_data <- pivot_longer(tree_ring_data, cols = all_of(species), names_to = 'species')
colnames(tree_ring_data) <- c('year','date','species','agb')
tree_ring_data$value_type <- 'Tree_Ring'

### load met data 

if (site == 'NRP'){
  site.alt = site
}


csv_clim_data <- paste0('/data/gallo/Met_updated/',site,'_met/weights/ensemble-weights-',site.alt,'-prism.csv')

clim_data <- read.csv(csv_clim_data)
clim_data <- clim_data[,1:3]

### create array for precip met data

precip_array <- array(0, c(116, 12, nrow(clim_data)))

### loading precip for each reconstruction

for (i in 1:nrow(clim_data)){
  clim_to_load <- clim_data[i,2]
  load(paste0('/data/gallo/Met_updated/',site,'_met/linkages/',clim_to_load,'.Rdata'))
  precip_array[,,i] <- tail(precip.mat, 116)
}

CI_precip_array <- array(0, c(116,12,3))

### confidence interval for all precip reconstructions

for (i in 1:116){
  for (j in 1:12){
    CI_precip_array[i,j,] <- CI(precip_array[i,j,])
  }
}

dimnames(CI_precip_array) <- list(c(1:116), c('January','February','March','April','May','June','July','August','September','October','November','December'), c('lower','mean','upper'))

### creating date vector for plotting and dataframes

dates <- seq(as.Date("01-01-1900",  format = "%d-%m-%Y"), length.out = 1392, by = "month")

### Collecting precip data into dataframe

CI_precip_df <- melt(CI_precip_array, varnames = names(dimnames(CI_precip_array)))
colnames(CI_precip_df) <- c('year','month','CI','precip')
CI_precip_df <- reshape(CI_precip_df, idvar = c('year','month'), timevar = 'CI', direction = 'wide')
CI_precip_df <- CI_precip_df[order(CI_precip_df$year, decreasing = FALSE),]
CI_precip_df$date <- dates

### scaling factor for secondary y axis in ggplot
scale = 100

### ggplot for comparing increment tree growth with precipitation 95% CI

ggplot(CI_precip_df, aes(x = date, y = precip.mean)) + geom_line(aes(color = "Precip")) + geom_ribbon(aes(ymin = precip.lower, ymax = precip.upper))+ 
   geom_line(data = tree_ring_data, aes(x = date, y = agb*scale, color = 'Species AGB')) + facet_grid(month~species)
  scale_y_continuous(sec.axis = sec_axis(~./scale, name = 'Yearly Diff in Agb')) + labs(x = 'Date', y = 'Precip (cm)', title = site) + scale_color_manual(values = c('red', 'black','orange','blue','darkgreen','grey'))
#ggplot(CI_precip_df, aes(x = year, y = precip.mean)) + geom_bar(stat = 'identity') + geom_errorbar(aes(ymin = precip.lower, ymax = precip.upper)) + facet_grid_paginate(~month, ncol=1, nrow = 1, page = 6) 




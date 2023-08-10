### Hayden Gallo  ###
### July 2023 ###
### MC runs vs tree ring ###

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
library(R.utils)

### Pick site 

site = 'HF'

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
setwd(dir)

### loading tree ring data and transforming to species difference in agb per year

load(paste0('/data/dbfiles/sda.obs.',site.alt,'.Rdata'))
tree_ring_data <- data.frame(obs.list$obs.mean)
rownames(tree_ring_data) <- species
final_year <- 1899+ncol(tree_ring_data)
tree_dates <- seq(as.Date("1900",  format = "%Y"), length.out = ncol(tree_ring_data), by = "year")
tree_ring_data <- data.frame(t(tree_ring_data))
tree_ring_data$year <- 1900:final_year
tree_ring_data$date <- tree_dates


tree_ring_data <- pivot_longer(tree_ring_data, cols = all_of(species), names_to = 'species')
colnames(tree_ring_data) <- c('year','date','species','agb')
tree_ring_data$value_type <- 'Tree_Ring'


### Load MC runs data

load(paste0('MC_runs/',site,'/params_all_mets_all.Rdata'))

### vector of parameters

all_parameters <- list('D3', 'MPLANT', 'DMAX', 'Frost', 'AGEMX', 'DMIN', 'G', 'SPRTND')


### Create df with CI of species specific biomass



CI_pft_biomass <- array(0, c(101,3,length(species),8))
dimnames(CI_pft_biomass) <- list(1915:2015, c('lower','mean','upper'), species, all_parameters)

spp_agb <- master_array_pft[,,,,,1]
spp_agb <- wrap(spp_agb, map = list(1,2,NA,4))
dimnames(spp_agb)[[1]] <- species

  for(j in 1:8){
    for (k in 1:101){
      for(l in 1:length(species)){
        
        CI_pft_biomass[k,,l,j] <-CI(spp_agb[l,k,,j])
      }
    }
  }

### Melting CI data to dataframe

CI_pft_biomass_melt <- melt(CI_pft_biomass, varnames = names(dimnames(CI_pft_biomass)))
colnames(CI_pft_biomass_melt) <- c('year','CI','species','param','agb')
CI_pft_biomass_melt <- reshape(CI_pft_biomass_melt, idvar = c('year','param','species'), timevar = 'CI', direction = 'wide')
CI_pft_biomass_melt$value_type <- "Linkages"

### Merge Tree Ring and MC runs data

p_diff_plot_df <- merge(CI_pft_biomass_melt, tree_ring_data, by  = c('year', 'species'))

colnames(CI_pft_biomass_melt)[5] <- 'agb'

### Plot Species specific LINKAGES vs Tree Rings

ggplot(CI_pft_biomass_melt, aes(x = year, y = agb)) + geom_line(aes(color = species, linetype = 'LINKAGES')) + facet_wrap(~param, ncol = 2) + 
  geom_line(data = tree_ring_data, aes(x = year, y = agb, color = species, linetype = 'Tree_Ring'))

### Plot Diff in Tree Rings vs LINKAGES

p_diff_plot_df$agb_diff <- (abs((p_diff_plot_df$agb - p_diff_plot_df$agb.mean)/p_diff_plot_df$agb.mean))

ggplot(p_diff_plot_df, aes(x = year, y = agb_diff)) + geom_line(aes(color = species)) + facet_wrap(~param, ncol =2)









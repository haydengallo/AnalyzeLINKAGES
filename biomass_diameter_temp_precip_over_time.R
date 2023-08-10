### Hayden Gallo
### June 2023

### Visualizations of biomass, diameter, temp, precip over time, MC runs

### need to load data formatting function

rm(list=ls())

library(ggplot2)
library(Rmisc)
library(reshape2)
library(tidyr)
library(dplyr)
library(pivottabler)

### setting working directory

dir = '/save/workflows/hgallo'
setwd(dir)

### choose site to investigate

site = 'Goose'

### sourcing function that collects MC runs data into arrays

source('/pecan/VM_scripts/HG/collating_data_MC_runs_function.R')

### use commented out row if you want to create new data product with function

#collating_MC_runs('NRP', 100, 101, 4, 'all','all','yes')

### this line loads already present data product 

load(paste0('MC_runs/',site,'/params_all_mets_all.Rdata'))

par(mfrow = c(1,1))

### vector of parameters

all_parameters <- list('D3', 'MPLANT', 'DMAX', 'Frost', 'AGEMX', 'DMIN', 'G', 'SPRTND')


### just grabs the unique mets used in these particular MC runs

clim <- dimnames((master_array_tot_biomass))[[4]]

### Parameter Distribution ###

### grabbing species specific agb from final year of runs

pft_agb_fin_yr_df <- melt(master_array_pft[,101,,,,1])
colnames(pft_agb_fin_yr_df) <- c('species','run','param','met','agb')

### gets parameter values from parameter that is varied

varied_param_df <- melt(master_array_varied_param)
colnames(varied_param_df) <- c('species','run','param','met','param_value')

param_dist_vs_agb_df <- merge(pft_agb_fin_yr_df, varied_param_df, by = c('species','run','param','met')) 


ggplot(param_dist_vs_agb_df, aes(x = param_value, y = agb, color = met)) + geom_point() + 
  facet_grid(vars(species),vars(param), scales = 'free') + labs(title = paste0('Parameter Distribution vs. Final Year AGB by Species at ',site))
ggsave(paste0(site,'_param_dist.png'), width = 12, height = 8)

### Biomass Overtime ###

### Non Species Specific Over Time ###

### creating confidence intervals for biomass over time

CI_tot_biomass <- array(0, c(101,3,8,10))

for (i in 1:10){
  for(j in 1:8){
    for (k in 1:101){

      CI_tot_biomass[k,,j,i] <-CI(master_array_tot_biomass[k,,j,i])
  
    }
  }
}

### melting CI data to df

dimnames(CI_tot_biomass) <- list(c(1:101), c('lower','mean','upper'), all_parameters, dimnames(master_array_tot_biomass)[[4]])
CI_tot_biomass_melt <- melt(CI_tot_biomass, varnames = names(dimnames(CI_tot_biomass)))
colnames(CI_tot_biomass_melt) <- c('year','CI','param','met','agb')
CI_tot_biomass_melt <- reshape(CI_tot_biomass_melt, idvar = c('year','param','met'), timevar = 'CI', direction = 'wide')
CI_tot_biomass_melt <- split(CI_tot_biomass_melt, f = CI_tot_biomass_melt$met)

### plotting for each individual climate scenario wrapped by parameter

for (i in 1:10){

temp <- ggplot(CI_tot_biomass_melt[[i]], aes(x = year, y = agb.mean)) + geom_line() +
geom_ribbon(aes(ymin = agb.lower, ymax = agb.upper, alpha = 0.1)) + facet_wrap(~ param, ncol = 2) +
  labs(title = clim[i])
print(temp)
#ggsave(paste0('/save/workflows/hgallo/MC_runs/',site,'/Plots/tot_agb_',clim[i],'.png'))
}

 
### Species Specific Over Time ###

### Confidence interval of species specific aboveground biomass overtime

CI_pft_biomass <- array(0, c(101,3,4,8,10))
dimnames(CI_pft_biomass) <- list(c(1:101), c('lower','mean','upper'), dimnames(master_array_pft)[[1]], all_parameters, dimnames(master_array_tot_biomass)[[4]])

spp_agb <- master_array_pft[,,,,,1]

for (i in 1:10){
  for(j in 1:8){
    for (k in 1:101){
      for(l in 1:4){
      
      CI_pft_biomass[k,,l,j,i] <-CI(spp_agb[l,k,,j,i])
      }
    }
  }
}

### Melting CI data to dataframe

CI_pft_biomass_melt <- melt(CI_pft_biomass, varnames = names(dimnames(CI_tot_biomass)))
colnames(CI_pft_biomass_melt) <- c('year','CI','species','param','met','agb')
CI_pft_biomass_melt <- reshape(CI_pft_biomass_melt, idvar = c('year','param','met','species'), timevar = 'CI', direction = 'wide')
CI_pft_biomass_melt <- split(CI_pft_biomass_melt, f = CI_pft_biomass_melt$met)

### Plots overtime for each climate scenario wrapped by parameter

for (i in 1:10){
  
  temp <- ggplot(CI_pft_biomass_melt[[i]], aes(x = year, y = agb.mean, color = species)) + geom_line() + 
    geom_ribbon(aes(ymin = agb.lower, ymax = agb.upper, alpha = 0.1), show.legend = FALSE) + facet_wrap(~ param, ncol = 2) +
    labs(title = clim[i]) + theme(legend.position = 'bottom')
  print(temp)
  #ggsave(paste0('/save/workflows/hgallo/MC_runs/',site,'/Plots/pft_agb_',clim[i],'.png'))
  
}


### Average Diameter Over Time ###

CI_pft_diam <- array(0, c(101,3,4,8,10))
dimnames(CI_pft_diam) <- list(c(1:101), c('lower','mean','upper'), dimnames(master_array_pft)[[1]], all_parameters, dimnames(master_array_tot_biomass)[[4]])

### grabbing diameter data from arrays

avg_diam <- master_array_pft[,,,,,5]

### creating confidence interval for diameter overtime

for (i in 1:10){
  for(j in 1:8){
    for (k in 1:101){
      for(l in 1:4){
        
        CI_pft_diam[k,,l,j,i] <-CI(avg_diam[l,k,,j,i])
      }
    }
  }
}

### melting diameter CI to dataframe

CI_pft_diam_melt <- melt(CI_pft_diam, varnames = names(dimnames(CI_pft_diam)))
colnames(CI_pft_diam_melt) <- c('year','CI','species','param','met','avg_diam')
CI_pft_diam_melt <- reshape(CI_pft_diam_melt, idvar = c('year','param','met','species'), timevar = 'CI', direction = 'wide')
CI_pft_diam_melt <- split(CI_pft_diam_melt, f = CI_pft_diam_melt$met)

### Plotting for each climate reconstruction wrapped by param

for (i in 1:10){
  
  temp <- ggplot(CI_pft_diam_melt[[i]], aes(x = year, y = avg_diam.mean, color = species)) + geom_line() + 
    geom_ribbon(aes(ymin = avg_diam.lower, ymax = avg_diam.upper, alpha = 0.1), show.legend = FALSE) + facet_wrap(~ param, ncol = 2) +
    labs(title = paste0('Average Diam Across Runs ', clim[i])) + theme(legend.position = 'bottom')
  print(temp)
  #ggsave(paste0('/save/workflows/hgallo/MC_runs/',site,'/Plots/avg_diam_',clim[i],'.png'))
  
}


### Plotting Binned Diameters Over Time

rm(master_array_algf)

### melting all dbh data to dataframe and then splitting by met

dbh_df <- melt(dbh_master_array, names(dimnames(dbh_master_array)))
colnames(dbh_df) <- c('individual','year','run','param','met','diameter')
dbh_df <- split(dbh_df, f = dbh_df$met)

### removing any zero entries in dbh for each met

for (i in 1:10){

dbh_df[[i]] <- dbh_df[[i]][dbh_df[[i]]$diameter != 0 ,]
print(i)
}

### binning diameter data for each met and then plotting 

for (i in 1:10){

dbh_df[[i]] <- dbh_df[[i]] %>% mutate(diam_bins = cut(diameter, breaks = c(0,10,20,30,40,50)))
test_plotting <- dbh_df[[i]] %>% group_by(year,param ,diam_bins) %>% summarise(n()/100)


temp_plot <- ggplot(test_plotting, aes(x = year, y = `n()/100`, color = diam_bins)) + geom_point() + facet_wrap(~ param, ncol = 2) +
  labs(title = paste0('Diam Across Runs ', clim[i])) + theme(legend.position = 'bottom') + geom_jitter()

print(temp_plot)
#ggsave(paste0('/save/workflows/hgallo/MC_runs/',site,'/Plots/bin_diam_',clim[i],'.png'))

}

### Plotting Precip and Temp Over Time 

### Precip ###

for (i in clim){
  
  ### Calling each climate scenario temporarily for plotting purposes
  
  met = i
  input = paste0('/data/gallo/Met_updated','/',site,'_met','/linkages/',met,'.Rdata')
  load(input)
  nyear <<- nyear_user
  end.year <<- 2015
  start.year <<- end.year-nyear_user+1
  precip.mat <<- tail(precip.mat, n = nyear_user)
  temp.mat <<- tail(temp.mat, n = nyear_user)
  
  ### Collecting precip data
  
  precip.mat <- data.frame(precip.mat)
  
  precip.mat[,13] <- 1915:2015
  
  colnames(precip.mat) <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Year')
  
  precip.mat <- pivot_longer(precip.mat, cols = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), names_to = 'month')
  
  
  temp_plot_precip <- ggplot(precip.mat, aes(x = Year, y = value)) + geom_point() + geom_smooth() + 
    facet_wrap( ~factor(month, levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))) + 
    labs(title = paste0('Precip Over Time ', met)) + xlab('Year') + ylab('Precipitation (cm)')
  
  print(temp_plot_precip)
  #ggsave(paste0('/save/workflows/hgallo/MC_runs/',site,'/Plots/precip_',met,'.png'))
  
  ### Collecting temp data 
  
  temp.mat <- data.frame(temp.mat)
  
  temp.mat[,13] <- 1915:2015
  
  colnames(temp.mat) <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Year')
  
  temp.mat <- pivot_longer(temp.mat, cols = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), names_to = 'month')
  
  
  temp_plot_precip <- ggplot(temp.mat, aes(x = Year, y = value)) + geom_point() + geom_smooth() + 
    facet_wrap( ~factor(month, levels = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))) + 
    labs(title = paste0('Temp Over Time ', met)) + xlab('Year') + ylab('Temp (C)')
  
  print(temp_plot_precip)
  #ggsave(paste0('/save/workflows/hgallo/MC_runs/',site,'/Plots/temp_',met,'.png'))
  
}

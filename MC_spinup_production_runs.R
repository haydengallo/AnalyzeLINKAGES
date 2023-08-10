### Hayden Gallo June 2023
### Spinup for Monte Carlo Simulations LINKAGES
### Script used for spinups for final production runs
### uses functions from script MC_helper_functions.R

### Monte Carlo Simulation Official Runs Simulation Routine
### Idea is to spinup LINKAGES for 1000 years prior to running MC sims
### To do this use highest weighted climate scenario for particular site
### Then use last 50 years and simulate using that 50 years x number of times until 1000 total years
### then save the state of the model there and use that state as starting place for MC sims

rm(list=ls())
setwd('/data/gallo')
#set correct directory for MC_helper_function script
source('~/pecan/VM_scripts/HG/MC_helper_functions.R')


library(linkages)
library(tidyverse)

### choose length of spinup

spinup_len = 1000

### Choose HF, Goose, NRP, Sylvania, Rooster
site  = 'Rooster'

if (site == 'HF'){
  site.alt = 'HARVARD'

  input = 'Start_Data/Harvard.input.updated.Rdata'
  load(input)

} else if (site == 'Goose'){
  site.alt = 'GOOSE'

  input = 'Start_Data/Goose.linkages.input.Rdata'
  load(input)

}  else if (site == 'Rooster'){
  site.alt = 'ROOSTER'

  input = 'Start_Data/Rooster.linkages.input.Rdata'
  load(input)

} else if (site == 'Sylvania'){
  site.alt = 'SYLVANIA'

  input = 'Start_Data/Sylvania.linkages.input.Rdata'
  load(input)

} else if (site == 'NRP'){
  site.alt = site

  input = 'Start_Data/NRP.linkages.input.Rdata'
  load(input)
}

met_weights <- read_csv(paste0('Met_updated/',site,'_met/weights/ensemble-weights-',site.alt,'-prism.csv'))
met_weights <- met_weights[,c(2:3)]
highest_wt_met <- met_weights[which.max(met_weights$wts),]


input = paste0('Met_updated','/',site,'_met','/linkages/',highest_wt_met$climate_model,'.Rdata')
load(input)

end.year = 2015
start.year = 1966
nyear <- 50
precip.mat <- tail(precip.mat, n = 50)
temp.mat <- tail(temp.mat, n = 50)

MC_sims_spinup(spinup_len = spinup_len)


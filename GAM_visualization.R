### GAM Visualization
### Hayden Gallo 
### July 2023

### clear global environment ###

rm(list=ls())

### load needed packages ###

library(mgcv)
library(reshape2)
library(tidyr)
library(abind)
library(gratia)
library(ggplot2)

### load run data ###

### select site to investigate

site = 'HF'

### setting working directory where GAM model outputs are stored, if your GAM outputs are not at this specific address, change working directory

dir = paste0('/save/workflows/hgallo/FROM CRC')
setwd(dir)

### loading gam output

load(paste0('GAM_',site,'_out.Rdata'))

### saving summary for easier acccess

summary_save <- summary(bam_higher_order_effects_by_species)

### test for viewing singular tensor plot

vis.gam(bam_higher_order_effects_by_species, view = c('ti(DMAX,precip_1,temp_1)'))

### Creating and saving pacf plot for understanding temporal lag

summary_save
png(paste0(site,'_pacf.png'))
pacf(residuals(bam_higher_order_effects_by_species))
dev.off()

### Creating and saving model diagnostic plots using gratia package

png(paste0(site,'_diagnostic_plots.png'))
#gam.check(bam_higher_order_effects_by_species)
gratia::appraise(bam_higher_order_effects_by_species)
dev.off()

### Creating and saving example tensor plot

png(paste0(site,'_example_tensor.png'))
gratia::draw(bam_higher_order_effects_by_species, select = c('ti(DMAX,precip_6,temp_6)'))
dev.off()

### Creating vector of all tensor plots for GAM

gam_tensor_plots <-                       c('ti(DMAX,precip_1,temp_1)', 
                                            'ti(DMAX,precip_2,temp_2)',
                                            'ti(DMAX,precip_3,temp_3)',
                                            'ti(DMAX,precip_4,temp_4)', 
                                            'ti(DMAX,precip_5,temp_5)', 
                                            'ti(DMAX,precip_6,temp_6)',
                                            'ti(DMAX,precip_7,temp_7)', 
                                            'ti(DMAX,precip_8,temp_8)', 
                                            'ti(DMAX,precip_9,temp_9)', 
                                            'ti(DMAX,precip_10,temp_10)', 
                                            'ti(DMAX,precip_11,temp_11)', 
                                            'ti(DMAX,precip_12,temp_12)', 
                                            'ti(AGEMX,precip_1,temp_1)', 
                                            'ti(AGEMX,precip_2,temp_2)', 
                                            'ti(AGEMX,precip_3,temp_3)', 
                                            'ti(AGEMX,precip_4,temp_4)', 
                                            'ti(AGEMX,precip_5,temp_5)', 
                                            'ti(AGEMX,precip_6,temp_6)', 
                                            'ti(AGEMX,precip_7,temp_7)', 
                                            'ti(AGEMX,precip_8,temp_8)', 
                                            'ti(AGEMX,precip_9,temp_9)', 
                                            'ti(AGEMX,precip_10,temp_10)',
                                            'ti(AGEMX,precip_11,temp_11)', 
                                            'ti(AGEMX,precip_12,temp_12)', 
                                            'ti(SPRTND,precip_1,temp_1)',  
                                            'ti(SPRTND,precip_2,temp_2)', 
                                            'ti(SPRTND,precip_3,temp_3)', 
                                            'ti(SPRTND,precip_4,temp_4)', 
                                            'ti(SPRTND,precip_5,temp_5)', 
                                            'ti(SPRTND,precip_6,temp_6)', 
                                            'ti(SPRTND,precip_7,temp_7)', 
                                            'ti(SPRTND,precip_8,temp_8)', 
                                            'ti(SPRTND,precip_9,temp_9)', 
                                            'ti(SPRTND,precip_10,temp_10)', 
                                            'ti(SPRTND,precip_11,temp_11)', 
                                            'ti(SPRTND,precip_12,temp_12)', 
                                            'ti(G,precip_1,temp_1)',  
                                            'ti(G,precip_2,temp_2)', 
                                            'ti(G,precip_3,temp_3)', 
                                            'ti(G,precip_4,temp_4)', 
                                            'ti(G,precip_5,temp_5)', 
                                            'ti(G,precip_6,temp_6)', 
                                            'ti(G,precip_7,temp_7)', 
                                            'ti(G,precip_8,temp_8)', 
                                            'ti(G,precip_9,temp_9)', 
                                            'ti(G,precip_10,temp_10)', 
                                            'ti(G,precip_11,temp_11)', 
                                            'ti(G,precip_12,temp_12)', 
                                            'ti(D3,precip_1,temp_1)',  
                                            'ti(D3,precip_2,temp_2)', 
                                            'ti(D3,precip_3,temp_3)', 
                                            'ti(D3,precip_4,temp_4)', 
                                            'ti(D3,precip_5,temp_5)', 
                                            'ti(D3,precip_6,temp_6)', 
                                            'ti(D3,precip_7,temp_7)', 
                                            'ti(D3,precip_8,temp_8)', 
                                            'ti(D3,precip_9,temp_9)', 
                                            'ti(D3,precip_10,temp_10)', 
                                            'ti(D3,precip_12,temp_12)',  
                                            'ti(D3,precip_11,temp_11)', 
                                            'ti(MPLANT,precip_1,temp_1)',  
                                            'ti(MPLANT,precip_2,temp_2)', 
                                            'ti(MPLANT,precip_3,temp_3)',
                                            'ti(MPLANT,precip_4,temp_4)', 
                                            'ti(MPLANT,precip_5,temp_5)', 
                                            'ti(MPLANT,precip_6,temp_6)', 
                                            'ti(MPLANT,precip_7,temp_7)', 
                                            'ti(MPLANT,precip_8,temp_8)', 
                                            'ti(MPLANT,precip_9,temp_9)', 
                                            'ti(MPLANT,precip_10,temp_10)', 
                                            'ti(MPLANT,precip_11,temp_11)', 
                                            'ti(MPLANT,precip_12,temp_12)', 
                                            'ti(FROST,precip_1,temp_1)',  
                                            'ti(FROST,precip_2,temp_2)', 
                                            'ti(FROST,precip_3,temp_3)', 
                                            'ti(FROST,precip_4,temp_4)', 
                                            'ti(FROST,precip_5,temp_5)', 
                                            'ti(FROST,precip_6,temp_6)', 
                                            'ti(FROST,precip_7,temp_7)', 
                                            'ti(FROST,precip_8,temp_8)', 
                                            'ti(FROST,precip_9,temp_9)',
                                            'ti(FROST,precip_10,temp_10)', 
                                            'ti(FROST,precip_12,temp_12)',  
                                            'ti(FROST,precip_11,temp_11)', 
                                            'ti(DMIN,precip_1,temp_1)',  
                                            'ti(DMIN,precip_2,temp_2)', 
                                            'ti(DMIN,precip_3,temp_3)', 
                                            'ti(DMIN,precip_4,temp_4)', 
                                            'ti(DMIN,precip_5,temp_5)', 
                                            'ti(DMIN,precip_6,temp_6)',
                                            'ti(DMIN,precip_7,temp_7)', 
                                            'ti(DMIN,precip_8,temp_8)', 
                                            'ti(DMIN,precip_9,temp_9)', 
                                            'ti(DMIN,precip_10,temp_10)', 
                                            'ti(DMIN,precip_11,temp_11)', 
                                            'ti(DMIN,precip_12,temp_12)')

### prints each individual tensor plot, this is the easiest way to print these plots in an effective and clear way that I found

for (i in 1:length(gam_tensor_plots)){
  print(gratia::draw(bam_higher_order_effects_by_species, select = gam_tensor_plots[i]))
}


#################################

### Gam Testing ###
### This script was part of testing phase of GAM development for MC runs ###
### Not intended for production quality models, used strictly for testing GAM formulations ###

#################################

### Hayden Gallo ###
### June 2023 ###

#################################

### clear global environment ###

rm(list=ls())

### load needed packages ###

library(mgcv)
library(reshape2)
library(tidyr)
library(abind)
library(gratia)

### load run data ###

sites <- 'HF'

for (i in sites){
  site = i

dir = paste0('/save/workflows/hgallo/MC_runs/HF/Gams')
setwd(dir)

nyear = 101
#load('merged_df_gam_testing.Rdata')
#load('gam_data_full_HF.Rdata')
data_to_load <- paste0('params_all_mets_all.Rdata')
load(data_to_load)
ctrl <- gam.control(trace = TRUE)

### Melting data to usable format ###

dimnames(master_array_pft)[[2]] <- 1:101
testing_data_agb <- master_array_pft[,,,,,1]
testing_data_params <- master_params[,,,,]

testing_data_params <- testing_data_params[,c('DMAX', "MPLANT", 'D3', 'AGEMX', 'FROST', 'DMIN', 'G', 'SPRTND'),,,]
test_df <- melt(testing_data_agb)
colnames(test_df) <- c('species', 'year', 'run', 'param','met', 'agb')

pft_df <- melt(testing_data_agb)
colnames(pft_df) <- c('species', 'year', 'run', 'param','met', 'agb')
pft_df <- pft_df[order(pft_df$species),]
pft_df$row_num <- 1:nrow(pft_df)

param_df <- melt(testing_data_params)
colnames(param_df) <- c('species', 'year', 'run', 'param','met', 'param_value')
merged_df <- merge(pft_df, param_df, by = c('species','run','param','met'))

merged_df <- merged_df %>% pivot_wider(names_from = year.y, values_from = param_value)
merged_df <- merged_df[order(merged_df$row_num),]


ar1_vector <- c(ifelse(merged_df$year.x == 1, TRUE, FALSE))

colnames(merged_df)[[5]] <- 'year'

### Incorporate Met Data ###

mets_used <- dimnames(testing_data_agb)[[5]]
mets_master_matrix <- array(0, c(nyear,24,length(mets_used)))
dimnames(mets_master_matrix) <- list(1:101, c(paste0('precip_',1:12),paste0('temp_',1:12)), mets_used)


count = 0

for (i in mets_used){
  count = count + 1
  met = i
  input = paste0('/data/gallo/Met_updated','/',site,'_met','/linkages/',met,'.Rdata')
  load(input)
  end.year <<- 2015
  start.year <<- end.year-nyear+1
  precip.mat <<- tail(precip.mat, n = nyear)
  temp.mat <<- tail(temp.mat, n = nyear)
  mets_master_matrix[,1:12,count] <- precip.mat
  mets_master_matrix[,13:24,count] <- temp.mat

}

mets_data <- melt(mets_master_matrix)
colnames(mets_data) <- c('year','month','met','value')
mets_data <- mets_data %>% pivot_wider(names_from = month, values_from = value)

merged_df <- merge(merged_df, mets_data, by = c('met', 'year'))

merged_df <- merged_df[order(merged_df$row_num),]

save(merged_df, file = paste0(site,'_full_data_for_gams.Rdata'))
}

### Building GAM models simplest to most complex ###
gamGS <- gam(agb ~ s(met, k = 3, bs = 're') + s(met, species, k = c(3,4), bs = 're'), control = ctrl, data = merged_df)
summary(gamGS)

# including gam with met smoothing, met by species smoothing, and smoothed interactions between
### This includes scaling up gam to test how long it will take to fit gam with all smoothers and also eventually all data per site


start_time <- Sys.time()

### Smoothers of all parameters by species

bam_all_params <- bam(agb ~ te(met, k = 3, bs = 're') + te(met, species, k = c(3,4), bs ='re') + s(DMAX, by = species) + s(DMIN, by = species) + s(MPLANT, by = species) + 
             s(D3, by = species) + s(FROST, by = species) + s(G,by = species) + s(AGEMX, by = species) + s(SPRTND, by = species), control = ctrl, data = merged_df)

end_time <- Sys.time()
total_run <- end_time-start_time
print(total_run)

summary(bam_all_params)

### Smoothers of all months of clim by species

bam_all_clim <- bam(agb ~ te(met, k = 3, bs = 're') + te(met, species, k = c(3,4), bs ='re') + s(precip_1, by = species) + s(precip_2, by = species) + s(precip_3, by = species) + s(precip_4, by = species) + s(precip_5, by = species) + s(precip_6, by = species) +
                                s(precip_7, by = species)+ s(precip_8, by = species) + s(precip_9, by = species) + s(precip_10, by = species) + s(precip_11, by = species) + s(precip_12, by = species) + s(temp_1, by = species) + s(temp_2, by = species) + s(temp_3, by = species) +
                                s(temp_4, by = species) + s(temp_5, by = species) + s(temp_6, by = species) + s(temp_7, by = species) + s(temp_8, by = species) + s(temp_9, by = species) + s(temp_10, by = species) + s(temp_11, by = species) + s(temp_12, by = species), discrete = TRUE, method = 'fREML',control = ctrl, data = merged_df)

summary(bam_all_clim)

start_time <- Sys.time()

bam_higher_order_effects <- bam(agb ~ s(met, by = species, k = c(3,4), bs = 're') + 
                                  # Higher order effects
                                  ti(DMAX, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(DMAX, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_12, temp_12, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(AGEMX, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(AGEMX, precip_12, temp_12, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(SPRTND, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(SPRTND, precip_12, temp_12, bs = c('cr','cr','cr')) +
                                  ti(G, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(G, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(G, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(G, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(G, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(G, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(G, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(G, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(G, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(G, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(G, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(G, precip_12, temp_12, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(D3, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(D3, precip_12, temp_12, bs = c('cr','cr','cr')) + 
                                  ti(D3, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(MPLANT, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(MPLANT, precip_12, temp_12, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(FROST, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(FROST, precip_12, temp_12, bs = c('cr','cr','cr')) + 
                                  ti(FROST, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(DMIN, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_4, temp_4, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_5, temp_5, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_6, temp_6, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_7, temp_7, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_8, temp_8, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_9, temp_9, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(DMIN, precip_12, temp_12, bs = c('cr','cr','cr')), method = 'fREML', data = merged_df, nthreads = 4, control = ctrl, discrete = TRUE)

end_time <- Sys.time()
total_run <- end_time-start_time
print(total_run)

### need to add in discrete == True and probably try to parallelize 

load('HF_full_data_for_gams.Rdata')

summary(bam_higher_order_effects, re.test = FALSE)
gam.check(bam_higher_order_effects)
qq_plot(bam_higher_order_effects)
pacf(residuals(bam_higher_order_effects))
bam_higher_order_effects$residuals
acf(residuals(bam_higher_order_effects))

plot(bam_higher_order_effects, shade = TRUE)



#### gamm with AR1 autocorrelation

### add ID column per individual time series

#merged_df$ID <- as.numeric(as.factor(with(merged_df,paste(species,run,param,met,sep='_'))))

gammGS_AR1 <- gamm(agb~ s(met, bs = 're') + te(met, species, k = c(3,4), bs ='re'), data = merged_df, correlation = corAR1(form = ~ year|ID))

### adding AR1 structure to model


#load(smaller_df)

merged_df$ID <- as.numeric(as.factor(with(merged_df,paste(species,run,param,met,sep='_'))))
ar1_vector <- c(ifelse(merged_df$year == 1, TRUE, FALSE))


bam_test_AR1 <- bam(agb ~ s(met, by = species, k = c(3,4), bs = 're') + 
                                  # Higher order effects
                                  ti(DMAX, precip_1, temp_1, bs = c('cr','cr','cr')) + 
                                  ti(DMAX, precip_2, temp_2, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_3, temp_3, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_10, temp_10, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_11, temp_11, bs = c('cr','cr','cr')) +
                                  ti(DMAX, precip_12, temp_12, bs = c('cr','cr','cr')), 
                                  method = 'fREML', data = merged_df, control = ctrl, discrete = TRUE, rho = .5, AR.start = ar1_vector)

summary(bam_test_AR1)
pacf(residuals(bam_test_AR1))
acf(residuals(bam_test_AR1))
gam.check(bam_test_AR1)


### Bam with higher order effects, met is global smoother and group specific smoother is species, there is also correlation strucutre for AR1 in this model


summary(test_species)
gam.check(test_species)

year_gam <- gam(agb ~ s(year) + s(year, met, k = c(5,3), bs = 're') + s(year, species, k =c(5,4), bs ='re'), data = merged_df)
summary(year_gam)
gam.check(year_gam)

species_gam <- gam(agb ~ s(species , k =4, bs = 're'), control = ctrl, data = merged_df)
gam.check(species_gam)


bam_I <- bam(agb ~ s(met, k =3, bs = 're') + ti(species, DMAX, precip_1, temp_1, bs = c('re','cr','cr','cr')), control = ctrl, data = merged_df)
summary(bam_I)
gam.check(bam_I)


bam_species <- bam(agb ~ s(met, by = species, k =3, bs = 're') + s(species, k = 4, bs = 're' ), data = merged_df)


ctrl <- gam.control(trace = TRUE)
ar1_vector <- c(ifelse(merged_df$year == 1, TRUE, FALSE))
start_time <- Sys.time()

bam_higher_order_effects_by_species <- bam(agb ~ species + s(met, k = 3, bs = 're') + 
                                    s(DMAX) + s(G) + s(D3) + s(FROST) + s(DMIN) + s(MPLANT) + s(SPRTND) + s(AGEMX) + 
                                    s(precip_1) + s(precip_2) + s(precip_3) + s(precip_4) + s(precip_5) + s(precip_6) + 
                                    s(precip_7) + s(precip_8) + s(precip_9) + s(precip_10) + s(precip_11) + s(precip_12) + 
                                    s(temp_1) + s(temp_2) + s(temp_3) + s(temp_4) + s(temp_5) + s(temp_6) +
                                    s(temp_7) + s(temp_8) + s(temp_9) + s(temp_10) + s(temp_11) + s(temp_12) +
                                    # Higher order effects
                                    ti(DMAX, temp_1, species, bs = c('cr','cr', 're')) + 
                                    ti(DMAX, temp_2, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_3, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_4, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_5, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_6, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX,temp_7, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_8, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_9, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_10, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_11, species, bs = c('cr','cr', 're')) +
                                    ti(DMAX, temp_12, species, bs = c('cr','cr', 're')) +
                                    ti(G, precip_1, temp_1, species, bs = c('cr','cr','cr','re')) + 
                                    ti(G, precip_2, temp_2, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_3, temp_3, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_4, temp_4, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_5, temp_5, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_6, temp_6, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_7, temp_7, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_8, temp_8, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_9, temp_9, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_10, temp_10, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_11, temp_11, species, bs = c('cr','cr','cr','re')) +
                                    ti(G, precip_12, temp_12, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_1, temp_1, species, bs = c('cr','cr','cr','re')) + 
                                    ti(D3, precip_2, temp_2, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_3, temp_3, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_4, temp_4, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_5, temp_5, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_6, temp_6, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_7, temp_7, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_8, temp_8, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_9, temp_9, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_10, temp_10, species, bs = c('cr','cr','cr','re')) +
                                    ti(D3, precip_12, temp_12, species, bs = c('cr','cr','cr','re')) + 
                                    ti(D3, precip_11, temp_11, species, bs = c('cr','cr','cr','re')) +
                                    ti(FROST, temp_1, species, bs = c('cr','cr', 're')) + 
                                    ti(FROST, temp_2, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_3, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_4, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_5, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_6, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_7, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_8, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_9, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_10, species, bs = c('cr','cr', 're')) +
                                    ti(FROST, temp_12, species, bs = c('cr','cr', 're')) + 
                                    ti(FROST, temp_11, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_1, species, bs = c('cr','cr', 're')) + 
                                    ti(DMIN, temp_2, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_3, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_4, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_5, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_6, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_7, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_8, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_9, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_10, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_11, species, bs = c('cr','cr', 're')) +
                                    ti(DMIN, temp_12, species, bs = c('cr','cr', 're')), method = 'fREML', data = merged_df, nthreads = 16, control = ctrl, discrete = TRUE, rho = .75, AR.start = ar1_vector)

end_time <- Sys.time()
run_time <- end_time - start_time
print(run_time)

save(bam_higher_order_effects_by_species, file = 'full_gam_HF.Rdata')




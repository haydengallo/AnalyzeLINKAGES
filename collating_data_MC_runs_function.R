### June 2023

##' @title Function for transforming final LINKAGES data into usable product
##' @author Hayden Gallo
##'
##' @param site                 specify site (HF, Rooster, Goose, Sylvania, NRP)
##' @param numruns              number of runs
##' @param nyear_user           number of years
##' @param nspec                number of species
##' @param parameters_wanted    list of parameters to populate arrays with, if all parameters wanted use 'all' as argument
##' @param mets_wanted          list of mets to populate arrays with, if all mets wanted use 'all' as argument
##' @param save_data            'yes' if user wants to save data as .Rdata product, 'no' otherwise

##' @description Collects and collates MC runs data and produces much more user friendly data product

##' @return master_array_pft
##' @return master_array_tot_biomass
##' @return master_array_varied_param
##' @return site
##' @return numruns
##' @return nyear_user
##' @return nspec


### This function is to be used when performing a monte carlo sensitivity analysis on LINKAGES using script called: MC_sims_production_runs.R

collating_MC_runs <- function(site, numruns, nyear_user, nspec, parameters_wanted, mets_wanted, save_data){
  
  nyear <- nyear_user
  
  all_parameters <- list('D3', 'MPLANT', 'DMAX', 'Frost', 'AGEMX', 'DMIN', 'G', 'SPRTND')
  
  load(paste0('MC_runs/',site,'/mets_used.Rdata'))
  met_runs_sample <- clim_use
  
  if(parameters_wanted == 'all'){params_to_use <- all_parameters}
  else{params_to_use <- match(parameters_wanted, all_parameters)}
  
  if(mets_wanted == 'all'){mets_to_use <- met_runs_sample}
  else{mets_to_use <- match(mets_wanted, met_runs_sample)}
  
  num_mets <- length(mets_to_use)
  
  # create save matrices
  avg_age <- array(0,c(nspec,nyear,numruns))
  avg_diam <- array(0,c(nspec,nyear,numruns))
  tot_agb <- matrix(0,nyear,numruns)
  spp_agb <- array(0,c(nspec,nyear,numruns))
  spp_trees <- array(0,c(nspec,nyear,numruns))
  spp_birth <- array(0,c(nspec,nyear,numruns))
  spp_death <- array(0,c(nspec,nyear,numruns))
  params <- array(0,c(nspec,26,numruns))
  varied_param <- matrix(0,nspec,numruns)
  spp_lgfval <- array(NA,c(nspec,nyear,numruns))
  spp_lgfname <- array(NA,c(nspec,nyear,numruns))
  spp_algf <- array(NA,c(200,nspec,nyear,numruns))
  master_array_tot_biomass <- array(0, c(nyear_user, numruns, length(all_parameters), num_mets))
  master_array_pft <-  array(0, c(nspec, nyear_user, numruns, length(all_parameters), num_mets, 6))
  master_array_varied_param <- array(0, c(nspec,numruns, length(all_parameters), num_mets))
  master_array_lgf <- array(NA, c(nspec, nyear_user, numruns, length(all_parameters), num_mets, 2))
  master_array_algf <- array(NA,dim=c(200,nspec,nyear,numruns,length(all_parameters),num_mets))
  dbh_master_array <- array(0, c(200, nyear,numruns, length(all_parameters), num_mets))
  master_params <- array(0, c(nspec, 26, numruns, length(all_parameters), num_mets))
  
  
  ### this loop, runs through all of the save folders for each parameter at the selected site and aggregates
  
  count_met = 0
  
  
  for (m in mets_to_use){
    met = m
    count_param = 0
    count_met = count_met + 1
    print(count_met)
    
    for (l in params_to_use){
      if (l == 'D3'){position = 14}
      else if (l == 'MPLANT'){position = 13}
      else if (l == 'DMAX'){position = 3}
      else if (l == 'Frost'){position = 15}
      else if (l == 'AGEMX'){position = 8}
      else if (l == 'DMIN'){position = 4}
      else if (l == 'G'){position = 9}
      else{position = 10}
      
      count_param = count_param + 1
      
      
      for (i in 1:numruns){
        
        # find directory with input and output file
        # this directory for use with changing met
        thisdir <- paste0('MC_runs/',site,'/',met,'/',l,'/',toString(i))
        load(paste0(thisdir,'/','linkages.out.Rdata'))
        load(paste0(thisdir,'/','linkages.input.Rdata'))
        
        # collect parameter information
        
        params[,,i] <- as.double(as.matrix(spp.params))
        varied_param[,i] <- spp.params[,position]
        dbh_master_array[,,i,count_param,count_met] <- dbh.save[,,1]
        master_params[,,i,count_param,count_met] <- as.double(as.matrix(spp.params))
        
        # collect annual data
        for (j in 1:nyear_user){
          st <- 1
          tot_agb[,i] <- ag.biomass[,1]
          
          
          # collect species-specific annual data
          for (k in 1:nspec){
            
            # determine indices of correct trees
            if (ntrees.kill[k,j,1]==0) next
            end <- st + ntrees.kill[k,j,1] - 1
            
            # store data
            spp_trees[k,j,i] <- ntrees.kill[k,j,1]
            spp_agb[k,j,i] <- agb.pft[k,j,1]
            spp_lgfval[k,j,i] <- min(gf.vec.save[k,,j,1],na.rm=TRUE)
            spp_lgfname[k,j,i] <- which.min(gf.vec.save[k,,j,1])
            spp_algf[,k,j,i] <- algf.save.keep[,k,j,1]
            if (j != 1){spp_birth[k,j,i] <- ntrees.birth[k,j,1]-ntrees.kill[k,(j-1),1]} else {spp_birth[k,j,i] = ntrees.birth[k,j,1]}
            spp_death[k,j,i] <- ntrees.birth[k,j,1]-ntrees.kill[k,j,1]
            avg_age[k,j,i] <- mean(iage.save[(st:end),j,1])
            avg_diam[k,j,i] <- mean(dbh.save[(st:end),j,1])
            st <- end + 1
          }
        }
      }
      
      
      master_array_pft[,,,count_param,count_met,1] <- spp_agb[,,]
      master_array_pft[,,,count_param,count_met,2] <- spp_birth[,,]
      master_array_pft[,,,count_param,count_met,3] <- spp_death[,,]
      master_array_pft[,,,count_param,count_met,4] <- avg_age[,,]
      master_array_pft[,,,count_param,count_met,5] <- avg_diam[,,]
      master_array_pft[,,,count_param,count_met,6] <- spp_trees[,,]
      
      master_array_tot_biomass[,,count_param,count_met] <- tot_agb[,]
      
      master_array_varied_param[,,count_param,count_met] <- varied_param[,]
      
      master_array_lgf[,,,count_param,count_met,1] <- spp_lgfval[,,]
      master_array_lgf[,,,count_param,count_met,2] <- spp_lgfname[,,]
      
      master_array_algf[,,,,count_param,count_met] <- spp_algf[,,,]
      
    }
    
    species_names <- spp.params$Spp_Name
    within_master_pft <- c('spp_agb', 'spp_birth', 'spp_death', 'avg_age', 'avg_diam', 'spp_trees')
    
    dimnames(master_array_pft) <- list(species_names, paste0('yr_',1:nyear_user), paste0('run_',1:numruns), params_to_use, mets_to_use, within_master_pft)
    dimnames(master_array_tot_biomass) <- list(paste0('yr_',1:nyear_user), paste0('run_',1:numruns), params_to_use, mets_to_use)
    dimnames(master_array_varied_param) <- list(species_names, paste0('run_',1:numruns), params_to_use, mets_to_use)
    dimnames(master_array_lgf) <-list(species_names, paste0('yr_',1:nyear_user), paste0('run_',1:numruns), params_to_use, mets_to_use, c('spp_lgfval','spp_lgfname'))
    dimnames(master_array_algf) <- list(paste0('tree_',1:200), species_names, paste0('yr_',1:nyear_user), paste0('run_',1:numruns), params_to_use, mets_to_use)
    dimnames(dbh_master_array) <- list(paste0('tree_',1:200), paste0('yr_',1:nyear_user),paste0('run_',1:numruns), params_to_use, mets_to_use)
    dimnames(master_params) <- list(species_names, colnames(spp.params), paste0('run_',1:numruns), params_to_use, mets_to_use)
    }
  
  concat_params <- paste(parameters_wanted, collapse = '_')
  concat_mets <- paste(mets_wanted, collapse = '_')
  
  if(save_data == 'yes'){
    
    save(master_array_pft = master_array_pft, master_array_tot_biomass = master_array_tot_biomass, master_array_varied_param= master_array_varied_param, master_array_lgf = master_array_lgf, master_array_algf = master_array_algf, dbh_master_array = dbh_master_array,
         site = site, numruns = numruns, nyear_user = nyear_user, nspec = nspec, parameters_wanted = parameters_wanted, mets_wanted = mets_wanted, master_params = master_params,
         file = paste0('MC_runs/',site,'/params_',concat_params,'_mets_',concat_mets,'.Rdata'))
    
  }
  
  return_list <- (list('master_array_pft' = master_array_pft, 'master_array_tot_biomass' = master_array_tot_biomass , 'master_array_varied_param' = master_array_varied_param, 'master_array_lgf' = master_array_lgf,
                       'master_array_algf' = master_array_algf, 'dbh_master_array' = dbh_master_array, 'master_params' = master_params, 'site' = site, 'numruns' = numruns, 'nyear_user' = nyear_user, 'nspec' = nspec))
  list2env(return_list, .GlobalEnv)
}



#collating_MC_runs(site = site, numruns = numruns, nyear_user = nyear_user, nspec = nspec, parameters_wanted = c('all'), mets_wanted = c('all'), save_data = 'yes')




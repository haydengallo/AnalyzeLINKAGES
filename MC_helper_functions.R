### Hayden Gallo
### May 2023
### Place to collect all helper functions for MC sims, spinup, testing etc.

running_linkages_MC_sims <- function(spp.params, parameter, runnum, met){

  # create new run folder
  # choose correct new_dir based on whether or not you have chosen to vary met
  if (vary_met == 'True'){new_dir = paste0('MC_runs/',site,'/',met,'/',parameter,'/',runnum)}
  else{new_dir = paste0('MC_runs/',site,'/',parameter,'/',runnum)}
  if (!dir.exists(new_dir)){
    dir.create(new_dir, recursive = T)
  }

  # save new input file with new parameters
  save(clat=clat,fdat=fdat,precip.mat=precip.mat,spp.params=spp.params,switch.mat=switch.mat,
       temp.mat=temp.mat,basesc=basesc,basesn=basesn,bgs=bgs,dry=dry,egs=egs,end.year=end.year,fc=fc,
       iplot=iplot,max.ind=max.ind,nspec=nspec,nyear=nyear,plat=plat,start.year=start.year,
       file=paste0(new_dir,'/','linkages.input.Rdata'))

  # run model and save output in ensemble run directory
  linkages(linkages.input = paste0(new_dir,'/','linkages.input.Rdata'),
           outdir = new_dir, spinup_input = paste0('MC_Sims_Spinup/',site,'/1/linkages.out.Rdata'))
}


monte_carlo_simulation <- function(parameter, numruns){
  start.time <<- Sys.time()
  count = 0
  count_mets = 0
  for (k in clim_use){
    met = k
    count_mets = count_mets + 1
    input = paste0('Met_updated','/',site,'_met','/linkages/',met,'.Rdata')
    load(input)
    nyear <<- nyear_user
    end.year <<- 2015
    start.year <<- end.year-nyear_user+1
    precip.mat <<- tail(precip.mat, n = nyear_user)
    temp.mat <<- tail(temp.mat, n = nyear_user)
    for (l in all_parameters){
      parameter = l
      count = count + 1
      if (l == 'D3'){position = 14}
      else if (l == 'MPLANT'){position = 13}
      else if (l == 'DMAX'){position = 3}
      else if (l == 'Frost'){position = 15}
      else if (l == 'AGEMX'){position = 8}
      else if (l == 'DMIN'){position = 4}
      else if (l == 'G'){position = 9}
      else if (l == 'SPRTND'){position = 10}
      temp_param <<- get(paste0(l,'_samples'))
      for (j in 1:numruns){
        if (vary_met == 'True'){new_dir = paste0('MC_runs/',site,'/',met,'/',parameter)}
        else{new_dir = paste0('MC_runs/',site,'/',parameter)}
        if (!dir.exists(new_dir)){dir.create(new_dir, recursive = T)}
        spp.params[,position] <<- temp_param[,j,count_mets]
        if (vary_met == 'True'){running_linkages_MC_sims(spp.params, parameter = parameter, runnum = j, met = met)}
        else {running_linkages_MC_sims(spp.params, parameter = parameter, runnum = j)}
      }
      spp.params <<- baseline_spp.params
      count = 0
    }
  }
  save(ens_wts = wts_use_1, file = paste0('MC_runs/',site,'/reweighted_weights.Rdata'))
  save(ens_wts = wts_use, file = paste0('MC_runs/',site,'/original_weights.Rdata'))
  save(clim_use, file = paste0('MC_runs/',site,'/mets_used.Rdata'))
  end.time <<- Sys.time()
}


running_linkages_spinup <- function(spp.params, runnum){

  new_dir = paste0('MC_Sims_Spinup/',site,'/',runnum)
  if (!dir.exists(new_dir)){
    dir.create(new_dir, recursive = T)
  }

  # save new input file with new parameters
  save(clat=clat,fdat=fdat,precip.mat=precip.mat,spp.params=spp.params,switch.mat=switch.mat,
       temp.mat=temp.mat,basesc=basesc,basesn=basesn,bgs=bgs,dry=dry,egs=egs,end.year=end.year,fc=fc,
       iplot=iplot,max.ind=max.ind,nspec=nspec,nyear=nyear,plat=plat,start.year=start.year,
       file=paste0(new_dir,'/','linkages.input.Rdata'))

  # run model and save output in ensemble run directory
  linkages(linkages.input = paste0(new_dir,'/','linkages.input.Rdata'),
           outdir = new_dir, spinup_input = NULL)
}


MC_sims_spinup <- function(spinup_len){

  temp.mat <<- rbind(temp.mat[rep(1:50, length = spinup_len-nyear), ], temp.mat)
  precip.mat <<- rbind(precip.mat[rep(1:50, length = spinup_len-nyear), ], precip.mat)
  nyear <<- nrow(temp.mat)

  new_dir = paste0('MC_Sims_Spinup/',site)
  if (!dir.exists(new_dir)){
    dir.create(new_dir, recursive = T)
  }
  running_linkages_spinup(spp.params, runnum = 1)

}







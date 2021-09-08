## Run the ZooMizer coupled model over multiple cores
##
## Adapted from Jason's ZooMSS code for UNSWs Katana
##
## Updated by Irene Richards Wed 8-Sep-2021

library(mizer)
library(tidyverse)
library(assertthat)

# Set Job Specifics

HPC <- FALSE # Select setup options for non-HPC runs
if(HPC == FALSE){
  ID <- 703 # Select a specific enviro row for local runs
  enviro <- readRDS("data/enviro_grid20210705.RDS")[ID,]
  # Generate unique job name
  jobname <- Sys.time() %>% format("%Y%m%d%H%M") %>% paste0("_grid",
                                                            "_sst-",enviro$sst,
                                                            "_chlo100-",
                                                            as.integer(enviro$chlo*100))
} else {
  jobname <- '20210823_grid' #job name used on queue
  ID <- as.integer(Sys.getenv('PBS_ARRAY_INDEX')) # Get the array run number on HPC
  enviro <- readRDS("data/enviro_grid20210705.RDS")[ID,]
}

Groups <- read.csv("data/TestGroups_mizer.csv") # Load in functional group information
ID_char <- sprintf("%04d",ID) # Sets the ID as a 4-digit character for sorting

# Choose environmental data to use (moved into if temporarily)
#enviro <- readRDS("data/enviro_grid20210705.RDS")[ID,]

source("uncoupledmodel.R")
source("ZooMizerResourceFunctions.R")

environment(new_project_simple) <- asNamespace('mizer') 
assignInNamespace("project_simple", new_project_simple, ns = "mizer")

# Run the ZooMSS model inside Mizer:
phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, kappa = 10^enviro$phyto_int[i], lambda= 1-enviro$phyto_slope[i], ... ) {
  npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
  npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  return(npp)
}

zoomss <- fZooMizer_run(groups = Groups, input = enviro)
saveRDS(zoomss, file = paste0("Output/", jobname, "_ZooMSS_", ID_char,".RDS"))



environment(project_simple) <- asNamespace('mizer')
assignInNamespace("project_simple", project_simple, ns = "mizer")

ZooMizer_coupled <- function(ID, tmax = 1000, effort = 0) {
  groups <- read.csv("data/TestGroups_mizer.csv") #load in ZooMSS groups
  zoo_groups <- groups[which(groups$Type == "Zooplankton"),] #just keep the zooplankton
  
  input <- readRDS("data/enviro_grid20210705.RDS")[ID,]
  
  stable_zoomizer <- readRDS(paste0("Output/", jobname, "_ZooMSS_", ID_char,".RDS"))
  
  times <- length(getTimes(stable_zoomizer))
  stable_zoo <- colSums(stable_zoomizer@n[ceiling(times/2+1):times,,])
  
  stable_zoomizer@n <- sweep(stable_zoomizer@n, 2, 2.5 * stable_zoomizer@params@species_params$Carbon, "*") #convert to carbon
  intercept <- mean(getCommunitySlope(stable_zoomizer,
                                      species = which(groups$Type == "Zooplankton"),
                                      max_w = 1e-5,
  )$intercept[ceiling(times/2+1):times])
  slope <- mean(getCommunitySlope(stable_zoomizer,
                                  species = which(groups$Type == "Zooplankton"),
                                  max_w = 1e-5,
  )$slope[ceiling(times/2+1):times])
  
  
  new <- stable_zoo[which(groups$Type == "Zooplankton"),]
  
  temp_eff <- 2.^((input$sst - 30)/10)
  
  fish_params <- newTraitParams(no_sp = 5,
                                min_w = 10^(-3),
                                min_w_inf = 10^3,
                                max_w_inf = 10^7,
                                no_w = (7-(-3))*10+1, # sets dx = 0.1 basically
                                min_w_pp = 10^(-14.5),
                                alpha = 0.25, # * temp_eff, #takes care of assimilation efficiency & temp effect on Encounter
                                kappa = exp(intercept), #10^(input$phyto_int) * 10^-1,
                                lambda = 1 - slope, # 1 - input$phyto_slope,
                                gamma = 1280 * temp_eff, #takes care of temp effect on PredRate and Encounter
                                # f0 = 0.6,
                                # h = 10^50,
                                R_factor = 1.01, #RMax = RFactor * RDI, default 4. Note RDI 
                                ks = 0, #set metabolism to zero
                                ext_mort_prop = 0, #currently zeroed since this is fitted. No need for temp effect here, calculates from PredMort
                                knife_edge_size = 10  # min size for fishing
  )
  
  # fish_params@species_params$gamma <- fish_params@species_params$gamma * temp_eff
  # fish_params <- setSearchVolume(fish_params)
  
  zoo_params <- newZooMizerParams(groups = zoo_groups, input = input, fish_params = fish_params)
  
  zoo_params@initial_n <- new[which(groups$Type == "Zooplankton"),]
  
  fish_params <- fish_params %>%
    setComponent(component = "zoo",
                 initial_value = zoo_params@initial_n,
                 dynamics_fun = "zoo_dynamics",
                 component_params = list(params = zoo_params)
    ) %>%
    setResource(resource_dynamics = "resource_zooMizer")
  
  # this line not needed if gamma is set to be temperature-dependent
  # fish_params <- setRateFunction(fish_params, "PredRate", "PredRate_temp")
  
  fish_params@other_params$temp_eff <- 2.^((input$sst - 30)/10)
  
  initialNResource(fish_params) <- resource_zooMizer(fish_params, fish_params@initial_n_other) #TODO: set this to be stable state of a regular ZooMizer run
  
  plotSpectra(fish_params, wlim = c(10^-14.5, NA), total = FALSE)
  
  initialN(fish_params) <- get_initial_n(fish_params)
  fish_params <- setParams(fish_params)
  
  return(project(fish_params, t_max = tmax, dt = 0.1, effort = effort))
}  

out <- ZooMizer_coupled(ID, tmax = 1000, effort = 0)
saveRDS(out, file = paste0("Output/", jobname, "_ZooMizer_", ID_char,".RDS"))

# Existing Resource Functions

# fZooMizer_run (uncoupledmodel.R[299]) is used to run ZooMSS within Mizer, 
# and it uses the phyto_fixed resource dynamics

# in fZooMizer_run, phytoplankton is classified thusly:
# min_w_pp = 10^(-14.5), #minimum phyto size. Note: use -14.4, not -14.5, otherwise it makes an extra size class
# w_pp_cutoff = 10^(input$phyto_max)* (1 + 1e-06), #maximum phyto size
# cc_phyto <- 0.1 #carbon content of phytoplankton



# uncoupledmodel & ZooMizerResourceFunctions
phyto_fixed <- function(params,
                        n,
                        n_pp, 
                        n_other, 
                        rates, 
                        dt, 
                        kappa=params@resource_params$kappa, 
                        lambda=params@resource_params$lambda, ... ) {
  #returns the fixed spectrum at every time step
  npp <- kappa*params@w_full^(1-lambda) / params@dw_full 
  npp[params@w_full > params@resource_params$w_pp_cutoff* (1 - 1e-06)] <- 0
  return(npp)
}
# ZooMizer setup file
phyto_fixed <- function(params, 
                        n, 
                        n_pp, 
                        n_other, 
                        rates, 
                        dt, 
                        kappa = 10^enviro$phyto_int[i], 
                        lambda= 1-enviro$phyto_slope[i], ... ) {
  npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
  npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  return(npp)
}

# Mizer
setResource() # Modify with resource_dynamics=newFunc for whatever I try
resource_semichemostat(params, n, n_pp, n_other, rates, t, dt, ...)

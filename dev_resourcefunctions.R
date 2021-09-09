# Existing Resource Functions

# phyto_fixed is used for phytoplankton, then zoo_dynamics for zooplankton and
# resource_zooMizer for the fish

# fZooMizer_run (uncoupledmodel.R[299]) is used to run ZooMSS within Mizer, 
# and it uses the phyto_fixed resource dynamics

# in fZooMizer_run, phytoplankton is classified thusly:
# min_w_pp = 10^(-14.5), #minimum phyto size. Note: use -14.4, not -14.5, otherwise it makes an extra size class
# w_pp_cutoff = 10^(input$phyto_max)* (1 + 1e-06), #maximum phyto size
# cc_phyto <- 0.1 #carbon content of phytoplankton



# Experimental Functions --------------------------------------------------

phyto_decrement <- function(params, n, n_pp, n_other, rates, dt, 
                        kappa = 10^enviro$phyto_int[i], 
                        lambda= 1-enviro$phyto_slope[i], ... )
{
  npp <- 0.99*kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
  npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  return(npp)
}
resource_dynamics = "phyto_decrement"

# ZooMizer Function -------------------------------------------------------

# Version found in uncoupledmodel.R & ZooMizerResourceFunctions.R
phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, 
                        kappa=params@resource_params$kappa, 
                        lambda=params@resource_params$lambda, ... )
{
  #returns the fixed spectrum at every time step
  npp <- kappa*params@w_full^(1-lambda) / params@dw_full 
  npp[params@w_full > params@resource_params$w_pp_cutoff* (1 - 1e-06)] <- 0
  return(npp)
}

# Version found in Setup_ZooMizer_grid.R and used within coupled ZooMSS
phyto_fixed <- function(params, n, n_pp, n_other, rates, dt, 
                        kappa = 10^enviro$phyto_int[i], 
                        lambda= 1-enviro$phyto_slope[i], ... )
{
  npp <- kappa*params@w_full^(-lambda) #returns the fixed spectrum at every time step
  npp[params@w_full > params@resource_params$w_pp_cutoff] <- 0
  return(npp)
}

# Function that actually runs ZooMSS with phyto_fixed
mf.params <- new_newMultispeciesParams(species_params=groups,
                                       interaction=NULL, #NULL sets all to 1, no strict herbivores
                                       min_w = 10^(-10.7),
                                       max_w = 10^7* (1 + 1e-06),
                                       no_w = no_w, #number of zoo+fish size classes;
                                       # w_full = 10^seq(from = -14.5, to = (log10(max(groups$w_inf)) + 0.1), by = 0.1),
                                       min_w_pp = 10^(-14.5), #minimum phyto size. Note: use -14.4, not -14.5, otherwise it makes an extra size class
                                       w_pp_cutoff = 10^(input$phyto_max)* (1 + 1e-06), #maximum phyto size
                                       n = 0.7, #The allometric growth exponent used in ZooMSS
                                       z0pre = 1, #external mortality (senescence)
                                       z0exp = 0.3,
                                       resource_dynamics = "phyto_fixed",
                                       kappa = kappa,
                                       lambda = lambda,
                                       RDD = constantRDD(species_params = groups), #first go at this
                                       #pred_kernel = ... #probably easiest to just import this/pre-calculate it, once dimensions are worked out
)
new_newMultispeciesParams <- function(
  species_params,
  interaction = NULL,
  no_w = 100,
  w_full = NA,         #added this line
  min_w = 0.001,
  max_w = NA,
  min_w_pp = NA,
  # setPredKernel()
  pred_kernel = NULL,
  # setSearchVolume()
  search_vol = NULL,
  # setMaxIntakeRate()
  intake_max = NULL,
  # setMetabolicRate()
  metab = NULL,
  p = 0.7,
  # setExtMort
  z0 = NULL,
  z0pre = 0.6,
  z0exp = n - 1,
  # setReproduction
  maturity = NULL,
  repro_prop = NULL,
  RDD = "BevertonHoltRDD",
  # setResource
  resource_rate = NULL,
  resource_capacity = NULL,
  n = 2 / 3,
  r_pp = 10,
  kappa = 1e11,
  lambda = 2.05,
  w_pp_cutoff = 10,
  resource_dynamics = "resource_semichemostat",
  # setFishing
  gear_params = data.frame(),
  selectivity = NULL,
  catchability = NULL,
  initial_effort = NULL) {
  no_sp <- nrow(species_params)
  
  ## For backwards compatibility, allow r_max instead of R_max
  if (!("R_max" %in% names(species_params)) &&
      "r_max" %in% names(species_params)) {
    names(species_params)[names(species_params) == "r_max"] <- "R_max"
  }
  
  ## Create MizerParams object ----
  params <- new_emptyParams(species_params,
                            gear_params,
                            w_full = w_full,
                            no_w = no_w,
                            min_w = min_w,
                            max_w = max_w,
                            min_w_pp = min_w_pp)
  
  ## Fill the slots ----
  params <- params %>%
    set_species_param_default("n", n) %>%
    set_species_param_default("p", p)
  params <- set_species_param_default(params, "q",
                                      lambda - 2 + params@species_params$n)
  if (is.null(interaction)) {
    interaction <- matrix(1, nrow = no_sp, ncol = no_sp)
  }
  params <-
    setParams(params,
              # setInteraction
              interaction = interaction,
              # setPredKernel()
              pred_kernel = pred_kernel,
              # setSearchVolume()
              search_vol = search_vol,
              # setMaxIntakeRate()
              intake_max = intake_max,
              # setMetabolicRate()
              metab = metab,
              # setExtMort
              z0 = z0,
              z0pre = z0pre,
              z0exp = z0exp,
              # setReproduction
              maturity = maturity,
              repro_prop = repro_prop,
              RDD = RDD,
              # setResource
              resource_rate = resource_rate,
              resource_capacity = resource_capacity,
              r_pp = r_pp,
              kappa = kappa,
              lambda = lambda,
              n = n,
              w_pp_cutoff = w_pp_cutoff,
              resource_dynamics = resource_dynamics,
              # setFishing
              gear_params = gear_params,
              selectivity = selectivity,
              catchability = catchability,
              initial_effort = initial_effort)
  
  params@initial_n <- get_initial_n(params)
  params@initial_n_pp <- params@cc_pp
  params@A <- rep(1, nrow(species_params))
  
  return(params)
}

# Mizer Functions ---------------------------------------------------------

# Modify with resource_dynamics=newFunc for whatever I try

setResource <- function (params, resource_rate = NULL, comment_rate = "set manually", 
                         resource_capacity = NULL, comment_capacity = "set manually", 
                         r_pp = resource_params(params)[["r_pp"]], kappa = resource_params(params)[["kappa"]], 
                         lambda = resource_params(params)[["lambda"]], n = resource_params(params)[["n"]], 
                         w_pp_cutoff = resource_params(params)[["w_pp_cutoff"]], 
                         resource_dynamics = NULL, ...) 
{
  assert_that(is(params, "MizerParams"), is.number(kappa), 
              kappa > 0, is.number(lambda), is.number(r_pp), r_pp > 
                0, is.number(w_pp_cutoff), is.number(n))
  params@resource_params[["kappa"]] <- kappa
  params@resource_params[["lambda"]] <- lambda
  params@resource_params[["r_pp"]] <- r_pp
  params@resource_params[["n"]] <- n
  params@resource_params[["w_pp_cutoff"]] <- w_pp_cutoff
  if (!is.null(resource_rate)) {
    if (is.null(comment(resource_rate))) {
      comment(resource_rate) <- comment_rate
    }
    assert_that(is.numeric(resource_rate), identical(length(resource_rate), 
                                                     length(params@rr_pp)))
    params@rr_pp[] <- resource_rate
    comment(params@rr_pp) <- comment(resource_rate)
  }
  else {
    rr_pp <- r_pp * params@w_full^(n - 1)
    if (!is.null(comment(params@rr_pp)) && different(params@rr_pp, 
                                                     rr_pp)) {
      message("The resource intrinsic growth rate has been commented and therefore will ", 
              "not be recalculated from the resource parameters.")
    }
    else {
      params@rr_pp[] <- rr_pp
    }
  }
  if (!is.null(resource_capacity)) {
    if (is.null(comment(resource_capacity))) {
      comment(resource_capacity) <- comment_capacity
    }
    assert_that(is.numeric(resource_capacity), identical(length(resource_capacity), 
                                                         length(params@cc_pp)))
    params@cc_pp[] <- resource_capacity
    comment(params@cc_pp) <- comment(resource_capacity)
  }
  else {
    cc_pp <- kappa * params@w_full^(-lambda)
    cc_pp[params@w_full > w_pp_cutoff] <- 0
    if (!is.null(comment(params@cc_pp)) && different(params@cc_pp, 
                                                     cc_pp)) {
      message("The resource carrying capacity has been commented and therefore will ", 
              "not be recalculated from the resource parameters.")
    }
    else {
      params@cc_pp[] <- cc_pp
    }
  }
  if (!is.null(resource_dynamics)) {
    assert_that(is.character(resource_dynamics))
    if (!is.function(get0(resource_dynamics))) {
      stop("The resource dynamics function \"", resource_dynamics, 
           "\" is not defined.")
    }
    params@resource_dynamics <- resource_dynamics
  }
  return(params)
}

resource_semichemostat <- function(params, n, n_pp, n_other, rates, t, dt, ...) 
{
  mur <- params@rr_pp + rates$resource_mort
  n_steady <- params@rr_pp * params@cc_pp/mur
  n_pp_new <- n_steady + (n_pp - n_steady) * exp(-mur * dt)
  sel <- mur == 0
  n_pp_new[sel] <- n_pp[sel]
  n_pp_new
}


# Modified fZooMizer_run for differing resource dynamics ------------------

fZooMizer_run <- function(groups, input, no_w = 178){
  
  kappa = 10^(input$phyto_int)
  lambda = 1-input$phyto_slope
  chlo = input$chlo
  sst = input$sst
  dt = input$dt
  tmax = input$tmax
  
  #data
  groups$w_min <- 10^groups$w_min #convert from log10 values
  groups$w_inf <- 10^groups$w_inf
  groups$w_mat <- 10^groups$w_mat
  groups$h <- 1e50 # should be Inf, but that breaks the calculations. Massive value still works out to effectively unlimited feeding as allowed in ZooMSS
  groups$ks <- 0 #turn off standard metabolism
  if (is.null(groups$interaction_resource)) {
    groups$interaction_resource <- 1
    groups$interaction_resource[which(groups$FeedType == "Carnivore")] <- 0
  }
  
  #todo - ramp up constant repro for coexistence
  
  mf.params <- new_newMultispeciesParams(species_params=groups,
                                         interaction=NULL, #NULL sets all to 1, no strict herbivores
                                         min_w = 10^(-10.7),
                                         max_w = 10^7* (1 + 1e-06),
                                         no_w = no_w, #number of zoo+fish size classes;
                                         # w_full = 10^seq(from = -14.5, to = (log10(max(groups$w_inf)) + 0.1), by = 0.1),
                                         min_w_pp = 10^(-14.5), #minimum phyto size. Note: use -14.4, not -14.5, otherwise it makes an extra size class
                                         w_pp_cutoff = 10^(input$phyto_max)* (1 + 1e-06), #maximum phyto size
                                         n = 0.7, #The allometric growth exponent used in ZooMSS
                                         z0pre = 1, #external mortality (senescence)
                                         z0exp = 0.3,
                                         resource_dynamics = resource_dynamics,
                                         kappa = kappa,
                                         lambda = lambda,
                                         RDD = constantRDD(species_params = groups), #first go at this
                                         #pred_kernel = ... #probably easiest to just import this/pre-calculate it, once dimensions are worked out
  )
  
  mf.params@species_params$w_min <- groups$w_min  #fix Mizer setting the egg weight to be one size larger for some groups.
  #mf.params@initial_n[] <- readRDS("data/initialn.RDS")
  
  temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))
  
  mf.params@other_params$assim_eff <- setassim_eff(groups)
  cc_phyto <- 0.1 #carbon content of phytoplankton
  mf.params@other_params$assim_phyto <- mf.params@species_params$GrossGEscale * cc_phyto #assimilation efficiency when eating phytoplankton
  
  mf.params@other_params$temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))
  
  mf.params <- setZooMizerConstants(params = mf.params, Groups = groups, sst = input$sst)
  #mf.params@initial_n[] <- readRDS("data/initialn.RDS")
  
  #mf.params <- setParams(mf.params)
  
  # mf.params <- setRateFunction(mf.params, "PredRate", "new_PredRate")
  mf.params <- setRateFunction(mf.params, "EReproAndGrowth", "new_EReproAndGrowth")
  mf.params <- setRateFunction(mf.params, "FeedingLevel", "newFeedingLevel")
  mf.params <- setRateFunction(mf.params, "Encounter", "new_Encounter")
  mf.params <- setRateFunction(mf.params, "PredRate", "new_PredRate")
  mf.params <- setReproduction(mf.params, repro_prop = matrix(0, nrow = nrow(mf.params@psi), ncol = ncol(mf.params@psi)))
  
  #mf.params@search_vol[] <- readRDS("data/SearchVol.rds")
  #mf.params <- setSearchVolume(mf.params)
  #mf.params <- setmort_test(mf.params, sst)
  # M_sb <- getExtMort(mf.params)
  # M_sb[] <- readRDS("data/mu_b.RDS")
  # temp_eff <-  matrix(2.^((sst - 30)/10), nrow = length(mf.params@species_params$species), ncol = length(mf.params@w))
  # M_sb <- temp_eff * M_sb *10 # Incorporate temp effect on senscence mortality
  
  # mf.params <- setExtMort(mf.params, z0 = M_sb)
  
  sim <- project(mf.params, dt = dt, t_max = tmax, t_save = 1) #TODO: make t_save an input to fZooMizer_run
  
  return(sim)
}

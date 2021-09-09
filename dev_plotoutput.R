source("ZooMizerPlots.R")
source("ZooMizerSummaryFunctions.R")

# Default phyto_fixed resource spectrum on grid cell 703
jobname <- "20210823_grid_sst-16_chlo100-398"
# Using phyto_decrement for 1% resource spectrum decrease each time step
jobname <- "202109091130_grid_sst-16_chlo100-398"

# Load RDS files and remove fish from ZooMSS output
zoomss <- readRDS(paste0("Output/", jobname, "_ZooMSS_", ID_char,".RDS"))
zoo_params <- removeSpecies(zoomss@params, c("Fish_Small","Fish_Med","Fish_Large"))
out <- readRDS(paste0("Output/", jobname, "_ZooMizer_", ID_char,".RDS"))



plot(out)
plot(zoomss)
plot(zoo_params)

# Plot the abundance spectra including ZooMizer resource
plotSpectra_ZooMizer(out, zoo_params) #check to see if there's a way to only plot Resource
plotSpectra_ZooMizer(out, zoo_params, wlim = c(1e-15, NA))
# Plot the abundance of species through time
plotAbundance_ZooMizer(out, zoo_params)
# Plot the biomass of species through time
plotBiomass_ZooMizer(out, zoo_params)


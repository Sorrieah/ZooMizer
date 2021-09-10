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

# Plot the abundance spectra including ZooMizer resource
plotSpectra_ZooMizer(out, zoo_params) #check to see if there's a way to only plot Resource
plotSpectra_ZooMizer(out, zoo_params, wlim = c(1e-15, NA))
# Plot the abundance of species through time
plotAbundance_ZooMizer(out, zoo_params)
# Plot the biomass of species through time
plotBiomass_ZooMizer(out, zoo_params)

# Attempt at making an auto-plot function
# dev_plot <- function(fish_object, 
#                      zoo_object, 
#                      jobname = jobname, 
#                      resource_function = resource_dynamics){
#   save_path <- paste0("Output/",jobname,"/")
#   zoo_params <- removeSpecies(zoo_object@params, c("Fish_Small","Fish_Med","Fish_Large"))
#   plot(fish_object)
#   ggsave(paste0(save_path,jobname,"_fish.png"), width = 862, height = 550, units = "px")
#   plot(zoo_object)
#   ggsave(paste0(save_path,jobname,"_zoo.png"), width = 862, height = 550, units = "px")
#   plotSpectra_ZooMizer(fish_object, zoo_params, wlim = c(1e-15, NA))
#   ggsave(paste0(save_path,jobname,"_spectra.png"), width = 862, height = 550, units = "px")
#   plotAbundance_ZooMizer(fish_object, zoo_params)
#   ggsave(paste0(save_path,jobname,"_abundance.png"), width = 862, height = 550, units = "px")
#   plotBiomass_ZooMizer(fish_object, zoo_params)
#   ggsave(paste0(save_path,jobname,"_biomass.png"), width = 862, height = 550, units = "px")
# }
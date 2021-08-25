#' Plot the abundance spectra including ZooMizer resource
#' 
#' Plots the number density multiplied by a power of the weight, with the power
#' specified by the `power` argument. This function plots the zooplankton
#' species from the ZooMizer resource.
#'
#' When called with a \linkS4class{MizerSim} object, the abundance is averaged
#' over the specified time range (a single value for the time range can be used
#' to plot a single time step). When called with a \linkS4class{MizerParams}
#' object the initial abundance is plotted.
#' 
#' @param object An object of class \linkS4class{MizerSim} or 
#'   \linkS4class{MizerParams}.
#' @inheritParams valid_species_arg
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step. Ignored when called with a \linkS4class{MizerParams}
#'   object.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the w axis. Use NA to refer to the existing minimum or maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param power The abundance is plotted as the number density times the weight
#' raised to `power`. The default \code{power = 1} gives the biomass 
#' density, whereas \code{power = 2} gives the biomass density with respect
#' to logarithmic size bins.
#' @param biomass `r lifecycle::badge("deprecated")`
#'  Only used if `power` argument is missing. Then
#'   \code{biomass = TRUE} is equivalent to \code{power=1} and 
#'   \code{biomass = FALSE} is equivalent to \code{power=0}
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Default is FALSE
#' @param resource A boolean value that determines whether resource is included.
#'   Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param ... Other arguments (currently unused)
#'   
#' @return A ggplot2 object
#' @export
#' @family plotting functions
#' @seealso [plotting_functions]
#' @examples
#' \donttest{
#' params <- suppressMessages(newMultispeciesParams(NS_species_params_gears, inter))
#' sim <- project(params, effort=1, t_max=20, t_save = 2, progress_bar = FALSE)
#' plotSpectra(sim)
#' plotSpectra(sim, wlim = c(1e-6, NA))
#' plotSpectra(sim, time_range = 10:20)
#' plotSpectra(sim, time_range = 10:20, power = 0)
#' plotSpectra(sim, species = c("Cod", "Herring"), power = 1)
#' }
plotSpectra_ZooMizer <- function(fish_object, zoo_object,
                        species = NULL,
                        time_range,
                        wlim = c(NA, NA), ylim = c(NA, NA),
                        power = 1, biomass = TRUE,
                        total = FALSE, resource = TRUE, 
                        background = TRUE,
                        highlight = NULL, ...) {
  # to deal with old-type biomass argument
  if (missing(power)) {
    power <- as.numeric(biomass)
  }
  species <- valid_species_arg(fish_object, species)
  if (is(fish_object, "MizerSim")) {
    if (missing(time_range)) {
      time_range  <- max(as.numeric(dimnames(fish_object@n)$time))
    }
    time_elements <- get_time_elements(fish_object, time_range)
    n <- apply(fish_object@n[time_elements, , , drop = FALSE], c(2, 3), mean)
    n_pp <- zoo_object@initial_n_pp
    n_zoo <- aperm(simplify2array(fish_object@n_other), c(3,1,2))
    n_zoo <- apply(n_zoo[time_elements, , , drop = FALSE], c(2, 3), mean)
    ps <- plot_spectra_zoomizer(fish_object@params, zoo_object,
                       n = n, n_pp = n_pp, n_zoo = n_zoo,
                       species = species, wlim = wlim,
                       ylim = ylim, power = power,
                       total = total, resource = resource,
                       background = background, highlight = highlight)
    return(ps)
  } else {
    ps <- plot_spectra_zoomizer(fish_object, zoo_object,
                       n = fish_object@initial_n,
                       n_zoo = zoo_object@initial_n,
                       n_pp = zoo_object@initial_n_pp,
                       species = species, wlim = wlim, ylim = ylim,
                       power = power,
                       total = total, resource = resource,
                       background = background, highlight = highlight)
    return(ps)
  }
}


plot_spectra_zoomizer <- function(fish_params, zoo_params, n, n_pp, n_zoo,
                         species, wlim, ylim, power,
                         total, resource, background, highlight) {
  fish_params <- validParams(fish_params)
  if (is.na(wlim[1])) {
    wlim[1] <- min(fish_params@w) / 100
  }
  if (is.na(wlim[2])) {
    wlim[2] <- max(fish_params@w_full)
  }
  # Need to keep species in order for legend
  species_levels <- c(dimnames(fish_params@initial_n)$sp,
                      "Background", dimnames(zoo_params@initial_n)$sp, "Total")
  fish_idx <- (length(fish_params@w_full) - length(fish_params@w) + 1):length(fish_params@w_full)
  zoo_idx <- (length(zoo_params@w_full) - length(zoo_params@w) + 1):length(zoo_params@w_full)
  if (total) {
    # Calculate total community abundance
    total_n <- n_pp
    total_n[zoo_idx] <- total_n[zoo_idx] + colSums(n_zoo)
    total_n[fish_idx] <- total_n[fish_idx] + colSums(n)
    total_n <- total_n * zoo_params@w_full^power
  }
  species <- c(valid_species_arg(zoo_params), valid_species_arg(fish_params, species))
  # Deal with power argument
  if (power %in% c(0, 1, 2)) {
    y_label <- c("Number density [1/g]", "Biomass density",
                 "Biomass density [g]")[power + 1]
  } else {
    y_label <- paste0("Number density * w^", power)
  }
  n <- sweep(n, 2, fish_params@w^power, "*")
  n_zoo <- sweep(n_zoo, 2, zoo_params@w^power, "*")
  n_mat <- matrix(0, nrow = length(species), ncol = length(zoo_params@w_full))
  n_mat[1:nrow(n_zoo), zoo_idx] <- n_zoo
  n_mat[(nrow(n_zoo)+1):length(species), fish_idx] <- n
  rownames(n_mat) <- species
  colnames(n_mat) <- zoo_params@w_full
  # Select only the desired species
  spec_n <- n_mat[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
  
  # Make data.frame for plot
  plot_dat <- data.frame(value = c(spec_n),
                         # ordering of factor is important for legend
                         Species = factor(dimnames(spec_n)[[1]],
                                          levels = species_levels),
                         w = rep(zoo_params@w_full,
                                 each = dim(spec_n)[[1]]))
  
  if (resource) {
    resource_sel <- (zoo_params@w_full >= wlim[1]) & 
      (zoo_params@w_full <= wlim[2])
    # Do we have any resource to plot?
    if (sum(resource_sel) > 0) {
      w_resource <- zoo_params@w_full[resource_sel]
      plank_n <- n_pp[resource_sel] * w_resource^power
      plot_dat <- rbind(plot_dat,
                        data.frame(value = c(plank_n),
                                   Species = "Resource",
                                   w = w_resource))
    }
  }
  if (total) {
    plot_dat <- rbind(plot_dat,
                      data.frame(value = c(total_n),
                                 Species = "Total",
                                 w = zoo_params@w_full))
  }
  # lop off 0s and apply wlim
  plot_dat <- plot_dat[(plot_dat$value > 0) & 
                         (plot_dat$w >= wlim[1]) &
                         (plot_dat$w <= wlim[2]), ]
  # Impose ylim
  if (!is.na(ylim[2])) {
    plot_dat <- plot_dat[plot_dat$value <= ylim[2], ]
  }
  if (is.na(ylim[1])) {
    ylim[1] <- 1e-20
  }
  plot_dat <- plot_dat[plot_dat$value > ylim[1], ]
  
  linecolours <- zoo_params@species_params$PlotColour
  names(linecolours) <- zoo_params@species_params$species
  linecolours <- c(linecolours, fish_params@linecolour)
  
  # Create plot
  

  p <- ggplot(plot_dat, aes(x = w, y = value)) +
    scale_x_continuous(name = "Size [g]", trans = "log10",
                       breaks = log_breaks()) +
    scale_y_continuous(name = y_label, trans = "log10",
                       breaks = log_breaks()) +
    scale_colour_manual(values = linecolours) +
    scale_linetype_manual(values = c(zoo_params@linetype[1:nrow(zoo_params@species_params)], fish_params@linetype))
  if (background) {
    back_n <- n[is.na(fish_params@A), , drop = FALSE]
    plot_back <- data.frame(value = c(back_n),
                            Species = as.factor(dimnames(back_n)[[1]]),
                            w = rep(zoo_params@w_full,
                                    each = dim(back_n)[[1]]))
    # lop off 0s and apply wlim
    plot_back <- plot_back[(plot_back$value > 0) & 
                             (plot_back$w >= wlim[1]) &
                             (plot_back$w <= wlim[2]), ]
    # Impose ylim
    if (!is.na(ylim[2])) {
      plot_back <- plot_back[plot_back$value <= ylim[2], ]
    }
    plot_back <- plot_back[plot_back$value > ylim[1], ]
    if (nrow(plot_back) > 0) {
      # Add background species
      p <- p +
        geom_line(aes(group = Species),
                  colour = fish_params@linecolour["Background"],
                  linetype = fish_params@linetype["Background"],
                  data = plot_back)
    }
  }
  linesize <- rep(0.8, length(c(zoo_params@linetype, fish_params@linetype))-4)
  names(linesize) <- c(zoo_params@species_params$species, names(fish_params@linetype))
  linesize[highlight] <- 1.6
  p <- p + scale_size_manual(values = linesize) + 
    geom_line(aes(colour = Species, linetype = Species, size = Species))
  return(p)
}

#' @rdname plotSpectra_ZooMizer
#' @export
plotlySpectra_ZooMizer <- function(fish_object, zoo_object, species = NULL,
                          time_range,
                          wlim = c(NA, NA), ylim = c(NA, NA),
                          power = 1, biomass = TRUE,
                          total = FALSE, resource = TRUE, 
                          background = TRUE,
                          highlight = NULL, ...) {
  argg <- as.list(environment())
  ggplotly(do.call("plotSpectra_ZooMizer", argg))
}


getBiomassFrame_ZooMizer <- function (sim, zoo_params, species = NULL, start_time = as.numeric(dimnames(sim@n)[[1]][1]), 
          end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]), 
          ylim = c(NA, NA), total = FALSE, ...) 
{
  species <- valid_species_arg(sim@params, species)
  b <- getBiomass_ZooMizer(sim, zoo_params, ...)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & (as.numeric(dimnames(b)[[1]]) <= 
                                                           end_time), , drop = FALSE]
  b_total <- rowSums(b)
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  bm <- reshape2::melt(b)
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value & (is.na(ylim[1]) | bm$value >= 
                                      ylim[1]) & (is.na(ylim[2]) | bm$value <= ylim[2]), ]
  names(bm) <- c("Year", "Species", "Biomass")
  species_levels <- c(zoo_params@species_params$species, dimnames(sim@n)$sp, "Background", 
                      "Resource", "Total")
  bm$Species <- factor(bm$Species, levels = species_levels)
  # bm <- bm[bm$Species %in% species, ]
  return(bm)
}

#' Plot the biomass of species through time
#'
#' @param sim An object of class MizerSim contianing the fish model
#' @param zoo_params An object of class MizerParams containing the zooplankton model
#' @param species A list of fish species to include in the plot
#' @param start_time 
#' @param end_time 
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param ylim 
#' @param total 
#' @param background 
#' @param highlight 
#' @param ... 
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
plotBiomass_ZooMizer <- function (sim, zoo_params, species = NULL, start_time, end_time, y_ticks = 6, 
          ylim = c(NA, NA), total = FALSE, background = FALSE, highlight = NULL, 
          ...) 
{
  species <- valid_species_arg(sim, species)
  if (missing(start_time)) 
    start_time <- as.numeric(dimnames(sim@n)[[1]][1])
  if (missing(end_time)) 
    end_time <- as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
  bm <- getBiomassFrame_ZooMizer(sim, zoo_params, species = dimnames(sim@n)$sp, 
                        start_time = start_time, end_time = end_time, ylim = ylim, 
                        total = total)#, ...)
  spec_bm <- bm[bm$Species %in% c("Total", zoo_params@species_params$species, species), ]
  x_label <- "Year"
  y_label <- "Biomass [g]"
  linecolours <- zoo_params@species_params$PlotColour
  names(linecolours) <- zoo_params@species_params$species
  linecolours <- c(linecolours, sim@params@linecolour)
  p <- ggplot(spec_bm, aes(x = Year, y = Biomass)) +
    scale_y_continuous(trans = "log10", 
                       breaks = mizer:::log_breaks(n = y_ticks),
                       labels = prettyNum, 
                       name = y_label) +
    scale_x_continuous(name = x_label) + 
    scale_colour_manual(values = linecolours)
    # scale_linetype_manual(values = sim@params@linetype)
  if (background) {
    back_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
    back_bm <- bm[bm$Species %in% back_sp, ]
    if (nrow(back_bm) > 0) {
      p <- p + geom_line(aes(group = Species), data = back_bm, 
                         colour = sim@params@linecolour["Background"], 
                         linetype = sim@params@linetype["Background"])
    }
  }
  linesize <- rep(0.8, length(sim@params@linetype))
  names(linesize) <- names(sim@params@linetype)
  linesize[highlight] <- 1.6
  p <- p + # scale_size_manual(values = linesize) +
         geom_line(aes(colour = Species)) #, linetype = Species)) #, size = Species))
  return(p)
}

plotlyBiomass_ZooMizer <- function (sim, zoo_params, species = NULL, start_time, end_time, y_ticks = 6, 
          ylim = c(NA, NA), total = FALSE, background = TRUE, highlight = NULL, 
          ...) 
{
  argg <- c(as.list(environment()), list(...))
  ggplotly(do.call("plotBiomass_ZooMizer", argg))
}

getAbundanceFrame_ZooMizer <- function (sim, zoo_params, species = NULL, start_time = as.numeric(dimnames(sim@n)[[1]][1]), 
                                      end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]), 
                                      ylim = c(NA, NA), total = FALSE, ...) 
{
  species <- valid_species_arg(sim@params, species)
  b <- getAbundance_ZooMizer(sim, zoo_params, ...)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) & (as.numeric(dimnames(b)[[1]]) <= 
                                                           end_time), , drop = FALSE]
  b_total <- rowSums(b)
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  bm <- reshape2::melt(b)
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value & (is.na(ylim[1]) | bm$value >= 
                                      ylim[1]) & (is.na(ylim[2]) | bm$value <= ylim[2]), ]
  names(bm) <- c("Year", "Species", "Biomass")
  species_levels <- c(zoo_params@species_params$species, dimnames(sim@n)$sp, "Background", 
                      "Resource", "Total")
  bm$Species <- factor(bm$Species, levels = species_levels)
  # bm <- bm[bm$Species %in% species, ]
  return(bm)
}

#' Plot the abundance of species through time
#'
#' @param sim An object of class MizerSim containing the fish model
#' @param zoo_params An object of class MizerParams containing the zooplankton model
#' @param species A list of fish species to include in the plot
#' @param start_time 
#' @param end_time 
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param ylim 
#' @param total 
#' @param background 
#' @param highlight 
#' @param ... 
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
plotAbundance_ZooMizer <- function (sim, zoo_params, species = NULL, start_time, end_time, y_ticks = 6, 
                                  ylim = c(NA, NA), total = FALSE, background = FALSE, highlight = NULL, 
                                  ...) 
{
  species <- valid_species_arg(sim, species)
  if (missing(start_time)) 
    start_time <- as.numeric(dimnames(sim@n)[[1]][1])
  if (missing(end_time)) 
    end_time <- as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
  bm <- getAbundanceFrame_ZooMizer(sim, zoo_params, species = dimnames(sim@n)$sp, 
                                 start_time = start_time, end_time = end_time, ylim = ylim, 
                                 total = total)#, ...)
  spec_bm <- bm[bm$Species %in% c("Total", zoo_params@species_params$species, species), ]
  x_label <- "Year"
  y_label <- "Abundance"
  linecolours <- zoo_params@species_params$PlotColour
  names(linecolours) <- zoo_params@species_params$species
  linecolours <- c(linecolours, sim@params@linecolour)
  p <- ggplot(spec_bm, aes(x = Year, y = Biomass)) +
    scale_y_continuous(trans = "log10", 
                       breaks = mizer:::log_breaks(n = y_ticks),
                       labels = prettyNum, 
                       name = y_label) +
    scale_x_continuous(name = x_label) + 
    scale_colour_manual(values = linecolours)
  # scale_linetype_manual(values = sim@params@linetype)
  if (background) {
    back_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
    back_bm <- bm[bm$Species %in% back_sp, ]
    if (nrow(back_bm) > 0) {
      p <- p + geom_line(aes(group = Species), data = back_bm, 
                         colour = sim@params@linecolour["Background"], 
                         linetype = sim@params@linetype["Background"])
    }
  }
  linesize <- rep(0.8, length(sim@params@linetype))
  names(linesize) <- names(sim@params@linetype)
  linesize[highlight] <- 1.6
  p <- p + # scale_size_manual(values = linesize) +
    geom_line(aes(colour = Species)) #, linetype = Species)) #, size = Species))
  return(p)
}

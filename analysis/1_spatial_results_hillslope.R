# Date: March 2023
# Author: Selina Baldauf
# Purpose: Reproduce Figures 3, S5-8 with the spatial hillslope results

library(data.table)
library(tidyverse)
library(ggspatial)
library(terra)
library(patchwork)
source(here::here("analysis/helper_functions/read_spatial.R"))

# Reading and preparing the spatial results -------------------------------

# Mask to cut out the hillslope cells from the square model results
the_mask <- terra::rast(
  as.matrix(
    data.table::fread(here::here("model/Parameters/outputMask_74.txt"))
  )
)

spatial_results <- read_spatial(
  type = "soil", time = "Months",
  result_path = here::here("data/model_results"),
  result_folder = "Results_tabernas",
  xsize = 74
)

spatial_results <- lapply(spatial_results, function(x) {
  r <- lapply(x, function(y) {
    # turn result into raster
    to_return <- terra::rast(as.matrix(y))
    # crop raster by the mask
    to_return <- mask(to_return, the_mask, maskvalue = 0)
  })
  names(r) <- month.abb[as.integer(stringr::str_replace(
    stringr::str_replace(names(r), "soil_m", ""),
    ".txt", ""
  ))]
  return(r)
})


# Creating all spatial plots  ---------------------------------------------

# Create plots for all processes and months in a loop, combine them into one
# plot and save the result in a plot list

# Define colors and create empty list to hold the plots
plot_list <- list()
moisture_col <- rev(paletteer::paletteer_c("grDevices::Blues 3", 30))
evap_col <- paletteer::paletteer_c("ggthemes::Classic Green", 30)

# Sequentially create plots
for (j in c("Jan", "Feb", "Mar", "Dec")) {
  L1 <- spatial_results[["L1"]][[j]] |>
    df_spatial() |>
    mutate(band1 = band1 * 100) |>
    filter(!(is.na(band1))) |>
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = band1)) +
    labs(fill = "Soil moisture upper (%)") +
    scale_fill_gradientn(colors = moisture_col, na.value = NA)

  L2 <- spatial_results[["L2"]][[j]] |>
    df_spatial() |>
    mutate(band1 = band1 * 100) |>
    filter(!(is.na(band1))) |>
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = band1)) +
    labs(fill = "Soil moisture lower (%)") +
    scale_fill_gradientn(colors = moisture_col, na.value = NA)

  qd <- spatial_results[["QD"]][[j]] |>
    df_spatial() |>
    filter(!(is.na(band1))) |>
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = band1)) +
    labs(fill = "Runoff (mm)") +
    scale_fill_gradientn(
      colors = moisture_col, na.value = NA, trans = "pseudo_log",
      breaks = scales::pretty_breaks(n = 3)
    )

  drain <- spatial_results[["deepdrain"]][[j]] |>
    df_spatial() |>
    filter(!(is.na(band1))) |>
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = band1)) +
    labs(fill = "Deep drainage (mm)") +
    scale_fill_gradientn(
      colors = moisture_col, na.value = NA, trans = "pseudo_log",
      breaks = scales::pretty_breaks(n = 2)
    )

  waterGain <- spatial_results[["waterGain"]][[j]] |>
    df_spatial() |>
    filter(!(is.na(band1))) |>
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = band1)) +
    labs(fill = "Water gain (mm)") +
    scale_fill_gradient2(
      low = "#C84528",
      high = "#077F82", na.value = NA, trans = "pseudo_log",
      breaks = scales::pretty_breaks(n = 2)
    )
  # scale_fill_gradient2(na.value = NA, trans = "pseudo_log")

  eptot <- spatial_results[["EPtot"]][[j]] |>
    df_spatial() |>
    filter(!(is.na(band1))) |>
    ggplot() +
    geom_raster(aes(x = x, y = y, fill = band1)) +
    labs(fill = "Evapotranspiration (mm)") +
    scale_fill_gradientn(colors = evap_col, na.value = NA, trans = "pseudo_log")

  plot <-
    (L1 +
      L2 +
      qd +
      waterGain +
      drain +
      eptot) &
      guides(
        fill = guide_colorbar(
          direction = "horizontal",
          title.position = "top"
        )
      ) &
      theme_void() &
      theme(
        text = element_text(size = 10),
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm")
      )
  plot_list[[j]] <- plot
}

# Look at the plots -------------------------------------------------------

# Plots can be accessed like this plot_list[[month]], e.g.

plot_list[["Jan"]]

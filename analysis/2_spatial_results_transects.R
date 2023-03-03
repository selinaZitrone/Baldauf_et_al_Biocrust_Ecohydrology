# Date: March 2023
# Author: Selina Baldauf
# Purpose: Reproduce Figures 4, S10 with the spatial transect results

library(data.table)
library(lubridate)
library(tidyverse)
library(ggspatial)
library(terra)
library(patchwork)

path_results <- here::here("data/model_results")
source("analysis/helper_functions/read_spatial.R")

# Custom theme for the plots ----------------------------------------------

theme_transect <- function() {
  theme_bw(base_size = 10) %+replace%
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      plot.margin = margin(),
      panel.border = element_rect(color = "black", size = 0.1, fill = NA),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.key.height = unit(0.3, "cm"),
      legend.position = "bottom",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9),
      axis.text.x = element_text(size = 7),
      rect = element_rect(fill = "transparent")
    )
}

# Read and summarize rainfall data --------------------------------------------
dates <- seq(
  from = parse_date_time("2010-1-1 0", orders = "%Y%m%d %H"),
  to = parse_date_time("2010-12-31 23", orders = "%Y%m%d %H"),
  by = "hours"
)
rainfall_h <- fread(paste0(
  path_results,
  "/Results_tabernas/tabernas_1years_hourly_tabernas.txt"
))[
  order(year, day, hour), datetime := dates
][, .(rain, datetime)]
rainfall_m <- rainfall_h[, .(rain = sum(rain)), by = month(datetime)][, round(rain, 0)]
rainfall_m <- tibble(month_int = 1:12, month = month.abb, rain = rainfall_m)

# Read and prepare the spatial results ------------------------------------

location <- rast("data/spatial/location.tif")
transect <- vect("data/spatial/transect.shp")
# Transform transect and hillslope into correct coordinate system
transect <- project(transect, location)

scenarios <- c("tabernas", "tabernas_halfrain", "zeroCrust", "zeroCrust_halfrain")
cross_section_list <- list()
for (i in scenarios) {
  soil_m <- read_spatial(
    type = "soil", time = "Months",
    result_path = path_results,
    result_folder = paste0("Results_", i),
    xsize = 74
  )

  # make soil a raster
  soil_m <- lapply(soil_m, function(y) {
    r <- lapply(y, function(z) {
      rast(as.matrix(z), crs = crs(location), extent = ext(location))
    })
    names(r) <- month.abb[as.integer(str_replace(
      str_replace(names(r), "soil_m", ""),
      ".txt", ""
    ))]
    return(r)
  })
  # extract only the transect from the hillslope cells
  soil_crossection <- lapply(soil_m, function(x) {
    lapply(x, terra::extract, transect, cells = TRUE, touches = TRUE)
  })

  cross_section_list[[i]] <- soil_crossection
}

# Get tables to plot ------------------------------------------------------
# Take spatial transect results and summarize the results in tables that can
# be used for plotting

soil_transect <- unlist(cross_section_list, recursive = FALSE) |>
  unlist(recursive = FALSE)
soil_transect <- rbindlist(soil_transect, idcol = "ident")
soil_transect <- soil_transect |> tidyr::separate(ident, c("scenario", "type", "month"),
  sep = "\\.",
  fill = "right"
)
# subset the soil transect by id
soil_transect <- soil_transect |> filter(ID == 1)

transect_to_plot <- soil_transect |>
  group_by(scenario, type, month) |>
  mutate(
    order = as.factor(n():1),
    colID = 1
  ) |>
  ungroup() |>
  rowwise() |>
  mutate(month_int = which(month.abb == month)) |>
  ungroup() |>
  left_join(rainfall_m, by = c("month", "month_int")) |>
  mutate(type = case_when(
    type == "L1" ~ "Moisture upper (%)",
    type == "L2" ~ "Moisture lower (%)",
    type == "waterGain" ~ "Water gain (mm)",
    type == "deepdrain" ~ "Deep drainage (mm)",
    type == "EPtot" ~ "Evapotranspiration (mm)",
    type == "QD" ~ "Runoff (mm)",
    TRUE ~ type
  ))

# calculate delta of crust and no crust for the supplement
transect_to_plot <- transect_to_plot |>
  separate(scenario, into = c("scenario", "rain_sc")) |>
  mutate(rain_sc = ifelse(is.na(rain_sc), "default", rain_sc)) |>
  pivot_wider(names_from = scenario, values_from = lyr.1) |>
  mutate(delta = tabernas - zeroCrust) |>
  mutate(type = factor(type, levels = c(
    "Moisture upper (%)", "Moisture lower (%)", "FL1", "Runoff (mm)",
    "runon", "Water gain (mm)", "Deep drainage (mm)", "Evapotranspiration (mm)",
    "crust_moisture"
  )))

# start plot in June
transect_to_plot <- mutate(transect_to_plot,
  month_int2 =
    ifelse(month_int >= 6, month_int - 5, month_int + 7)
)

# Create single transect plots ------------------------------------------------

# Plots for absolute values and delta between crust-no crust comparison

moisture_col <- rev(paletteer::paletteer_c("grDevices::Blues 3", 30))
evap_col <- paletteer::paletteer_c("ggthemes::Classic Green", 30)
rain_multiplier <- 5

rain_scen <- "halfrain"
rain_color <- "deeppink"
rain_width <- 0.4

moisture_upper <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Moisture upper (%)") %>%
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = tabernas),
    width = 0.8
  ) +
  scale_fill_gradientn(colors = moisture_col, na.value = NA, guide = guide_colorbar(
    title = "Moisture upper\n(%)",
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

moisture_lower <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Moisture lower (%)") %>%
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = tabernas),
    width = 0.8
  ) +
  scale_fill_gradientn(
    colors = moisture_col, na.value = NA,
    guide = guide_colorbar(
      title = "Moisture lower\n(%)",
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()


moisture_upper_delta <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Moisture upper (%)") %>%
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = delta),
    width = 0.8
  ) +
  scale_fill_gradient2(
    guide = guide_colorbar(
      title = expression("\u0394 Moisture upper\n(%)"),
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

moisture_lower_delta <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Moisture lower (%)") %>%
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = delta),
    width = 0.8
  ) +
  scale_fill_gradient2(
    guide = guide_colorbar(
      title = expression("\u0394 Moisture lower\n(%)"),
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()


moisture <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type %in% c("Moisture upper (%)", "Moisture lower (%)")) |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = tabernas),
    width = 0.8
  ) +
  facet_wrap(~type, ncol = 1) +
  labs(fill = "soil moisture\n[%]") +
  scale_fill_gradientn(colors = moisture_col, na.value = NA) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

moisture_delta <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type %in% c("Moisture upper (%)", "Moisture lower (%)")) |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = delta),
    width = 0.8
  ) +
  facet_wrap(~type, ncol = 1) +
  labs(fill = "soil moisture\n[%]") +
  scale_fill_gradient2() +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

waterGain <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Water gain (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = tabernas),
    width = 0.8
  ) +
  scale_fill_gradient2(guide = guide_colorbar(
    title = "Water gain\n(mm)",
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

waterGain_delta <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Water gain (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = delta),
    width = 0.8
  ) +
  scale_fill_gradient2(guide = guide_colorbar(
    title = expression("\u0394 Water gain\n(mm)"),
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

runoff <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Runoff (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = tabernas),
    width = 0.8
  ) +
  scale_fill_gradient2(guide = guide_colorbar(
    title = "Runoff\n(mm)",
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

runoff_delta <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Runoff (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = delta),
    width = 0.8
  ) +
  scale_fill_gradient2(guide = guide_colorbar(
    title = expression("\u0394 Runoff\n(mm)"),
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

deepdrain <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Deep drainage (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = tabernas),
    width = 0.8
  ) +
  scale_fill_gradientn(colors = moisture_col, na.value = NA, guide = guide_colorbar(
    title = "Deep drainage\n(mm)",
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain [mm]",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

deepdrain_delta <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Deep drainage (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = delta),
    width = 0.8
  ) +
  scale_fill_gradient2(guide = guide_colorbar(
    title = expression("\u0394 Deep drainage\n(mm)"),
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier), color = rain_color, size = rain_width) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect()

evaporation <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Evapotranspiration (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = tabernas),
    width = 0.8
  ) +
  scale_fill_gradientn(
    colors = evap_col,
    na.value = NA,
    guide = guide_colorbar(
      title = "Evapotranspiration\n(mm)",
      direction = "horizontal",
      title.position = "top"
    )
  ) +
  geom_line(aes(y = rain / rain_multiplier),
    color = rain_color,
    linewidth = rain_width
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect() +
  theme(
    axis.title.y.right = element_text(color = rain_color, size = 9),
    axis.text.y.right = element_text(color = rain_color, size = 7),
    axis.ticks.y.right = element_line(color = rain_color)
  )

evaporation_delta <- transect_to_plot |>
  filter(rain_sc == rain_scen) |>
  filter(type == "Evapotranspiration (mm)") |>
  ggplot(aes(x = month_int2)) +
  geom_col(aes(y = colID, group = order, fill = delta),
    width = 0.8
  ) +
  scale_fill_gradient2(guide = guide_colorbar(
    title = expression("\u0394 Evapotranspiration\n(mm)"),
    direction = "horizontal",
    title.position = "top"
  )) +
  geom_line(aes(y = rain / rain_multiplier),
    color = rain_color,
    linewidth = rain_width
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * rain_multiplier,
      name = "Rain (mm)",
      breaks = scales::breaks_pretty()
    )
  ) +
  scale_x_continuous(breaks = c(1, 4, 7, 10), labels = month.abb[c(6:12, 1:5)][c(1, 4, 7, 10)]) +
  coord_cartesian(xlim = c(1, 12), clip = "off") +
  theme_transect() +
  theme(
    axis.title.y.right = element_text(color = rain_color, size = 9),
    axis.text.y.right = element_text(color = rain_color, size = 7),
    axis.ticks.y.right = element_line(color = rain_color)
  )

# Combine the plots -------------------------------------------------------

proc_combined <-
  moisture_upper + moisture_lower + waterGain + deepdrain + evaporation +
  plot_layout(nrow = 1)
# Figure 4
proc_combined

proc_combined_delta <-
  moisture_upper_delta + moisture_lower_delta + deepdrain_delta + runoff_delta +
  evaporation_delta +
  plot_layout(nrow = 1)
# Figure S10
proc_combined_delta
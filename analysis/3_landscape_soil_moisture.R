# Date: March 2023
# Author: Selina Baldauf
# Purpose: Reproduce Figures 4, S9 for hillslope scale soil moisture
# with and without biocrusts

library(data.table)
library(tidyverse)
library(terra)
source("analysis/helper_functions/read_spatial.R")

# Read and format results -------------------------------------------------
path_results <- here::here("data/model_results")
# Mask to cut out the hillslope cells from the square model results
the_mask <- terra::rast(
  as.matrix(
    data.table::fread(here::here("model/Parameters/outputMask_74.txt"))
  )
)
scenarios <- c(
  "tabernas", "tabernas_halfrain", "zeroCrust",
  "zeroCrust_halfrain"
)
pl <- list()
for (i in scenarios) {
  # ---------------------------------------------------------
  # Monthly moisture
  # ---------------------------------------------------------

  soil_m <- read_spatial(
    type = "soil", time = "Months",
    result_path = path_results,
    result_folder = paste0("Results_", i),
    xsize = 74
  )

  # make soil a raster
  soil_m <- lapply(soil_m, function(y) {
    r <- lapply(y, function(z) { # turn result into raster
      to_return <- terra::rast(as.matrix(z))
      # crop raster by the mask
      to_return <- mask(to_return, the_mask, maskvalue = 0)
    })
    names(r) <- month.abb[as.integer(str_replace(
      str_replace(names(r), "soil_m", ""),
      ".txt", ""
    ))]
    return(r)
  })

  soil_m <- lapply(soil_m, function(y) {
    lapply(y, ggspatial::df_spatial)
  })

  soil_m <- lapply(soil_m, function(x) {
    dt <- rbindlist(x, idcol = "month")
  })


  pl[[i]] <- rbindlist(soil_m, idcol = "variable")
}

# Brint all the data of the scenarios into one table
plot_list <- rbindlist(pl, idcol = "scenario")
plot_list <- plot_list[!is.na(band1)]

plot_list <- mutate(plot_list,
  scenario = case_when(
    scenario == "tabernas_halfrain" ~ "baseline drought",
    scenario == "zeroCrust_halfrain" ~ "no-biocrust drought",
    scenario == "zeroCrust" ~ "no-biocrust",
    scenario == "tabernas" ~ "baseline"
  )
) |>
  mutate(scenario = factor(scenario, levels = c(
    "baseline", "baseline drought",
    "no-biocrust", "no-biocrust drought"
  )))


# Summarize the data to plot hillslope moisture
to_plot <- plot_list |>
  mutate(rain = case_when(
    str_detect(scenario, "drought") ~ "dry",
    TRUE ~ "default"
  )) |>
  mutate(scenario = str_replace(scenario, " drought", "")) |>
  filter(variable %in% c("L1", "L2")) |>
  mutate(month = factor(month, levels = month.abb[c(6:12, 1:5)])) |>
  mutate(variable = ifelse(variable == "L1", "upper", "lower")) |>
  mutate(variable = factor(variable, levels = c("upper", "lower")))

to_plot_summary <- to_plot |>
  group_by(scenario, month, rain, variable) |>
  mutate(band1 = band1 * 100) |>
  summarize(
    mean_moisture = mean(band1),
    median_moisture = median(band1),
    min_moisture = min(band1),
    max_moisture = max(band1),
    sd_moisture = sd(band1)
  ) |>
  mutate(scenario = factor(scenario, levels = c("baseline", "no-biocrust")))

# Make the plot -----------------------------------------------------------

to_plot_summary |>
  filter(!(scenario %in% c("IC", "L", "C", "L drought", "C drought"))) |>
  mutate(variable = case_when(
    variable == "upper" ~ "upper layer",
    variable == "lower" ~ "lower layer"
  )) %>%
  mutate(variable = factor(variable, levels = c("upper layer", "lower layer"))) %>%
  mutate(rain = case_when(
    rain == "default" ~ "wet year",
    rain == "dry" ~ "dry year"
  )) %>%
  mutate(rain = factor(rain, levels = c("wet year", "dry year"))) %>%
  mutate(scenario = fct_recode(scenario, "baseline-biocrust" = "baseline")) %>%
  ggplot(aes(color = scenario, x = month)) +
  geom_pointrange(
    aes(
      y = median_moisture,
      ymin = median_moisture - sd_moisture,
      ymax = median_moisture + sd_moisture
    ),
    position = position_dodge(width = .95), size = .3
  ) +
  facet_grid(variable ~ rain, scales = "free") +
  scale_color_manual(values = c("#529985", "#767676")) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
  theme_bw() +
  labs(x = "Month", y = "Soil moisture (%)", color = "Scenario") +
  theme(
    legend.position = "bottom"
  )


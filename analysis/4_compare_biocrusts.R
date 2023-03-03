# Date: March 2023
# Author: Selina Baldauf
# Purpose: Reproduce Figure 5 to compare the effects of different biocrusts

library(data.table)
library(lubridate)
library(tidyverse)
library(patchwork)

path_results <- here::here("data/model_results")
crust_colors <- c("#A489A4", "#6F556F", "#DC8861", "#B3AC3C", "#5A561E")


# Read the results --------------------------------------------------------
scenarios <- c("zeroCrust", "cyano", "physical", "lichen")
dates <- seq(
  from = parse_date_time("2010-1-1 0", orders = "%Y%m%d %H"),
  to = parse_date_time("2010-12-31 23", orders = "%Y%m%d %H"),
  by = "hours"
)

temporal_cs <- lapply(scenarios, function(x) {
  temporal_cs <- fread(
    paste0(path_results, "/Results_", x, "/", "/tabernas_1years_CS_", x, ".txt")
  )[
    order(year, day, hour), datetime := dates
  ]
  temporal_cs <- dplyr::select(temporal_cs, !c(year, hour, day))
  return(temporal_cs)
})

names(temporal_cs) <- scenarios
temporal_cs <- rbindlist(temporal_cs, idcol = "scenario")

to_plot_hour <- temporal_cs %>%
  pivot_longer(!c(datetime, scenario)) %>%
  separate(name,
    into = c("variable", "layer", "cell_type"),
    sep = "_", fill = "right"
  )

mean_cols <- names(temporal_cs)[str_detect(
  names(temporal_cs),
  "temperature|moisture"
)]

sum_cols <- names(temporal_cs)[str_detect(
  names(temporal_cs),
  paste(
    c(
      "rain", "water",
      "infiltration",
      "drainage", "QD", "EP"
    ),
    collapse = "|"
  )
)]

day <- temporal_cs %>%
  group_by(date = date(datetime), scenario) %>%
  summarize(
    across(.cols = all_of(mean_cols), .fns = mean),
    across(.cols = all_of(sum_cols), .fns = sum)
  ) %>%
  pivot_longer(!c(date, scenario)) %>%
  separate(name, into = c("name", "layer", "cell_type"), sep = "_", fill = "right") %>%
  unite("name", name:layer, sep = "_") %>%
  filter(!name %in% c(
    "water_L0", "rain_NA", "temperature_NA",
    "infiltration_L2", "drainage_L1"
  )) %>%
  mutate(name = factor(name, levels = c(
    "moisture_L1", "moisture_L2",
    "infiltration_L1", "EP_L1",
    "drainage_L2", "QD_L0"
  ))) %>%
  filter(
    !((scenario == "cyano" & cell_type %in% c("crust0", "crust1", "crust3")) |
      (scenario == "physical" & cell_type %in% c("crust0", "crust2", "crust3")) |
      (scenario == "lichen" & cell_type %in% c("crust0", "crust1", "crust2")) |
      (scenario == "zeroCrust" & cell_type %in% c("crust1", "crust2", "crust3")))
  )

day <- mutate(day, value = ifelse(
  name %in% c("moisture_L1", "moisture_L2"), value * 100, value
))

crust_labeller <- c(
  crust0 = "no-biocrust",
  crust1 = "IC",
  crust2 = "C",
  crust3 = "L"
)

scenario_labeller <- c(
  cyano = "C",
  lichen = "L",
  physical = "IC",
  zeroCrust = "no-biocrust"
)


# Make the plots of the crusted cells -------------------------------------
to_plot_day <- day %>%
  filter(month(date) %in% c(12, 1:3)) |>
  mutate(date = case_when(
    month(date) == 12 ~ date - years(1),
    TRUE ~ date
  )) |>
  filter(cell_type != "crust0") |>
  mutate(cell_type = factor(cell_type, levels = c("crust0", "crust1", "crust2", "crust3")))

moist_1_d <- to_plot_day %>%
  filter(cell_type != "veg" & name == "moisture_L1") %>%
  ggplot() +
  labs(y = "Moisture\nupper (%)") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

moist_2_d <- to_plot_day %>%
  filter(cell_type != "veg" & name == "moisture_L2") %>%
  ggplot() +
  labs(y = "Moisture\nlower (%)") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

infiltration_d <- to_plot_day %>%
  filter(cell_type != "veg" & name == "infiltration_L1") %>%
  ggplot(base_size = 10) +
  labs(y = "Infiltration upper (mm)") +
  theme(axis.text.x = element_blank())

evap_d <- to_plot_day %>%
  filter(cell_type != "veg" & name == "EP_L1") %>%
  ggplot() +
  labs(y = "Evaporation\nupper (mm)") +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank())

deepdrain_d <- to_plot_day %>%
  filter(cell_type != "veg" & name == "drainage_L2") %>%
  ggplot() +
  labs(y = "Deep\ndrainage (mm)") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

runoff_d <- to_plot_day %>%
  filter(cell_type != "veg" & name == "QD_L0") %>%
  ggplot() +
  labs(y = "Runoff (mm)") +
  theme_bw(base_size = 10) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Combine the plots to Fig. 5
(moist_1_d / moist_2_d / deepdrain_d / runoff_d / evap_d) +
  plot_layout(ncol = 1, guides = "collect", ) &
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  geom_line(aes(x = date, y = value, color = cell_type)) &
  scale_x_date(expand = c(0, 0)) &
  scale_color_manual(values = crust_colors, labels = crust_labeller) &
  labs(color = "Biocrust type") &
  theme(
    plot.tag.position = c(0.01, 0.98),
    plot.tag = element_text(face = "bold")
  )

# Make the plots of the vegetated cells (Fig. S11) ------------------------
to_plot_day <- day %>%
  filter(month(date) %in% c(12, 1:3)) |>
  mutate(date = case_when(
    month(date) == 12 ~ date - years(1),
    TRUE ~ date
  )) |>
  filter(scenario != "zeroCrust") |>
  mutate(scenario = factor(scenario, levels = c("physical", "cyano", "lichen", "zeroCrust")))

moist_1_d_v <- to_plot_day %>%
  filter(cell_type == "veg" & name == "moisture_L1") %>%
  ggplot() +
  labs(y = "Moisture upper (%)") +
  theme(axis.text.x = element_blank())

moist_2_d_v <- to_plot_day %>%
  filter(cell_type == "veg" & name == "moisture_L2") %>%
  ggplot() +
  labs(y = "Moisture lower (%)")

# Combine to get actual plot (S10)
moist_1_d_v + moist_2_d_v +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  geom_line(aes(x = date, y = value, color = scenario)) &
  scale_color_manual(values = crust_colors, labels = scenario_labeller) &
  labs(color = "Biocrust type") &
  theme_bw() &
  theme(
    axis.title.x = element_blank(),
    plot.margin = margin(t = 6, r = 4, b = 0, l = 2),
    plot.tag.position = c(0.05, 1.05)
  )

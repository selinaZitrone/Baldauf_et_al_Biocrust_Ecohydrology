# Date: January 2022
# Author: Selina Baldauf
# Purpose: control running ecohyd scenarios

# get command line argumenst
args <- commandArgs(trailingOnly = TRUE)

elevation <- args[1]
soil <- args[2]
vegetation <- args[3]
parameter <- args[4]
name <- args[5]

# read in parameter file to get xsize
parameter_file <- paste0("modelparameters_", parameter, ".txt")

if (!file.exists(paste0("Parameters/", parameter_file))) {
    stop(paste0("File ", parameter_file, " does not exist in /Parameters_all!"))
}
# read the content
parameter_file_content <- readLines(
    con = paste0("./Parameters/", parameter_file),
    skipNul = TRUE
)
# find the xsize from there
xsize_string <- parameter_file_content[grepl(
    pattern = "xsize",
    x = parameter_file_content
)]
xsize <- sub(".* ([0-9]+).*", "\\1", xsize_string)

# check if the other input files that are needed exist
elevation_file <- paste0("elevation_", xsize, "_", elevation, ".txt")
soil_file <- paste0("soilparameters_", soil, ".txt")
vegetation_file <- paste0("vegetationparameters_", vegetation, ".txt")

for (i in c(elevation_file, soil_file, vegetation_file)) {
    if (!file.exists(paste0("Parameters/", i))) {
        stop(paste0("File ", i, "does not exist in /Parameters_all!"))
    }
}

# write the correct model scenarios file
scenarios <- c(
    "#General scenario parameters",
    "Scenarios: 1",
    "ClimateRepetitions: 1",
    "ModelRepetitions: 1",
    "",
    "###################",
    "#   Scenario 1    #",
    "##################*",
    paste0("OutputID: ", name),
    paste0("ModelParametersID: ", parameter),
    paste0("WeatherFileID: ", "tabernas"),
    paste0("ElevationFileID: ", elevation),
    paste0("SoilParametersID: ", soil),
    paste0("VegetationParametersID: ", vegetation)
)

writeLines(scenarios, con = "./Parameters/modelscenarios.txt")

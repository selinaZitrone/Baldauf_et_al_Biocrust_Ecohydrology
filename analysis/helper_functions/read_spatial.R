# Author: Selina Baldauf
# Date: December 2021
# Purpose: Functions to read spatial output for soil and crust
#    Functions can be included in any file analysing spatial Ecohyd Output


#' @title Read spatial biocrust output
#' @param time Type of output to be read: Can be "day", "month", or "hour"
#' @param type Type of output to be read: Can be "crust" or "soil"
#' @param spatial_folder name of folder of spatial output
#' @param result_folder name of folder of results to read
#' @param xsize size of grid to read
#' @return list with data.tables containing the matrix to be plotted.
#'    In case of soil output, moisture from both layers will be included in the
#'    list
read_spatial <- function(type, time,
                         spatial_folder = "Spatial",
                         result_folder = "Results",
                         result_path = here::here(),
                         xsize = 30) {
  if (!(time %in% c("Days", "Months", "Hours"))) {
    stop(paste0(
      "The time argument should be either
            \"Days\", \"Months\", or \"Hours\" but is ", time, " instead"
    ))
  }
  if (!(type %in% c("soil", "crust"))) {
    stop(paste0("The ype argument should be either
        \"soil\" or \"crust\" but is ", type, " instead"))
  }
  # Filenames without full path
  filenames <- list.files(
    path = paste0(result_path, "/", result_folder, "/", spatial_folder, "/", time),
    pattern = type
  )
  
  # Files to read
  file_list <- list.files(
    path = paste0(result_path, "/", result_folder, "/", spatial_folder, "/", time),
    pattern = type,
    full.names = TRUE
  )
  
  # Read files
  
  files <- lapply(file_list, data.table::fread, skip = 1, nrow = xsize)
  # name the list
  names(files) <- filenames
  if (type == "soil") {
    files_L2 <- lapply(file_list, data.table::fread,
                       skip = xsize + 4,
                       nrows = xsize
    )
    # name the list
    names(files_L2) <- filenames
    files_QD <- lapply(file_list, data.table::fread,
                       skip = 2 * xsize + 7,
                       nrow = xsize
    )
    names(files_QD) <- filenames
    files_infiltration <- lapply(file_list, data.table::fread,
                                 skip = 3 * xsize + 10, nrow = xsize
    )
    names(files_infiltration) <- filenames
    files_runon <- lapply(file_list, data.table::fread,
                          skip = 4 * xsize + 13, nrow = xsize
    )
    names(files_runon) <- filenames
    files_deepdrain <- lapply(file_list, data.table::fread,
                              skip = 5 * xsize + 16, nrow = xsize
    )
    names(files_deepdrain) <- filenames
    
    if (time != "Hours") {
      files_EPtot <- lapply(file_list, data.table::fread,
                            skip = 6 * xsize + 19, nrow = xsize
      )
      names(files_EPtot) <- filenames
    }
    
    # calculate water gain by subtracting runon and runoff
    files_waterGain <- lapply(names(files_runon), function(x) {
      files_runon[[x]] - files_QD[[x]]
    })
    names(files_waterGain) <- filenames
    if (time != "Hours") {
      files <- list(
        L1 = files,
        L2 = files_L2,
        QD = files_QD,
        FL1 = files_infiltration,
        runon = files_runon,
        deepdrain = files_deepdrain,
        waterGain = files_waterGain,
        EPtot = files_EPtot # ,#
        # EPL1 = files_EPL1,
        # EP2 = files_EPL2
      )
    } else {
      files <- list(
        L1 = files,
        L2 = files_L2,
        QD = files_QD,
        FL1 = files_infiltration,
        runon = files_runon,
        deepdrain = files_deepdrain,
        waterGain = files_waterGain
      )
    }
  } else {
    if (time == "Hours") {
      files_surface <- lapply(file_list, data.table::fread,
                              skip = xsize + 2, nrow = xsize
      )
      names(files_surface) <- filenames
      
      files <- list(
        crust_moisture = files,
        waterL0 = files_surface
      )
    } else {
      files <- list(
        crust_moisture = files
      )
    }
  }
  return(files)
}

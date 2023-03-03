#pragma once

#include <map>
#include <string>

#include "BiocrustLandscape.h"
#include "Parameters.h"
#include "WaterLandscape.h"

class modelOutput {
 private:
  // number of days for each month
  map<const int, const int> months{{1, 31}, {2, 28},  {3, 31},  {4, 30},
                                   {5, 31}, {6, 30},  {7, 31},  {8, 31},
                                   {9, 30}, {10, 31}, {11, 30}, {12, 31}};

  unsigned int current_month = 1;
  unsigned int day_in_month = 1;
  bool output_month =
      false;  // If day is the last day of the month, output is written

  // Names of spatial output directories
  const string output_dir_months =
      "Results" + delimiter + "Spatial" + delimiter + "Months" + delimiter;
  const string output_dir_days =
      "Results" + delimiter + "Spatial" + delimiter + "Days" + delimiter;

  const string output_dir_hours =
      "Results" + delimiter + "Spatial" + delimiter + "Hours" + delimiter;

  vector<string> spatial_output_times;  // spatial output times (days and hours)

  /*************************************************************************************
    \author Selina Baldauf
    \date December 2021
    \brief Read in the days and hours for which to output spatial moisture data,
    reads in the file Parameters/spatial_output.yml. Is called in
    modelOutput->spatialOutput
  ****************************************************************************************/
  void read_spatial_output_times();

  /*************************************************************************************
    \author Selina Baldauf
    \date December 2021
    \brief Write mean crust moisture for each month of the year for each
    cell in spatial output file. Is called in modelOutput->spatialOutput
    \param current_month current month that is
    written
    \param bl BiocrustLandscape object to extract monthly moisture sums from
    \param gridSize dimension of the grid
  ****************************************************************************************/
  void spatialOutput_month_crust(unsigned int current_month,
                                 BiocrustLandscape* bl, unsigned int gridSize);

  /*************************************************************************************
    \author Selina Baldauf
    \date December 2021
    \brief Write mean soil moisture for each month of the year for each cell
    in spatial output file. Is called in modelOutput->spatialOutput
    \param current_month current month
    \param wl WaterLandscape object to extract monthly moisture sums from
    \param gridSize dimension of the grid
  ****************************************************************************************/
  void spatialOutput_month_soil(unsigned int current_month, WaterLandscape* wl,
                                unsigned int gridSize);

  /*************************************************************************************
    \author Selina Baldauf
    \date March 2022
    \brief Reset accumulators for mean/sum monthly values per cell
    \param wl WaterLandscape object to extract monthly moisture sums from
    \param gridSize dimension of the grid
  ****************************************************************************************/
  void resetSpatialOutput_month_soil(WaterLandscape* wl, unsigned int gridSize);

  /*************************************************************************************
    \author Selina Baldauf
    \date December 2021
    \brief Write mean crust moisture for each day that is specified by user and
    for each cell in spatial output file. Is called in
    modelOutput->spatialOutput
    \param current_month current month
    \param current_day current day
    \param bl BiocrustLandscape object to extract daily moisture sums from
    \param gridSize dimension of the grid
  ****************************************************************************************/
  void spatialOutput_day_crust(unsigned int current_month,
                               unsigned int current_day, BiocrustLandscape* bl,
                               unsigned int gridSize);

  /*************************************************************************************
    \author Selina Baldauf
    \date December 2021
    \brief Write mean soil moisture for each day that is specified by user and
    for each cell in spatial output file. Is called in
    modelOutput->spatialOutput
    \param current_month current month
    \param current_day current day
    \param wl WaterLandscape object to extract monthly moisture sums from
    \param gridSize dimension of the grid
  ****************************************************************************************/
  void spatialOutput_day_soil(unsigned int month, unsigned int current_day,
                              WaterLandscape* wl, unsigned int gridSize);

  /*************************************************************************************
    \author Selina Baldauf
    \date December 2021
    \brief Write mean crust moisture for each hour that is specified by user and
    for each cell in spatial output file. Is called in
    modelOutput->spatialOutput
    \param current_month current month
    \param current_day current day
    \param hour current hour
    \param bl BiocrustLandscape object to extract daily moisture sums from
    \param gridSize dimension of the grid
  ****************************************************************************************/
  void spatialOutput_hour_crust(unsigned int current_month,
                                unsigned int current_day, unsigned int hour,
                                BiocrustLandscape* bl, unsigned int gridSize);

  /*************************************************************************************
    \author Selina Baldauf
    \date December 2021
    \brief Write mean soil moisture for each hour that is specified by user and
    for each cell in spatial output file. Is called in
    modelOutput->spatialOutput
    \param current_month current month
    \param current_day current day
    \param hour current hour
    \param wl WaterLandscape object to extract monthly moisture sums from
    \param gridSize dimension of the grid
  ****************************************************************************************/
  void spatialOutput_hour_soil(unsigned int current_month,
                               unsigned int current_day, unsigned int hour,
                               WaterLandscape* wl, unsigned int gridSize);

 public:
  /*************************************************************************************
   \author Selina Baldauf
   \date December 2021
   \brief Function to write all spatial output of crust moisture and soil
   moisture (L1 & L2) for: monthly means, daily and hourly means as specified
   by user input. Is called in WaterLandscape->CalculateProcesses
   \param current_month current month
   \param current_day current day
   \param bl BiocrustLandscape object to extract daily moisture sums from
 ****************************************************************************************/
  void spatialOutput(unsigned int day, unsigned int hour, BiocrustLandscape* bl,
                     WaterLandscape* wl);
};
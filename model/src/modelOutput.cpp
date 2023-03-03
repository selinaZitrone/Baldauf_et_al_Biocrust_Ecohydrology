#include "modelOutput.h"

#include <vector>

#include "json.hpp"

// for convenience
using json = nlohmann::json;

/*************************************************************************************
Function to write all spatial output of crust moisture and soil moisture (L1 &
L2) for: monthly means, daily and hourly means as specified by user input. Is
called in WaterLandscape->CalculateProcesses
****************************************************************************************/
void modelOutput::spatialOutput(unsigned int day, unsigned int hour,
                                BiocrustLandscape* bl, WaterLandscape* wl) {
  const int xsize = bl->getGridSize();
  const int ysize = bl->getGridSize();

  // if beginning of year, read spatial output times
  if (day == 0 & hour == 0) {
    read_spatial_output_times();
  }

  string temp_day = to_string(day_in_month) + "-" + to_string(current_month);
  string temp_day_hour = to_string(day_in_month) + "-" +
                         to_string(current_month) + "-" + to_string(hour);

  // write hourly output if its time
  if (std::count(spatial_output_times.begin(), spatial_output_times.end(),
                 temp_day_hour) == 1) {
    spatialOutput_hour_crust(current_month, day_in_month, hour, bl, xsize);
    spatialOutput_hour_soil(current_month, day_in_month, hour, wl, xsize);
  }

  // write daily output and monthly output if it'S time
  if (hour == 23) {
    if (std::count(spatial_output_times.begin(), spatial_output_times.end(),
                   temp_day) == 1) {
      spatialOutput_day_crust(current_month, day_in_month, bl, xsize);
      spatialOutput_day_soil(current_month, day_in_month, wl, xsize);
    }

    // reset the accumulation of daily moisture
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        wl->daily_soil_moistureL1[x][y] = 0.0;
        wl->daily_soil_moistureL2[x][y] = 0.0;
        bl->daily_crust_moisture[x][y] = 0.0;
      }
    }
  }

  // write monthly output if it'S the time
  if (hour == 23 & day_in_month == months[current_month]) {
    spatialOutput_month_crust(current_month, bl, xsize);
    spatialOutput_month_soil(current_month, wl, xsize);

    // reset accumulation of soil variables
    resetSpatialOutput_month_soil(wl, xsize);

    // reset accumulation of monthly moisture
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        bl->monthly_crust_moisture[x][y] = 0;
        wl->monthly_soil_moistureL1[x][y] = 0;
        wl->monthly_soil_moistureL2[x][y] = 0;
      }
    }
  }

  // update day and month for next day
  if (hour == 23) {
    if (day_in_month < months[current_month]) {
      day_in_month++;
    } else {
      day_in_month = 1;
      current_month++;
    }
  }
}

/*******************************************************************************************
Write mean crust moisture for each month of the year for each cell in spatial
output file. Is called in modelOutput->spatialOutput
******************************************************************************************/
void modelOutput::spatialOutput_month_crust(unsigned int current_month,
                                            BiocrustLandscape* bl,
                                            unsigned int gridSize) {
  string outputFileName =
      output_dir_months + "crust_m" + to_string(current_month) + ".txt";
  fstream monthlyOutFile;
  monthlyOutFile.open(outputFileName.c_str(), ios::out | ios::trunc);
  unsigned int hours_per_month = months[current_month] * 24;

  if (!monthlyOutFile) {
    cerr << outputFileName << " could not be opened";
    exit(-1);
  }
  monthlyOutFile << "Crust moisture [-]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << bl->monthly_crust_moisture[x][y] / hours_per_month;
    }
    monthlyOutFile << "\n";
  }
}

/*******************************************************************************************
Write mean soil moisture for each month of the year for each cell in spatial
output file. Is called in modelOutput->spatialOutput
******************************************************************************************/
void modelOutput::spatialOutput_month_soil(unsigned int current_month,
                                           WaterLandscape* wl,
                                           unsigned int gridSize) {
  string outputFileName =
      output_dir_months + "soil_m" + to_string(current_month) + ".txt";
  fstream monthlyOutFile;
  monthlyOutFile.open(outputFileName.c_str(), ios::out | ios::trunc);
  unsigned int hours_per_month = months[current_month] * 24;

  if (!monthlyOutFile) {
    cerr << outputFileName << " could not be opened";
    exit(-1);
  }
  monthlyOutFile << "Soil moisture L1 [-]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_soil_moistureL1[x][y] / hours_per_month;
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "Soil moisture L2 [-]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_soil_moistureL2[x][y] / hours_per_month;
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "Runoff sum [mm]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_QD_cell[x][y];
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "Infiltration sum L1 [mm]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_infiltrationL1_cell[x][y];
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "Runon sum [mm]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_runon_cell[x][y];
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "Drainage deep [mm]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_deepdrain_cell[x][y];
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "EPtot [mm]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_evapotranspiration_tot_cell[x][y];
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "EPL1 [mm]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_evapotranspiration_L1_cell[x][y];
    }
    monthlyOutFile << "\n";
  }
  monthlyOutFile << "\n \n";
  monthlyOutFile << "EPL2 [mm]"
                 << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) monthlyOutFile << "\t";
      monthlyOutFile << wl->monthly_evapotranspiration_L2_cell[x][y];
    }
    monthlyOutFile << "\n";
  }
}

void modelOutput::resetSpatialOutput_month_soil(WaterLandscape* wl,
                                                unsigned int gridSize) {
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      wl->monthly_soil_moistureL1[x][y] = 0;
      wl->monthly_soil_moistureL2[x][y] = 0;
      wl->monthly_deepdrain_cell[x][y] = 0;
      wl->monthly_QD_cell[x][y] = 0;
      wl->monthly_infiltrationL1_cell[x][y] = 0;
      wl->monthly_runon_cell[x][y] = 0;
      wl->monthly_evapotranspiration_L1_cell[x][y] =
          wl->monthly_evapotranspiration_L2_cell[x][y] =
              wl->monthly_evapotranspiration_tot_cell[x][y] = 0;
    }
  }
}

/*******************************************************************************************
Write mean crust moisture for each day that is specified by user and for each
cell in spatial output file. Is called in modelOutput->spatialOutput
******************************************************************************************/
void modelOutput::spatialOutput_day_crust(unsigned int current_month,
                                          unsigned int current_day,
                                          BiocrustLandscape* bl,
                                          unsigned int gridSize) {
  string outputFileName = output_dir_days + "crust_d" + to_string(current_day) +
                          "_m" + to_string(current_month) + ".txt";
  fstream dailyOutFile;
  dailyOutFile.open(outputFileName.c_str(), ios::out | ios::trunc);

  if (!dailyOutFile) {
    cerr << outputFileName << " could not be opened";
    exit(-1);
  }
  dailyOutFile << "Crust moisture [-]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << bl->daily_crust_moisture[x][y] / 24;
    }
    dailyOutFile << "\n";
  }
}

/*******************************************************************************************
Write mean soil moisture for each day that is specified by user and for each
cell in spatial output file. Is called in modelOutput->spatialOutput
******************************************************************************************/
void modelOutput::spatialOutput_day_soil(unsigned int current_month,
                                         unsigned int current_day,
                                         WaterLandscape* wl,
                                         unsigned int gridSize) {
  string outputFileName = output_dir_days + "soil_d" + to_string(current_day) +
                          "_m" + to_string(current_month) + ".txt";
  fstream dailyOutFile;
  dailyOutFile.open(outputFileName.c_str(), ios::out | ios::trunc);

  if (!dailyOutFile) {
    cerr << outputFileName << " could not be opened";
    exit(-1);
  }
  dailyOutFile << "Soil moisture L1 [-]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_soil_moistureL1[x][y] / 24;
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "Soil moisture L2 [-]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_soil_moistureL2[x][y] / 24;
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "Runoff sum [mm]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_QD_cell[x][y];
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "Infiltration L1 sum [mm]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_infiltrationL1_cell[x][y];
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "Runon sum [mm]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_runon_cell[x][y];
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "Drainage deep [mm]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_deepdrain_cell[x][y];
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "EPtot [mm]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_evapotranspiration_tot_cell[x][y];
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "EPL1 [mm]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_evapotranspiration_tot_cell[x][y];
    }
    dailyOutFile << "\n";
  }
  dailyOutFile << "\n \n";
  dailyOutFile << "EPL2 [mm]"
               << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) dailyOutFile << "\t";
      dailyOutFile << wl->daily_evapotranspiration_L2_cell[x][y];
    }
    dailyOutFile << "\n";
  }
}

/*******************************************************************************************
Write mean crust moisture for each hour that is specified by user and for each
cell in spatial output file. Is called in modelOutput->spatialOutput
******************************************************************************************/
void modelOutput::spatialOutput_hour_crust(unsigned int current_month,
                                           unsigned int current_day,
                                           unsigned int hour,
                                           BiocrustLandscape* bl,
                                           unsigned int gridSize) {
  string outputFileName =
      output_dir_hours + "crust_d" + to_string(current_day) + "_m" +
      to_string(current_month) + "_h" + to_string(hour) + ".txt";
  fstream hourlyOutFile;
  hourlyOutFile.open(outputFileName.c_str(), ios::out | ios::trunc);

  if (!hourlyOutFile) {
    cerr << outputFileName << " could not be opened";
    exit(-1);
  }
  hourlyOutFile << "Crust moisture [-]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << bl->getCellMoisture(x, y);
    }
    hourlyOutFile << "\n";
  }
  hourlyOutFile << "Surface water [mm]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << bl->getCellSurfaceWater(x, y);
    }
    hourlyOutFile << "\n";
  }
}

/*******************************************************************************************
Write mean soil moisture and runoff for each hour that is specified by user
and for each cell in spatial output file. Is called in
modelOutput->spatialOutput
******************************************************************************************/
void modelOutput::spatialOutput_hour_soil(unsigned int current_month,
                                          unsigned int current_day,
                                          unsigned int hour, WaterLandscape* wl,
                                          unsigned int gridSize) {
  string outputFileName = output_dir_hours + "soil_d" + to_string(current_day) +
                          "_m" + to_string(current_month) + "_h" +
                          to_string(hour) + ".txt";
  fstream hourlyOutFile;
  hourlyOutFile.open(outputFileName.c_str(), ios::out | ios::trunc);

  if (!hourlyOutFile) {
    cerr << outputFileName << " could not be opened";
    exit(-1);
  }
  hourlyOutFile << "Soil moisture L1 [-]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << wl->soilgrid[x][y]->waterL1_rel;
    }
    hourlyOutFile << "\n";
  }
  hourlyOutFile << "\n \n";
  hourlyOutFile << "Soil moisture L2 [-]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << wl->soilgrid[x][y]->waterL2_rel;
    }
    hourlyOutFile << "\n";
  }
  hourlyOutFile << "\n \n";
  hourlyOutFile << "Runoff [mm]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << wl->soilgrid[x][y]->QD;
    }
    hourlyOutFile << "\n";
  }
  hourlyOutFile << "\n \n";
  hourlyOutFile << "Infiltration L1 sum [mm]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << wl->soilgrid[x][y]->FL1;
    }
    hourlyOutFile << "\n";
  }
  hourlyOutFile << "\n \n";
  hourlyOutFile << "Runon sum [mm]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << wl->runon[x][y];
    }
    hourlyOutFile << "\n";
  }
  hourlyOutFile << "\n \n";
  hourlyOutFile << "Drainage deep [mm]"
                << "\n";
  for (int x = 0; x < gridSize; x++) {
    for (int y = 0; y < gridSize; y++) {
      if (y != 0) hourlyOutFile << "\t";
      hourlyOutFile << wl->soilgrid[x][y]->draindeep;
    }
    hourlyOutFile << "\n";
  }
}

/*******************************************************************************************
Read in the days and hours for which to output spatial moisture data, reads in
the file Parameters/spatial_output.yml. Is called in
modelOutput->spatialOutput
*******************************************************************************************/
void modelOutput::read_spatial_output_times() {
  string filename = "Parameters" + delimiter + "spatial_output.json";
  ifstream times_file(filename);
  json j;
  times_file >> j;

  for (auto& element : j["days"]) {
    spatial_output_times.push_back(element.get<string>());
  }
  for (auto& element : j["hours"]) {
    spatial_output_times.push_back(element.get<string>());
  }
}
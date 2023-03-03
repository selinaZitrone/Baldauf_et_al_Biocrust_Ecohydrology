#pragma once
#include <map>
#include <string>

#include "CrustParameters.h"

using namespace std;

class BiocrustCell {
 public:
  map<string, double> crust_params_cell;  // crust parameters of this cell
  bool with_hydrophobicity =
      false;  // should the cell be run with or without hydrophobicity

  // infiltration variables
  // assuming a vegetation cover of 0% in cells occupied by biocrusts
  const double veg_function = 0.013 / 0.067;
  bool crust_wettingFront = false;  // is there wettingfront in timestep

  unsigned int crust_type = 0;  // type of crust that occupies this cell (0
                                // means no crust) process variables

  // process variables
  double moisture_crust = 0.0;      // relative crust moisture (-)
  double surface_water = 0.0;       // surface water (mm)
  double infiltration_crust = 0.0;  // infiltration into crust (mm)
  double leakage_crust = 0.0;       // leakage from crust into soil (mm)
  double evaporation_crust = 0.0;   // evaporation from crust (mm)
  double runoff_crust = 0.0;        // runoff from the crust (mm)
  double ponding_crust = 0.0;       // ponding on the crust surface (mm)
  double ponded_evaporation = 0.0;  // evaporation of ponded surface water (mm)

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Calculate crust infiltration. Is called in
  BiocrustCell->calculateCellProcessHour
  \param surface_water surface water (mm) available for infiltration
  \param moisture_crust current relative moisture content of the biocrust
  \param crust_wettingFront boolean holding whether there is already a
  wettingfront from the last time step
  \param with_hydrophobicity boolean saying
  whether infiltration should consider hydrophobicity or not
  ****************************************************************************************/
  double crustInfiltration(double surface_water, double moisture_crust,
                           bool crust_wettingFront, bool with_hydrophobicity);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Calculate hydrophobicity factor. This factor is multiplied with
  infiltration to reduce infiltration if crust is dry. Is called in
  BiocrustCell->crustInfiltration
  \param moisture_crust current relative moisture content of the biocrust
  ****************************************************************************************/
  double calculateHydrophobicityFactor(double moisture_crust);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Update crust moisture with current (positive or negative) water input.
  Is called in BiocrustCell->calculateCellProcessHour and
  BiocrustCell->calculateCellProcessesDay
  \param current_moisture current relative moisture content of the biocrust
  \param water_input water input (mm)
  ****************************************************************************************/
  double updateCrustMoisture(double current_moisture, double water_input);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Calculate leakage of water from biocrust to soil. Is called in
  BiocrustCell->calculateCellProcessHour
  \param current_moisture current relative moisture content of the biocrust
  \param surface_water water on the biocrust surface to be leaked
  ****************************************************************************************/
  double crustLeakageIntoSoil(double current_moisture, double surface_water);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Calculate evaporation from the crust layer. Is called in
  BiocrustCell->calculateCellProcessDay
  \param current_moisture current relative moisture content of the biocrust
  \param potEP potential Evaporation (from WaterLandscape.cpp)
  \param openDemand Evaporative demand that is still open after ponded
  evaporation
  ****************************************************************************************/
  double crustEvaporation(double current_moisture, double potEP,
                          double openDemand);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Calculate evaporation from the ponded water. All surface water can
  freely evaporate if it is less than the evaporative demand. Evaporative demand
  is calculated as potential evaporation. Is called in
  BiocrustCell->calculateCellProcessDay
  \param surface_water ponded water [mm]
  \param openDemand Evaporative demand that is open (pot EP)
  ****************************************************************************************/
  double pondedEvaporation(double openDemand, double surface_water);

  /*************************************************************************************
   \author Selina Baldauf
   \date October 2021
   \brief Calculate all biocrust cell processes that run on hourly timescale. Is
   called in BiocrustLandscape->calculateProcessHour
   ****************************************************************************************/
  void calculateCellProcessesHour();

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Calculate all biocrust cell processes that run on daily timescale. Is
   called in BiocrustLandscape->calculateProcessDay
  \param potEP potential evaporation
  ****************************************************************************************/
  void calculateCellProcessesDay(double potEP);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Biocrust Cell constructor
  \param parameters a map of the parameters for the crust that is located in the
  current cell
  \param moisture_init initial crust moisture in the cell (defaults
  to 0.08)
  \param crust_type type of biocrust present in the cell (default is 0 (no
  crust))
  ****************************************************************************************/
  BiocrustCell(vector<map<string, double>> crust_parameters,
               unsigned int crust_type, double moisture_init = 0.08);
};
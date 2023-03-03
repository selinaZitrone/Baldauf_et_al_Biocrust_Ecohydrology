#pragma once

#include <fstream>
#include <vector>

#include "BiocrustCell.h"
#include "Parameters.h"
#include "Weather.h"

class BiocrustLandscape {
 private:
  unsigned int xsize, ysize;
  std::ofstream resultFileCrust;            // daily output file
  std::ofstream resultFileCrustHour;        // hourly output file
  std::ofstream resultFileCrustHour_CS;     // hourly output file
  vector<vector<BiocrustCell*>> crustgrid;  // Biocrust grid

  /*************************************************************************************
    \author Selina Baldauf
    \date November 2021
    \brief Read crust cover file from /Parameters/crustcover_xsize.txt and
       initialize crustgrid of Biocrust cells with the respective crust
       parameters
    \param c Pointer to CrustParameters object from which to take the
       crust parameters for the respective crust types
    ****************************************************************************************/
  void initaliseCrustGridFromFile(CrustParameters* c);

  /*************************************************************************************
  \author Selina Baldauf
  \date November 2021
  \brief Write daily crust output file. Is called in
  BiocrustLandscape->writeBiocrustOutput
  \param year current year
  \param day current day
  \param hour current hour
  ****************************************************************************************/
  void writeCrustOutputFile(unsigned int year, unsigned int day,
                            unsigned int hour);

  /*************************************************************************************
  \author Selina Baldauf
  \date November 2021
  \brief Initialise daily crust output file. Is called in BiocrustLandscape
  constructor.
  ****************************************************************************************/
  void initialiseCrustResultFile();
  void initialiseCrustResultFile_CS();

  /*************************************************************************************
  \author Selina Baldauf
  \date November 2021
  \brief Write hourly crust output file. Is called in
  BiocrustLandscape->writeBiocrustOutput
  \param year current year
  \param day current day
  \param hour current hour
  ****************************************************************************************/
  void writeHourlyCrustOutputFile(unsigned int year, unsigned int day,
                                  unsigned int hour);

  void writeHourlyCrustOutputFile_CS(unsigned int year, unsigned int day,
                                     unsigned int hour);

  /*************************************************************************************
  \author Selina Baldauf
  \date November 2021
  \brief Initialise hourly crust output file. Is called in BiocrustLandscape
  constructor.
  ****************************************************************************************/
  void initialiseHourlyCrustResultFile();

  /*************************************************************************************
  \author Selina Baldauf
  \date November 2021
  \brief Calculate grid means and standard deviation of all hydrological
  processes for every hour, is called in BiocrustLandscape->writeBiocrustOutput
  \param outputCells output mask for the cells to include in the output
  ****************************************************************************************/
  void calculateGridMeans(Int2D outputCells);
  void calculateGridMeans_CS(Int2D outputCells);

  /*************************************************************************************
   \author Selina Baldauf
   \date December 2021
   \brief Calculate cell mean moisture for each day and each month. Is needed
   for spatial output. Is called in BiocrustLandscape->writeBiocrustOutput
   ****************************************************************************************/
  void calculateCellMeanMoisture();

  // variables for accumulating daily values
  double crustMoisture_dayMean = 0.0;
  double crustInfiltration_daySum = 0.0;
  double crustLeakage_daySum = 0.0;
  double crustEvaporation_daySum = 0.0;
  double crustRunoff_daySum = 0.0;
  double crustPonding_daySum = 0.0;
  double crustWaterL0_daySum = 0.0;

  // variables for daily standard deviation over grid
  double crustMoisture_daySD = 0.0;
  double crustInfiltration_daySD = 0.0;
  double crustLeakage_daySD = 0.0;
  double crustEvaporation_daySD = 0.0;
  double crustRunoff_daySD = 0.0;
  double crustPonding_daySD = 0.0;
  double crustWaterL0_daySD = 0.0;

  // variables for hourly grid mean values
  double crustMoisture_gridMean = 0.0;
  double crustInfiltration_gridMean = 0.0;
  double crustLeakage_gridMean = 0.0;
  double crustEvaporation_gridMean = 0.0;
  double crustRunoff_gridMean = 0.0;
  double crustPonding_gridMean = 0.0;
  double crustWaterL0_gridMean = 0.0;
  double pondedEvaporation_gridMean = 0.0;

  double crustWaterL0_gridMean_veg = 0.0;
  double crustWaterL0_gridMean_crust = 0.0;
  double crustLeakage_gridMean_veg = 0.0;
  double crustLeakage_gridMean_crust = 0.0;

  // variables for hourly grid standard deviation
  double crustMoisture_gridSD = 0.0;
  double crustInfiltration_gridSD = 0.0;
  double crustLeakage_gridSD = 0.0;
  double crustEvaporation_gridSD = 0.0;
  double crustRunoff_gridSD = 0.0;
  double crustPonding_gridSD = 0.0;
  double crustWaterL0_gridSD = 0.0;
  double pondedEvaporation_gridSD = 0.0;

 public:
  // vectors for monthly and daily output per cell
  // is used to produce monthly and daily spatial output
  vector<vector<double>> monthly_crust_moisture;
  vector<vector<double>> daily_crust_moisture;

  /*************************************************************************************
   \author Selina Baldauf
   \date October 2021
   \brief Helper function to print the biocrust map that is read in.
   Just to check if map is read correctly. Is called in the
   BiocrustLandscape constructor
  ****************************************************************************************/
  void printCrustMap();

  /*************************************************************************************
   \author Selina Baldauf
   \date October 2021
   \brief Loop through biocrust grid and calculate biocrust cell processes on
   hourly timestep (infiltration, leakage). Is called in
   WaterLandscape->calculateProcesses
   \param prec precipitation (mm)
  ****************************************************************************************/
  void calculateProcessesHour(double prec);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Loop through biocrust grid and calculate biocrust cell processes on
  daily timestep (ponded and biocrust evaporation. Is called in
  WaterLandscape->calculateProcesses
  \param pot_ep potential evaporation (mm) on that day
  ****************************************************************************************/
  void calculateProcessesDay(double pot_ep);

  /*************************************************************************************
   \author Selina Baldauf
   \date October 2021
   \brief Coordinate biocrust output: calculate daily grid means,
   mean moisture for each cell, call functions to write hourly and daily output.
   Is called in WaterLandscape->calculateProcesses
   \param year current year
   \param day current day
   \param hour current hour
   \param outputCells output mask for the hillslope
  ****************************************************************************************/
  void writeBiocrustOutput(unsigned int year, unsigned int day,
                           unsigned int hour, Int2D outputCells);

  /*************************************************************************************
   \author Selina Baldauf
   \date October 2021
   \brief Close crust result files. Is called in controller->runSimulation
  ****************************************************************************************/
  void finishSimulation();

  /*************************************************************************************
  \author Selina Baldauf
  \date December 2021
  \brief Public getter for crust moisture in cell. Is called in
  modelOutput->spatialOutput_hour_crust
  \param x x-coordinatee of cell
  \param y y-coordinate of cell
  ****************************************************************************************/
  double getCellMoisture(int x, int y);

  /*************************************************************************************
  \author Selina Baldauf
  \date January 2022
  \brief Public getter for crust cell surface water. Is called in
  WaterLandscape->calculateProcesses
  \param x x-coordinatee of cell
  \param y y-coordinate of cell
  ****************************************************************************************/
  double getCellSurfaceWater(int x, int y);

  /*************************************************************************************
  \author Selina Baldauf
  \date January 2022
  \brief Public getter for crust leakage. Is called in
  WaterLandscape->calculateProcesses
  \param x x-coordinatee of cell
  \param y y-coordinate of cell
  ****************************************************************************************/
  double getCellLeakage(int x, int y);

  /*************************************************************************************
  \author Selina Baldauf
  \date January 2022
  \brief Public getter for xsize to be used in model output. Is called in
  modelOutput->
  \todo this should be solved more elegantly e.g. using a globals.h with
  global const variables
  ****************************************************************************************/
  int getGridSize();

  /*************************************************************************************
 \author Selina Baldauf
 \date January 2022
 \brief Public getter for crust type Is called in
 WaterLandscape->CalculateProcesses
 \param x x-coordinatee of cell
  \param y y-coordinate of cell
 ****************************************************************************************/
  int getCrusttype(int x, int y);

  /*************************************************************************************
   \author Selina Baldauf
   \date January 2022
   \brief Public getter for EP factor of crust Is called in
   WaterLandscape->CalculateProcesses
   \param x x-coordinatee of cell
    \param y y-coordinate of cell
   ****************************************************************************************/
  double getEPfactor(int x, int y);

  /*************************************************************************************
   \author Selina Baldauf
   \date January 2022
   \brief Public getter for mannings n of crust Is called in
   WaterLandscape->CalculateRunoff
   \param x x-coordinatee of cell
    \param y y-coordinate of cell
   ****************************************************************************************/
  double getManningn(int x, int y);

  /*************************************************************************************
  \author Selina Baldauf
  \date January 2022
  \brief Public setter to update crust cell surface water from Waterlandscape.
  Is called in WaterLandscape->CalculateProcesses
  \param x x-coordinatee of cell
  \param y y-coordinate of cell
  \param amount amount of water to add from soil ponding
  ****************************************************************************************/
  void addSurfaceWater(int x, int y, double amount);

  /*************************************************************************************
  \author Selina Baldauf
  \date October 2021
  \brief Constructor
  \param p Pointer to Parameter object for general parameters
  \param c Pointer to crust parameters
  ****************************************************************************************/
  BiocrustLandscape(Parameter* p, CrustParameters* c);
};

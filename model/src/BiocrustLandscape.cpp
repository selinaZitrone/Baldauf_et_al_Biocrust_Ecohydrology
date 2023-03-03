#include "BiocrustLandscape.h"

#include <math.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

/*************************************************************************************
Loop through biocrust grid and calculate biocrust cell processes on hourly
timestep (infiltration, leakage). Is called in
WaterLandscape->calculateProcesses
****************************************************************************************/
void BiocrustLandscape::calculateProcessesHour(double prec) {
  for (size_t i = 0; i < xsize; i++) {
    for (size_t j = 0; j < ysize; j++) {
      // add rain to surface water
      crustgrid[i][j]->surface_water += prec;
      if (crustgrid[i][j]->crust_type == 0) {
        crustgrid[i][j]->leakage_crust = crustgrid[i][j]->surface_water;
        crustgrid[i][j]->moisture_crust = 0;
        crustgrid[i][j]->infiltration_crust = 0;
        crustgrid[i][j]->surface_water = 0.0;
        continue;
      }
      crustgrid[i][j]->calculateCellProcessesHour();
    }
  }
}

/*************************************************************************************
Loop through biocrust grid and calculate biocrust cell processes on
daily timestep (ponded and biocrust evaporation.
Is called in WaterLandscape->calculateProcesses
****************************************************************************************/
void BiocrustLandscape::calculateProcessesDay(double pot_ep) {
  // return the water that is leaked throught each biocrust cell to the soil
  for (size_t i = 0; i < xsize; i++) {
    for (size_t j = 0; j < ysize; j++) {
      if (crustgrid[i][j]->crust_type == 0) {
        continue;
      }
      crustgrid[i][j]->calculateCellProcessesDay(pot_ep);
    }
  }
}

/*************************************************************************************
Initialise daily crust output file. Is called in BiocrustLandscape constructor.
****************************************************************************************/
void BiocrustLandscape::initaliseCrustGridFromFile(CrustParameters* c) {
  ifstream CrustCoverFile;
  const std::string CrustcoverFileName =
      "Parameters" + delimiter + "crustcover_" + std::to_string(xsize) + ".txt";
  cout << "Read crust cover from: " << CrustcoverFileName << endl;
  CrustCoverFile.open(CrustcoverFileName.c_str(),
                      std::ios::binary | std::ios::in);

  if (!CrustCoverFile) {
    std::cerr << CrustcoverFileName << " could not be opened" << std::endl;
    exit(-1);
  }
  int testsum = 0;
  unsigned int temp_type = 0;
  // int test[xsize][ysize];
  for (size_t i = 0; i < xsize; i++) {
    for (size_t j = 0; j < ysize; j++) {
      CrustCoverFile >> temp_type;
      crustgrid[i][j] = new BiocrustCell(c->crust_params, temp_type);
      testsum += this->getCrusttype(i, j);
    }
  }
  cout << "crust types test sum at reading: " << testsum << endl;
}

/*************************************************************************************
Helper function to print the biocrust map that is read in. Just to check if map
is read correctly. Is called in the BiocrustLandscape constructor
****************************************************************************************/
void BiocrustLandscape::printCrustMap() {
  for (size_t i = 0; i < xsize; i++) {
    for (size_t j = 0; j < ysize; j++) {
      cout << this->crustgrid[i][j]->crust_type << " ";
    }
    cout << endl;
  }
}

/*************************************************************************************
Initialise daily crust output file. Is called in BiocrustLandscape constructor.
****************************************************************************************/
void BiocrustLandscape::initialiseCrustResultFile() {
  // create and open result file
  string resultFileNameCrust = "Results" + delimiter + "Result.txt";

  resultFileCrust.open(resultFileNameCrust.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!resultFileCrust) {
    cerr << resultFileNameCrust << " could not be opened in main!\n";
    exit(-1);
  }

  // initialise headers
  resultFileCrust << "year"
                  << "\t"
                  << "day"
                  << "\t"
                  << "crustMoisture"
                  << "\t"
                  << "SD_crustMoisture"
                  << "\t"
                  << "infiltration"
                  << "\t"
                  << "SD_infiltration"
                  << "\t"
                  << "leakage"
                  << "\t"
                  << "SD_leakage"
                  << "\t"
                  << "crust_evaporation"
                  << "\t"
                  << "SD_evaporation"
                  << "\t"
                  << "runoff"
                  << "\t"
                  << "SD_runoff"
                  << "\t"
                  << "ponding"
                  << "\t"
                  << "SD_ponding"
                  << "\t"
                  << "waterL0" << endl;
}

void BiocrustLandscape::initialiseCrustResultFile_CS() {
  // create and open result file
  string resultFileNameCrust = "Results" + delimiter + "Result_hour_CS.txt";

  resultFileCrustHour_CS.open(resultFileNameCrust.c_str(),
                              ios::out | ios::trunc);

  // if something goes wrong...
  if (!resultFileCrustHour_CS) {
    cerr << resultFileNameCrust << " could not be opened in main!\n";
    exit(-1);
  }

  // initialise headers
  resultFileCrustHour_CS << "year"
                         << "\t"
                         << "day"
                         << "\t"
                         << "hour"
                         << "\t"
                         << "waterL0_crust"
                         << "\t"
                         << "waterL0_veg"
                         << "\t"
                         << "leakage_crust"
                         << "\t"
                         << "leakage_veg" << endl;
}

/*************************************************************************************
Initialise hourly crust output file. Is called in BiocrustLandscape constructor.
****************************************************************************************/
void BiocrustLandscape::initialiseHourlyCrustResultFile() {
  // create and open result file
  string resultFileNameCrust = "Results" + delimiter + "Result_hour.txt";

  resultFileCrustHour.open(resultFileNameCrust.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!resultFileCrustHour) {
    cerr << resultFileNameCrust << " could not be opened in main!\n";
    exit(-1);
  }

  // initialise headers
  resultFileCrustHour << "year"
                      << "\t"
                      << "day"
                      << "\t"
                      << "hour"
                      << "\t"
                      << "crustMoisture"
                      << "\t"
                      << "crustMoisture_SD"
                      << "\t"
                      << "infiltration"
                      << "\t"
                      << "infiltration_SD"
                      << "\t"
                      << "leakage"
                      << "\t"
                      << "leakage_SD"
                      << "\t"
                      << "ponded_evaporation"
                      << "\t"
                      << "ponded_evaporation_SD"
                      << "\t"
                      << "crust_evaporation"
                      << "\t"
                      << "crust_evaporation_SD"
                      << "\t"
                      << "runoff"
                      << "\t"
                      << "runoff_SD"
                      << "\t"
                      << "ponding"
                      << "\t"
                      << "ponding_SD"
                      << "\t"
                      << "waterL0" << endl;
}

/*************************************************************************************
Calculate grid means and standard deviation of all hydrological processes for
every hour, is called in BiocrustLandscape->writeBiocrustOutput
*************************************************************************************/
void BiocrustLandscape::calculateGridMeans_CS(Int2D outputCells) {
  // calculate the sums (removing the border cells)
  int count_vegCells = 0;
  int count_crustCells = 0;
  crustWaterL0_gridMean_veg = 0.0;
  crustWaterL0_gridMean_crust = 0.0;
  crustLeakage_gridMean_crust = 0.0;
  crustLeakage_gridMean_veg = 0.0;

  // grid means without border cells
  for (size_t i = 0 + 1; i < (xsize - 1); i++) {
    for (size_t j = 0 + 1; j < (ysize - 1); j++) {
      if (outputCells[i][j]) {
        if (crustgrid[i][j]->crust_type == 0) {
          count_vegCells++;
          crustWaterL0_gridMean_veg += crustgrid[i][j]->surface_water;
          crustLeakage_gridMean_veg += crustgrid[i][j]->leakage_crust;
        } else {
          count_crustCells++;
          crustWaterL0_gridMean_crust += crustgrid[i][j]->surface_water;
          crustLeakage_gridMean_crust += crustgrid[i][j]->leakage_crust;
        }
      }
    }
  }
  // number of cells without the border cells
  // unsigned int numcells = (xsize - 2) * (ysize - 2);

  crustWaterL0_gridMean_veg /= count_vegCells;
  crustWaterL0_gridMean_crust /= count_crustCells;
  crustLeakage_gridMean_crust /= count_crustCells;
  crustLeakage_gridMean_veg /= count_vegCells;
}

void BiocrustLandscape::calculateGridMeans(Int2D outputCells) {
  // reset mean variables

  crustMoisture_gridMean = 0.0;
  crustInfiltration_gridMean = 0.0;
  crustLeakage_gridMean = 0.0;
  crustEvaporation_gridMean = 0.0;
  crustRunoff_gridMean = 0.0;
  crustPonding_gridMean = 0.0;
  crustWaterL0_gridMean = 0.0;
  pondedEvaporation_gridMean = 0.0;

  // reset standard deviations
  crustMoisture_gridSD = 0.0;
  crustInfiltration_gridSD = 0.0;
  crustLeakage_gridSD = 0.0;
  crustEvaporation_gridSD = 0.0;
  crustRunoff_gridSD = 0.0;
  crustPonding_gridSD = 0.0;
  crustWaterL0_gridSD = 0.0;
  pondedEvaporation_gridSD = 0.0;

  // grid means without border cells
  for (size_t i = 0 + 1; i < (xsize - 1); i++) {
    for (size_t j = 0 + 1; j < (ysize - 1); j++) {
      if (outputCells[i][j]) {
        crustMoisture_gridMean += this->crustgrid[i][j]->moisture_crust;
        crustInfiltration_gridMean += this->crustgrid[i][j]->infiltration_crust;
        crustLeakage_gridMean += this->crustgrid[i][j]->leakage_crust;
        crustEvaporation_gridMean += this->crustgrid[i][j]->evaporation_crust;
        crustRunoff_gridMean += this->crustgrid[i][j]->runoff_crust;
        crustPonding_gridMean += this->crustgrid[i][j]->ponding_crust;
        crustWaterL0_gridMean += this->crustgrid[i][j]->surface_water;
        pondedEvaporation_gridMean += this->crustgrid[i][j]->ponded_evaporation;
      }
    }
  }
  // number of cells without the border cells
  // unsigned int numcells = (xsize - 2) * (ysize - 2);
  unsigned int numcells = 0;
  for (auto& n : outputCells) {
    for (auto& n2 : n) {
      numcells += n2;
    }
  }

  crustMoisture_gridMean /= numcells;
  crustInfiltration_gridMean /= numcells;
  crustLeakage_gridMean /= numcells;
  crustEvaporation_gridMean /= numcells;
  crustRunoff_gridMean /= numcells;
  crustPonding_gridMean /= numcells;
  crustWaterL0_gridMean /= numcells;
  pondedEvaporation_gridMean /= numcells;

  // calculate grid standard deviations without boarder cells
  for (size_t i = 0 + 1; i < (xsize - 1); i++) {
    for (size_t j = 0 + 1; j < (ysize - 1); j++) {
      if (outputCells[i][j]) {
        crustMoisture_gridSD += pow(
            this->crustgrid[i][j]->moisture_crust - crustMoisture_gridMean, 2);
        crustInfiltration_gridSD +=
            pow(this->crustgrid[i][j]->infiltration_crust -
                    crustInfiltration_gridMean,
                2);
        crustLeakage_gridSD += pow(
            this->crustgrid[i][j]->leakage_crust - crustLeakage_gridMean, 2);
        crustEvaporation_gridSD +=
            pow(this->crustgrid[i][j]->evaporation_crust -
                    crustEvaporation_gridMean,
                2);
        crustRunoff_gridSD +=
            pow(this->crustgrid[i][j]->runoff_crust - crustRunoff_gridMean, 2);
        crustPonding_gridSD += pow(
            this->crustgrid[i][j]->ponding_crust - crustPonding_gridMean, 2);
        crustWaterL0_gridSD += pow(
            this->crustgrid[i][j]->surface_water - crustWaterL0_gridMean, 2);
        pondedEvaporation_gridSD +=
            pow(this->crustgrid[i][j]->ponded_evaporation -
                    pondedEvaporation_gridMean,
                2);
      }
    }
  }

  crustMoisture_gridSD = sqrt(crustMoisture_gridSD / numcells);
  crustInfiltration_gridSD = sqrt(crustInfiltration_gridSD / numcells);
  crustLeakage_gridSD = sqrt(crustLeakage_gridSD / numcells);
  crustEvaporation_gridSD = sqrt(crustEvaporation_gridSD / numcells);
  crustRunoff_gridSD = sqrt(crustRunoff_gridSD / numcells);
  crustPonding_gridSD = sqrt(crustPonding_gridSD / numcells);
  crustWaterL0_gridSD = sqrt(crustWaterL0_gridSD / numcells);
  pondedEvaporation_gridSD = sqrt(pondedEvaporation_gridSD / numcells);
}

/*************************************************************************************
Write daily crust output file. Is called in
BiocrustLandscape->writeBiocrustOutput
*************************************************************************************/
void BiocrustLandscape::writeCrustOutputFile(unsigned int year,
                                             unsigned int day,
                                             unsigned int hour) {
  unsigned int daynumber = year * daysPerYear + day;

  if (hour == 0) {
    // set all counters to 0
    crustMoisture_dayMean = 0.0;
    crustInfiltration_daySum = 0.0;
    crustLeakage_daySum = 0.0;
    crustEvaporation_daySum = 0.0;
    crustRunoff_daySum = 0.0;
    crustPonding_daySum = 0.0;
    crustWaterL0_daySum = 0.0;

    // variables for daily standard deviation over grid
    crustMoisture_daySD = 0.0;
    crustInfiltration_daySD = 0.0;
    crustLeakage_daySD = 0.0;
    crustEvaporation_daySD = 0.0;
    crustRunoff_daySD = 0.0;
    crustPonding_daySD = 0.0;
    crustWaterL0_daySD = 0.0;
  }

  crustMoisture_dayMean += crustMoisture_gridMean;
  crustInfiltration_daySum += crustInfiltration_gridMean;
  crustLeakage_daySum += crustLeakage_gridMean;
  crustEvaporation_daySum += crustEvaporation_gridMean;
  crustRunoff_daySum += crustRunoff_gridMean;
  crustPonding_daySum += crustPonding_gridMean;
  crustWaterL0_daySum += crustWaterL0_gridMean;

  // variables for daily standard deviation over grid
  crustMoisture_daySD += crustMoisture_gridSD;
  crustInfiltration_daySD += crustInfiltration_gridSD;
  crustLeakage_daySD += crustLeakage_gridSD;
  crustEvaporation_daySD += crustEvaporation_gridSD;
  crustRunoff_daySD += crustRunoff_gridSD;
  crustPonding_daySD += crustPonding_gridSD;
  crustWaterL0_daySD += crustWaterL0_gridSD;

  if (hour == 23) {
    // write results
    crustMoisture_dayMean /= 24;
    crustMoisture_daySD /= 24;
    resultFileCrust << year << "\t" << daynumber << "\t"
                    << crustMoisture_dayMean << "\t" << crustMoisture_daySD
                    << "\t" << crustInfiltration_daySum << "\t"
                    << crustInfiltration_daySD << "\t" << crustLeakage_daySum
                    << "\t" << crustLeakage_daySD << "\t"
                    << crustEvaporation_daySum << "\t" << crustEvaporation_daySD
                    << "\t" << crustRunoff_daySum << "\t" << crustRunoff_daySD
                    << "\t" << crustPonding_daySum << "\t" << crustPonding_daySD
                    << "\t" << crustWaterL0_daySum << endl;
  }
}

/*************************************************************************************
Write hourly crust output file. Is called in
BiocrustLandscape->writeBiocrustOutput
*************************************************************************************/
void BiocrustLandscape::writeHourlyCrustOutputFile(unsigned int year,
                                                   unsigned int day,
                                                   unsigned int hour) {
  unsigned int daynumber = year * daysPerYear + day;

  // write results
  resultFileCrustHour << year << "\t" << daynumber << "\t" << hour << "\t"
                      << crustMoisture_gridMean << "\t" << crustMoisture_gridSD
                      << "\t" << crustInfiltration_gridMean << "\t"
                      << crustInfiltration_gridSD << "\t"
                      << crustLeakage_gridMean << "\t" << crustLeakage_gridSD
                      << "\t" << pondedEvaporation_gridMean << "\t"
                      << pondedEvaporation_gridSD << "\t"
                      << crustEvaporation_gridMean << "\t"
                      << crustEvaporation_gridSD << "\t" << crustRunoff_gridMean
                      << "\t" << crustRunoff_gridSD << "\t"
                      << crustPonding_gridMean << "\t" << crustWaterL0_gridSD
                      << "\t" << crustWaterL0_gridMean << endl;
}

void BiocrustLandscape::calculateCellMeanMoisture() {
  for (size_t i = 0; i < xsize; i++) {
    for (size_t j = 0; j < ysize; j++) {
      monthly_crust_moisture[i][j] += crustgrid[i][j]->moisture_crust;
      daily_crust_moisture[i][j] += crustgrid[i][j]->moisture_crust;
    }
  }
}

void BiocrustLandscape::writeHourlyCrustOutputFile_CS(unsigned int year,
                                                      unsigned int day,
                                                      unsigned int hour) {
  resultFileCrustHour_CS << year << "\t" << day << "\t" << hour << "\t"
                         << crustWaterL0_gridMean_crust << "\t"
                         << crustWaterL0_gridMean_veg << "\t"
                         << crustLeakage_gridMean_crust << "\t"
                         << crustLeakage_gridMean_veg << endl;
}

/*************************************************************************************
Coordinate biocrust output: calculate daily grid means,
mean moisture for each cell, call functions to write hourly and daily output.
Is called in WaterLandscape->calculateProcesses
****************************************************************************************/
void BiocrustLandscape::writeBiocrustOutput(unsigned int year, unsigned int day,
                                            unsigned int hour,
                                            Int2D outputCells) {
  calculateGridMeans(outputCells);
  calculateGridMeans_CS(outputCells);
  calculateCellMeanMoisture();
  writeHourlyCrustOutputFile(year, day, hour);
  writeHourlyCrustOutputFile_CS(year, day, hour);
  writeCrustOutputFile(year, day, hour);
}

/*************************************************************************************
Close crust result files. Is called in controller->runSimulation
****************************************************************************************/
void BiocrustLandscape::finishSimulation() {
  resultFileCrust.close();
  resultFileCrustHour.close();
}

/*************************************************************************************
Public getter for crust moisture in cell. Is called in
modelOutput->spatialOutput_hour_crust
****************************************************************************************/
double BiocrustLandscape::getCellMoisture(int x, int y) {
  return crustgrid[x][y]->moisture_crust;
}

/*************************************************************************************
Public getter for crust cell surface water. Is called in
waterLandscape->CalculateProcesses
****************************************************************************************/
double BiocrustLandscape::getCellSurfaceWater(int x, int y) {
  return crustgrid[x][y]->surface_water;
}

/*************************************************************************************
Public setter to update crust cell surface water from Waterlandscape.
Is called in WaterLandscape->CalculateProcesses
****************************************************************************************/
void BiocrustLandscape::addSurfaceWater(int x, int y, double amount) {
  crustgrid[x][y]->surface_water += amount;
}

/*************************************************************************************
Public getter for crust leakage. Is called in
waterLandscape->CalculateProcesses
****************************************************************************************/
double BiocrustLandscape::getCellLeakage(int x, int y) {
  return crustgrid[x][y]->leakage_crust;
}

/*************************************************************************************
 Public getter for crust type Is called in
 WaterLandscape->CalculateProcesses
****************************************************************************************/
int BiocrustLandscape::getCrusttype(int x, int y) {
  return crustgrid[x][y]->crust_type;
}

/*************************************************************************************
 Public getter for EP factor Is called in
 WaterLandscape->CalculateProcesses
****************************************************************************************/
double BiocrustLandscape::getEPfactor(int x, int y) {
  return crustgrid[x][y]->crust_params_cell.find("EP_factor")->second;
}

/*************************************************************************************
 Public getter for Manning's n Is called in
 WaterLandscape->CalculateProcesses
****************************************************************************************/
double BiocrustLandscape::getManningn(int x, int y) {
  return crustgrid[x][y]->crust_params_cell.find("manning_n")->second;
}

/*************************************************************************************
Public getter for xsize to be used in model output. Is called in
  modelOutput->
****************************************************************************************/
int BiocrustLandscape::getGridSize() { return xsize; }

/*************************************************************************************
Constructor
****************************************************************************************/
BiocrustLandscape::BiocrustLandscape(Parameter* p, CrustParameters* c) {
  xsize = p->xsize;
  ysize = p->ysize;
  crustgrid.resize(xsize, vector<BiocrustCell*>(ysize));

  monthly_crust_moisture.resize(xsize, vector<double>(ysize, 0.0));
  daily_crust_moisture.resize(xsize, vector<double>(ysize, 0.0));

  initaliseCrustGridFromFile(c);
  // printCrustMap();
  // cout << "-----------------------------------------" << endl;

  int test_sum = 0;
  for (int i = 0; i < 74; i++) {
    for (int j = 0; j < 74; j++) {
      test_sum += this->getCrusttype(i, j);
    }
  }

  cout << "crust sum initialization crust landscape: " << test_sum << endl;

  // printCrustMap();
  initialiseCrustResultFile();
  initialiseCrustResultFile_CS();
  initialiseHourlyCrustResultFile();
}
/********************************************************************************************
 * WaterLandscape.cpp
 * This file describes the water dynamics in the landscape.
 * It contains a water landscape which has lots of grid cells (soilCell).
 ********************************************************************************************/

#include "WaterLandscape.h"

#include <math.h>  //to use mathematical functions

#include <algorithm>   //to use min/max
#include <cstdlib>     //to use exit()
#include <filesystem>  // to check if file exists
#include <fstream>     //to read and write files
#include <iostream>    //to use cin / cout
#include <string>      //to use strings

#include "SortedList.h"
#include "VegetationLandscape.h"
#include "modelOutput.h"

modelOutput output;

using namespace std;

/*******************************************************************************************
 * calculate soil moisture processes
 * hourly: infiltration, runoff
 * daily: evaporation
 * is called in --> controller::runSimulation
 *******************************************************************************************/
void WaterLandscape::calculateProcesses(int year, int day, Weather* w,
                                        BiocrustLandscape* bl) {
  double rain;
  bool surfaceWater;
  bool leak = false;
  int daynumber = year * daysPerYear + day;

  for (int hour = 0; hour < hoursPerDay; hour++) {
    rain = w->prec[daynumber][hour];
    //---------------------------------------------
    /// PRECIPITATION.
    /// Run precipitation routine only if there is rain.
    /// In this case the surface water flag is set to true.
    //---------------------------------------------
    // if (rain > 0) {
    //   /// \implements \fn WaterLandscape::precipitation(double actualPrec)
    //   precipitation(rain);
    // }

    //--------------------------------------------
    // CRUST Processes
    // Before soil infiltration water has to pass through the biocrust layer if
    // present
    bl->calculateProcessesHour(rain);

    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        soilgrid[x][y]->waterL0_abs = bl->getCellLeakage(x, y);
        // if (soilgrid[x][y]->waterL0_abs > 0) {
        //   leak = true;
        // }
      }
    }

    // ------------------------------------------
    // only if there is surface water, calculate infiltration
    // if (leak) {
    double meanWaterL1 = 0;
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        meanWaterL1 += soilgrid[x][y]->waterL1_rel;
      }
    }
    meanWaterL1 = meanWaterL1 / (xsize * ysize * 1.0);

    //---------------------------------------------
    /// INFILTRATION.
    /// Returns true, if there is still surface water.
    //---------------------------------------------
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        // If cell has surface water, calculate infiltration
        if (soilgrid[x][y]->waterL0_abs > 0) {
          /// \implements \fn WaterCell::twoLayerInfiltration(double
          /// meanWaterL1)
          surfaceWater = soilgrid[x][y]->twoLayerInfiltration(meanWaterL1);
          if (surfaceWater) {
            // if there is still surface water after infiltration:
            // Add the ponded water back to the surface water layer of the crust
            bl->addSurfaceWater(x, y, soilgrid[x][y]->waterL0_abs);
            soilgrid[x][y]->waterL0_abs = 0;
          }
        } else {
          soilgrid[x][y]->FL1 = 0;
          soilgrid[x][y]->FL2 = 0;
          soilgrid[x][y]->drainL1L2 = 0;
          soilgrid[x][y]->draindeep = 0;
        }
      }
    }
    //   leak = false;
    // //} else {
    //   for (int x = 0; x < xsize; x++) {
    //     for (int y = 0; y < ysize; y++) {
    //       soilgrid[x][y]->FL1 = 0;
    //       soilgrid[x][y]->FL2 = 0;
    //       soilgrid[x][y]->drainL1L2 = 0;
    //       soilgrid[x][y]->draindeep = 0;
    //     }
    //   }
    //}

    // check if there is still surface water
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        if (bl->getCellSurfaceWater(x, y) > 0) surfaceWaterFlag = true;
      }
    }

    //---------------------------------------------
    /// CALCULATE RUNOFF.
    /// Only if there is still surface water after infiltration calculate
    /// runoff.
    //---------------------------------------------
    if (surfaceWaterFlag) {
      /// \implements \fn WaterLandscape::calculateRunoff()
      calculateRunoff(bl);
      surfaceWaterFlag = false;
    } else {
      for (int x = 0; x < xsize; x++) {
        for (int y = 0; y < ysize; y++) {
          soilgrid[x][y]->QD = 0;
        }
      }
    }

    //---------------------------------------------
    /// LAYER DIFFUSION.
    /// Vertical diffusion of soil water between layers.
    //---------------------------------------------
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        /// \implements \fn WaterCell::layerDiffusion()
        soilgrid[x][y]->layerDiffusion();
      }
    }

    //---------------------------------------------
    /// POT EVAPOTRANSPIRATION.
    /// Evapotranspiration is calculated daily.
    //---------------------------------------------
    if (hour == 23) {
      /// \implements \fn WaterLandscape::potEvapotranspiration(int year, int
      /// day, Weather* w)
      potEvapotranspiration(year, day, w);
      meanEP = 0;
      // first calculate crust evaporation for all cells with crusts in
      bl->calculateProcessesDay(potEP);

      //---------------------------------------------
      /// ACT LAYER EVAPOTRANSPIRATION AND SEPERATION OF EVAPOTRANSPIRATION.
      /// Is the potential evapotranspiration noteworthy?
      //---------------------------------------------
      // is the current cell a crust cell or not
      bool crust_cell = false;
      // ep factor of the current cell. if 1 then it's a soil cell
      double cell_ep_factor = 1;
      if (potEP >= 0.00001) {
        for (int x = 0; x < xsize; x++) {
          for (int y = 0; y < ysize; y++) {
            if (bl->getCrusttype(x, y) == 0) {
              crust_cell = false;
              cell_ep_factor = soilgrid[x][y]->EPfactor;
            } else {
              crust_cell = true;
              cell_ep_factor = bl->getEPfactor(x, y);
            }
            /// \implements \fn
            /// WaterCell::actLayerEvapotranspiration(double
            /// rainTot, double potEP)
            soilgrid[x][y]->actLayerEvapotranspiration(
                w->precSum[daynumber], potEP, crust_cell, cell_ep_factor);
            /// \implements \fn WaterCell::transpirationupper()
            soilgrid[x][y]->transpirationupper(); /*Tong*/
            /// \implements \fn WaterCell::evaporationupper()
            soilgrid[x][y]->evaporationupper(); /*Tong*/
            /// \implements \fn WaterCell::transpirationlower()
            soilgrid[x][y]->transpirationlower(); /*Tong*/
            /// \implements \fn WaterCell::transpirationPFTU()
            soilgrid[x][y]->transpirationPFTU(); /*Tong*/
            /// \implements \fn WaterCell::transpirationPFTL()
            soilgrid[x][y]->transpirationPFTL(); /*Tong*/
          }
        }
      }
    }

    //---------------------------------------------
    /// TOTAL WATER.
    /// Calculate the total water amount.
    //---------------------------------------------
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        /// \implements \fn WaterCell::totalWater()
        soilgrid[x][y]->totalWater();
      }
    }

    //---------------------------------------------
    /// MEAN MOISTURE.
    /// Calculate mean values for soil moisture.
    /// \implements \fn WaterLandscape::meanMoisture()
    //---------------------------------------------
    meanMoisture();

    //---------------------------------------------
    /// SPATIAL ANALYSIS.
    /// Calculate soil moisture dependent on vegetation cover and on the number
    /// of cells that provide runoff to the respective cell (log-scale).
    /// \implements \fn WaterLandscape::spatialAnalysis(int year, int day, int
    /// hour)
    //---------------------------------------------
    // spatialAnalysis(year, day, hour);

    //---------------------------------------------
    /// MEAN PROCESSES.
    /// Calculate the spatial mean values for the calculated processes
    /// infiltration, drainage, evapotranspiration and runoff. \implements \fn
    /// WaterLandscape::meanProcesses()
    //---------------------------------------------;
    meanProcesses();
    meanProcesses_CS(bl);

    //---------------------------------------------
    /// WRITE OUTPUT FILE.
    /// \implements \fn WaterLandscape::writeOutputFile(int year, int day, int
    /// hour, Weather* weather)
    //---------------------------------------------
    // writeOutputFile(year, day, hour, w);
    writeOutputFileSelina(year, day, hour, w);
    bl->writeBiocrustOutput(year, day, hour, this->outputCells);

    output.spatialOutput(day, hour, bl, this);

    // reset runon for the next time step
    for (int i = 0; i < xsize; i++) {
      for (int j = 0; j < ysize; j++) {
        runon[i][j] = 0;
      }
    }

    //---------------------------------------------
    /// EXPORT LANDSCAPE FILE.
    /// \implements \fn WaterLandscape::exportCurrentMoisture(int year,
    /// Parameter* p)
    //---------------------------------------------
  }
  //---------------------------------------------
  /// EXPORT LANDSCAPE FILE.
  /// \implements \fn WaterLandscape::exportCurrentMoisture(int year, Parameter*
  /// p) KI
  //---------------------------------------------
  //    if (day == 364) {
  //        exportCurrentMoisture(<#int year#>, <#Parameter *p#>);
  //    }
}

/*******************************************************************************************
 * routine that is called from outside at a specific time step:
 * the current moisture values for each grid cell are exported
 * to be able to look at patterns
 * \todo Btodo: shift this to output
 * is called in --> NONE
 *******************************************************************************************/
void WaterLandscape::exportCurrentMoisture(int year, Parameter* p) {
  // create and open result file
  string y = std::to_string(year + 1);
  string currentMoistureFileName =
      "Results" + delimiter + "CurrentSoilMoisture_" + y + ".txt";
  fstream currentMoistureFile;
  currentMoistureFile.open(currentMoistureFileName.c_str(),
                           ios::out | ios::trunc);

  // if something goes wrong...
  if (!currentMoistureFile) {
    cerr << currentMoistureFileName
         << " could not be opened in waterLandscape->exportCurrentMoisture!\n";
    exit(-1);
  }
  //    currentMoistureFile << year << "Year" << "\t" ;             // add time
  //    component KI
  currentMoistureFile << "Surface Water [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentMoistureFile << soilgrid[x][y]->waterL0_abs << "\t";
    }
    currentMoistureFile << "\n";
  }
  currentMoistureFile << "\n \n";
  currentMoistureFile << "Soil Moisture in First Layer [%]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentMoistureFile << soilgrid[x][y]->waterL1_rel << "\t";
    }
    currentMoistureFile << "\n";
  }
  currentMoistureFile << "\n \n";
  currentMoistureFile << "Soil Moisture in Second Layer [%]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentMoistureFile << soilgrid[x][y]->waterL2_rel << "\t";
    }
    currentMoistureFile << "\n";
  }
  currentMoistureFile << "\n \n";
  currentMoistureFile << "Soil Moisture in Both Layers [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentMoistureFile << soilgrid[x][y]->totWater << "\t";
    }
    currentMoistureFile << "\n";
  }
  currentMoistureFile.close();
}

/******************************************************************************************
 * routine that is called from outside at a specific time step:
 * the current process rates / values for each grid cell are exported
 * to be able to look at patterns
 * \todo Btodo: shift this to output
 * is called in --> NONE
 ******************************************************************************************/
void WaterLandscape::exportCurrentProcesses(int year, Parameter* p) {
  // create and open result file
  string y = std::to_string(year + 1);
  string currentProcessFileName =
      "Results" + delimiter + "CurrentProcesses" + y + ".txt";
  fstream currentProcessFile;
  currentProcessFile.open(currentProcessFileName.c_str(),
                          ios::out | ios::trunc);
  double waterGain;

  // if something goes wrong...
  if (!currentProcessFile) {
    cerr << currentProcessFileName
         << " could not be opened in waterLandscape->exportCurrentProcesses!\n";
    exit(-1);
  }

  currentProcessFile << "ACTUAL RUNOFF\n \n";
  currentProcessFile << "Runoff [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentProcessFile << soilgrid[x][y]->QD << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "Runon [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentProcessFile << runon[x][y] << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "Runon - Runoff (Total Gain) [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      waterGain = runon[x][y] - soilgrid[x][y]->QD;
      currentProcessFile << waterGain << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "EVAPOTRANSPIRATION\n \n";
  currentProcessFile << "Evapotranspiration from Surface [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      // divided by 24 since it is daily data
      currentProcessFile << (soilgrid[x][y]->EPL0) / 24.0 << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "Evapotranspiration from Upper Soil Layer [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      // divided by 24 since it is daily data
      currentProcessFile << soilgrid[x][y]->EPL1 / 24.0 << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "Evapotranspiration from Lower Soil Layer [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      // divided by 24 since it is daily data
      currentProcessFile << soilgrid[x][y]->EPL2 / 24.0 << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "Total Evapotranspiration [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      // divided by 24 since it is daily data
      currentProcessFile << soilgrid[x][y]->EPtot / 24.0 << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "Infiltration L1 [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentProcessFile << soilgrid[x][y]->FL1 << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile << "Infiltration L2 [mm]\n \n";

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      currentProcessFile << soilgrid[x][y]->FL2 << "\t";
    }
    currentProcessFile << "\n";
  }
  currentProcessFile << "\n \n";
  currentProcessFile.close();
}

/*******************************************************************************************
 * calculate the spatial mean values for the calculated processes infiltration,
 * drainage, evapotranspiration and runoff
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterLandscape::meanProcesses() { /*Tong, based on PFT-oriented*/

  // set mean values to 0
  meanFL1 = meanFL2 = meanFtot = 0;
  meanDrainL1L2 = meanDrainDeep = meanLayerDiff = 0;
  meanQD = QDlost = 0;
  meanEPL0 = meanEPL1 = meanEPL2 = meanEP = 0;
  meanTL1 = meanTL2 = meanEL1 = 0; /*Tong*/

  sdFL1 = sdFL2 = sdQD = sdEP = sdEPL1 = sdEPL2 = 0.0;

  for (int m = 0; m < PFTs; m++) { /*Tong*/
    meanPftTL1[m] = 0;
    meanPftTL2[m] = 0;
    meanPftCover[m] = 0;
  }

  // calculate the sums (removing the border cells and only the cells that
  // are covered by the output mask)

  for (int i = 0 + 1; i < (xsize ); i++) {
    for (int j = 0 + 1; j < (ysize ); j++) {
      if (outputCells[i][j]) {
        meanFL1 += soilgrid[i][j]->FL1;
        meanFL2 += soilgrid[i][j]->FL2;
        meanFtot += (soilgrid[i][j]->FL1 + soilgrid[i][j]->FL2);
        meanDrainL1L2 += soilgrid[i][j]->drainL1L2;
        meanDrainDeep += soilgrid[i][j]->draindeep;
        meanQD += soilgrid[i][j]->QD;
        meanEPL0 += (soilgrid[i][j]->EPL0);
        meanEPL1 += (soilgrid[i][j]->EPL1);
        meanEPL2 += (soilgrid[i][j]->EPL2);
        meanEP += (soilgrid[i][j]->EPL0 + soilgrid[i][j]->EPL1 +
                   soilgrid[i][j]->EPL2);
        meanTL1 += (soilgrid[i][j]->TL1); /*Tong*/
        meanTL2 += (soilgrid[i][j]->TL2); /*Tong*/
        meanEL1 += (soilgrid[i][j]->EL1); /*Tong*/

        for (int m = 0; m < PFTs; m++) { /*Tong*/
          meanPftTL1[m] += (soilgrid[i][j]->PftTL1[m]);
          meanPftTL2[m] += (soilgrid[i][j]->PftTL2[m]);
          meanPftCover[m] += (soilgrid[i][j]->coverPFT[m]);
        }
        meanLayerDiff += (soilgrid[i][j]->layerDiff);
      }
    }
  }

  // lost runoff at border cells
  for (int i = 0; i < xsize; i++) {
    if (runoffdirection[i][0] == -1) QDlost += soilgrid[i][0]->QD;
    if (runoffdirection[i][ysize - 1] == -1)
      QDlost += soilgrid[i][ysize - 1]->QD;
  }
  for (int j = 1; j < ysize - 1; j++) {
    if (runoffdirection[0][j] == -1) QDlost += soilgrid[0][j]->QD;
    if (runoffdirection[xsize - 1][j] == -1)
      QDlost += soilgrid[xsize - 1][j]->QD;
  }

  // number of cells without the border cells
  // double numcells = (xsize - 2) * (ysize - 2) * 1.0;
  int numcells = 0;
  for (auto& n : outputCells) {
    for (auto& n2 : n) {
      numcells += n2;
    }
  }

  meanFL1 = meanFL1 / numcells;
  meanFL2 = meanFL2 / numcells;
  meanFtot = meanFtot / numcells;
  meanDrainL1L2 = meanDrainL1L2 / numcells;
  meanDrainDeep = meanDrainDeep / numcells;
  meanQD = meanQD / numcells;
  meanEPL0 = meanEPL0 / numcells;
  meanEPL1 = meanEPL1 / numcells;
  meanEPL2 = meanEPL2 / numcells;
  meanEP = meanEP / numcells;   /*Tong*/
  meanTL1 = meanTL1 / numcells; /*Tong*/
  meanTL2 = meanTL2 / numcells; /*Tong*/

  for (int m = 0; m < PFTs; m++) { /*Tong*/
    meanPftTL1[m] = meanPftTL1[m] / numcells;
    meanPftTL2[m] = meanPftTL2[m] / numcells;
    meanPftCover[m] = meanPftCover[m] / numcells;
  }

  meanEL1 = meanEL1 / numcells; /*Tong*/
  meanLayerDiff = meanLayerDiff / numcells;

  // calculate sd
  for (int i = 0 + 1; i < (xsize); i++) {
    for (int j = 0 + 1; j < (ysize); j++) {
      if (outputCells[i][j]) {
        sdFL1 = pow(soilgrid[i][j]->FL1 - meanFL1, 2);
        sdFL2 = pow(soilgrid[i][j]->FL2 - meanFL2, 2);
        sdQD = pow(soilgrid[i][j]->QD - meanQD, 2);
        sdEP = pow(soilgrid[i][j]->EPtot - meanEP, 2);
        sdEPL1 = pow(soilgrid[i][j]->EPL1 - meanEPL1, 2);
        sdEPL2 = pow(soilgrid[i][j]->EPL2 - meanEPL2, 2);
      }
    }
  }

  sdFL1 = sqrt(sdFL1 / numcells);
  sdFL2 = sqrt(sdFL2 / numcells);
  sdQD = sqrt(sdQD / numcells);
  sdEP = sqrt(sdEP / numcells);
  sdEPL1 = sqrt(sdEPL1 / numcells);
  sdEPL2 = sqrt(sdEPL2 / numcells);
}

/***************************************************************************************
   Calculate the spatial mean values for the calculated processes
   *infiltration, drainage, evapotranspiration and runoff.
   *******************************************************************************************/
void WaterLandscape::meanProcesses_CS(BiocrustLandscape* bl) {
  // set mean values to 0
  meanWaterL0_abs_crust0 = meanWaterL0_abs_crust1 = meanWaterL0_abs_crust2 =
      meanWaterL0_abs_crust3 = meanWaterL0_abs_veg = 0;
  meanWaterL1_rel_crust0 = meanWaterL1_rel_crust1 = meanWaterL1_rel_crust2 =
      meanWaterL1_rel_crust3 = meanWaterL1_rel_veg = 0;
  meanWaterL2_rel_crust0 = meanWaterL2_rel_crust1 = meanWaterL2_rel_crust2 =
      meanWaterL2_rel_crust3 = meanWaterL2_rel_veg = 0;
  meanFL1_crust0 = meanFL1_crust1 = meanFL1_crust2 = meanFL1_crust3 =
      meanFL1_veg = 0;
  meanDrainL1L2_crust0 = meanDrainL1L2_crust1 = meanDrainL1L2_crust2 =
      meanDrainL1L2_crust3 = meanDrainL1L2_veg = 0;
  meanDrainDeep_crust0 = meanDrainDeep_crust1 = meanDrainDeep_crust2 =
      meanDrainDeep_crust3 = meanDrainDeep_veg = 0;
  meanQD_crust0 = meanQD_crust1 = meanQD_crust2 = meanQD_crust3 = meanQD_veg =
      0;
  meanEPL1_crust0 = meanEPL1_crust1 = meanEPL1_crust2 = meanEPL1_crust3 =
      meanEPL1_veg = 0;

  // calculate the sums (removing the border cells)
  int count_vegCells = 0;
  int count_crustCells1 = 0;
  int count_crustCells2 = 0;
  int count_crustCells3 = 0;
  int count_crustCells0 = 0;

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {  
      if (outputCells[i][j]) {
        // if cell is vegetated
        if (soilgrid[i][j]->coverPFT[0] + soilgrid[i][j]->coverPFT[1] > 0) {
          count_vegCells++;
          meanWaterL0_abs_veg += soilgrid[i][j]->waterL0_abs;
          meanWaterL1_rel_veg += soilgrid[i][j]->waterL1_rel;
          meanWaterL2_rel_veg += soilgrid[i][j]->waterL2_rel;
          meanFL1_veg += soilgrid[i][j]->FL1;
          meanDrainL1L2_veg += soilgrid[i][j]->drainL1L2;
          meanDrainDeep_veg += soilgrid[i][j]->draindeep;
          meanQD_veg += soilgrid[i][j]->QD;
          meanEPL1_veg += (soilgrid[i][j]->EPL1);
        } else {
          if (bl->getCrusttype(i, j) == 0) {
            count_crustCells0++;
            meanWaterL0_abs_crust0 += soilgrid[i][j]->waterL0_abs;
            meanWaterL1_rel_crust0 += soilgrid[i][j]->waterL1_rel;
            meanWaterL2_rel_crust0 += soilgrid[i][j]->waterL2_rel;
            meanFL1_crust0 += soilgrid[i][j]->FL1;
            meanDrainL1L2_crust0 += soilgrid[i][j]->drainL1L2;
            meanDrainDeep_crust0 += soilgrid[i][j]->draindeep;
            meanQD_crust0 += soilgrid[i][j]->QD;
            meanEPL1_crust0 += (soilgrid[i][j]->EPL1);
          } else if (bl->getCrusttype(i, j) == 1) {
            count_crustCells1++;
            meanWaterL0_abs_crust1 += soilgrid[i][j]->waterL0_abs;
            meanWaterL1_rel_crust1 += soilgrid[i][j]->waterL1_rel;
            meanWaterL2_rel_crust1 += soilgrid[i][j]->waterL2_rel;
            meanFL1_crust1 += soilgrid[i][j]->FL1;
            meanDrainL1L2_crust1 += soilgrid[i][j]->drainL1L2;
            meanDrainDeep_crust1 += soilgrid[i][j]->draindeep;
            meanQD_crust1 += soilgrid[i][j]->QD;
            meanEPL1_crust1 += (soilgrid[i][j]->EPL1);
          } else if (bl->getCrusttype(i, j) == 2) {
            count_crustCells2++;
            meanWaterL0_abs_crust2 += soilgrid[i][j]->waterL0_abs;
            meanWaterL1_rel_crust2 += soilgrid[i][j]->waterL1_rel;
            meanWaterL2_rel_crust2 += soilgrid[i][j]->waterL2_rel;
            meanFL1_crust2 += soilgrid[i][j]->FL1;
            meanDrainL1L2_crust2 += soilgrid[i][j]->drainL1L2;
            meanDrainDeep_crust2 += soilgrid[i][j]->draindeep;
            meanQD_crust2 += soilgrid[i][j]->QD;
            meanEPL1_crust2 += (soilgrid[i][j]->EPL1);
          } else if (bl->getCrusttype(i, j) == 3) {
            count_crustCells3++;
            meanWaterL0_abs_crust3 += soilgrid[i][j]->waterL0_abs;
            meanWaterL1_rel_crust3 += soilgrid[i][j]->waterL1_rel;
            meanWaterL2_rel_crust3 += soilgrid[i][j]->waterL2_rel;
            meanFL1_crust3 += soilgrid[i][j]->FL1;
            meanDrainL1L2_crust3 += soilgrid[i][j]->drainL1L2;
            meanDrainDeep_crust3 += soilgrid[i][j]->draindeep;
            meanQD_crust3 += soilgrid[i][j]->QD;
            meanEPL1_crust3 += (soilgrid[i][j]->EPL1);
          }
        }
      }
    }
  }

  // mean of vegetated cells
  if (count_vegCells != 0) {
    meanWaterL0_abs_veg = meanWaterL0_abs_veg / count_vegCells;
    meanWaterL1_rel_veg = meanWaterL1_rel_veg / count_vegCells;
    meanWaterL2_rel_veg = meanWaterL2_rel_veg / count_vegCells;
    meanFL1_veg = meanFL1_veg / count_vegCells;
    meanDrainL1L2_veg = meanDrainL1L2_veg / count_vegCells;
    meanDrainDeep_veg = meanDrainDeep_veg / count_vegCells;
    meanQD_veg = meanQD_veg / count_vegCells;
    meanEPL1_veg = meanEPL1_veg / count_vegCells;
  }
  // mean of crust 0 cells
  if (count_crustCells0 != 0) {
    meanQD_crust0 = meanQD_crust0 / count_crustCells0;
    meanWaterL0_abs_crust0 = meanWaterL0_abs_crust0 / count_crustCells0;
    meanWaterL1_rel_crust0 = meanWaterL1_rel_crust0 / count_crustCells0;
    meanWaterL2_rel_crust0 = meanWaterL2_rel_crust0 / count_crustCells0;
    meanFL1_crust0 = meanFL1_crust0 / count_crustCells0;
    meanDrainL1L2_crust0 = meanDrainL1L2_crust0 / count_crustCells0;
    meanDrainDeep_crust0 = meanDrainDeep_crust0 / count_crustCells0;
    meanEPL1_crust0 = meanEPL1_crust0 / count_crustCells0;
  }
  // mean of crust 1 cells
  if (count_crustCells1 != 0) {
    meanQD_crust1 = meanQD_crust1 / count_crustCells1;
    meanWaterL0_abs_crust1 = meanWaterL0_abs_crust1 / count_crustCells1;
    meanWaterL1_rel_crust1 = meanWaterL1_rel_crust1 / count_crustCells1;
    meanWaterL2_rel_crust1 = meanWaterL2_rel_crust1 / count_crustCells1;
    meanFL1_crust1 = meanFL1_crust1 / count_crustCells1;
    meanDrainL1L2_crust1 = meanDrainL1L2_crust1 / count_crustCells1;
    meanDrainDeep_crust1 = meanDrainDeep_crust1 / count_crustCells1;
    meanEPL1_crust1 = meanEPL1_crust1 / count_crustCells1;
  }
  // mean of crust 2 cells
  if (count_crustCells2 != 0) {
    meanQD_crust2 = meanQD_crust2 / count_crustCells2;
    meanWaterL0_abs_crust2 = meanWaterL0_abs_crust2 / count_crustCells2;
    meanWaterL1_rel_crust2 = meanWaterL1_rel_crust2 / count_crustCells2;
    meanWaterL2_rel_crust2 = meanWaterL2_rel_crust2 / count_crustCells2;
    meanFL1_crust2 = meanFL1_crust2 / count_crustCells2;
    meanDrainL1L2_crust2 = meanDrainL1L2_crust2 / count_crustCells2;
    meanDrainDeep_crust2 = meanDrainDeep_crust2 / count_crustCells2;
    meanEPL1_crust2 = meanEPL1_crust2 / count_crustCells2;
  }
  // mean of crust 3 cells
  if (count_crustCells3 != 0) {
    meanQD_crust3 = meanQD_crust3 / count_crustCells3;
    meanWaterL0_abs_crust3 = meanWaterL0_abs_crust3 / count_crustCells3;
    meanWaterL1_rel_crust3 = meanWaterL1_rel_crust3 / count_crustCells3;
    meanWaterL2_rel_crust3 = meanWaterL2_rel_crust3 / count_crustCells3;
    meanFL1_crust3 = meanFL1_crust3 / count_crustCells3;
    meanDrainL1L2_crust3 = meanDrainL1L2_crust3 / count_crustCells3;
    meanDrainDeep_crust3 = meanDrainDeep_crust3 / count_crustCells3;
    meanEPL1_crust3 = meanEPL1_crust3 / count_crustCells3;
  }
}
/*******************************************************************************************
 * calculate the mean soil moisture of the grid
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterLandscape::meanMoisture() {
  meanWaterL1_rel = meanWaterL2_rel = 0;
  meanWaterL0_abs = meanWaterL1_abs = meanWaterL2_abs = 0;
  meanWater_rel = meanWater_abs = 0;
  sdWaterL1_rel = sdWaterL2_rel = 0.0;

  // sum up all values for grid (without border cells)
  for (int i = 0 + 1; i < (xsize); i++) {
    for (int j = 0 + 1; j < (ysize); j++) {
      meanWaterL0_abs += soilgrid[i][j]->waterL0_abs;
      meanWaterL1_abs += soilgrid[i][j]->waterL1_abs;
      meanWaterL2_abs += soilgrid[i][j]->waterL2_abs;
      meanWaterL1_rel += soilgrid[i][j]->waterL1_rel;
      meanWaterL2_rel += soilgrid[i][j]->waterL2_rel;
      meanWater_rel +=
          (soilgrid[i][j]->waterL1_abs + soilgrid[i][j]->waterL2_abs) /
          (soilgrid[i][j]->depthL1 + soilgrid[i][j]->depthL2);
    }
  }

  // calculate the grid means

  // number of grid cells to average over (number of cells without border
  // cells)
  double numcells = (double)((xsize ) * (ysize ));

  meanWaterL0_abs = meanWaterL0_abs / (numcells);
  meanWaterL1_abs = meanWaterL1_abs / (numcells);
  meanWaterL2_abs = meanWaterL2_abs / (numcells);

  meanWaterL1_rel = meanWaterL1_rel / (numcells);
  meanWaterL2_rel = meanWaterL2_rel / (numcells);

  meanWater_abs = meanWaterL1_abs + meanWaterL2_abs;
  meanWater_rel = meanWater_rel / (numcells);

  // calculate the standard deviation

  for (int i = 0; i < (xsize ); i++) {
    for (int j = 0; j < (ysize ); j++) {
      sdWaterL1_rel += pow(soilgrid[i][j]->waterL1_rel - meanWaterL1_rel, 2);
      sdWaterL2_rel += pow(soilgrid[i][j]->waterL2_rel - meanWaterL2_rel, 2);
    }
  }

  sdWaterL1_rel = sqrt(sdWaterL1_rel / numcells);
  sdWaterL2_rel = sqrt(sdWaterL2_rel / numcells);
}

/*******************************************************************************************
 * sum up moisture and runoff for each day and month. Will be written as
 *output in modelOutput
 *******************************************************************************************/
void WaterLandscape::calculateCellMeans() {
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      // cout << monthly_soil_moistureL1[i][j] << endl;
      monthly_soil_moistureL1[i][j] += soilgrid[i][j]->waterL1_rel;
      monthly_soil_moistureL2[i][j] += soilgrid[i][j]->waterL2_rel;
      monthly_QD_cell[i][j] += soilgrid[i][j]->QD;
      monthly_infiltrationL1_cell[i][j] += soilgrid[i][j]->FL1;
      monthly_runon_cell[i][j] += runon[i][j];
      monthly_deepdrain_cell[i][j] += soilgrid[i][j]->draindeep;
      monthly_evapotranspiration_L1_cell[i][j] += soilgrid[i][j]->EPL1;
      monthly_evapotranspiration_L2_cell[i][j] += soilgrid[i][j]->EPL2;
      monthly_evapotranspiration_tot_cell[i][j] += soilgrid[i][j]->EPtot;
      daily_soil_moistureL1[i][j] += soilgrid[i][j]->waterL1_rel;
      daily_soil_moistureL2[i][j] += soilgrid[i][j]->waterL2_rel;
      daily_QD_cell[i][j] += soilgrid[i][j]->QD;
      daily_infiltrationL1_cell[i][j] += soilgrid[i][j]->FL1;
      daily_runon_cell[i][j] += runon[i][j];
      daily_deepdrain_cell[i][j] += soilgrid[i][j]->draindeep;
      daily_evapotranspiration_L1_cell[i][j] += soilgrid[i][j]->EPL1;
      daily_evapotranspiration_L2_cell[i][j] += soilgrid[i][j]->EPL2;
      daily_evapotranspiration_tot_cell[i][j] += soilgrid[i][j]->EPtot;
    }
  }
}

/*******************************************************************************************
 * run routine to calculate the potential evapotranspiration
 * according to Hargreaves (1974)
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterLandscape::potEvapotranspiration(int year, int day, Weather* w) {
  int daynumber = year * daysPerYear + day;

  // General values
  if (w->tempAvg[daynumber] != nodata) {
    // Hargreaves
    potEP = 0.0023 * (w->tempAvg[daynumber] + 17.8) *
            pow((w->tempMax[daynumber] - w->tempMin[daynumber]), 0.5) *
            w->extraterrRadiation[daynumber];
  }

  else
    potEP = 0;
}

/*******************************************************************************************
 * follow the actual runoff to the next "knot" and fill the list in the right
 *order this procedure is called at the beginning of the simulation to find
 *out the correct order of the cells is called in -->
 *WaterLandscape::calculateRunonCells
 *******************************************************************************************/
void WaterLandscape::followShortRunoffPath(int xpos, int ypos) {
  Element* el = new Element(xpos, ypos, numberOfFlows[xpos][ypos], -1);
  list2->insertEnd(el);
  calculated[xpos][ypos] = true;
  if ((runoffdirection[xpos][ypos] != -1) &&
      (numberOfFlows[(runoffdirection[xpos][ypos] -
                      runoffdirection[xpos][ypos] % ysize) /
                     ysize][runoffdirection[xpos][ypos] % ysize] ==
       numberOfFlows[xpos][ypos]) &&
      (calculated[(runoffdirection[xpos][ypos] -
                   runoffdirection[xpos][ypos] % ysize) /
                  ysize][runoffdirection[xpos][ypos] % ysize] == false)) {
    followShortRunoffPath(
        (runoffdirection[xpos][ypos] - runoffdirection[xpos][ypos] % ysize) /
            ysize,
        runoffdirection[xpos][ypos] % ysize);
  }
}

/*******************************************************************************************
 * add the actual precipitation to the surface water
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterLandscape::precipitation(double actualPrec) {
  // this case should not happen, but you never know...
  // no data is treated as zero
  if (actualPrec < 0) actualPrec = 0;

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      soilgrid[i][j]->waterL0_abs += actualPrec;
    }
  }

  // flag is set to true, if there is surface water
  if (actualPrec > 0) {
    surfaceWaterFlag = true;
  }
}

/*******************************************************************************************
 * distribution of actual runoff to each cell
 * in this routine, the order of the cells in the list is followed;
 * a cell does not get runoff until all cells contributing runoff to this
 * cell have been evaluated
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterLandscape::calculateRunoff(BiocrustLandscape* bl) {
  // the following route follows the order of the cells in the list,
  // first the runon to cells with only one flow contributing to the
  // runon are calculated, afterwards the ones with two flows etc
  // this guaranties, that cells are not double counted and that
  // water can already infiltrate into upslope cells before it reaches
  // downslope cells

  int xpos, ypos;

  // follow sorted list and calculate the infiltrated water in the
  // cell and the water that is passed to the next cell
  list2->reset();
  while (!list2->endpos()) {
    list2->advance();

    xpos = list2->position->xval;
    ypos = list2->position->yval;

    runon[xpos][ypos] = 0;

    // it is only necessary to evaluate the runoff that reaches a cell, if
    // there are cells, that distribute runoff to the actual cell
    if (runoncells[xpos][ypos] != 0) {
      // look at the surrounding cells and calculate the total runon
      for (int k = max(0, xpos - 1); k <= min(xpos + 1, xsize - 1); k++) {
        for (int l = max(0, ypos - 1); l <= min(ypos + 1, ysize - 1); l++) {
          // if the cell provides runoff to this cell and it is not the cell
          // itself and it has runoff
          if ((runoffdirection[k][l] == xpos * ysize + ypos) &&
              !(k == xpos && l == ypos) && (soilgrid[k][l]->QD > 0)) {
            // does this make sense?
            runon[xpos][ypos] += soilgrid[k][l]->QD;
          }
        }
      }
      bl->addSurfaceWater(xpos, ypos, runon[xpos][ypos]);
      // soilgrid[xpos][ypos]->waterL0_abs += runon[xpos][ypos];
    }

    double cStotal = 0; /*Tong*/
    double cGtotal = 0;
    double cAtotal = 0;

    for (int i = 0; i < PFTs; i++) { /*Tong*/
      if (i < Shrubs)
        cStotal += soilgrid[xpos][ypos]->coverPFT[i];
      else if (i >= Shrubs && i < Shrubs + Perennials)
        cGtotal += soilgrid[xpos][ypos]->coverPFT[i];
      else
        cAtotal += soilgrid[xpos][ypos]->coverPFT[i];
    }

    // Routine according to Manning-Strickler
    if (soilgrid[xpos][ypos]->inclination < 0) {
      cout << " Problem bei Manning-Strickler"
           << soilgrid[xpos][ypos]->inclination << "\n";
    }
    // if it's a crust cell, retrieve Manning's n of crust
    double n = 1;
    if (bl->getCrusttype(xpos, ypos) > 0 && cGtotal + cAtotal + cStotal == 0) {
      n = 0.28;
      // n = bl->getManningn(xpos, ypos);
    } else {
      n = 1;  // 0.15 + 0.25 * (cGtotal + cAtotal + cStotal);
    }
    // cout << "cGtotal " << cGtotal <<  "cAtotal" << cAtotal << " cStotal" <<
    // cStotal<<endl;
    soilgrid[xpos][ypos]->QD =
        1 / n * pow(bl->getCellSurfaceWater(xpos, ypos), (2.0 / 3.0)) *
        (1 - (cGtotal + cAtotal + cStotal) / 2.0) *
        sqrt(soilgrid[xpos][ypos]->inclination);
    soilgrid[xpos][ypos]->QD = max(0.0, min(bl->getCellSurfaceWater(xpos, ypos),
                                            soilgrid[xpos][ypos]->QD));
    bl->addSurfaceWater(xpos, ypos, -soilgrid[xpos][ypos]->QD);
    // soilgrid[xpos][ypos]->waterL0_abs -= soilgrid[xpos][ypos]->QD;
  }

  // set the flag for surface water:
  // false = no surface water left, true = still surface water
  // surfaceWaterFlag = false;
  // for (int i = 0; i < xsize; i++) {
  //   for (int j = 0; j < ysize; j++) {
  //     if (bl->getCellSurfaceWater(xpos, ypos) > 0) surfaceWaterFlag = true;
  //   }
  // }
}

/*******************************************************************************************
 * follow the runoff path
 * this is a recursive method, it stops, if a cell is reached which does not
 * distribute runoff to other cells;
 * this routine is only called at the beginning of the simulation to
 * find out the order of the cells
 * is called in --> WaterLandscape::calculateRunonCells
 *******************************************************************************************/
void WaterLandscape::followRunoff(int i, int j, int addrunon, bool first) {
  // add the additional runon to the cell
  runoncells[i][j] += addrunon;
  numberOfFlows[i][j] += 1;

  // stop condition: no runoff to other cells
  if (runoffdirection[i][j] == -1) return;
  // stop condition: the routine was called from outside with a new
  // cell, but the cell has already distributed its runoff
  if (first && calculated[i][j]) return;

  // if the runoff of this cell was not calculated yet,
  // the runoff of the cell itself is added
  if (!calculated[i][j]) {
    addrunon = runoncells[i][j] + 1;
    calculated[i][j] = true;
  }

  // call the method again for the cell where the runoff is distributed to
  int x = (runoffdirection[i][j] - runoffdirection[i][j] % ysize) / ysize;
  int y = runoffdirection[i][j] % ysize;

  followRunoff(x, y, addrunon, false);
}

/*******************************************************************************************
 * write the number of cells that contribute runoff to each cell in a file
 * this file is used to calculate the initial vegetation distribution:
 * water availability is assessed via amount of runoff
 * \todo Btodo: shift this to output
 * can be called in --> WaterLandscape::calculateRunonCells
 *******************************************************************************************/
void WaterLandscape::writeRunonCells() {
  // create and open result file
  string runonfn = "Results" + delimiter + "cellsContributingWater.txt";
  ofstream runonf;
  runonf.open(runonfn.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!runonf) {
    cerr << runonfn << " could not be opened in main!\n";
    exit(-1);
  }
  runonf << "Cells contributing runoff to each cell\n \n";
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      if (j != 0) runonf << "\t";
      runonf << runoncells[i][j];
    }
    runonf << "\n";
  }
  runonf << "\n"
         << "\n";
  runonf << "runoff direction \n";
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      if (j != 0) runonf << "\t";
      runonf << runoffdirection[i][j];
    }
    runonf << "\n";
  }
  runonf << "\n"
         << "\n";
runonf << "inclination \n";
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      if (j != 0) runonf << "\t";
      runonf << inclination[i][j];
    }
    runonf << "\n";
  }
  runonf.close();
}

/*******************************************************************************************
 * calculate the number of cells that produce runon for a specific cell
 * this routine is only called at the beginning of the simulation
 * is called in --> WaterLandscape::initialiseSoilParameters
 *******************************************************************************************/
void WaterLandscape::calculateRunonCells() {
  // help variable to store coordinates of so far lower cell
  int x_lower;
  int y_lower;

  // direction in which the runoff of a cell flows
  // integer number indicates absolute number of cell,
  // e.g. cell [x][y]->ysize*x + y = runoffdirection
  //		y = runoffdirection%ysize
  //		x = (runoffdirection - y)/ysize

  // initialise values and calculate runoff direction

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      calculated[i][j] = false;
      runoncells[i][j] = 0;
      numberOfFlows[i][j] = -1;

      x_lower = i;
      y_lower = j;

      // look at surrounding 8 cells (plus itself)
      for (int k = max(0, i - 1); k <= min(i + 1, xsize - 1); k++) {
        for (int l = max(0, j - 1); l <= min(j + 1, ysize - 1); l++) {
          if (elevation[k][l] < elevation[x_lower][y_lower]) {
            x_lower = k;
            y_lower = l;
          }
        }
      }

      if ((x_lower == i) && (y_lower == j)) {
        runoffdirection[i][j] =
            -1;  // this cell does not distribute water to other cells
      } else {
        runoffdirection[i][j] = x_lower * ysize + y_lower;
      }
    }
  }

  // calculate number of cells that produce runon to a cell
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      followRunoff(i, j, 0, true);
    }
  }

  list1->delAll();
  int helpFlow;
  // make new list in order of number of flows that contribute to the runon
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      helpFlow = 0;
      for (int k = max(0, i - 1); k <= min(i + 1, xsize - 1); k++) {
        for (int l = max(0, j - 1); l <= min(j + 1, ysize - 1); l++) {
          if ((k != i) && (l != j) &&
              ((runoffdirection[k][l] - runoffdirection[k][l] % ysize) /
                   ysize ==
               i) &&
              (runoffdirection[k][l] % ysize == j)) {
            helpFlow = numberOfFlows[k][l];
          }
        }
      }
      Element* el = new Element(i, j, numberOfFlows[i][j], helpFlow);
      list1->insertSorted(el);
    }
  }

  int xpos, ypos;
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      calculated[i][j] = false;
    }
  }

  list2->delAll();
  list1->reset();
  while (!list1->endpos()) {
    list1->advance();
    xpos = list1->position->xval;
    ypos = list1->position->yval;

    // if the runon to this cell is not calulated yet
    if (calculated[xpos][ypos] == false) {
      followShortRunoffPath(xpos, ypos);
    }
  }

}

/*******************************************************************************************
 * set elevation of each cell
 * this routine is called at the beginning of the simulation to
 * initialise the landscape.
 * is called in --> WaterLandscape::initialiseSoilParameters
 *******************************************************************************************/
void WaterLandscape::setElevation(Parameter* p) {
  string parameterFileName;
  ifstream parameterFile;

  parameterFileName = "Parameters" + delimiter + "elevation_" +
                      std::to_string(xsize) + "_" +
                      p->elevationFileID[p->scenario] + ".txt";

  cout << "Read elevation from: " << parameterFileName << endl;

  parameterFile.open(parameterFileName.c_str(), ios::binary | ios::in);

  // if something goes wrong...
  if (!parameterFile) {
    cerr << parameterFileName << " could not be opened!\n";
    exit(-1);
  }

  // read in parameters and close file
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      parameterFile >> elevation[i][j];
    }
  }
  parameterFile.close();
}

/*******************************************************************************************
 * set elevation of each cell
 * this routine is called at the beginning of the simulation to
 * initialise the landscape.
 * is called in --> WaterLandscape::initialiseSoilParameters
 *******************************************************************************************/
void WaterLandscape::setOutputCells() {
  string parameterFileName;
  ifstream parameterFile;

  parameterFileName =
      "Parameters" + delimiter + "outputMask_" + std::to_string(xsize) + ".txt";

  // check if the file exists
  std::filesystem::path parPath{parameterFileName};
  if (std::filesystem::exists(parPath)) {
    cout << "Read output Mask from: " << parameterFileName << endl;

    parameterFile.open(parameterFileName.c_str(), ios::binary | ios::in);

    // if something goes wrong...
    if (!parameterFile) {
      cerr << parameterFileName << " could not be opened!\n";
      exit(-1);
    }
    // read in parameters and close file
    for (int i = 0; i < xsize; i++) {
      for (int j = 0; j < ysize; j++) {
        parameterFile >> outputCells[i][j];
      }
    }
    parameterFile.close();

  } else {
    cout << "File " << parameterFileName
         << " does not exist. Outputting all cells without a mask" << endl;

    for (int i = 0; i < xsize; i++) {
      for (int j = 0; j < ysize; j++) {
        outputCells[i][j] = 1;
      }
    }
  }
  int count_in = 0;
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      count_in += outputCells[i][j];
    }
  }
  cout << "output cells: " << count_in << endl;
}

/*******************************************************************************************
 * set aspect of each cell
 * this routine is called at the beginning of the simulation to
 * initialise the landscape
 * is called in --> WaterLandscape::initialiseSoilParameters
 *******************************************************************************************/
void WaterLandscape::setAspectAndInclination() {
  int x = 0;
  int y = 0;

  // ofstream inclinationFile;
  // string inclinationFile_name = "Results" + delimiter + "inclination.txt";

  // inclinationFile.open(inclinationFile_name.c_str(), ios::out | ios::trunc);

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      if (runoffdirection[i][j] != -1) {
        // find out x and y coordinates of lowest neighbouring cell
        x = (runoffdirection[i][j] - runoffdirection[i][j] % ysize) / ysize;
        y = runoffdirection[i][j] % ysize;
        // calculate inclination: arctangens of opposite leg / adjacent leg
        //							= atan (
        //(elevation1-elevation2)/length ) in radian
        inclination[i][j] =
            atan((elevation[i][j] - elevation[x][y]) / (cellsize * 1.0));
      } else {  // the inclination is set equal to the inclination to the
                // neighbour cell
        if (i == 0)
          inclination[i][j] =
              atan((elevation[1][j] - elevation[0][j]) / (cellsize * 1.0));
        else if (i == xsize - 1)
          inclination[i][j] =
              atan((elevation[xsize - 2][j] - elevation[xsize - 1][j]) /
                   (cellsize * 1.0));
        else if (j == 0)
          inclination[i][j] =
              atan((elevation[i][1] - elevation[i][0]) / (cellsize * 1.0));
        else if (j == ysize - 1)
          inclination[i][j] =
              atan((elevation[i][ysize - 2] - elevation[i][ysize - 1]) /
                   (cellsize * 1.0));
        else
          inclination[i][j] = 0;
        if ((i == 0) && (j == 0))
          inclination[i][j] =
              atan((elevation[1][1] - elevation[i][j]) / (cellsize * 1.0));
        if ((i == 0) && (j == ysize - 1))
          inclination[i][j] = atan((elevation[1][ysize - 2] - elevation[i][j]) /
                                   (cellsize * 1.0));
        if ((i == xsize - 1) && (j == 0))
          inclination[i][j] = atan((elevation[xsize - 2][1] - elevation[i][j]) /
                                   (cellsize * 1.0));
        if ((i == xsize - 1) && (j == ysize - 1))
          inclination[i][j] =
              atan((elevation[xsize - 2][ysize - 2] - elevation[i][j]) /
                   (cellsize * 1.0));
      }

      // inclination in degree
      inclination_deg[i][j] = inclination[i][j] * 180.0 / PI;

      // if the inclinination is small, the cell is assumed to be flat
      if (inclination_deg[i][j] < 5) aspect[i][j] = FL;

      // it is assumed, that the upper boundary of the MAP is oriented north
      else {
        if (i - x == -1) {
          if (j - y == -1) aspect[i][j] = SE;
          if (j - y == 0) aspect[i][j] = SS;
          if (j - y == 1) aspect[i][j] = SW;
        } else if (i - x == 0) {
          if (j - y == -1) aspect[i][j] = EE;
          if (j - y == 1) aspect[i][j] = WW;
        } else {
          if (j - y == -1) aspect[i][j] = NE;
          if (j - y == 0) aspect[i][j] = NN;
          if (j - y == 1) aspect[i][j] = NW;
        }
      }
      // if(outputCells[i][j]){
      // inclinationFile << inclination[i][j] << endl;
      // }
    }
  }
}

/*******************************************************************************************
 * write results in result file
 * \todo Btodo: shift to output
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterLandscape::writeOutputFile(int year, int day, int hour,
                                     Weather* weather) {
  int daynumber = year * daysPerYear + day;

  if (hour == 0) {
    dailyMeanWaterL0_abs = 0;
    dailyMeanWaterL1_rel = 0;
    dailyMeanWaterL2_rel = 0;
    dailyMeanFL1 = 0;
    dailyMeanFL2 = 0;
    dailyMeanFtot = 0;
    dailyMeanDrainL1L2 = 0;
    dailyMeanDrainDeep = 0;
    dailyMeanQD = 0;
    dailyQDlost = 0;
    dailyMeanLayerDiff = 0;
  }

  dailyMeanWaterL0_abs += meanWaterL0_abs;
  dailyMeanWaterL1_rel += meanWaterL1_rel;
  dailyMeanWaterL2_rel += meanWaterL2_rel;
  dailyMeanFL1 += meanFL1;
  dailyMeanFL2 += meanFL2;
  dailyMeanFtot += meanFtot;
  dailyMeanDrainL1L2 += meanDrainL1L2;
  dailyMeanDrainDeep += meanDrainDeep;
  dailyMeanQD += meanQD;
  dailyQDlost += QDlost;
  dailyMeanLayerDiff += meanLayerDiff;

  // daily results
  if (hour == 23) {
    dailyMeanWaterL0_abs = dailyMeanWaterL0_abs / 24.0;
    dailyMeanWaterL1_rel = dailyMeanWaterL1_rel / 24.0;
    dailyMeanWaterL2_rel = dailyMeanWaterL2_rel / 24.0;

    resultFileDay << year << "\t" << day << "\t" << weather->precSum[daynumber]
                  << "\t" << weather->tempAvg[daynumber] << "\t"
                  << dailyMeanWaterL0_abs << "\t" << dailyMeanWaterL1_rel * 100
                  << "\t" << dailyMeanWaterL2_rel * 100 << "\t";
    resultFileDay << dailyMeanFL1 << "\t" << dailyMeanFL2 << "\t"
                  << dailyMeanFtot << "\t" << dailyMeanDrainL1L2 << "\t"
                  << dailyMeanDrainDeep << "\t" << dailyMeanQD << "\t"
                  << dailyQDlost << "\t";

    for (int i = 0; i < Shrubs; i++) {
      resultFileDay << meanPftCover[i] << "\t" << meanPftTL1[i] << "\t"
                    << meanPftTL2[i] << "\t"; /*Tong*/
    }

    for (int i = Shrubs; i < Shrubs + Perennials; i++) {
      resultFileDay << meanPftCover[i] << "\t" << meanPftTL1[i] << "\t"
                    << meanPftTL2[i] << "\t"; /*Tong*/
    }

    for (int i = Shrubs + Perennials; i < PFTs; i++) {
      resultFileDay << meanPftCover[i] << "\t" << meanPftTL1[i]
                    << "\t"; /*Tong*/
    }

    resultFileDay << potEP << "\t" << meanEPL0 << "\t" << meanEPL1 << "\t"
                  << meanEPL2 << "\t" << meanEP << "\t" /*Tong*/
                  << meanTL1 << "\t" << meanTL2 << "\t" << meanEL1 << "\t"
                  << dailyMeanLayerDiff << "\n";
    //
  }
}

/*******************************************************************************************
 * write results in result file
 * \todo Btodo: shift to output
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterLandscape::writeOutputFileSelina(int year, int day, int hour,
                                           Weather* weather) {
  int daynumber = year * daysPerYear + day;

  // write hourly output
  resultFileSelinaHour << year << "\t" << day << "\t" << hour << "\t"
                       << weather->prec[day][hour] << "\t"
                       << weather->temp[day][hour] << "\t" << meanWaterL0_abs
                       << "\t" << meanWaterL1_rel << "\t" << meanWaterL2_rel
                       << "\t" << meanFL1 << "\t" << meanFL2 << "\t" << meanFtot
                       << "\t" << meanDrainL1L2 << "\t" << meanDrainDeep << "\t"
                       << meanQD << "\t" << QDlost << "\t" << potEP << "\t"
                       << meanEPL0 << "\t" << meanEPL1 << "\t" << meanEPL2
                       << "\t" << meanEP << "\t" << meanTL1 << "\t" << meanTL2
                       << "\t" << meanEL1 << "\t" << meanLayerDiff << "\n";

  resultFileSelinaHour_CS
      << year << "\t" << day << "\t" << hour << "\t" << weather->prec[day][hour]
      << "\t" << weather->temp[day][hour] << "\t" << meanWaterL0_abs_crust0
      << "\t" << meanWaterL0_abs_crust1 << "\t" << meanWaterL0_abs_crust2
      << "\t" << meanWaterL0_abs_crust3 << "\t" << meanWaterL0_abs_veg << "\t"
      << meanWaterL1_rel_crust0 << "\t" << meanWaterL1_rel_crust1 << "\t"
      << meanWaterL1_rel_crust2 << "\t" << meanWaterL1_rel_crust3 << "\t"
      << meanWaterL1_rel_veg << "\t" << meanWaterL2_rel_crust0 << "\t"
      << meanWaterL2_rel_crust1 << "\t" << meanWaterL2_rel_crust2 << "\t"
      << meanWaterL2_rel_crust3 << "\t" << meanWaterL2_rel_veg << "\t"
      << meanFL1_crust0 << "\t" << meanFL1_crust1 << "\t" << meanFL1_crust2
      << "\t" << meanFL1_crust3 << "\t" << meanFL1_veg << "\t"
      << meanDrainL1L2_crust0 << "\t" << meanDrainL1L2_crust1 << "\t"
      << meanDrainL1L2_crust2 << "\t" << meanDrainL1L2_crust3 << "\t"
      << meanDrainL1L2_veg << "\t" << meanDrainDeep_crust0 << "\t"
      << meanDrainDeep_crust1 << "\t" << meanDrainDeep_crust2 << "\t"
      << meanDrainDeep_crust3 << "\t" << meanDrainDeep_veg << "\t"
      << meanQD_crust0 << "\t" << meanQD_crust1 << "\t" << meanQD_crust2 << "\t"
      << meanQD_crust3 << "\t" << meanQD_veg << "\t" << meanEPL1_crust0 << "\t"
      << meanEPL1_crust1 << "\t" << meanEPL1_crust2 << "\t" << meanEPL1_crust3
      << "\t" << meanEPL1_veg << "\n";

  // write daily output

  calculateCellMeans();

  if (hour == 0) {
    dailyMeanWaterL0_abs = 0;
    dailyMeanWaterL1_rel = 0;
    dailyMeanWaterL2_rel = 0;
    dailyMeanFL1 = 0;
    dailyMeanFL2 = 0;
    dailyMeanFtot = 0;
    dailyMeanDrainL1L2 = 0;
    dailyMeanDrainDeep = 0;
    dailyMeanQD = 0;
    dailyQDlost = 0;
    dailyMeanLayerDiff = 0;
    dailySD_WaterL1_rel = dailySD_WaterL2_rel = dailySD_FL1 = dailySD_FL2 =
        dailySD_QD = 0.0;
  }

  dailyMeanWaterL0_abs += meanWaterL0_abs;
  dailyMeanWaterL1_rel += meanWaterL1_rel;
  dailyMeanWaterL2_rel += meanWaterL2_rel;
  dailyMeanFL1 += meanFL1;
  dailyMeanFL2 += meanFL2;
  dailyMeanFtot += meanFtot;
  dailyMeanDrainL1L2 += meanDrainL1L2;
  dailyMeanDrainDeep += meanDrainDeep;
  dailyMeanQD += meanQD;
  dailyQDlost += QDlost;
  dailyMeanLayerDiff += meanLayerDiff;
  dailySD_WaterL1_rel += dailySD_WaterL1_rel;
  dailySD_WaterL2_rel += dailySD_WaterL2_rel;
  dailySD_FL1 += dailySD_FL1;
  dailySD_FL2 += dailySD_FL2;
  dailySD_QD += dailySD_QD;

  // daily results
  if (hour == 23) {
    dailyMeanWaterL0_abs = dailyMeanWaterL0_abs / 24.0;
    dailyMeanWaterL1_rel = dailyMeanWaterL1_rel / 24.0;
    dailyMeanWaterL2_rel = dailyMeanWaterL2_rel / 24.0;

    dailySD_WaterL1_rel = dailySD_WaterL1_rel / 24.0;
    dailySD_WaterL2_rel = dailySD_WaterL2_rel / 24.0;

    resultFileSelina << year << "\t" << day << "\t"
                     << weather->precSum[daynumber] << "\t"
                     << weather->tempAvg[daynumber] << "\t"
                     << dailyMeanWaterL0_abs << "\t"
                     << dailyMeanWaterL1_rel * 100 << "\t"
                     << dailySD_WaterL1_rel << "\t"
                     << dailyMeanWaterL2_rel * 100 << "\t"
                     << dailySD_WaterL2_rel << "\t" << dailyMeanFL1 << "\t"
                     << dailySD_FL1 << "\t" << dailyMeanFL2 << "\t"
                     << dailySD_FL2 << "\t" << dailyMeanFtot << "\t"
                     << dailyMeanDrainL1L2 << "\t" << dailyMeanDrainDeep << "\t"
                     << dailyMeanQD << "\t" << dailySD_QD << "\t" << dailyQDlost
                     << "\t";

    for (int i = 0; i < Shrubs; i++) {
      resultFileSelina << meanPftCover[i] << "\t" << meanPftTL1[i] << "\t"
                       << meanPftTL2[i] << "\t"; /*Tong*/
    }

    for (int i = Shrubs; i < Shrubs + Perennials; i++) {
      resultFileSelina << meanPftCover[i] << "\t" << meanPftTL1[i] << "\t"
                       << meanPftTL2[i] << "\t"; /*Tong*/
    }

    for (int i = Shrubs + Perennials; i < PFTs; i++) {
      resultFileSelina << meanPftCover[i] << "\t" << meanPftTL1[i]
                       << "\t"; /*Tong*/
    }

    resultFileSelina << potEP << "\t" << meanEPL0 << "\t" << meanEPL1 << "\t"
                     << sdEPL1 << "\t" << meanEPL2 << "\t" << sdEPL2 << "\t"
                     << meanEP << "\t" << sdEP << "\t" /*Tong*/
                     << meanTL1 << "\t" << meanTL2 << "\t" << meanEL1 << "\t"
                     << dailyMeanLayerDiff << "\n";
    //
  }
}

/*******************************************************************************************
 * mainly to close files, but also to write some last results
 * is called in --> controller::runSimulation
 *******************************************************************************************/
void WaterLandscape::finishSimulation() { resultFileDay.close(); }

/*******************************************************************************************
 * read soil parameters from file
 * initialise soil grid: elevation, inclination, runoff-path...
 * is called in --> WaterLandscape::WaterLandscape
 *******************************************************************************************/
void WaterLandscape::initialiseSoilParameters(
    Parameter* p) { /*adjusted by Tong*/

  double suction_, Ks_, waterL1_sat_, WP_, FC_, depthL1_, depthL2_,
      iniWaterL0abs_, iniWaterL1_, iniWaterL2_, EPfactor_, resWater_, maxFL2_,
      maxAmountFL2_, diffConst_, TE_, Vegcoef1_, Vegcoef2_, Veg_ETcoef1_,
      Veg_ETcoef2_;

  // open soil parameters file
  ifstream parameterFile;
  string soilParametersName;

  soilParametersName = "Parameters" + delimiter + "soilparameters_" +
                       p->soilFileID[p->scenario] + ".txt";

  cout << "Read soil parameters from: " << soilParametersName << endl;

  // open file
  parameterFile.open(soilParametersName.c_str(), ios::binary | ios::in);

  // if something goes wrong...
  if (!parameterFile) {
    cerr << soilParametersName << " could not be opened!\n";
    exit(-1);
  }

  int fsize = 1000;
  parameterFile.ignore(fsize, '#');
  parameterFile.ignore(fsize, '*');
  // read in parameters and close file
  parameterFile.ignore(fsize, '#');
  parameterFile.ignore(fsize, ':');
  parameterFile >> suction_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> Ks_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> WP_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> FC_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> waterL1_sat_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> depthL1_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> depthL2_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> iniWaterL0abs_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> iniWaterL1_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> iniWaterL2_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> EPfactor_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> resWater_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> maxFL2_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> maxAmountFL2_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> diffConst_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> TE_;

  parameterFile.ignore(fsize, '#');
  parameterFile.ignore(fsize, ':');
  parameterFile >> Vegcoef1_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> Vegcoef2_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> Veg_ETcoef1_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> Veg_ETcoef2_;

  parameterFile.close();

  //---------------------------------------------
  /// SET ELEVATION.
  /// Initialise the elevation of each cell from a file.
  /// \implements \fn  WaterLandscape::setElevation(Parameter* p)
  //---------------------------------------------
  setElevation(p);

  //---------------------------------------------
  /// CALCULATE RUNON CELLS.
  /// Calculate the resulting number of cells,
  /// that provide runon to each cell.
  /// \implements \fn  WaterLandscape::calculateRunonCells()
  //---------------------------------------------
  calculateRunonCells();

  //---------------------------------------------
  /// SET ASPECT AND INCLINATION.
  /// Calculate the aspect and inclination of each cell.
  /// \implements \fn  WaterLandscape::setAspectAndInclination()
  //---------------------------------------------
  setAspectAndInclination();

  writeRunonCells();

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      //---------------------------------------------
      /// NEW SOIL.
      /// Every grid cell gets the same values...
      /// This is not needed at the moment but it offers the opportunity to
      /// account for heterogeneities later on.
      /// \implements \fn WaterCell::WaterCell(int x, int y, double suction_,
      /// double Ks_, double waterL1_sat_, double WP_, double FC_, double
      /// depthL1_, double depthL2_, double iniWaterL0abs_, double
      /// iniWaterL1_, double iniWaterL2_, double EPfactor_, double
      /// resWater_,double maxFL2_,double maxAmountFL2_,
      ///                            double diffCnst_, double TE_, int PFTs_
      ///                            int Perennials_, int Shrubs_, int
      ///                            Annuals_, double Vegcoef1_, double
      ///                            Vegcoef2_, double Veg_ETcoef1_, double
      ///                            Veg_ETcoef2_)
      //---------------------------------------------
      soilgrid[i][j] = new WaterCell(
          i, j, suction_, Ks_, waterL1_sat_, WP_, FC_, depthL1_, depthL2_,
          iniWaterL0abs_, iniWaterL1_, iniWaterL2_, EPfactor_, resWater_,
          maxFL2_, maxAmountFL2_, diffConst_, TE_, PFTs, Perennials, Shrubs,
          Annuals, Vegcoef1_, Vegcoef2_, Veg_ETcoef1_, Veg_ETcoef2_);
      soilgrid[i][j]->aspect = aspect[i][j];
      soilgrid[i][j]->inclination = inclination[i][j];

      //---------------------------------------------
      /// SET ASPECT FACTOR.
      /// \implements \fn WaterCell::setAspectFactor()
      //---------------------------------------------
      soilgrid[i][j]->setAspectFactor();
    }
  }
}

/*******************************************************************************************
 * initialise file for soil hydrology time series
 * can be called in --> WaterLandscape::WaterLandscape
 *******************************************************************************************/
void WaterLandscape::initialiseResultFile(Parameter* p) {
  // create and open result file
  string resultFileNameDay =
      "Results" + delimiter + p->site + "_" + std::to_string(p->simYears) +
      "years_daily_" + p->outputFileID[p->scenario] + "_climrep-" +
      std::to_string(p->climateRepetition + 1) + "_modelrep-" +
      std::to_string(p->modelRepetition + 1) + ".txt";

  resultFileDay.open(resultFileNameDay.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!resultFileDay) {
    cerr << resultFileNameDay << " could not be opened in main!\n";
    exit(-1);
  }

  resultFileDay << "year"
                << "\t"
                << "day"
                << "\t"
                << "weather->precSum[daynumber]"
                << "\t"
                << "weather->tempAvg[daynumber]"
                << "\t"
                << "dailyMeanWaterL0_abs"
                << "\t"
                << "dailyMeanWaterL1_rel*100"
                << "\t"
                << "dailyMeanWaterL2_rel*100"
                << "\t"
                << "dailyMeanFL1"
                << "\t"
                << "dailyMeanFL2"
                << "\t"
                << "dailyMeanFtot"
                << "\t"
                << "dailyMeanDrainL1L2"
                << "\t"
                << "dailyMeanDrainDeep"
                << "\t"
                << "dailyMeanQD"
                << "\t"
                << "dailyQDlost"
                << "\t";

  for (int i = 0; i < Shrubs; i++) {
    resultFileDay << "meanScover"
                  << "\t"
                  << "meanSTranspirationL1"
                  << "\t"
                  << "meanSTranspirationL2"
                  << "\t";
  }

  for (int i = 0; i < Perennials; i++) {
    resultFileDay << "meanGcover"
                  << "\t"
                  << "meanGTranspirationL1"
                  << "\t"
                  << "meanGTranspirationL2"
                  << "\t";
  }

  for (int i = 0; i < Annuals; i++) {
    resultFileDay << "meanAcover"
                  << "\t"
                  << "meanATranspirationL1"
                  << "\t";
  }

  resultFileDay << "potEP"
                << "\t"
                << "meanEPL0"
                << "\t"
                << "meanEPL1"
                << "\t"
                << "meanEPL2"
                << "\t"
                << "meanEP"
                << "\t"
                << "meanTranspirationL1"
                << "\t"
                << "meanTranspirationL2"
                << "\t"
                << "meanEL1"
                << "\t"
                << "dailyMeanLayerDiff"
                << "\n";
}

/*******************************************************************************************
 * initialise file for soil hydrology time series
 * can be called in --> WaterLandscape::WaterLandscape
 *******************************************************************************************/
void WaterLandscape::initialiseResultFileSelina(Parameter* p) {
  // create and open result file
  string resultFileNameDay = "Results" + delimiter + p->site + "_" +
                             std::to_string(p->simYears) + "years_daily_" +
                             p->outputFileID[p->scenario] + ".txt";

  resultFileSelina.open(resultFileNameDay.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!resultFileSelina) {
    cerr << resultFileNameDay << " could not be opened in main!\n";
    exit(-1);
  }

  resultFileSelina << "year"
                   << "\t"
                   << "day"
                   << "\t"
                   << "rain"
                   << "\t"
                   << "temperature"
                   << "\t"
                   << "dailyMeanWaterL0_abs"
                   << "\t"
                   << "moisture_L1"
                   << "\t"
                   << "SD_moisture_L1"
                   << "\t"
                   << "moisture_L2"
                   << "\t"
                   << "SD_moisture_L2"
                   << "\t"
                   << "infiltration_L1"
                   << "\t"
                   << "SD_infiltration_L1"
                   << "\t"
                   << "infiltration_L2"
                   << "\t"
                   << "SD_infiltration_L2"
                   << "\t"
                   << "infiltration_tot"
                   << "\t"
                   << "drainage_L1"
                   << "\t"
                   << "drainage_L2"
                   << "\t"
                   << "QD"
                   << "\t"
                   << "SD_QD"
                   << "\t"
                   << "QDlost"
                   << "\t";

  for (int i = 0; i < Shrubs; i++) {
    resultFileSelina << "SCover"
                     << "\t"
                     << "STranspiration_L1"
                     << "\t"
                     << "STranspiration_L2"
                     << "\t";
  }

  for (int i = 0; i < Perennials; i++) {
    resultFileSelina << "Gcover"
                     << "\t"
                     << "GTranspiration_L1"
                     << "\t"
                     << "GTranspiration_L2"
                     << "\t";
  }

  for (int i = 0; i < Annuals; i++) {
    resultFileSelina << "Acover"
                     << "\t"
                     << "ATranspiration_L1"
                     << "\t";
  }

  resultFileSelina << "potEP"
                   << "\t"
                   << "EP_L0"
                   << "\t"
                   << "EP_L1"
                   << "\t"
                   << "SD_EP_L1"
                   << "\t"
                   << "EP_L2"
                   << "\t"
                   << "SD_EP_L2"
                   << "\t"
                   << "EP_tot"
                   << "\t"
                   << "SD_EP_tot"
                   << "\t"
                   << "transpiration_L1"
                   << "\t"
                   << "transpiration_L2"
                   << "\t"
                   << "E_L1"
                   << "\t"
                   << "layerDiff"
                   << "\n";
}

/*******************************************************************************************
 * initialise file for soil hydrology time series
 * can be called in --> WaterLandscape::WaterLandscape
 *******************************************************************************************/
void WaterLandscape::initialiseResultFileSelinaHour(Parameter* p) {
  // create and open result file
  string resultFileName = "Results" + delimiter + p->site + "_" +
                          std::to_string(p->simYears) + "years_hourly_" +
                          p->outputFileID[p->scenario] + ".txt";

  resultFileSelinaHour.open(resultFileName.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!resultFileSelinaHour) {
    cerr << resultFileName << " could not be opened in main!\n";
    exit(-1);
  }

  resultFileSelinaHour << "year"
                       << "\t"
                       << "day"
                       << "\t"
                       << "hour"
                       << "\t"
                       << "rain"
                       << "\t"
                       << "temperature"
                       << "\t"
                       << "MeanWaterL0_abs"
                       << "\t"
                       << "moisture_L1"
                       << "\t"
                       << "moisture_L2"
                       << "\t"
                       << "infiltration_L1"
                       << "\t"
                       << "infiltration_L2"
                       << "\t"
                       << "infiltration_tot"
                       << "\t"
                       << "drainage_L1"
                       << "\t"
                       << "drainage_L2"
                       << "\t"
                       << "QD"
                       << "\t"
                       << "QDlost"
                       << "\t";

  // for (int i = 0; i < Shrubs; i++) {
  //   resultFileSelinaHour << "SCover"
  //                        << "\t"
  //                        << "STranspiration_L1"
  //                        << "\t"
  //                        << "STranspiration_L2"
  //                        << "\t";
  // }

  // for (int i = 0; i < Perennials; i++) {
  //   resultFileSelinaHour << "Gcover"
  //                        << "\t"
  //                        << "GTranspiration_L1"
  //                        << "\t"
  //                        << "GTranspiration_L2"
  //                        << "\t";
  // }

  // for (int i = 0; i < Annuals; i++) {
  //   resultFileSelinaHour << "Acover"
  //                        << "\t"
  //                        << "ATranspiration_L1"
  //                        << "\t";
  // }

  resultFileSelinaHour << "potEP"
                       << "\t"
                       << "EP_L0"
                       << "\t"
                       << "EP_L1"
                       << "\t"
                       << "EP_L2"
                       << "\t"
                       << "EP_tot"
                       << "\t"
                       << "transpiration_L1"
                       << "\t"
                       << "transpiration_L2"
                       << "\t"
                       << "E_L1"
                       << "\t"
                       << "layerDiff"
                       << "\n";
}

/*******************************************************************************************
 * initialise file for soil hydrology time series
 * can be called in --> WaterLandscape::WaterLandscape
 *******************************************************************************************/
void WaterLandscape::initialiseResultFileSelinaHour_CS(Parameter* p) {
  // create and open result file
  string resultFileName = "Results" + delimiter + p->site + "_" +
                          std::to_string(p->simYears) + "years_CS_" +
                          p->outputFileID[p->scenario] + ".txt";

  resultFileSelinaHour_CS.open(resultFileName.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!resultFileSelinaHour_CS) {
    cerr << resultFileName << " could not be opened in main!\n";
    exit(-1);
  }

  resultFileSelinaHour_CS << "year"
                          << "\t"
                          << "day"
                          << "\t"
                          << "hour"
                          << "\t"
                          << "rain"
                          << "\t"
                          << "temperature"
                          << "\t"
                          << "water_L0_crust0"
                          << "\t"
                          << "water_L0_crust1"
                          << "\t"
                          << "water_L0_crust2"
                          << "\t"
                          << "water_L0_crust3"
                          << "\t"
                          << "water_L0_veg"
                          << "\t"
                          << "moisture_L1_crust0"
                          << "\t"
                          << "moisture_L1_crust1"
                          << "\t"
                          << "moisture_L1_crust2"
                          << "\t"
                          << "moisture_L1_crust3"
                          << "\t"
                          << "moisture_L1_veg"
                          << "\t"
                          << "moisture_L2_crust0"
                          << "\t"
                          << "moisture_L2_crust1"
                          << "\t"
                          << "moisture_L2_crust2"
                          << "\t"
                          << "moisture_L2_crust3"
                          << "\t"
                          << "moisture_L2_veg"
                          << "\t"
                          << "infiltration_L1_crust0"
                          << "\t"
                          << "infiltration_L1_crust1"
                          << "\t"
                          << "infiltration_L1_crust2"
                          << "\t"
                          << "infiltration_L1_crust3"
                          << "\t"
                          << "infiltration_L1_veg"
                          << "\t"
                          << "drainage_L1_crust0"
                          << "\t"
                          << "drainage_L1_crust1"
                          << "\t"
                          << "drainage_L1_crust2"
                          << "\t"
                          << "drainage_L1_crust3"
                          << "\t"
                          << "drainage_L1_veg"
                          << "\t"
                          << "drainage_L2_crust0"
                          << "\t"
                          << "drainage_L2_crust1"
                          << "\t"
                          << "drainage_L2_crust2"
                          << "\t"
                          << "drainage_L2_crust3"
                          << "\t"
                          << "drainage_L2_veg"
                          << "\t"
                          << "QD_L0_crust0"
                          << "\t"
                          << "QD_L0_crust1"
                          << "\t"
                          << "QD_L0_crust2"
                          << "\t"
                          << "QD_L0_crust3"
                          << "\t"
                          << "QD_L0_veg"
                          << "\t"
                          << "EP_L1_crust0"
                          << "\t"
                          << "EP_L1_crust1"
                          << "\t"
                          << "EP_L1_crust2"
                          << "\t"
                          << "EP_L1_crust3"
                          << "\t"
                          << "EP_L1_veg"
                          << "\n";
}

/*******************************************************************************************
 * default constructor for water landscape
 * for reasons of simplification the constructor also initialises the
 *parameters is called in --> controller::initialiseSimulation
 *******************************************************************************************/
WaterLandscape::WaterLandscape(Parameter* p,
                               VegetationLandscape* vegLandscape) {
  list1 = new CellList;
  list2 = new CellList;

  surfaceWaterFlag = false;

  xsize = p->xsize;
  ysize = p->ysize;
  vegTimeStep = p->vegTimeStep;
  cellsize = p->cellsize;
  PFTs = vegLandscape->PFTs;
  Perennials = vegLandscape->Perennials;
  Annuals = vegLandscape->Annuals;
  Shrubs = vegLandscape->Shrubs;

  elevation.resize(xsize, Double1D(ysize, 0));
  aspect.resize(xsize, Int1D(ysize, 0));
  inclination.resize(xsize, Double1D(ysize, 0));
  inclination_deg.resize(xsize, Double1D(ysize, 0));
  runoncells.resize(xsize, Int1D(ysize, 0));
  runoffdirection.resize(xsize, Int1D(ysize, 0));
  calculated.resize(xsize, Bool1D(ysize, 0));
  numberOfFlows.resize(xsize, Int1D(ysize, 0));
  runon.resize(xsize, Double1D(ysize, 0));
  soilgrid.resize(xsize, vector<WaterCell*>(ysize));
  meanPftCover.resize(PFTs, 0);
  meanPftTL1.resize(PFTs, 0);
  meanPftTL2.resize(PFTs, 0);

  // Resize output vectors
  monthly_soil_moistureL1.resize(xsize, vector<double>(ysize, 0.0));
  monthly_soil_moistureL2.resize(xsize, vector<double>(ysize, 0.0));
  monthly_QD_cell.resize(xsize, vector<double>(ysize, 0.0));
  monthly_infiltrationL1_cell.resize(xsize, vector<double>(ysize, 0.0));
  monthly_runon_cell.resize(xsize, vector<double>(ysize, 0.0));
  monthly_deepdrain_cell.resize(xsize, vector<double>(ysize, 0.0));
  monthly_evapotranspiration_L1_cell.resize(xsize, vector<double>(ysize, 0.0));
  monthly_evapotranspiration_L2_cell.resize(xsize, vector<double>(ysize, 0.0));
  monthly_evapotranspiration_tot_cell.resize(xsize, vector<double>(ysize, 0.0));
  daily_soil_moistureL1.resize(xsize, vector<double>(ysize, 0.0));
  daily_soil_moistureL2.resize(xsize, vector<double>(ysize, 0.0));
  daily_QD_cell.resize(xsize, vector<double>(ysize, 0.0));
  daily_infiltrationL1_cell.resize(xsize, vector<double>(ysize, 0.0));
  daily_runon_cell.resize(xsize, vector<double>(ysize, 0.0));
  daily_deepdrain_cell.resize(xsize, vector<double>(ysize, 0.0));
  daily_evapotranspiration_L1_cell.resize(xsize, vector<double>(ysize, 0.0));
  daily_evapotranspiration_L2_cell.resize(xsize, vector<double>(ysize, 0.0));
  daily_evapotranspiration_tot_cell.resize(xsize, vector<double>(ysize, 0.0));

  meanFL1 = 0;
  meanFL2 = 0;
  meanFtot = 0;
  meanDrainL1L2 = 0;
  meanDrainDeep = 0;
  meanQD = 0;
  QDlost = 0;
  potEP = 0;
  meanEPL0 = 0;
  meanEPL1 = 0;
  meanEPL2 = 0;
  meanEP = 0;
  meanTL1 = 0; /*Tong*/
  meanTL2 = 0; /*Tong*/
  meanEL1 = 0; /*Tong*/

  for (int i = 0; i < PFTs; i++) { /*Tong*/
    meanPftTL1[i] = 0;
    meanPftTL2[i] = 0;
  }

  meanWaterL0_abs = 0;
  meanWaterL1_rel = 0;
  meanWaterL2_rel = 0;
  dailyMeanFL1 = 0;
  dailyMeanFL2 = 0;
  dailyMeanFtot = 0;
  dailyMeanDrainL1L2 = 0;
  dailyMeanDrainDeep = 0;
  dailyMeanQD = 0;
  dailyQDlost = 0;
  dailyMeanWaterL0_abs = 0;
  dailyMeanWaterL1_rel = 0;
  dailyMeanWaterL2_rel = 0;

  //---------------------------------------------
  /// INITIALISE RESULT FILE.
  /// \implements \fn  WaterLandscape::initialiseResultFile(Parameter* p)
  //---------------------------------------------
  initialiseResultFileSelinaHour(p);
  initialiseResultFileSelina(p);
  initialiseResultFileSelinaHour_CS(p);

  // set cells from output mask (if such a file exists)
  outputCells.resize(xsize, Int1D(ysize, 0));
  setOutputCells();

  //---------------------------------------------
  /// INITIALISE SOIL PARAMETERS.
  /// Read soil parameters from file.
  /// \implements \fn  WaterLandscape::initialiseSoilParameters(Parameter* p)
  //---------------------------------------------
  initialiseSoilParameters(p);

  // vegetation model gets depth of soil layers
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      vegLandscape->veggrid[i][j]->soilDepthL1 = soilgrid[i][j]->depthL1;
      vegLandscape->veggrid[i][j]->soilDepthL2 = soilgrid[i][j]->depthL2;
    }
  }
}

/*******************************************************************************************
 * destructor for water landscape
 *******************************************************************************************/
WaterLandscape::~WaterLandscape() {
  for (unsigned i = 0; i < soilgrid.size(); i++)
    for (unsigned j = 0; j < soilgrid.at(i).size(); j++) {
      delete soilgrid[i][j];
      soilgrid[i][j] = NULL;
    }

  delete list1;
  delete list2;
  /// \todo Btodo: what does this line want to say us?
  // delete elevation;
}

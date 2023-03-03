/*******************************************************************************************
 * VegetationLandscape.cpp
 * In this file, all functions are described that are relevant for vegetation
 * dynamics.
 * It contains the vegetation landscape vegLandscape with lots of grid cells
 * of vegetation
 * Functions that occur on the landscape scale are defined for the whole
 * landscape, while functions for the single grid cell are defined for the
 * variable vegetation
 *******************************************************************************************/

#include "VegetationLandscape.h"

#include <random>

#include "Annual.h"
#include "Perennial.h"
#include "RandomNumberGenerator.h"
#include "Shrub.h"


using namespace std;

/*******************************************************************************************
 * all concerning vegetation change are controlled here
 * \todo Functions calculateUptakeRateCover and calculateWaterUptake will be
 *called every day, so why calling various time in several- if-statements and
 *not at the beginning once? is called in --> controller::runSimulation
 *******************************************************************************************/
void VegetationLandscape::calculateProcesses(
    int year, int day, WaterLandscape* waterlandscape, Weather* we,
    Parameter* p) { /*Tong, based on PFT-oriented*/

  ///----------------------------FIRST DAY OF THE
  ///YEAR-------------------------------------------
  if (day == 0) {
    ///----------------------------FIRST DAY OF THE FIRST
    ///YEAR---------------------------------
    if (year == 0) {
      for (int x = 0; x < xsize; x++) {
        for (int y = 0; y < ysize; y++) {
          //--------------
          /// SET TO ZERO.
          //--------------
          /// \todo Btodo: all PFTs should have the same parameters
          veggrid[x][y]->reserveBM = 0;
          veggrid[x][y]->resCover = 0;
          for (int i = 0; i < PFTs; i++) {
            veggrid[x][y]->pftList[i]->Biomass = 0;
            veggrid[x][y]->edibleBM_reserve[i] = 0;
            if (i < Shrubs) {
              veggrid[x][y]->pftList[i]->est = 0;
            } else {
              veggrid[x][y]->pftList[i]->resBM = 0;
            }
          }
        }
      }
    }
    ///----------------------------END: FIRST DAY OF THE FIRST
    ///YEAR---------------------------------

    for (int m = 0; m < Shrubs; m++) {
      estCount[m] = 0;
    }

    firemem = fire2;  // Übertrage fire2 Wert auf fire Mem, so dass man weiss ob
                      // im Vorjahr feuer war!
    fire2 = false;

    /// \todo Btodo: do we need these?
    for (int l = 0; l < 10;
         l++) {  // count of cells above 10 levels is set to zero
      for (int i = 0; i < Shrubs + Perennials; i++) {
        Cabove[l][i] = 0;
      }
    }

    /// \todo Btodo: all PFTs should have the same parameters
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        for (int i = Shrubs + Perennials; i < PFTs; i++) {
          veggrid[x][y]->pftList[i]->Cover = 0;
        }
        for (int i = 0; i < Shrubs + Perennials; i++) {
          veggrid[x][y]->pftList[i]->totalDryMort = 0;
          veggrid[x][y]->pftList[i]->totalGrowth = 0;
        }
        for (int i = 0; i < Shrubs; i++) {
          estCount[i] += veggrid[x][y]->pftList[i]->est;
        }
      }
    }

    // fire is initiated at the end of the dry season, before growth period
    // starts, hence plants have the possibility to grow hereafter, but not so,
    // shrubs..? fire_ini
    /// \todo fire routine needs to be reviewed! Fire parameters (e.g. which
    /// year?) should be a parameter read in from file
    if (fire == true && year % 5 == 0 && year > 3) {
      fire2 = true;
      /// \implements \fn VegetationLandscape::Fire()
      calculateFire();
    }
    //-------------------------
    /// CALCULATE WATER UPTAKE.
    //-------------------------
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        for (int i = 0; i < PFTs; i++) {
          /// \implements \fn
          /// PlantFunctionalType::calculateUptakeRateCover(VegetationCell*
          /// vegetation, Weather* weather, int year)
          veggrid[x][y]->pftList[i]->calculateUptakeRateCover(veggrid[x][y], we,
                                                              year);
        }
        // \implements \fn VegetationCell::calculateWaterUptake()
        veggrid[x][y]->calculateWaterUptake();
      }
    }
  }  // day==0 end
  ///----------------------------END: FIRST DAY OF THE
  ///YEAR-------------------------------------------

  ///----------------------------EVERY VEGTIMESTEP WITHIN WET
  ///SEASON-------------------------------------------
  if (day >= wetStart && day <= wetEnd &&
      day % vegTimeStep == vegTimeStep - 1) {
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        //-----------------------------------------------------------------
        /// CALCULATE AVERAGE MOISTURE IN BOTH LAYERS AND WATER UPTAKE RATE.
        //-----------------------------------------------------------------
        for (int i = 0; i < PFTs; i++) {
          /// \implements \fn
          /// PlantFunctionalType::calculateAvMoistL1(VegetationCell*
          /// vegetation)
          veggrid[x][y]->pftList[i]->calculateAvMoistL1(veggrid[x][y]);
          /// \implements \fn
          /// PlantFunctionalType::calculateAvMoistL2(VegetationCell*
          /// vegetation)
          veggrid[x][y]->pftList[i]->calculateAvMoistL2(veggrid[x][y]);
        }

        //---------------
        /// UPDATE COVER.
        //---------------
        /// \implements \fn VegetationCell::calculatePerennialCover()
        // veggrid[x][y]->gCover = veggrid[x][y]->calculatePerennialCover();
        // /// \implements \fn VegetationCell::calculateShrubCover()
        // veggrid[x][y]->sCover = veggrid[x][y]->calculateShrubCover();
        // /// \implements \fn VegetationCell::calculateAnnualCover()
        // veggrid[x][y]->aCover = veggrid[x][y]->calculateAnnualCover();
        // /// \implements \fn VegetationCell::calculateTotalCover()
        // veggrid[x][y]->totalPFTcover = veggrid[x][y]->calculateTotalCover();
        //--------------------------------------
        /// CALCULATE GROWTH, COVER AND BIOMASS.
        //--------------------------------------
        // for (int i = 0; i < PFTs; i++){
        //     /// \implements \fn
        //     PlantFunctionalType::calculateGrowthL1(VegetationCell*
        //     vegetation, int index)
        //     veggrid[x][y]->pftList[i]->calculateGrowthL1(veggrid[x][y], i);
        //     /// \implements \fn
        //     PlantFunctionalType::calculateGrowthL2(VegetationCell*
        //     vegetation, int index)
        //     veggrid[x][y]->pftList[i]->calculateGrowthL2(veggrid[x][y], i);
        //     /// \implements \fn PlantFunctionalType::calculateGrowth()
        //     veggrid[x][y]->pftList[i]->calculateGrowth();
        //     /// \implements \fn PlantFunctionalType::calculateTotalGrowth()
        //     veggrid[x][y]->pftList[i]->calculateTotalGrowth();
        //     /// \implements \fn PlantFunctionalType::addGrowthToCover()
        //     veggrid[x][y]->pftList[i]->addGrowthToCover();
        //     /// \implements \fn
        //     PlantFunctionalType::calculateBiomass(VegetationCell* vegetation,
        //     Weather* we2, int yr2)
        //     veggrid[x][y]->pftList[i]->calculateBiomass(veggrid[x][y], we,
        //     year);
        // }

        //-------------------------
        /// CALCULATE WATER UPTAKE.
        //-------------------------
        for (int i = 0; i < PFTs; i++) {
          /// \implements \fn
          /// PlantFunctionalType::calculateUptakeRateCover(VegetationCell*
          /// vegetation, Weather* weather, int year)
          veggrid[x][y]->pftList[i]->calculateUptakeRateCover(veggrid[x][y], we,
                                                              year);
        }
        /// \implements \fn VegetationCell::calculateWaterUptake()
        veggrid[x][y]->calculateWaterUptake();
      }
      ///----------------------------END: EVERY VEGTIMESTEP WITHIN WET
      ///SEASON-------------------------------------------
    }
  }

  ///----------------------------END: WITHIN WET
  ///SEASON-------------------------------------------

  ///----------------------------AT THE END OF THE GROWING
  ///SEASON-------------------------------------------
  // if the soil was relatively dry during the growth season: mortality due to
  // water shortage (vegMortality())
  //  a basic stochastic mortality  due to several possible but not explicitly
  //  simulated reasons (frost, disease, browse)
  if (day == growEnd) {
    // //--------------------------------------------------------------------------------------------
    // /// CALCULATE WATER UPTAKE RATE, VEGMORTALITY, BASIC MORTALITY, BASIC
    // MORTALITY RATE, BIOMASS.
    // //--------------------------------------------------------------------------------------------
    // for (int x = 0; x < xsize; x++) {
    //     for (int y = 0; y < ysize; y++) {
    //         for (int i = 0; i < PFTs; i++){
    //             if (i < Shrubs) {
    //                 if (veggrid[x][y]->pftList[i]->Cover > 0) {
    //                     veggrid[x][y]->cohortAge[i]++;  // increase age of
    //                     shrub cohort by 1 veggrid[x][y]->cohortAgeFire[i]++;
    //                     // increase age of shrub cohort by 1
    //                 } else {
    //                     veggrid[x][y]->cohortAge[i] = 0;
    //                     veggrid[x][y]->cohortAgeFire[i] = 0;
    //                 }
    //             }
    // }

    //         for (int i = 0; i < PFTs; i++){
    //             /// \implements \fn
    //             PlantFunctionalType::vegMortality(VegetationCell* vegetation,
    //             int day, int growStart, int growEnd, int year, Weather*
    //             weather)
    //             veggrid[x][y]->pftList[i]->vegMortality(veggrid[x][y], day,
    //             growStart, growEnd, year, we);
    //         }

    //         veggrid[x][y]->gBasicMortRate_min =
    //         veggrid[x][y]->pftList[Shrubs]->BasicMortRate; /// \todo why only
    //         for first Perennial? for (int i = 0; i < PFTs; i++){
    //             /// \implements \fn
    //             PlantFunctionalType::basicMortality(VegetationCell*
    //             vegetation, int index)
    //             veggrid[x][y]->pftList[i]->basicMortality(veggrid[x][y], i);

    //             if (i >= Shrubs && i < Perennials) {
    //                 if (veggrid[x][y]->pftList[i]->BasicMortRate <
    //                 veggrid[x][y]->gBasicMortRate_min ){
    //                     veggrid[x][y]->gBasicMortRate_min =
    //                     veggrid[x][y]->pftList[i]->BasicMortRate;  /// \todo
    //                     why here only for last Perennial?
    //                 }
    //             }
    //             /// \implements \fn
    //             PlantFunctionalType::calculateBiomass(VegetationCell*
    //             vegetation, Weather *we2, int yr2)
    //             veggrid[x][y]->pftList[i]->calculateBiomass(veggrid[x][y],
    //             we, year);
    //         }
    //         //-------------------------
    //         /// CALCULATE WATER UPTAKE.
    //         //-------------------------
    //         for (int i = 0; i < PFTs; i++){
    //             /// \implements \fn
    //             PlantFunctionalType::calculateUptakeRateCover(VegetationCell*
    //             vegetation, Weather* weather, int year)
    //             veggrid[x][y]->pftList[i]->calculateUptakeRateCover(veggrid[x][y],
    //             we, year);
    //         }
    //         /// \implements \fn VegetationCell::calculateWaterUptake()
    //    	    veggrid[x][y]->calculateWaterUptake();
    //     }
    // }
    /// \implements \fn VegetationLandscape::meanCoverBiomass()
    // meanCoverBiomass();
  }
  ///----------------------------END: AT THE END OF THE GROWING
  ///SEASON-------------------------------------------

  ///----------------------------LAST DAY OF THE
  ///YEAR------------------------------------------------------------
  if (day == 364) {
    //-------------------------------------
    /// CALCULATE GRAZING.
    //-------------------------------------
    /// \implements \fn void VegetationLandscape::Grazing(Weather* we, int yr2)
    // grazing(we, year);

    //-------------------------------------
    /// CALCULATE DISPERSAL.
    /// Dead plants won't disperse anymore.
    //-------------------------------------
    // oldstock = realstock; /// \todo why new variable and not just hand over
    // realstock?
    // /// \implements \fn VegetationLandscape::dispersal(int Stockrate)
    // dispersal(oldstock);

    // //--------------------
    // /// CALCULATE BIOMASS.
    // //--------------------
    // for (int x = 0; x < xsize; x++) {
    //     for (int y = 0; y < ysize; y++){
    //         for (int i = 0; i < PFTs; i++){
    //             /// \implements \fn
    //             PlantFunctionalType::calculateBiomass(VegetationCell*
    //             vegetation, Weather* we2, int yr2)
    //             veggrid[x][y]->pftList[i]->calculateBiomass(veggrid[x][y],
    //             we, year);
    //         }
    //     }
    // }

    //------------------------------
    /// ELIMINATE RUDIMENTARY COVER.
    /// \todo why only cover of shrubs and annuals?
    //------------------------------
    // for (int x = 0 ; x < xsize; x++) {
    //     for (int y = 0; y < ysize; y++){
    //         for (int i = 0; i < Shrubs; i++){
    //             if (veggrid[x][y]->pftList[i]->Cover < sclim) {
    //                 veggrid[x][y]->pftList[i]->Cover = 0;
    //             }
    //         }
    //         for (int i = Shrubs + Perennials; i < PFTs; i++){
    //             veggrid[x][y]->pftList[i]->Cover = 0;
    //         }
    //     }
    // }

    //------------------------------------------
    /// SET GROWTH AND AVERAGE MOISTURE TO ZERO.
    //------------------------------------------
    for (int x = 0; x < xsize; x++) {
      for (int y = 0; y < ysize; y++) {
        for (int i = 0; i < PFTs; i++) {
          veggrid[x][y]->pftList[i]->AvMoistL1 = 0;
          veggrid[x][y]->pftList[i]->DryMort = 0;
          veggrid[x][y]->pftList[i]->DryMortL1 = 0;
          veggrid[x][y]->pftList[i]->DryMortL2 = 0;
          veggrid[x][y]->pftList[i]->AvMoistSuppL1 = 0;
          veggrid[x][y]->pftList[i]->AvMoistSuppL2 = 0;

          if (i < Shrubs + Perennials) {
            veggrid[x][y]->pftList[i]->Growth = 0;
            veggrid[x][y]->pftList[i]->GrowthL1 = 0;
            veggrid[x][y]->pftList[i]->GrowthL2 = 0;
            veggrid[x][y]->pftList[i]->AvMoistL2 = 0;
          }
        }
      }
    }

    //----------------------------------------------------------------------------------------------------------
    /// WRITE YEARLY SUMMARY FILE.
    /// \implements \fn VegetationLandscape::writeYearlySummary(WaterLandscape
    /// *waterlandscape, Weather *we, int year)
    //----------------------------------------------------------------------------------------------------------
    //writeYearlySummary(waterlandscape, we, year);

    if (year == p->simYears - 1) {
      yearlySummaryFile.close();
    }
    fire = false;  /// \todo why?

  }  // Day 364 end
     ///----------------------------END: LAST DAY OF THE
     ///YEAR------------------------------------------------------------

  //-------------------------------------------------------
  /// CALCULATE MEAN.
  /// \implements \fn VegetationLandscape::calculateMean(int day)
  //-------------------------------------------------------
  calculateMean(day);

  //-------------------------------------------------------
  /// CALCULATE SHRUB.
  /// \implements \fn VegetationLandscape::calculateShrub()
  //-------------------------------------------------------
  calculateShrub();
}

/*******************************************************************************************
 * calculate total cover in each grid cell
 * is called in -->
 *******************************************************************************************/
void VegetationLandscape::calculateTotalCover() {
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      //-------------------------------------------------------
      /// CALCULATE TOTAL COVER.
      /// \implements \fn VegetationCell::calculateTotalCover()
      //-------------------------------------------------------
      veggrid[x][y]->calculateTotalCover();
    }
  }
}

/*******************************************************************************************
 * calculate "diffusive" dispersal for grasses and shrubs: each plant type has
 * a defined dispersal function (decreasing with distance). Within this distance
 * a small fraction of the plant (e.g. by seeds, clonal) is dispersed
 * \todo needs to be revised!
 * is called in --> VegetationLandscape::dispersal
 *******************************************************************************************/
void VegetationLandscape::calculateShrubDispersal(
    int x1, int y1, int x2, int y2) { /*Tong, based on PFT-oriented*/

  for (int i = 0; i < Shrubs; i++) {
    dispHelpS[i] = 0;
    zuf_dispS[i] = 0;
  }

  // dirk faktor von Stocking Rate erhoeht ausbreitungsdistanz und flacht Kurve
  // ab
  for (int i = 0; i < Shrubs; i++) {
    newdistconst[i] = veggrid[x1][y1]->pftList[i]->distConstS +
                      veggrid[x1][y1]->pftList[i]->distConstS2 * stockingRate;
  }
  // euklidian distance
  double help2 = pow(x1 - x2, 2) + pow(y1 - y2, 2);
  if (help2 < 0) {
    cout << "SQRT Problem in CalculateShrubDispersal" << help2 << "\n";
  }
  double dist = sqrt(help2);
  /// \implements \fn VegetationCell::calculateShrubCover()
  double scover = veggrid[x2][y2]->calculateShrubCover();
  /// \implements \fn VegetationCell::calculatePerennialCover()
  double pcover = veggrid[x2][y2]->calculatePerennialCover();

  double zuf = RandomNumberGenerator::generate_random_float(
      0.0, 1.0);  // random effect of dispersal (maybe zero, maybe double of
                  // calculated value)

  for (int i = 0; i < Shrubs; i++) {
    dispHelpS[i] =
        veggrid[x1][y1]->pftList[i]->disFrac *
        veggrid[x1][y1]->pftList[i]->dist0S * exp(-newdistconst[i] * dist) *
        veggrid[x1][y1]->pftList[i]->Cover * veggrid[x2][y2]->pftList[i]->est *
        max((1 - scover - pcover), 0.0);
    zuf_dispS[i] = zuf * 2.0 * dispHelpS[i];
  }

  for (int i = 0; i < Shrubs; i++) {
    veggrid[x2][y2]->pftList[i]->Disp += zuf_dispS[i];
  }

  for (int i = 0; i < Shrubs; i++) {
    veggrid[x2][y2]->pftList[i]->DispSend += zuf_dispS[i];
  }
}

/*******************************************************************************************
 * calculate "diffusive" dispersal:
 * each plant type (shrubs / grasses) has
 * a defined dispersal function (decreasing with distance). Within this distance
 * a small fraction of the plant (e.g. by seeds, clonal) is dispersed
 * \todo needs to be revised!
 * is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationLandscape::dispersal(
    int Stockrate) { /*Tong, based on PFT-oriented*/

  double scover;
  double gcover;
  double scoverbeforedis;
  double gcoverbeforedis;
  double sadded;
  double gadded;
  double coverleft;
  double coveradded;

  for (int i = Shrubs; i < Shrubs + Perennials; i++) {
    dispHelp[i] = 0;
  }

  // first set the dispersal to zero
  // and afterwards add the dispersal from each cell
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int i = Shrubs; i < Shrubs + Perennials; i++) {
        veggrid[x][y]->pftList[i]->Disp = 0;
        veggrid[x][y]->pftList[i]->DispSend = 0;
        dispHelp[i] += veggrid[x][y]->pftList[i]->disFrac *
                       veggrid[x][y]->pftList[i]->Cover;
      }
    }
  }

  for (int i = Shrubs; i < Shrubs + Perennials; i++) {
    dispHelp[i] /= (1.0 * xsize * ysize);
  }

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int i = 0; i < Shrubs; i++) {
        veggrid[x][y]->pftList[i]->Disp = 0;
        veggrid[x][y]->pftList[i]->DispSend = 0;
      }
    }
  }

  /*----------------------- dispersal of grasses ------------------------*/
  // assumption: grass is not limited in dispersal  //the amount per cell
  // depends on the mean cover of the total grid
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      /// \implements \fn VegetationCell::calculateShrubCover()
      scover = veggrid[x][y]->calculateShrubCover();
      /// \implements \fn VegetationCell::calculatePerennialCover()
      gcover = veggrid[x][y]->calculatePerennialCover();

      // sending part of dispersal seeds of perennial grass
      for (int i = Shrubs; i < Shrubs + Perennials; i++) {
        veggrid[x][y]->pftList[i]->DispSend =
            veggrid[x][y]->pftList[i]->Cover *
            veggrid[x][y]->pftList[i]->disFrac;
      }

      double zuf = RandomNumberGenerator::generate_random_float(0.0, 1.0);

      // receiving part of dispersal seeds of perennial grass
      for (int i = Shrubs; i < Shrubs + Perennials; i++) {
        veggrid[x][y]->pftList[i]->Disp = dispHelp[i] *
                                          veggrid[x][y]->pftList[i]->est *
                                          max(1 - gcover - scover, 0.0);
        if (zuf < 0.01) {
          // there has a problem, for multi-perennials and one perennial, it has
          // a difference, one perennial has one possibility of increasing 0.02,
          //  multi-perennials has several possibilities of increasing 0.02
          veggrid[x][y]->pftList[i]->Disp +=
              0.02 / Perennials;  /// \todo Btodo: why is this random process
                                  /// added this way, shouldn't it be neutral,
                                  /// i.e. between -0.01 and +0.01?
        }
      }
    }
  }

  /*----------------------- dispersal of shrubs ------------------------*/
  // neue "distablishment rule vom 15.04.2010 mit neuer SR Anpassung
  // vom 5.Mai2010
  // dirk faktor von Stocking Rate erhoeht ausbreitungsdistanz und flacht Kurve
  // ab

  // if (Stockrate > 30)  newdistconst = distConstS;
  // evaluate the dispersal of each shrub cell and where it disperses
  for (int x1 = 0; x1 < xsize; x1++) {
    for (int y1 = 0; y1 < ysize; y1++) {
      for (int i = 0; i < Shrubs; i++) {
        newdistconst[i] = veggrid[x1][y1]->pftList[i]->distConstS +
                          veggrid[x1][y1]->pftList[i]->distConstS2 * Stockrate;
        double ExitParam = veggrid[x1][y1]->pftList[i]->dist0S *
                           maxDistS;  // maxdistS% von VMax
        int helpDist = 0;
        double helpMaxDist;

        do {
          helpMaxDist = veggrid[x1][y1]->pftList[i]->dist0S *
                        exp(-newdistconst[i] * helpDist);
          helpDist++;
        } while (helpMaxDist > ExitParam);

        // square around x1 and y1 (within distance), with boundary  only for
        // cells with shrubs that are "adult"
        if (veggrid[x1][y1]->pftList[i]->Cover >
            scrlim) {  //_dirk if cover below scrlim, we assume shrubs to be
                       //juvenile and non-reproducing
          for (int x2 = max(x1 - int((helpDist) / (1.0 * cellsize)), 0);
               x2 < min(x1 + int((helpDist) / (1.0 * cellsize)), xsize - 1);
               x2++) {
            for (int y2 = max(y1 - int((helpDist) / (1.0 * cellsize)), 0);
                 y2 < min(y1 + int((helpDist) / (1.0 * cellsize)), ysize - 1);
                 y2++) {
              //---------------------------------------------
              /// CALCULATE SHRUB DISPERSAL.
              /// \implements \fn
              /// VegetationLandscape::calculateShrubDispersal(int x1, int y1,
              /// int x2, int y2)
              //---------------------------------------------
              calculateShrubDispersal(x1, y1, x2, y2);
            }
          }
        }
      }
    }  // y1 end
  }    // x1 end

  /*----------------------- add dispersed amount to cover
   * ------------------------*/
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      /// \implements \fn VegetationCell::calculateShrubCover()
      scoverbeforedis = veggrid[x][y]->calculateShrubCover();
      /// \implements \fn VegetationCell::calculatePerennialCover()
      gcoverbeforedis = veggrid[x][y]->calculatePerennialCover();

      sadded = 0;
      gadded = 0;

      for (int i = 0; i < Shrubs; i++) {
        sadded += veggrid[x][y]->pftList[i]->Disp;
      }

      for (int i = Shrubs; i < Shrubs + Perennials; i++) {
        gadded += veggrid[x][y]->pftList[i]->Disp;
      }

      // limit the total cover of PFTs
      // add the present cover of all PFTs to see how much space is left
      coverleft = 1 - scoverbeforedis - gcoverbeforedis;
      if (coverleft < 0) {
        coverleft = 0;
      }

      //	add the potential additional cover of all PFTs due to dispersal
      coveradded = sadded + gadded;
      if (coveradded > 1) {
        coveradded = 1;
      }

      //  limit dispersal of each PFT if potential dispersal is bigger than left
      //  space
      if (coveradded > coverleft) {
        for (int i = 0; i < Shrubs + Perennials; i++) {
          veggrid[x][y]->pftList[i]->Cover +=
              (veggrid[x][y]->pftList[i]->Disp / coveradded) * coverleft;
        }
      } else {  // operate the original algorithm of dispersal for each PFT
        for (int i = Shrubs; i < Shrubs + Perennials; i++) {
          veggrid[x][y]->pftList[i]->Cover += veggrid[x][y]->pftList[i]->Disp;
        }

        double zuf = RandomNumberGenerator::generate_random_float(0.0, 1.0);

        for (int i = 0; i < Shrubs; i++) {
          veggrid[x][y]->pftList[i]->Cover += veggrid[x][y]->pftList[i]->Disp;

          if (zuf < 0.001) {
            // there has a problem, for multi-shrubs and one shrub, it has a
            // difference, one shrub has one possibility of increasing 0.02,
            //  multi-shrubs has several possibilities of increasing 0.02
            veggrid[x][y]->pftList[i]->Cover += 2.0 * sclim / Shrubs;
          }
        }
      }
    }
  }
}

/*******************************************************************************************
 * calculate mean cover in the grid and biomass per ha
 * \todo Btodo: do we need this function?
 * is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationLandscape::meanCoverBiomass() { /*Tong, based on PFT-oriented*/

  for (int i = 0; i < PFTs; i++) {
    meanC[i] = 0;
    meanB[i] = 0;
    netMortR[i] = 0;
  }

  meanRB = 0;

  // to define cover thresholds to compare with
  double thresh[10];  //_dirk
  thresh[0] = 0.001;
  thresh[1] = 0.005;
  thresh[2] = 0.01;
  thresh[3] = 0.05;
  thresh[4] = 0.25;
  thresh[5] = 0.3;
  thresh[6] = 0.4;
  thresh[7] = 0.5;
  thresh[8] = 0.6;
  thresh[9] = 0.1;

  for (int i = 0; i < 10; i++) {
    for (int m = 0; m < PFTs; m++) {
      Cabove[i][m] = 0;
    }
    /// \todo what is sCsumabove? sCsumabove was not correctly initilised, only
    /// Shrubs but here until 10
    sCsumabove[i] = 0;
  }

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int i = 0; i < PFTs; i++) {
        meanC[i] += veggrid[x][y]->pftList[i]->Cover;
        netMortR[i] += (veggrid[x][y]->pftList[i]->totalDryMort +
                        veggrid[x][y]->pftList[i]->basicM);  //  _dirk
      }

      for (int i = 0; i < 10; i++) {
        veggrid[x][y]->sCover = 0;
        for (int m = 0; m < PFTs; m++) {
          if (veggrid[x][y]->pftList[m]->Cover > thresh[i]) {
            Cabove[i][m]++;
          }
          if (m < Shrubs) {
            veggrid[x][y]->sCover += veggrid[x][y]->pftList[m]->Cover;
          }
        }
        if (veggrid[x][y]->sCover > thresh[i]) {
          sCsumabove[i]++;
        }
      }
    }
  }

  for (int i = 0; i < PFTs; i++) {
    meanC[i] /= (xsize * ysize * 1.0);
    netMortR[i] /= (xsize * ysize * 1.0);  // mean mortality per cell
  }

  // mean biomass per ha!
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int i = 0; i < PFTs; i++) {
        meanB[i] += veggrid[x][y]->pftList[i]->Biomass;
      }
      meanRB += veggrid[x][y]->reserveBM;
    }
  }
  for (int i = 0; i < PFTs; i++) {
    meanB[i] =
        meanB[i] / ((xsize * ysize * cellsize * cellsize) / (100.0 * 100));
  }
  meanRB = meanRB / ((xsize * ysize * cellsize * cellsize) / (100.0 * 100));
}

/*******************************************************************************************
 * organizes fire, calls calculateFire and stores fire results
 * only called, if a fire actually occurs
 * is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationLandscape::calculateFire() {
  fires++;
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      /// \implements \fn VegetationCell::calculatePFTfire()
      veggrid[x][y]->calculatePFTfire();
    }
  }
}

/*******************************************************************************************
 * \todo include brief description here
 * is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationLandscape::grazing(
    Weather* we, int yr2) {  // grazing algorithm has been rewritten by Tong

  // one purpose: remove the facilitation relation among PFTs (such as
  // competition and coexistence relation in one function) calculate the
  // relative preference for each plant type, this parameter will be used to
  // adjust the grazed amount
  // if we have a type that is preferred above average, relprefer[PFT] is more
  // than 1; if it is below average, relprefer[PFT] is less than 1.
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      // calculate the mean preference
      meanprefer = 0;
      for (int m = 0; m < PFTs; m++) {
        meanprefer += veggrid[i][j]->pftList[m]->grazeprefer;
      }

      meanprefer /= (1.0 * PFTs);

      // calculate the relative preference
      for (int m = 0; m < PFTs; m++) {
        veggrid[i][j]->relprefer[m] =
            veggrid[i][j]->pftList[m]->grazeprefer / meanprefer;
      }
    }
  }

  // calculate the needed biomass dependent on livestock density
  realstock = stockingRate;  // realized stocking rate might change if not
                             // enough fodder is available
  double BMneed = 0;
  BMneed = 365 * percapability * BiomassNeed *
           (double)(1.0 /
                    realstock);  // Total biomass need per year and ha (what one
                                 // LSU needs per per year is 365*0.020*450000 )

  // calculate the average needed biomass per grid cell
  double avBMneed = 0.0;
  double areagrazed = (double)(xsize * cellsize * ysize * cellsize) /
                      (100.0 * 100);            // area grazed in ha
  BMneed = areagrazed * BMneed;                 // total absolute biomass demand
  avBMneed = (double)BMneed / (xsize * ysize);  // average need per cell

  // calculate the available biomass in each grid cell
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      for (int m = 0; m < PFTs; m++) {
        // calculate alive biomass
        /// \implements \fn
        /// PlantFunctionalType::calculateBiomass(VegetationCell* vegetation,
        /// Weather *we2, int yr2)
        veggrid[i][j]->pftList[m]->calculateBiomass(veggrid[i][j], we, yr2);
        // store the value as old biomass
        veggrid[i][j]->pftList[m]->Biomassold =
            veggrid[i][j]->pftList[m]->Biomass;
        // calculate the safe part of alive biomass
        veggrid[i][j]->safeBM_alive[m] = veggrid[i][j]->pftList[m]->grazeLim *
                                         veggrid[i][j]->pftList[m]->Biomass;
        // calculate the edible part of alive biomass
        veggrid[i][j]->edibleBM_alive[m] =
            veggrid[i][j]->pftList[m]->Biomass - veggrid[i][j]->safeBM_alive[m];
      }
    }
  }

  // total edible biomass in a cell (annual grass, perennial grass and shrub)
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      veggrid[i][j]->AnnualresBM = 0;
      veggrid[i][j]->PerennialresBM = 0;
      veggrid[i][j]->AnnualaliveBM = 0;
      veggrid[i][j]->PerennialaliveBM = 0;
      veggrid[i][j]->ShrubaliveBM = 0;

      for (int m = 0; m < PFTs; m++) {
        // Grazing is different for shrub, shrub has no reserve biomass.
        if (m < Shrubs) {
          // calculate the edible alive biomass for all shrubs in one cell
          veggrid[i][j]->ShrubaliveBM += veggrid[i][j]->edibleBM_alive[m];
          // calculate the total edible biomass in a cell
          veggrid[i][j]->edibleBM_total[m] = veggrid[i][j]->edibleBM_alive[m];
        } else if (m >= Shrubs && m < Shrubs + Perennials) {
          // sum of edible reserve biomass for all perennial grass in a cell
          veggrid[i][j]->PerennialresBM += veggrid[i][j]->edibleBM_reserve[m];
          // sum of edible alive biomass for all perennial grass in a cell
          veggrid[i][j]->PerennialaliveBM += veggrid[i][j]->edibleBM_alive[m];
          // calculate the total edible biomass in a cell for for one annual
          // grass
          veggrid[i][j]->edibleBM_total[m] =
              veggrid[i][j]->edibleBM_reserve[m] +
              veggrid[i][j]->edibleBM_alive[m];
        } else {
          // sum of edible reserve biomass for all annual grass in a cell
          veggrid[i][j]->AnnualresBM += veggrid[i][j]->edibleBM_reserve[m];
          // sum of edible alive biomass for all annual grass in a cell
          veggrid[i][j]->AnnualaliveBM += veggrid[i][j]->edibleBM_alive[m];
          // calculate the total edible biomass in a cell for for one annual
          // grass
          veggrid[i][j]->edibleBM_total[m] =
              veggrid[i][j]->edibleBM_reserve[m] +
              veggrid[i][j]->edibleBM_alive[m];
        }
      }

      // calculate the edible alive biomass for all shrubs in a cell, there is
      // no reserve biomass for shrub
      veggrid[i][j]->edibleShrub = veggrid[i][j]->ShrubaliveBM;
      // calculate the total edible biomass in a cell for all annual grasses
      veggrid[i][j]->edibleAnnual =
          veggrid[i][j]->AnnualresBM + veggrid[i][j]->AnnualaliveBM;
      // calculate the total edible biomass in a cell for all perennial grasses
      veggrid[i][j]->ediblePerennial =
          veggrid[i][j]->PerennialresBM + veggrid[i][j]->PerennialaliveBM;

      // calculate the edible reserve biomass of annual grass and perennial
      // grass in a cell
      veggrid[i][j]->reserveBM =
          veggrid[i][j]->AnnualresBM + veggrid[i][j]->PerennialresBM;
      // calculate the edible alive biomass of annual grass, perennial grass and
      // shrub in a cell
      veggrid[i][j]->aliveBM = veggrid[i][j]->AnnualaliveBM +
                               veggrid[i][j]->PerennialaliveBM +
                               veggrid[i][j]->ShrubaliveBM;
      // calculate the total edible biomass in a cell
      veggrid[i][j]->totalBM =
          veggrid[i][j]->reserveBM + veggrid[i][j]->aliveBM;
    }
  }
  // total edible biomass in whole grid
  double BMreserve = 0;
  double BMalive = 0;
  double BMtotal = 0;
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      BMtotal += veggrid[i][j]->totalBM;
      BMreserve += veggrid[i][j]->reserveBM;
      BMalive += veggrid[i][j]->aliveBM;
    }
  }
  // if biomass needed is less than available biomass, then determine how much
  // is taken from each grid cell and each type
  grazedBM = 0;  // total amount of Biomass that has been grazed so far
  for (int m = 0; m < PFTs; m++) {
    grazedPftBM[m] = 0;  // grazed biomass of PFTs
  }

  // compare demand with avaible biomass
  // if biomass needed is more than total edible biomass, livestock density
  // should be adjusted
  int st = stockingRate;
  if (BMtotal < BMneed) {
    double need = 0;
    do {
      need = 365 * percapability * BiomassNeed *
             (double)(1.0 / st);  // what they need per year and ha
      need = areagrazed * need;
      st++;
    } while (need > BMtotal);
    realstock = st;

    // remove edible part of reserve biomass and alive biomass except for safe
    // biomass in all cells
    for (int i = 0; i < xsize; i++) {
      for (int j = 0; j < ysize; j++) {
        for (int m = 0; m < PFTs; m++) {
          if (m < Shrubs) {
            veggrid[i][j]->pftList[m]->Biomass = veggrid[i][j]->safeBM_alive[m];
            veggrid[i][j]->pftList[m]->resBM = 0;
          } else {
            veggrid[i][j]->pftList[m]->Biomass = veggrid[i][j]->safeBM_alive[m];
            veggrid[i][j]->pftList[m]->resBM = veggrid[i][j]->safeBM_reserve[m];
          }
          grazedPftBM[m] += veggrid[i][j]->edibleBM_total[m];
          grazedBM += veggrid[i][j]->edibleBM_total[m];
        }
      }
    }
  } else {
    double demand_help;  // heterogeneity factor
    int zufx, zufy;      // random iterate cells until demand is fulfilled
    int count = 0;
    do {
      count++;
      // individual cells are chosen randomly and from every cell biomass is
      // removed interval of uniform int random number generator is [a,b], i.e.
      // includes upper bound! therefore range for zufx, zufy has to be
      // [0,xsize-1/ysize-1]
      zufx = RandomNumberGenerator::generate_random_int(0, xsize - 1);
      zufy = RandomNumberGenerator::generate_random_int(0, ysize - 1);
      // \todo why is this demand_help needed as additional variable, and not
      // avBMneed used instead?
      demand_help = avBMneed;  // calculate demand of specific cell
      // first determine how much is eaten by each type and then how much of
      // this comes from the reserve and alive biomass. Reasoning: if cattle
      // does not like type 1, then it will also not like the reserve biomass of
      // type 1

      double edibleBMG = 0;
      double edibleBMA = 0;
      double edibleBMS = 0;

      // after grazing once, the reserve biomass and alive biomass of PFTs in
      // one cell have changed, so should recalculate the edible part of PFTs
      for (int m = 0; m < PFTs; m++) {
        if (m < Shrubs) {
          veggrid[zufx][zufy]->edibleBM_alive[m] =
              veggrid[zufx][zufy]->pftList[m]->Biomass *
              (1 - veggrid[zufx][zufy]->pftList[m]->grazeLim);
          edibleBMS += veggrid[zufx][zufy]->edibleBM_alive[m];
        } else if (m >= Shrubs && m < Shrubs + Perennials) {
          veggrid[zufx][zufy]->edibleBM_total[m] =
              (veggrid[zufx][zufy]->pftList[m]->resBM +
               veggrid[zufx][zufy]->pftList[m]->Biomass) *
              (1 - veggrid[zufx][zufy]->pftList[m]->grazeLim);
          veggrid[zufx][zufy]->edibleBM_reserve[m] =
              veggrid[zufx][zufy]->pftList[m]->resBM *
              (1 - veggrid[zufx][zufy]->pftList[m]->grazeLim);
          ratio_Reserve[m] = veggrid[zufx][zufy]->edibleBM_reserve[m] /
                             veggrid[zufx][zufy]->edibleBM_total[m];
          edibleBMG += veggrid[zufx][zufy]->edibleBM_total[m];
        } else {
          veggrid[zufx][zufy]->edibleBM_total[m] =
              (veggrid[zufx][zufy]->pftList[m]->resBM +
               veggrid[zufx][zufy]->pftList[m]->Biomass) *
              (1 - veggrid[zufx][zufy]->pftList[m]->grazeLim);
          veggrid[zufx][zufy]->edibleBM_reserve[m] =
              veggrid[zufx][zufy]->pftList[m]->resBM *
              (1 - veggrid[zufx][zufy]->pftList[m]->grazeLim);
          ratio_Reserve[m] = veggrid[zufx][zufy]->edibleBM_reserve[m] /
                             veggrid[zufx][zufy]->edibleBM_total[m];
          edibleBMA += veggrid[zufx][zufy]->edibleBM_total[m];
        }
      }

      veggrid[zufx][zufy]->totalBM = edibleBMG + edibleBMA + edibleBMS;

      for (int m = 0; m < PFTs; m++) {
        if (veggrid[zufx][zufy]->ediblePerennial > 0 ||
            veggrid[zufx][zufy]->edibleAnnual > 0 ||
            veggrid[zufx][zufy]->edibleShrub > 0) {
          if (m < Shrubs) {
            reldemand[m] = veggrid[zufx][zufy]->relprefer[m] *
                           (veggrid[zufx][zufy]->edibleBM_alive[m] /
                            veggrid[zufx][zufy]->totalBM);
          } else {
            reldemand[m] = veggrid[zufx][zufy]->relprefer[m] *
                           (veggrid[zufx][zufy]->edibleBM_total[m] /
                            veggrid[zufx][zufy]->totalBM);
          }
        } else
          reldemand[m] = 0;
      }

      // another algorithm of grazing
      for (int m = 0; m < PFTs; m++) {
        if (m < Shrubs) {
          if (veggrid[zufx][zufy]->edibleBM_alive[m] >
              (reldemand[m] * demand_help)) {
            grazedPftBM[m] += (reldemand[m] * demand_help);
            grazedBM += (reldemand[m] * demand_help);
            veggrid[zufx][zufy]->edibleBM_alive[m] -=
                (reldemand[m] * demand_help);
            veggrid[zufx][zufy]->pftList[m]->Biomass =
                (veggrid[zufx][zufy]->pftList[m]->Biomass *
                     veggrid[zufx][zufy]->pftList[m]->grazeLim +
                 veggrid[zufx][zufy]->edibleBM_alive[m]);
          } else {  // in this cell ,shrub should not be grazed any more
            grazedPftBM[m] += veggrid[zufx][zufy]->edibleBM_alive[m];
            grazedBM += veggrid[zufx][zufy]->edibleBM_alive[m];
            veggrid[zufx][zufy]->pftList[m]->Biomass =
                veggrid[zufx][zufy]->pftList[m]->Biomass *
                veggrid[zufx][zufy]->pftList[m]->grazeLim;
          }
        } else {
          // if the sum of edible reserve biomass and edible alive biomass is
          // more than biomass need, biomass need will be distributed based on
          // the ratio of edible reserve biomass and edible alive biomass,
          // applicable on perennial grass and annual grass
          if (veggrid[zufx][zufy]->edibleBM_total[m] >
              (reldemand[m] * demand_help)) {
            grazedPftBM[m] += (reldemand[m] * demand_help);
            grazedBM += (reldemand[m] * demand_help);
            veggrid[zufx][zufy]->pftList[m]->resBM -=
                (reldemand[m] * demand_help * ratio_Reserve[m]);
            veggrid[zufx][zufy]->pftList[m]->Biomass -=
                (reldemand[m] * demand_help * (1 - ratio_Reserve[m]));
          }
          // if not, edible reserve biomass and edible alive biomass will be
          // removed, and left safe reserve biomass and safe alive biomass,
          // applicable to annual grass and perennial grass
          else {  // in this cell, after one grazing, annual grass should not be
                  // grazed any more
            grazedPftBM[m] += veggrid[zufx][zufy]->edibleBM_total[m];
            grazedBM += veggrid[zufx][zufy]->edibleBM_total[m];
            veggrid[zufx][zufy]->pftList[m]->resBM =
                veggrid[zufx][zufy]->pftList[m]->resBM *
                veggrid[zufx][zufy]->pftList[m]->grazeLim;
            veggrid[zufx][zufy]->pftList[m]->Biomass =
                veggrid[zufx][zufy]->pftList[m]->Biomass *
                veggrid[zufx][zufy]->pftList[m]->grazeLim;
          }
        }
      }

    } while (grazedBM < BMneed && count < (xsize * ysize * 5000));

    for (int i = Shrubs; i < Shrubs + Perennials; i++) {
      if (isnan(veggrid[zufx][zufy]->pftList[i]->Biomass)) {
        cout << "Grazing4"
             << "\t" << veggrid[zufx][zufy]->pftList[i]->Biomass << "\n";
      }
    }
  }
  // calculate grazing biomass for each PFT per ha
  grazedBM = (double)grazedBM /
             ((double)(xsize * ysize * cellsize * cellsize) / (100.0 * 100));
  for (int i = 0; i < PFTs; i++) {
    grazedPftBM[i] =
        (double)grazedPftBM[i] /
        ((double)(xsize * ysize * cellsize * cellsize) / (100.0 * 100));
  }

  // calculate the new vegetation cover and account for the destroyed biomass
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      for (int m = 0; m < PFTs; m++) {
        if (veggrid[i][j]->pftList[m]->Biomassold == 0)
          fraction[m] = 0;
        else
          fraction[m] = (double)(veggrid[i][j]->pftList[m]->Biomassold -
                                 veggrid[i][j]->pftList[m]->Biomass) /
                        veggrid[i][j]->pftList[m]->Biomassold;
        if (fraction[m] < 0)
          cout << "Grazing fraction zu klein"
               << "\t" << veggrid[i][j]->pftList[m]->Biomassold << "\t"
               << veggrid[i][j]->pftList[m]->Biomass << "\n";
        /// \implements \fn PlantFunctionalType::calculateCover(VegetationCell*
        /// vegetation, Weather* we2, int yr2, double factor)
        veggrid[i][j]->pftList[m]->calculateCover(veggrid[i][j], we, yr2,
                                                  fraction[m]);
      }
      // calculate the cover of reserve part for annual grass and perennial
      // grass
      /// \implements \fn VegetationCell::calculateResAnnualcover()
      veggrid[i][j]->calculateResAnnualcover();
      /// \implements \fn VegetationCell::calculateResPerennialcover()
      veggrid[i][j]->calculateResPerennialcover();
      /// \implements \fn VegetationCell::calculateResCover()
      veggrid[i][j]->calculateResCover();
    }
  }

  // plant reserve biomass is calculated after grazing, their biomass
  // contributes to a certain fraction to the fodder(Reserve biomass)
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      for (int m = 0; m < PFTs; m++) {
        if (m < Shrubs) {
          // There is no reserve biomass for shrub
          veggrid[i][j]->pftList[m]->resBM = 0;
        } else {
          // reserve biomass is decrease from the reserve biomass of last year
          veggrid[i][j]->pftList[m]->resBM *=
              veggrid[i][j]->pftList[m]->reserveLast;
          // reserve biomass receives from the biomass of this year
          veggrid[i][j]->pftList[m]->resBM +=
              veggrid[i][j]->pftList[m]->reserveThis *
              veggrid[i][j]->pftList[m]->Biomass;
          // calculate the safe part of reserve biomass
          veggrid[i][j]->safeBM_reserve[m] =
              veggrid[i][j]->pftList[m]->grazeLim *
              veggrid[i][j]->pftList[m]->resBM;
          // calculate the edible part of reserve biomass
          veggrid[i][j]->edibleBM_reserve[m] =
              veggrid[i][j]->pftList[m]->resBM -
              veggrid[i][j]->safeBM_reserve[m];
        }
      }
    }
  }
}

/*******************************************************************************************
 * calculate biodiversity indices (Richness, ...) including total PFTs present
 *in the landscape will add one to Richness if 0.5% cover are found in at least
 *1 cell is called in --> VegetationLandscape::calculateBiodiv
 *******************************************************************************************/

int VegetationLandscape::calculateRichness() {
  Double1D totalPFTcover;
  totalPFTcover.resize(PFTs, 0);
  int Richness{0};

  for (int i = 0; i < PFTs; i++) {
    [&] {  // I used a lambda function here
      for (int x = 0; x < xsize; x++) {
        for (int y = 0; y < ysize; y++) {
          if (veggrid[x][y]->pftList[i]->Cover > 0.05) {
            Richness++;
            return;
          }
        }
      }
    }();
  }
  return Richness;
}

/*******************************************************************************************
 * calculate biodiversity indices Shannon Diversity including total PFTs cover
 * is called in --> VegetationLandscape::calculateShannon
 *******************************************************************************************/

double VegetationLandscape::calculateShannon() {
  Double1D relativePFTcover;
  relativePFTcover.resize(PFTs, 0);
  double ShannonDiv{0};
  double totalPFTcover{0};

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int i = 0; i < PFTs; i++) {
        relativePFTcover[i] += veggrid[x][y]->pftList[i]->Cover;
      }
    }
  }

  for (int i = 0; i < PFTs; i++) {
    relativePFTcover[i] /= (1.0 * xsize * ysize);
    totalPFTcover += relativePFTcover[i];
  }

  for (int i = 0; i < PFTs; i++) {
    if (relativePFTcover[i] > 0)  //
      ShannonDiv += relativePFTcover[i] / totalPFTcover *
                    log(relativePFTcover[i] / totalPFTcover);
  }

  ShannonDiv *= -1;

  return ShannonDiv;
}

/*******************************************************************************************
 * calculate biodiversity index Evenness
 * is called in --> VegetationLandscape::calculateEvenness
 *******************************************************************************************/

double VegetationLandscape::calculateEvenness() {
  double Evenness = calculateShannon() / calculateRichness();

  return Evenness;
}

/*******************************************************************************************
 * calculate maximum, minimum and standard deviation of shrub cover for the
 *whole grid is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationLandscape::calculateShrub() { /*Tong, calculate the maximum,
                                                minimum of shrub cover in each
                                                cell*/

  for (int m = 0; m < Shrubs; m++) {
    MaxSC[m] = 0;
    MinSC[m] = 0;
  }
  // calculate the maximum shrub cover
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int m = 0; m < Shrubs; m++) {
        if (veggrid[x][y]->pftList[m]->Cover > MaxSC[m]) {
          MaxSC[m] = veggrid[x][y]->pftList[m]->Cover;
        }
      }
    }
  }

  // calculate the minimum shrub cover
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int m = 0; m < Shrubs; m++) {
        if (veggrid[x][y]->pftList[m]->Cover < MinSC[m]) {
          MinSC[m] = veggrid[x][y]->pftList[m]->Cover;
        }
      }
    }
  }
}

/*******************************************************************************************
 * calculate mean vegetation cover for the whole grid
 * is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationLandscape::calculateMean(
    int day) { /*Tong, based on PFT-oriented*/

  for (int i = 0; i < PFTs; i++) {
    meanCover[i] = 0;
    meanML1[i] = 0;
    meanML2[i] = 0;
    meandML1[i] = 0;
    meandML2[i] = 0;
    meanGrowthL1[i] = 0;
    meanGrowthL2[i] = 0;
    meanBasicMort[i] = 0;
    meanDryMort[i] = 0;
  }

  meanRCover = 0;
  meanGtotalcover = 0;
  meanStotalcover = 0;
  meanAtotalcover = 0;

  // meanG1BM = 0; meanS1BM = =; meanA1BM = 0;
  if (day == 0) {
    for (int i = 0; i < PFTs; i++) {
      meanAnnualGrowthL1[i] = 0;
      meanAnnualGrowthL2[i] = 0;
      meanAnnualBasicMort[i] = 0;
      meanGrowth[i] = 0;
      sumGrowth[i] = 0;
    }
  }

  if (day == growEnd) {
    for (int i = 0; i < PFTs; i++) {
      meanDryMortL1[i] = 0;
      meanDryMortL2[i] = 0;
      meanDisp[i] = 0;
      meanDispSend[i] = 0;
    }
  }

  meanMoistSeasonL1 = 0;
  meanMoistSeasonL2 = 0;
  meanMoistWetSeasonL1 = 0;
  meanMoistWetSeasonL2 = 0;

  for (int i = 0; i < PFTs; i++) {
    meanT1[i] = 0;
    meanT2[i] = 0;
  }

  meanTL1 = 0;
  meanTL2 = 0;
  meanE = 0;

  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int i = 0; i < PFTs; i++) {
        meanCover[i] += veggrid[x][y]->pftList[i]->Cover;
        meanML1[i] += veggrid[x][y]->pftList[i]->AvMoistL1;
        meanML2[i] += veggrid[x][y]->pftList[i]->AvMoistL2;
        meandML1[i] += veggrid[x][y]->pftList[i]->AvMoistSuppL1;
        meandML2[i] += veggrid[x][y]->pftList[i]->AvMoistSuppL2;
        meanGrowthL1[i] += veggrid[x][y]->pftList[i]->GrowthL1;
        meanGrowthL2[i] += veggrid[x][y]->pftList[i]->GrowthL2;
        meanBasicMort[i] += veggrid[x][y]->pftList[i]->BasicMort;
        sumGrowth[i] += (veggrid[x][y]->pftList[i]->GrowthL1 +
                         veggrid[x][y]->pftList[i]->GrowthL2);
      }

      meanRCover += veggrid[x][y]->resCover;

      if (day == growEnd) {
        for (int i = 0; i < PFTs; i++) {
          meanDryMort[i] += veggrid[x][y]->pftList[i]->DryMort;
          meanDryMortL1[i] += veggrid[x][y]->pftList[i]->DryMortL1;
          meanDryMortL2[i] += veggrid[x][y]->pftList[i]->DryMortL2;
          meanDisp[i] += veggrid[x][y]->pftList[i]->Disp;
          meanDispSend[i] += veggrid[x][y]->pftList[i]->DispSend;
        }
      }

      meanMoistSeasonL1 += veggrid[x][y]->moistSeasonL1;
      meanMoistSeasonL2 += veggrid[x][y]->moistSeasonL2;
      meanMoistWetSeasonL1 += veggrid[x][y]->moistWetSeasonL1;
      meanMoistWetSeasonL2 += veggrid[x][y]->moistWetSeasonL2;

      for (int i = 0; i < PFTs; i++) {
        meanT1[i] += veggrid[x][y]->yearlyPftTL1[i];
        meanT2[i] += veggrid[x][y]->yearlyPftTL2[i];
      }
      meanTL1 += veggrid[x][y]->yearlyTL1;
      meanTL2 += veggrid[x][y]->yearlyTL2;
      meanE += veggrid[x][y]->yearlyE;
    }
  }

  // number of grid cells to average over
  double numcells = (double)(xsize * ysize);

  if (meanRCover < 0.01) {
    meanRCover = 0;
  } else
    meanRCover = meanRCover / (numcells);

  for (int i = 0; i < PFTs; i++) {
    meanCover[i] = meanCover[i] / (numcells);
    meanML1[i] = meanML1[i] / (numcells);
    meanML2[i] = meanML2[i] / (numcells);
    meandML1[i] = meandML1[i] / (numcells);
    meandML2[i] = meandML2[i] / (numcells);
    meanGrowthL1[i] = meanGrowthL1[i] / (numcells);
    meanGrowthL2[i] = meanGrowthL2[i] / (numcells);
    meanBasicMort[i] = meanBasicMort[i] / (numcells);

    if (i < Shrubs)
      meanStotalcover += meanCover[i];
    else if (i >= Shrubs && i < Shrubs + Perennials)
      meanGtotalcover += meanCover[i];
    else
      meanAtotalcover += meanCover[i];
  }

  // calculate the standard deviation of shrub cover in all cells  /*Tong*/
  for (int i = 0; i < Shrubs; i++) {
    sumshrub[i] = 0;
  }
  for (int x = 0; x < xsize; x++) {
    for (int y = 0; y < ysize; y++) {
      for (int i = 0; i < Shrubs; i++) {
        sumshrub[i] +=
            pow(fabs(veggrid[x][y]->pftList[i]->Cover - meanCover[i]), 2.0);
      }
    }
  }
  for (int i = 0; i < Shrubs; i++) {
    SDSC[i] = sqrt(sumshrub[i] / (numcells - 1.0));
  }

  if (day == growEnd) {
    for (int i = 0; i < PFTs; i++) {
      meanDryMort[i] = meanDryMort[i] / (numcells);
      meanDryMortL1[i] = meanDryMortL1[i] / (numcells);
      meanDryMortL2[i] = meanDryMortL2[i] / (numcells);
      meanDisp[i] = meanDisp[i] / (numcells);
      meanDispSend[i] = meanDispSend[i] / (numcells);
      sumGrowth[i] =
          sumGrowth[i] / (numcells);  //_dirk summiert Wachstum uebers ganze
                                      //Jahr und mittelt auf die Zelle
    }
  }

  meanMoistSeasonL1 = meanMoistSeasonL1 / (numcells);
  meanMoistSeasonL2 = meanMoistSeasonL2 / (numcells);
  meanMoistWetSeasonL1 = meanMoistWetSeasonL1 / (numcells);
  meanMoistWetSeasonL2 = meanMoistWetSeasonL2 / (numcells);

  for (int i = 0; i < PFTs; i++) {
    meanT1[i] = meanT1[i] / (numcells);
    meanT2[i] = meanT2[i] / (numcells);
  }

  meanTL1 = meanTL1 / (numcells);
  meanTL2 = meanTL2 / (numcells);
  meanE = meanE / (numcells);

  if (day % vegTimeStep == vegTimeStep - 1) {
    for (int i = 0; i < PFTs; i++) {
      meanAnnualBasicMort[i] += meanBasicMort[i];
    }
  }

  if ((day % vegTimeStep == vegTimeStep - 1) && (day >= wetStart) &&
      (day <= wetEnd)) {
    for (int i = 0; i < PFTs; i++) {
      if (i < Shrubs + Perennials) {
        meanAnnualGrowthL1[i] += meanGrowthL1[i];
        meanAnnualGrowthL2[i] += meanGrowthL2[i];
        meanGrowth[i] += (meanAnnualGrowthL1[i] + meanAnnualGrowthL2[i]);
      } else {
        meanAnnualGrowthL1[i] += meanGrowthL1[i];
      }
    }
  }
}

/*******************************************************************************************
 * write a file with yearly results
 * is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationLandscape::writeYearlySummary(
    WaterLandscape* waterlandscape, Weather* we,
    int year) { /*Tong, based on PFT-oriented*/

  for (int m = 0; m < PFTs; m++) {
    grest[m] = 0;
  }

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      for (int m = 0; m < PFTs; m++) {
        grest[m] += veggrid[i][j]->pftList[m]->Biomass;
      }
    }
  }

  for (int m = 0; m < PFTs; m++) {
    grest[m] /= ((xsize * ysize * cellsize * cellsize) / (100.0 * 100));
  }

  int fir;
  if (fire2) {
    fir = 1;
  } else {
    fir = 0;
  }

  yearlySummaryFile << year << "\t" << fir << "\t";
  for (int i = Shrubs; i < Shrubs + Perennials; i++) {
    yearlySummaryFile << meanCover[i] << "\t" << meanB[i] << "\t" << grest[i]
                      << "\t" << grazedPftBM[i] << "\t" << meanT1[i] << "\t"
                      << meanT2[i] << "\t";
  }
  for (int i = 0; i < Shrubs; i++) {
    yearlySummaryFile << meanCover[i] << "\t" << MaxSC[i] << "\t" << MinSC[i]
                      << "\t" << SDSC[i] << "\t" << meanB[i] << "\t" << grest[i]
                      << "\t" << grazedPftBM[i] << "\t" << ScoverNr[i] << "\t"
                      << estCount[i] << "\t" << Cabove[0][i] << "\t"
                      << Cabove[2][i] << "\t" << Cabove[9][i] << "\t"
                      << meanT1[i] << "\t" << meanT2[i] << "\t";
  }
  for (int i = Shrubs + Perennials; i < PFTs; i++) {
    yearlySummaryFile << meanCover[i] << "\t" << meanB[i] << "\t" << grest[i]
                      << "\t" << grazedPftBM[i] << "\t" << meanT1[i] << "\t";
  }

  yearlySummaryFile << meanGtotalcover << "\t" << meanStotalcover << "\t"
                    << meanAtotalcover << "\t" << sCsumabove[0] << "\t"
                    << meanRCover << "\t" << realstock << "\t" << meanRB << "\t"
                    << grazedBM << "\t" << we->precSumYear[year] << "\t"
                    << meanE << "\t" << meanTL1 << "\t" << meanTL2 << "\t"
                    << meanMoistSeasonL1 << "\t" << meanMoistSeasonL2 << "\t"
                    << calculateRichness() << "\t" << calculateShannon() << "\t"
                    << calculateEvenness() << "\n";
}

/*******************************************************************************************
 * initialize vegetation cover
 * is called in --> VegetationLandscape::VegetationLandscape
 *******************************************************************************************/
void VegetationLandscape::initialiseVegCover() { /*Tong, based on PFT-oriented*/

  for (int i = 0; i < Shrubs; i++) {
    ScoverNr[i] = 0;
  }

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      double zuf = RandomNumberGenerator::generate_random_float(0.0, 1.0);
      for (int m = 0; m < Shrubs; m++) {
        if (zuf < SCoverbase[m]) {
          double zuf2 = RandomNumberGenerator::generate_random_float(0.0, 1.0);
          veggrid[i][j]->pftList[m]->Cover =
              zuf2 * (MaxSCover[m] - MinSCover[m]) + MinSCover[m];
          ScoverNr[m]++;

          /// \todo why these both for-loops within m-for-loop?
          for (int v = Shrubs; v < Shrubs + Perennials; v++) {
            veggrid[i][j]->pftList[v]->Cover = PcoverS[v - Shrubs];
          }
          for (int v = Shrubs + Perennials; v < PFTs; v++) {
            veggrid[i][j]->pftList[v]->Cover = 0;
          }
          veggrid[i][j]->cohortAge[m] =
              RandomNumberGenerator::generate_random_int(Mincohort[m],
                                                         Maxcohort[m]);
          veggrid[i][j]->cohortAgeFire[m] = 0;
        } else {
          veggrid[i][j]->pftList[m]->Cover = 0;

          /// \todo why these both for-loops within m-for-loop?
          for (int v = Shrubs; v < Shrubs + Perennials; v++) {
            veggrid[i][j]->pftList[v]->Cover =
                PerennialNS * veggrid[i][j]->pftList[v]->cmax;
            //                        veggrid[i][j]->pftList[v]->Cover =
            //                        PcoverNS[v-Shrubs] *
            //                        veggrid[i][j]->pftList[v]->cmax;
          }
          for (int v = Shrubs + Perennials; v < PFTs; v++) {
            veggrid[i][j]->pftList[v]->Cover = 0;
          }
          veggrid[i][j]->cohortAge[m] = 0;
          veggrid[i][j]->cohortAgeFire[m] = 0;
        }
      }
    }
  }
}

/*******************************************************************************************
 * initialise vegetation cover from files
 * is called in --> VegetationLandscape::VegetationLandscape
 * At the moment only for 3 meta pfts
 *******************************************************************************************/
void VegetationLandscape::initialiseVegCoverFromFile() {
  // read the initial vegetation cover matrix of shrub, perennial grass and
  // crust cover
  ifstream VegshrubcoverFile, VegperennialcoverFile;
  string VegshrubcoverName = "Parameters" + delimiter + "vegshrubcover_" +
                             std::to_string(xsize) + ".txt";
  string VegperennialcoverName = "Parameters" + delimiter +
                                 "vegperennialcover_" + std::to_string(xsize) +
                                 ".txt";
  // CrustcoverName = "Parameters" + delimiter + "crustcover" + ".txt";
  // firstly read the shrub cover
  cout << "Read shrub cover from: " << VegshrubcoverName << endl;
  VegshrubcoverFile.open(VegshrubcoverName.c_str(), ios::binary | ios::in);

  // if something goes wrong...
  if (!VegshrubcoverFile) {
    cerr << VegshrubcoverName << " could not be opened shrub cover file!\n";
    exit(-1);
  }

  // read in cover and close file
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      // VegshrubcoverFile >> vegshrubcover[i][j];
      for (int m = 0; m < Shrubs; m++) {
        VegshrubcoverFile >> veggrid[i][j]->pftList[m]->Cover;
        if (veggrid[i][j]->pftList[m]->Cover > 0) {
          veggrid[i][j]->cohortAge[m] =
              rand() % (Maxcohort[m] - Mincohort[m]);  // random(99)
          veggrid[i][j]->cohortAgeFire[m] = 0;
        } else {
          veggrid[i][j]->cohortAge[m] = 0;
          veggrid[i][j]->cohortAgeFire[m] = 0;
        }
      }
    }
  }
  VegshrubcoverFile.close();

  // second read the perennial cover
  cout << "Read perennial cover from: " << VegperennialcoverName << endl;
  VegperennialcoverFile.open(VegperennialcoverName.c_str(),
                             ios::binary | ios::in);

  // if something goes wrong...
  if (!VegperennialcoverFile) {
    cerr << VegperennialcoverName
         << " could not be opened perennial cover file!\n";
    exit(-1);
  }

  // read in cover and close file
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      for (int v = Shrubs; v < Shrubs + Perennials; v++) {
        VegperennialcoverFile >> veggrid[i][j]->pftList[v]->Cover;
      }
    }
  }
  VegperennialcoverFile.close();

  // initialize Annual Vegetation cover with 0
  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      for (int v = Shrubs + Perennials; v < PFTs; v++) {
        veggrid[i][j]->pftList[v]->Cover = 0;
      }
    }
  }
}

/*******************************************************************************************
 * initialise vegetation parameters
 * is called in --> VegetationLandscape::VegetationLandscape
 *******************************************************************************************/
void VegetationLandscape::initialiseVegParams(
    Parameter* p) { /*Tong, based on PFT-oriented*/

  // constant of the vegetation model
  double aGrowthConst_, resWater_, moistL1_, moistL2_, overlap_, grazeParam_,
      MAP_, bm_c_rain_, grazeparA_, grazeparB_, firecoeff_, fuelLim_, EnScover_,
      initialshare_;

  // open vegetation parameters file
  ifstream parameterFile;
  string vegetationParametersName;

  vegetationParametersName = "Parameters" + delimiter +
                             "vegetationparameters_" +
                             p->vegetationFileID[p->scenario] + ".txt";

  cout << "Read vegetation parameters from: " << vegetationParametersName
       << endl;

  // open file
  parameterFile.open(vegetationParametersName.c_str(), ios::binary | ios::in);

  // if something goes wrong...
  if (!parameterFile) {
    cerr << vegetationParametersName << " could not be opened!\n";
    exit(-1);
  }

  int fsize = 1000;
  parameterFile.ignore(fsize, '#');
  parameterFile.ignore(fsize, '#');
  // read in constant parameters of the model
  parameterFile.ignore(fsize, ':');
  parameterFile >> resWater_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> moistL1_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> moistL2_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> maxDistS;
  parameterFile.ignore(fsize, ':');
  parameterFile >> growStart;
  parameterFile.ignore(fsize, ':');
  parameterFile >> growEnd;
  parameterFile.ignore(fsize, ':');
  parameterFile >> aGrowthConst_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> sclim;
  parameterFile.ignore(fsize, ':');
  parameterFile >> scrlim;
  parameterFile.ignore(fsize, ':');
  parameterFile >> pgshare;
  parameterFile.ignore(fsize, ':');
  parameterFile >> inishare;
  parameterFile.ignore(fsize, ':');
  parameterFile >> overlap_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> grazeParam_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> bm_c_rain_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> MAP_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> grazeparA_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> grazeparB_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> firecoeff_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> fuelLim_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> EnScover_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> initialshare_;
  parameterFile.ignore(fsize, ':');
  parameterFile >> stockingRate;
  parameterFile.ignore(fsize, ':');
  parameterFile >> percapability;
  parameterFile.ignore(fsize, ':');
  parameterFile >> BiomassNeed;
  parameterFile.ignore(fsize, ':');
  parameterFile >> fire;
  parameterFile.ignore(fsize, ':');
  parameterFile >> Shrubs;
  parameterFile.ignore(fsize, ':');
  parameterFile >> Perennials;
  parameterFile.ignore(fsize, ':');
  parameterFile >> Annuals;
  parameterFile.ignore(fsize, ':');
  parameterFile >> PerennialNS;

  PerennialNS = PerennialNS / Perennials;

  PFTs = Annuals + Perennials + Shrubs;

  // initialize vectors
  yearlyPTL1.resize(p->simYears, Double1D(Perennials, 0));
  veggrid.resize(xsize, vector<VegetationCell*>(ysize));
  meanML1.resize(PFTs, 0);
  meanML2.resize(PFTs, 0);
  meandML1.resize(PFTs, 0);
  meandML2.resize(PFTs, 0);
  meanGrowthL1.resize(PFTs, 0);
  meanGrowthL2.resize(PFTs, 0);
  meanBasicMort.resize(PFTs, 0);
  meanDryMort.resize(PFTs, 0);
  meanDryMortL1.resize(PFTs, 0);
  meanDryMortL2.resize(PFTs, 0);
  meanAGrowthL1.resize(Annuals, 0);
  meanDisp.resize(PFTs, 0);
  meanDispSend.resize(PFTs, 0);
  meanC.resize(PFTs, 0);
  meanB.resize(PFTs, 0);
  meanBurnedG.resize(Perennials, 0);
  meanBurnedS.resize(Shrubs, 0);
  meanAnnualGrowthL1.resize(PFTs, 0);
  meanAnnualGrowthL2.resize(PFTs, 0);
  meanAnnualBasicMort.resize(PFTs, 0);
  meanEndSeasonCoverG.resize(Perennials, 0);
  meanEndSeasonCoverS.resize(Shrubs, 0);

  estCount.resize(Shrubs, 0);
  meanCover.resize(PFTs, 0);
  Maxcohort.resize(Shrubs, 0);
  Mincohort.resize(Shrubs, 0);
  MaxSC.resize(Shrubs, 0);
  MinSC.resize(Shrubs, 0);
  SDSC.resize(Shrubs, 0);
  SCoverbase.resize(Shrubs, 0);
  MaxSCover.resize(Shrubs, 0);
  MinSCover.resize(Shrubs, 0);
  PcoverS.resize(Perennials, 0);
  //    PcoverNS.resize(Perennials, 0);
  meanT1.resize(PFTs, 0);
  meanT2.resize(PFTs, 0);
  Cabove.resize(10, Int1D(PFTs, 0));
  sCsumabove.resize(10, 0);  /// \todo is 10 correct?
  ScoverNr.resize(Shrubs, 0);
  dispHelpS.resize(Shrubs, 0);
  zuf_dispS.resize(Shrubs, 0);
  newdistconst.resize(Shrubs, 0);
  dispHelp.resize(PFTs, 0);
  grazedPftBM.resize(PFTs, 0);
  fraction.resize(PFTs, 0);
  grest.resize(PFTs, 0);
  reldemand.resize(PFTs, 0);
  ratio_Reserve.resize(PFTs, 0);
  netMortR.resize(PFTs, 0);
  meanGrowth.resize(PFTs, 0);
  sumGrowth.resize(PFTs, 0);
  sumshrub.resize(Shrubs, 0);

  name_.resize(PFTs);
  type_.resize(PFTs);
  UptakeRate_.resize(PFTs, 0);
  RootL1_.resize(PFTs, 0);
  RootL2_.resize(PFTs, 0);
  GrowthR_.resize(PFTs, 0);
  Transcoef_.resize(PFTs, 0);
  Kcoef_.resize(PFTs, 0);
  ecoef_.resize(PFTs, 0);
  temcoef_.resize(PFTs, 0);
  LAI_C_.resize(PFTs, 0);
  MortR_.resize(PFTs, 0);
  MortBasic_.resize(PFTs, 0);
  Grazed_.resize(PFTs, 0);
  WPWaterL1_.resize(PFTs, 0);
  WPWaterL2_.resize(PFTs, 0);
  dispFrac_.resize(PFTs, 0);
  cmax_.resize(PFTs, 0);
  conversionC_BM_.resize(PFTs, 0);
  grazeLim_.resize(PFTs, 0);
  grazeprefer_.resize(PFTs, 0);
  grazepreferintra_.resize(PFTs, 0);
  estW_.resize(PFTs, 0);
  specLoss_.resize(PFTs, 0);
  dist0_.resize(PFTs, 0);
  distConst_.resize(PFTs, 0);
  distConst2_.resize(PFTs, 0);
  reserveLast_.resize(PFTs, 0);
  reserveThis_.resize(PFTs, 0);

  for (int k = 0; k < Shrubs; k++) {
    parameterFile.ignore(fsize, '*');  // read parameters from first shrub
    parameterFile.ignore(fsize, ':');
    parameterFile >> name_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> type_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> UptakeRate_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> RootL1_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> RootL2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> GrowthR_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Transcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Kcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> ecoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> temcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> LAI_C_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MortR_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MortBasic_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Grazed_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> WPWaterL1_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> WPWaterL2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> cmax_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> conversionC_BM_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazeLim_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazeprefer_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazepreferintra_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> dispFrac_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> estW_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> specLoss_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> dist0_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> distConst_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> distConst2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> reserveLast_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> reserveThis_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Maxcohort[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Mincohort[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> SCoverbase[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MaxSCover[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MinSCover[k];
  }

  for (int k = Shrubs; k < Shrubs + Perennials; k++) {
    parameterFile.ignore(fsize,
                         '*');  // read parameters from first perennial grass
    parameterFile.ignore(fsize, ':');
    parameterFile >> name_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> type_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> UptakeRate_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> RootL1_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> RootL2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> GrowthR_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Transcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Kcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> ecoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> temcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> LAI_C_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MortR_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MortBasic_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Grazed_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> WPWaterL1_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> WPWaterL2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> cmax_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> conversionC_BM_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazeLim_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazeprefer_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazepreferintra_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> dispFrac_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> estW_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> specLoss_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> dist0_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> distConst_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> distConst2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> reserveLast_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> reserveThis_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> PcoverS[k - Shrubs];
    //        parameterFile.ignore(fsize, ':');
    //        parameterFile >> PcoverNS[k-Shrubs];
  }

  for (int k = Shrubs + Perennials; k < PFTs; k++) {
    parameterFile.ignore(fsize, '*');  // read parameters from first annual
                                       // grass
    parameterFile.ignore(fsize, ':');
    parameterFile >> name_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> type_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> UptakeRate_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> RootL1_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> RootL2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> GrowthR_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Transcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Kcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> ecoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> temcoef_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> LAI_C_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MortR_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> MortBasic_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> Grazed_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> WPWaterL1_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> WPWaterL2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> cmax_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> conversionC_BM_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazeLim_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazeprefer_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> grazepreferintra_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> dispFrac_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> estW_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> specLoss_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> dist0_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> distConst_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> distConst2_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> reserveLast_[k];
    parameterFile.ignore(fsize, ':');
    parameterFile >> reserveThis_[k];
  }

  parameterFile.close();

  for (int i = 0; i < xsize; i++) {
    for (int j = 0; j < ysize; j++) {
      veggrid[i][j] = new VegetationCell(
          PFTs, Annuals, Perennials, Shrubs, resWater_, moistL1_, moistL2_,
          overlap_, grazeParam_, aGrowthConst_, bm_c_rain_, grazeparA_,
          grazeparB_, MAP_, firecoeff_, fuelLim_, EnScover_, initialshare_,
          vegTimeStep, cellsize);

      for (int k = 0; k < Shrubs; k++) {
        veggrid[i][j]->pftList[k] = new Shrub(
            name_[k], type_[k], UptakeRate_[k], RootL1_[k], RootL2_[k],
            GrowthR_[k], Transcoef_[k], Kcoef_[k], ecoef_[k], temcoef_[k],
            LAI_C_[k], MortR_[k], MortBasic_[k], Grazed_[k], WPWaterL1_[k],
            WPWaterL2_[k], cmax_[k], overlap_, conversionC_BM_[k], grazeParam_,
            grazeLim_[k], grazeprefer_[k], grazepreferintra_[k], dispFrac_[k],
            estW_[k], bm_c_rain_, MAP_, grazeparA_, grazeparB_, specLoss_[k],
            firecoeff_, fuelLim_, dist0_[k], distConst_[k], distConst2_[k],
            reserveLast_[k], reserveThis_[k], vegTimeStep, cellsize);
      }

      for (int k = Shrubs; k < Shrubs + Perennials; k++) {
        veggrid[i][j]->pftList[k] = new Perennial(
            name_[k], type_[k], UptakeRate_[k], RootL1_[k], RootL2_[k],
            GrowthR_[k], Transcoef_[k], Kcoef_[k], ecoef_[k], temcoef_[k],
            LAI_C_[k], MortR_[k], MortBasic_[k], Grazed_[k], WPWaterL1_[k],
            WPWaterL2_[k], cmax_[k], overlap_, conversionC_BM_[k], grazeParam_,
            grazeLim_[k], grazeprefer_[k], grazepreferintra_[k], dispFrac_[k],
            estW_[k], bm_c_rain_, MAP_, grazeparA_, grazeparB_, specLoss_[k],
            firecoeff_, fuelLim_, dist0_[k], distConst_[k], distConst2_[k],
            reserveLast_[k], reserveThis_[k], vegTimeStep, cellsize);
      }

      for (int k = Shrubs + Perennials; k < PFTs; k++) {
        veggrid[i][j]->pftList[k] = new Annual(
            name_[k], type_[k], UptakeRate_[k], RootL1_[k], RootL2_[k],
            GrowthR_[k], Transcoef_[k], Kcoef_[k], ecoef_[k], temcoef_[k],
            LAI_C_[k], MortR_[k], MortBasic_[k], Grazed_[k], WPWaterL1_[k],
            WPWaterL2_[k], cmax_[k], overlap_, conversionC_BM_[k], grazeParam_,
            grazeLim_[k], grazeprefer_[k], grazepreferintra_[k], dispFrac_[k],
            estW_[k], bm_c_rain_, MAP_, grazeparA_, grazeparB_, specLoss_[k],
            firecoeff_, fuelLim_, dist0_[k], distConst_[k], distConst2_[k],
            reserveLast_[k], reserveThis_[k], vegTimeStep, cellsize);
      }
    }
  }
}

/*******************************************************************************************
 * initialise file for vegetation cover time series
 * is called in --> VegetationLandscape::VegetationLandscape
 *******************************************************************************************/
void VegetationLandscape::initialiseResultFile(
    Parameter* p) { /*Tong, based on PFT-oriented*/

  yearlySummaryFileName =
      "Results" + delimiter + p->site + "_" + std::to_string(p->simYears) +
      "years_yearly_" + p->outputFileID[p->scenario] + "_climrep-" +
      std::to_string(p->climateRepetition + 1) + "_modelrep-" +
      std::to_string(p->modelRepetition + 1) + ".txt";

  if (!yearlySummaryFile) {
    cout << yearlySummaryFileName << " could not be opened!\n";
  }
  yearlySummaryFile.open(yearlySummaryFileName.c_str(), ios::out | ios::trunc);

  // if something goes wrong...
  if (!yearlySummaryFile) {
    cout << yearlySummaryFileName
         << " could not be opened in vegLandscape->initialiseResultFile!\n";
  }

  yearlySummaryFile << "year"
                    << "\t"
                    << "fire"
                    << "\t";

  for (int i = 0; i < Perennials; i++) {
    yearlySummaryFile << "meanGCover" << i << "\t"
                      << "meanGB" << i << "\t"
                      << "pgrest" << i << "\t"
                      << "grazedGBM" << i << "\t"
                      << "PTranspirationL1" << i << "\t"
                      << "PTranspirationL2" << i << "\t";
  }

  for (int i = 0; i < Shrubs; i++) {
    yearlySummaryFile << "meanSCover" << i << "\t"
                      << "MaxSCover" << i << "\t"
                      << "MinSCover" << i << "\t"
                      << "SDSCover" << i << "\t"
                      << "meanSB" << i << "\t"
                      << "srest" << i << "\t"
                      << "grazedSBM" << i << "\t"
                      << "shrubCells" << i << "\t"
                      << "estCells" << i << "\t"
                      << "sc01" << i << "\t"
                      << "sc1" << i << "\t"
                      << "sc10" << i << "\t"
                      << "STranspirationL1" << i << "\t"
                      << "STranspirationL2" << i << "\t";
  }

  for (int i = 0; i < Annuals; i++) {
    yearlySummaryFile << "meanACover" << i << "\t"
                      << "meanAB" << i << "\t"
                      << "agrest" << i << "\t"
                      << "grazedABM" << i << "\t"
                      << "ATranspirationL1" << i << "\t";
  }

  yearlySummaryFile << "meanGtotalcover"
                    << "\t"
                    << "meanStotalcover"
                    << "\t"
                    << "meanAtotalcover"
                    << "\t"
                    << "scsum01"
                    << "\t"
                    << "meanRCover"
                    << "\t"
                    << "effSR"
                    << "\t"
                    << "meanRB"
                    << "\t"
                    << "BMGrazed"
                    << "\t"
                    << "AnnualRain"
                    << "\t"
                    << "Annualevaporation"
                    << "\t"
                    << "AnnualtranspirationL1"
                    << "\t"
                    << "AnnualtranspirationL2"
                    << "\t"
                    << "ML1"
                    << "\t"
                    << "ML2"
                    << "\t"
                    << "Richness"
                    << "\t"
                    << "ShannonDiv"
                    << "\t"
                    << "Evenness"
                    << "\n";
}

/*******************************************************************************************
 * \todo nothing happens here!?
 * is called in --> controller::runSimulation
 *******************************************************************************************/
void VegetationLandscape::finishSimulation() { //yearlySummaryFile.close(); 
}

/*******************************************************************************************
 * constructor for class VegetationLandscape
 * is called in --> controller::initialiseSimulation
 *******************************************************************************************/
VegetationLandscape::VegetationLandscape(
    Parameter* p,
    WaterLandscape* waterlandscape) {  // constructor has been changed by Tong

  // initialize parameters and result files
  xsize = p->xsize;
  ysize = p->ysize;
  vegTimeStep = p->vegTimeStep;
  cellsize = p->cellsize;

  restings = 0;
  fires = 0;
  oldstock = stockingRate;

  /// \implements \fn VegetationLandscape::initialiseVegParams(Parameter* p)
  initialiseVegParams(p);
  /// \implements \fn VegetationLandscape::initialiseResultFile(Parameter* p)
  //initialiseResultFile(p);
  /// \implements \fn VegetationLandscape::initialiseVegCover()
  initialiseVegCoverFromFile();
}

/*******************************************************************************************
 * destructor for class VegetationLandscape
 * is called in --> controller::deleteObjects
 *******************************************************************************************/
VegetationLandscape::~VegetationLandscape() {
  for (unsigned int i = 0; i < veggrid.size(); i++)
    for (unsigned int j = 0; j < veggrid.at(i).size(); j++) {
      delete veggrid[i][j];
      veggrid[i][j] = NULL;
    }
}

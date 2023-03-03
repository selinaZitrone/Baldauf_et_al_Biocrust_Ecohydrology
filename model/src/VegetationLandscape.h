/***************************************************************************************
 * \class   VegetationLandscape VegetationLandscape.h
 * \brief   In this file, all functions are described that are relevant for
 *vegetation dynamics in a landscape. \author  Britta Tietjen
 *******************************************************************************************/

#ifndef VEGETATIONLANDSCAPE_H_
#define VEGETATIONLANDSCAPE_H_

#include <math.h>  //to use mathematical functions

#include <algorithm>  //to use min/max
#include <cstdlib>    //to use exit() and for srand() and rand()
#include <fstream>    //to read and write files
#include <iostream>   //to use cin / cout
#include <numeric>
#include <string>  //to use strings
#include <vector>

#include "Parameters.h"
#include "VegetationCell.h"
#include "Weather.h"


using namespace std;

// Forward declaration
class PlantFunctionalType;
class WaterLandscape;

class VegetationLandscape {
 private:
  // Private parameters
  bool fire;
  bool fire2;
  bool firemem;
  fstream yearlySummaryFile;
  string yearlySummaryFileName;

  Double1D meanML1;  //!< mean available moisture in layer 1
  Double1D meanML2;  //!< mean available moisture in layer 2
  Double1D
      meandML1;  //!< mean moisture to support life in layer 1 (dryMortality)
  Double1D
      meandML2;  //!< mean moisture to support life in layer 2 (dryMortality)

  Double1D meanGrowthL1;  //!< mean growth and mortality
  Double1D meanGrowthL2;
  Double1D meanBasicMort;
  Double1D meanDryMort;

  Double1D meanDryMortL1;
  Double1D meanDryMortL2;

  Double1D meanAGrowthL1;  //_dirk
  Double1D meanDisp;       //!< dispersal that is reached by a cell
  Double1D meanDispSend;   //!< dispersal that is dispersed from a cell

  Double1D meanC;  //!< mean cover
  Double1D meanB;  //!< mean total biomass per ha (g/ha)
  double meanRB;   // _dirk

  Double1D meanBurnedG;  //!< mean cover per cell that was burned
  Double1D meanBurnedS;  //!< mean cover per cell that was burned

  double meanMoistWetSeasonL1;  //!< mean soil moistur during wet season
  double meanMoistWetSeasonL2;
  Double1D meanAnnualGrowthL1;  //!< total growth during wet season
  Double1D meanAnnualGrowthL2;
  Double1D meanAnnualBasicMort;  //!< annual basic mortality
  Double1D meanEndSeasonCoverG;  //!< vegetation cover at the end of the growing
                                 //!< or wet season
  Double1D meanEndSeasonCoverS;

  int xsize;
  int ysize;
  int vegTimeStep;
  int cellsize;

  // Private member functions
  /***************************************************************************************
   * \brief Calculate total cover in each grid cell.
   *******************************************************************************************/
  void calculateTotalCover();

  /***************************************************************************************
   * \brief Calculate mean vegetation cover and other things for the whole grid.
   * \param p pointer to parameter object;
   *******************************************************************************************/
  void calculateMean(int day);

  /***************************************************************************************
   * \brief Calculate maximum, minimum and standard deviation of shrub cover for
   *the whole grid.
   *******************************************************************************************/
  void calculateShrub();

  /***************************************************************************************
   * \brief Calculate mean cover in the grid and biomass per ha.
   *******************************************************************************************/
  void meanCoverBiomass();

  /***************************************************************************************
   * \brief  Organizes fire, calls calculatePFTfire and stores fire results
   * only called, if a fire actually occurs
   *******************************************************************************************/
  void calculateFire();

  /***************************************************************************************
   * \brief Initialise file for veg cover time series.
   * \todo Explain parameters and include brief description here
   * \param p pointer to parameter object;
   * \param fire fire yes or not;
   *******************************************************************************************/
  void initialiseResultFile(Parameter* p);

  /***************************************************************************************
   * \brief Write a file with yearly results.
   * \todo Explain parameters and include brief description here
   * \param waterlandscape pointer to WaterLandscape object;
   * \param we pointer to weather object;
   * \param year current year;
   *******************************************************************************************/
  void writeYearlySummary(WaterLandscape* waterlandscape, Weather* we,
                          int year);

  /***************************************************************************************
   * \brief Calculate diversity indices Richness, Shannon & Evenness
   *******************************************************************************************/
  int calculateRichness();
  double calculateShannon();
  double calculateEvenness();

 public:
  // Constructor and Destructor
  /***************************************************************************************
   * \brief Constructor
   * \param p pointer to Parameter object;
   * \param waterlandscape pointer to WaterLandscape object;
   *******************************************************************************************/
  VegetationLandscape(Parameter* p, WaterLandscape* waterlandscape);

  /***************************************************************************************
   * \brief Destructor
   * \param p pointer to parameter object;
   *******************************************************************************************/
  ~VegetationLandscape();

  // Public member functions
  /***************************************************************************************
   * \brief Initialize vegetation cover.
   *******************************************************************************************/
  void initialiseVegCover();

  void initialiseVegCoverFromFile();
  /***************************************************************************************
   * \brief Initialise vegetation parameters.
   * \todo Explain parameters here.
   * \param p pointer to parameter object;
   *******************************************************************************************/
  void initialiseVegParams(Parameter* p);

  /***************************************************************************************
   * \brief All concerning vegetation change are controlled here.
   * \param year current year;
   * \param day current day;
   * \param waterlandscape pointer to WaterLandscape object;
   * \param we pointer to weather object;
   * \param p pointer to parameter object;
   *******************************************************************************************/
  void calculateProcesses(int year, int day, WaterLandscape* waterlandscape,
                          Weather* we, Parameter* p);

  /***************************************************************************************
   * \brief Calculate "diffusive" dispersal:
   * Each plant type (shrubs / grasses) has
   * a defined dispersal function (decreasing with distance). Within this
   *distance a small fraction of the plant (e.g. by seeds, clonal) is dispersed.
   * \todo Explain parameters here.
   * \param Stockrate ???;
   *******************************************************************************************/
  void dispersal(int Stockrate);

  /***************************************************************************************
   * \brief Calculate "diffusive" dispersal for grasses and shrubs: Each plant
   *type has a defined dispersal function (decreasing with distance). Within
   *this distance a small fraction of the plant (e.g. by seeds, clonal) is
   *dispersed. \todo Explain parameters here. \param x1 ???; \param y1 ???;
   * \param x2 ???;
   * \param y2 ???;
   *******************************************************************************************/
  void calculateShrubDispersal(int x1, int y1, int x2, int y2);

  /***************************************************************************************
   * \brief ???
   * \todo Explain parameters and include description here.
   * \param we pointer to weather object;
   * \param yr2 current year;
   *******************************************************************************************/
  void grazing(Weather* we, int yr2); /*Tong*/

  /***************************************************************************************
   * \brief ???
   * \todo Include description here.
   *******************************************************************************************/
  void finishSimulation();

  // Public parameters
  Int1D estCount;  //!< counter for cells in which establishments of seeds was
                   //!< possible
  Double1D meanCover;  //!< mean cover for the whole grid
  double meanRCover;
  Double1D MaxSC;  //!< maximum shrub cover for all cells /*Tong*/
  Double1D MinSC;  //!< minimum shrub cover for all cells /*Tong*/
  Double1D SDSC;   //!< standard deviation of shrub cover for all cells /*Tong*/

  Double1D SCoverbase;  //!< cells for generated shrub cover /*Tong*/
  Double1D MaxSCover;   //!< Maximum shrub cover generated /*Tong*/
  Double1D MinSCover;   //!< Minimum shrub cover generated /*Tong*/
  Double1D PcoverS;     //!< perennial cover for cells having shrub /*Tong*/
  Int1D Maxcohort;      //!< Maximum cohort age for shrub /*Tong*/
  Int1D Mincohort;      //!< Minimum cohort age for shrub /*Tong*/
  double PerennialNS;   // max perennial cover divided by number of perennilas
                        // when having no shrub // Katja
  Double1D PcoverNS;    //!< perennial cover for cells having no shrub /*Tong*/

  int PFTs;
  int Perennials; /*Tong*/
  int Shrubs;     /*Tong*/
  int Annuals;    /*Tong*/

  Double2D yearlyPTL1;

  Double1D meanT1;  //!< mean yearly transpiration in the upper layer /*Tong*/
  Double1D meanT2;  //!< mean yearly transpiration in the lower layer /*Tong*/
  double meanTL1;   //!< mean yearly transpiration in the upper layer /*Tong*/
  double meanTL2;   //!< mean yearly transpiration in the lower layer /*Tong*/
  double meanE;     //!< mean yearly evaporation in the upper layer /*Tong*/

  double meanMoistSeasonL1;  //!< mean soil moisture during 'official' growing
                             //!< season
  double meanMoistSeasonL2;
  Int2D Cabove;     //!< _dirk to count number of cells above 10 levels of cover
                    //!< (10%,20%...90%)
  Int1D sCsumabove; /*Tong*/
  Int1D ScoverNr;   /*Tong*/

  int growStart;  //!< when does the 'official' growing season start
  int growEnd;    //!< and when does it end
  int wetStart;   //!< first rain > 5mm
  int wetEnd;     //!< last rain > 5mm + 14

  // equation: dispFracS*dist0*exp(-frac*distance)
  double maxDistS;
  Double1D dispHelpS;
  Double1D dispHelp;
  Double1D zuf_dispS;
  Double1D newdistconst;

  double sclim;   //_dirk
  double scrlim;  //_dirk

  double pgshare;   //_dirk
  double inishare;  //_dirk

  double meanGtotalcover;  //!< calculate the mean total cover of all the
                           //!< perennial grass for the whole cells
  double meanStotalcover;  //!< calculate the mean total cover of all the shrub
                           //!< for the whole cells
  double meanAtotalcover;  //!< calculate the mean total cover of all the annual
                           //!< grass for the whole cells

  int BMlack;
  double avBMneed;
  int BMneed;
  double grazedBM;
  Double1D grazedPftBM;
  Double1D fraction;
  int stockingRate;
  int realstock;   //!< realized stocking rate
  Double1D grest;  //!< rest biomass after grazing

  // new parameter of grazing process
  double meanprefer;  //!< mean of grazing preference               /*Tong*/
  Double1D
      reldemand;  //!< relative grazing biomass demand for each PFT /*Tong*/
  Double1D ratio_Reserve; /*Tong*/
  double percapability;   //!< the amount of the livestock body mass per day for
                          //!< biomass demand
  double
      BiomassNeed;  //!< the related parameter to calculate the biomass per year

  int oldstock;
  int restings;
  int fires;
  Double1D netMortR;    //!< _dirk mean mortality per cell
  Double1D meanGrowth;  //!< mean growth per cell
  Double1D sumGrowth;
  Double1D sumshrub; /*Tong*/

  // parameter about plant functional types                    /*Tong*/
  vector<string> name_, type_;
  Double1D UptakeRate_, RootL1_, RootL2_, GrowthR_, Transcoef_, Kcoef_, ecoef_,
      temcoef_, LAI_C_, MortR_, MortBasic_, Grazed_, WPWaterL1_, WPWaterL2_,
      dispFrac_, cmax_, conversionC_BM_, grazeLim_, grazeprefer_,
      grazepreferintra_, estW_, specLoss_, dist0_, distConst_, distConst2_,
      reserveLast_, reserveThis_;

  // Fields
  vector<vector<VegetationCell*> > veggrid;
};

#endif  // VEGETATIONLANDSCAPE_H

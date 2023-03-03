/****************************************************************************************
 * \class   WaterLandscape WaterLandscape.h
 * \brief   This file describes the water dynamics in the landscape.
 *          It contains a water landscape which has lots of grid cells
 *(soilCell). \author  Britta Tietjen
 *******************************************************************************************/

#ifndef WATERLANDSCAPE_H_
#define WATERLANDSCAPE_H_

#include <fstream>  //to read and write files
#include <string>   //to use strings

#include "BiocrustLandscape.h"
#include "Parameters.h"
#include "VegetationLandscape.h"
#include "WaterCell.h"
#include "Weather.h"

//! An enum.
/*! Constants to describe aspect. */
enum Aspects {
  FL = 0,  //!< enum value 1 for flat, no aspect = 0
  NN,      //!< enum value 2 for north = 1
  NE,      //!< enum value 3 for north-east = 2
  EE,      //!< enum value 4 for east = 3
  SE,      //!< enum value 5 for south-east = 4
  SS,      //!< enum value 6 for south = 5
  SW,      //!< enum value 7 for south-west = 6
  WW,      //!< enum value 8 for west = 7
  NW       //!< enum value 9 for north-west = 8
};

const double PI = 3.14159265358979323846264338327950288;

// Forward declaration
class CellList;

class WaterLandscape {
 private:
  // Private parameters
  Double2D elevation;  //!< elevation of each cell
  Int2D outputCells;   //!< elevation of each cell
  Int2D aspect;        //!< aspect of each cell --> direction of inclination
  Double2D
      inclination;  //!< inclination of cell in radians --> to lowest neighbor
  Double2D
      inclination_deg;     //!< inclination of cell to lowest neighbor (degree)
  Int2D runoncells;        //!< number of cells that provide runon
  Int2D runoffdirection;   //!< which is the lowest neighbor
  ofstream resultFile;     //!< file to store the results...
  ofstream resultFileDay;  //!< daily mean values (or sums, whatever makes more
                           //!< sense)
  ofstream resultFileSelina;  // Daily hydrology output for Selina' analysis
  ofstream resultFileSelinaHour;
  ofstream resultFileSelinaHour_CS;

  double potEP;            //!< potential evaporation [mm/h]
  double meanWaterL1_rel;  //!< mean relative soil moisture (of the landscape)
                           //!< in the upper layer
  double meanWaterL2_rel;  //!< mean relative soil moisture (of the landscape)
                           //!< in the lower layer
  double sdWaterL1_rel = 0.0;  // standard deviation of landscape soil moisture
  double sdWaterL2_rel = 0.0;

  double meanWaterL0_abs;  //!< mean surface water (of the landscape)
  double meanWaterL1_abs;  //!< mean absolute soil moisture (of the landscape)
                           //!< in the upper layer
  double meanWaterL2_abs;  //!< mean absolute soil moisture (of the landscape)
                           //!< in the lower layer
  double meanWater_rel;  //!< mean relative soil moisture (of the landscape) in
                         //!< both layers
  double meanWater_abs;  //!< mean absolute soil moisture (of the landscape) in
                         //!< both layers
  double meanFL1;        //!< mean infiltration into layer 1
  double meanFL2;        //!< mean infiltration into layer 2

  double sdFL1 = 0.0;  // standard deviation of landscape infiltration
  double sdFL2 = 0.0;

  double meanFtot;       //!< mean infiltration into both layers
  double meanDrainL1L2;  //!< mean drainage from layer 1 into layer 2
  double meanDrainDeep;  //!< mean deep drainage (lost from layer 2
  double meanQD;         //!< mean runoff per cell
  double sdQD = 0.0;     // standard deviation of landscape runoff

  double QDlost;    //!< total runoff losses at border cells
  double meanEPL0;  //!< mean evapotranspiration from surface
  double meanEPL1;  //!< mean evapotranspiration from layer 1
  double meanEPL2;  //!< mean evapotranspiration from layer 2
  double meanEP;    //!< mean total evapotranspiration

  double sdEPL1 = 0.0;  // standard deviation of evapotranspiration
  double sdEPL2 = 0.0;
  double sdEP = 0.0;

  double meanTL1;  //!< mean transpiration(mm) of all cells in the upper layer
                   //!< /*Tong*/
  double meanTL2;  //!< mean transpiration(mm) of all cells in the lower layer
                   //!< /*Tong*/
  double meanEL1;  //!< mean evaporation(mm) of all cells in the upper layer
                   //!< /*Tong*/

  double meanLayerDiff;

  // separate means for vegetation and crust cells
  double meanWaterL0_abs_crust1, meanWaterL0_abs_crust2, meanWaterL0_abs_crust3,
      meanWaterL0_abs_crust0, meanWaterL0_abs_veg = 0;
  double meanWaterL1_rel_crust1, meanWaterL1_rel_crust2, meanWaterL1_rel_crust3,
      meanWaterL1_rel_crust0, meanWaterL1_rel_veg = 0;
  double meanWaterL2_rel_crust1, meanWaterL2_rel_crust2, meanWaterL2_rel_crust3,
      meanWaterL2_rel_crust0, meanWaterL2_rel_veg = 0;
  double meanFL1_crust1, meanFL1_crust2, meanFL1_crust3, meanFL1_crust0,
      meanFL1_veg = 0;
  double meanFL2_crust1, meanFL2_crust2, meanFL2_crust3, meanFL2_crust0,
      meanFL2_veg = 0;
  double meanDrainL1L2_crust1, meanDrainL1L2_crust2, meanDrainL1L2_crust3,
      meanDrainL1L2_crust0, meanDrainL1L2_veg = 0;
  double meanDrainDeep_crust1, meanDrainDeep_crust2, meanDrainDeep_crust3,
      meanDrainDeep_crust0, meanDrainDeep_veg = 0;
  double meanQD_crust1, meanQD_crust2, meanQD_crust3, meanQD_crust0,
      meanQD_veg = 0;
  double meanEPL1_crust1, meanEPL1_crust2, meanEPL1_crust3, meanEPL1_crust0,
      meanEPL1_veg = 0;

  // daily means

  double dailyMeanWaterL0_abs;  //!< daily mean of surface water height [mm]
                                //!< (average)
  double dailyMeanWaterL1_rel;  //!< daily mean of soil moisture in upper layer
                                //!< [vol%] (average)
  double dailyMeanWaterL2_rel;  //!< daily mean of soil moisture in lower layer
                                //!< [vol%] (average)
  double dailyMeanFL1;          //!< daily mean infiltration into layer 1 (sum)
  double dailyMeanFL2;          //!< daily mean infiltration into layer 2 (sum)
  double dailyMeanFtot;  //!< daily mean infiltration into both layers (sum)
  double dailyMeanDrainL1L2;  //!< daily mean drainage from layer 1 into layer 2
                              //!< (sum)
  double dailyMeanDrainDeep;  //!< daily mean deep drainage (lost from layer 2)
                              //!< (sum)
  double dailyMeanQD;         //!< daily mean runoff per cell (sum)
  double dailyQDlost;  //!< total daily runoff losses at border cells	(sum)
  double dailyMeanLayerDiff;

  // daily mean standard deviations

  double dailySD_WaterL1_rel, dailySD_WaterL2_rel, dailySD_FL1, dailySD_FL2,
      dailySD_QD = 0.0;

  Bool2D calculated;  //!< indicates if the water of this cell has already been
                      //!< distributed
  bool surfaceWaterFlag;  //!< flag to indicate whether there is a cell with
                          //!< surface water
  Int2D numberOfFlows;

  int PFTs;
  int Perennials; /*Tong*/
  int Shrubs;     /*Tong*/
  int Annuals;    /*Tong*/
  int xsize;
  int ysize;
  int vegTimeStep;
  int cellsize;
  Double1D meanPftCover;  //!< mean cover of all cells for each PFT

  // Fields
  CellList* list1;
  CellList* list2;

  // Private member functions
  /***************************************************************************************
   * \brief Set elevation of each cell.
   * This routine is called at the beginning of the simulation to initialise the
   *landscape. \param p pointer to parameter object;
   *******************************************************************************************/
  void setElevation(Parameter* p);

  /***************************************************************************************
   * \brief Set output Mask of each cell.
   * This routine is called at the beginning of the simulation to initialise the
   *landscape.
   *******************************************************************************************/
  void setOutputCells();

  /***************************************************************************************
   * \brief Set aspect of each cell.
   * This routine is called at the beginning of the simulation to initialise the
   *landscape.
   *******************************************************************************************/
  void setAspectAndInclination();

  /***************************************************************************************
   * \brief Write the number of cells that contribute runoff to each cell in a
   *file. This file is used to calculate the initial vegetation distribution:
   * water availability is assessed via amount of runoff.
   * \todo Btodo: shift this to output
   *******************************************************************************************/
  void writeRunonCells();

  /***************************************************************************************
   * \brief Calculate the number of cells that produce runon for a specific
   *cell. This routine is only called at the beginning of the simulation.
   *******************************************************************************************/
  void calculateRunonCells();

  /***************************************************************************************
   * \brief Follow the runoff path.
   * This is a recursive method, it stops, if a cell is reached, which does not
   *distribute runoff to other cells. This routine is only called at the
   *beginning of the simulation to find out the order of the cells; \param i
   *current position x; \param j current position y; \param addrunon number of
   *cells that provide runon already; \param first is this the first cell in the
   *runoff-row?;
   *******************************************************************************************/
  void followRunoff(int i, int j, int addrunon, bool first);

  /***************************************************************************************
   * \brief Follow the actual runoff to the next "knot" and fill the list in the
   *right order this procedure is called at the beginning of the simulation to
   *find out the correct order of the cells \param i current position x; \param
   *j current position y;
   *******************************************************************************************/
  void followShortRunoffPath(int i, int j);

  /***************************************************************************************
   * \brief Distribution of actual runoff to each cell
   * in this routine, the order of the cells in the list is followed;
   * a cell does not get runoff until all cells contributing runoff to this cell
   *have been evaluated
   \param bl Biocrustlandscape object to get current surface water from
   *******************************************************************************************/
  void calculateRunoff(BiocrustLandscape* bl);

  /***************************************************************************************
   * \brief add the actual precipitation to the surface water
   * \param rain current amount of precipitation [mm];
   *******************************************************************************************/
  void precipitation(double rain);

  /***************************************************************************************
   * \brief Write results in result file.
   * \param year current year;
   * \param day current day;
   * \param hour current hour;
   * \param weather pointer to weather object;
   *******************************************************************************************/
  void writeOutputFile(int year, int day, int hour, Weather* weather);

  /***************************************************************************************
   * \brief Read soil parameters from file to initialise soil grid.
   * Initialise soil grid: elevation, inclination, runoff-path...
   * \param p pointer to parameter object;
   *******************************************************************************************/
  void initialiseSoilParameters(Parameter* p);

  /***************************************************************************************
   * \brief Initialise result file for soil hydrology time series.
   * \param p pointer to parameter object;
   *******************************************************************************************/
  void initialiseResultFile(Parameter* p);

  /***************************************************************************************
   * \brief Initialise result file for soil hydrology time series with variables
   *that Selina needs
    \param p pointer to parameter object;
   *******************************************************************************************/
  void initialiseResultFileSelina(Parameter* p);
  void initialiseResultFileSelinaHour(Parameter* p);
  void initialiseResultFileSelinaHour_CS(Parameter* p);
  void writeOutputFileSelina(int year, int day, int hour, Weather* weather);
  /***************************************************************************************
   * \brief Run routine to calculate the potential evapotranspiration
   * according to Hargreaves (1974).
   * \param year current year;
   * \param day current day;
   * \param w pointer to weather object;
   *******************************************************************************************/
  void potEvapotranspiration(int year, int day, Weather* w);

  /***************************************************************************************
   * \brief Calculate the spatial mean values for the calculated processes
   *infiltration, drainage, evapotranspiration and runoff.
   *******************************************************************************************/
  void meanProcesses();

  /***************************************************************************************
   * \brief Calculate the spatial mean values for the calculated processes
   separate for crust and non-crust cells
   *******************************************************************************************/
  void meanProcesses_CS(BiocrustLandscape* bl);

  /***************************************************************************************
   * \brief Calculate the mean soil moisture of the grid.
   *******************************************************************************************/
  void meanMoisture();

  /***************************************************************************************
   * \brief Sum up moisture values for each cell to output monthly and daily
   *means
   *******************************************************************************************/
  void calculateCellMeans();

 public:
  // Constructor and Destructor
  /***************************************************************************************
   * \brief Constructor
   * \param p pointer to parameter object;
   * \param vegLandscape pointer to VegetationLandscape object;
   *******************************************************************************************/
  WaterLandscape(Parameter* p, VegetationLandscape* vegLandscape);

  /***************************************************************************************
   * \brief Destructor
   *******************************************************************************************/
  ~WaterLandscape();

  // Public parameters
  Double2D runon;      //!< runon of each cell --> simpler to have it here
  Double1D meanPftTL1; /*Tong*/
  Double1D meanPftTL2; /*Tong*/

  // Fields
  vector<vector<WaterCell*>>
      soilgrid;  //!< grid of water cells --> these are objects

  // vectors for monthly and daily output per cell
  vector<vector<double>> monthly_soil_moistureL1;
  vector<vector<double>> monthly_soil_moistureL2;
  vector<vector<double>> monthly_QD_cell;
  vector<vector<double>> monthly_infiltrationL1_cell;
  vector<vector<double>> monthly_runon_cell;
  vector<vector<double>> monthly_deepdrain_cell;
  vector<vector<double>> monthly_evapotranspiration_L1_cell;
  vector<vector<double>> monthly_evapotranspiration_L2_cell;
  vector<vector<double>> monthly_evapotranspiration_tot_cell;
  vector<vector<double>> daily_soil_moistureL1;
  vector<vector<double>> daily_soil_moistureL2;
  vector<vector<double>> daily_QD_cell;
  vector<vector<double>> daily_infiltrationL1_cell;
  vector<vector<double>> daily_runon_cell;
  vector<vector<double>> daily_deepdrain_cell;
  vector<vector<double>> daily_evapotranspiration_L1_cell;
  vector<vector<double>> daily_evapotranspiration_L2_cell;
  vector<vector<double>> daily_evapotranspiration_tot_cell;

  // Public member functions
  /***************************************************************************************
   * \brief Calculate soil moisture processes,
   * hourly: infiltration, runoff, daily: evaporation.
   * \param year current year;
   * \param day current day;
   * \param w pointer to weather object;
   *******************************************************************************************/
  void calculateProcesses(int year, int day, Weather* w, BiocrustLandscape* bl);

  /***************************************************************************************
   * \brief Mainly to close files, but also to write some last results.
   *******************************************************************************************/
  void finishSimulation();

  /***************************************************************************************
   * \brief Routine that is called from outside at a specific time step
   * the current moisture values for each grid cell are exported.
   * to be able to look at patterns
   * \param year current year;
   * \param p pointer to parameter object;
   *******************************************************************************************/
  void exportCurrentMoisture(int year, Parameter* p);

  /***************************************************************************************
   * \brief Routine that is called from outside at a specific time step:
   * the current process rates / values for each grid cell are exported
   * to be able to look at patterns.
   * \param year current year;
   * \param p pointer to parameter object;
   *******************************************************************************************/
  void exportCurrentProcesses(int year, Parameter* p);
};

#endif  // WATERLANDSCAPE_H

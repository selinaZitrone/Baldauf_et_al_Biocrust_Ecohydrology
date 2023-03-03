/****************************************************************************************
 * Controller.cpp
 * \author  Britta Tietjen
 * \brief   This file contains the main routine to start the simulation.
 *          Objects are constructed, the communication between objects is
 *          handled, here is also the main time loop.
 *******************************************************************************************/
/*------------------------------------------------------------------------------------------
    INSTRUCTIONS FOR CODING IN ECOHYD
    ---------------------------------
    (1) Please be aware of keeping the code clear! Therefore, check whether your
  changes could be moved to a new class or at least to a new member function!

    (2) Stick to the coding style used in this program!

    (3) Please make sure that you comment everything (explain new parameters
  etc.) using the doxygen commands. Follow the documentation style in EcoHyD.
        http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html
        http://www.stack.nl/~dimitri/doxygen/manual/commands.html

    (4) Please avoid comments like "changed by ..." and don't comment out that
  much! Instead use the \author-command for new functions or classes you create,
  and upload your versions in a repository to allow for tracking changes you
  have made.

    (5) If you find a bug or something else which might be of everyone's
  interest, please contact all those who are currently working on EcoHyD!
  ------------------------------------------------------------------------------------------*/

#include <ctime>

#include "BiocrustLandscape.h"
#include "CrustParameters.h"
#include "Parameters.h"
#include "RandomNumberGenerator.h"
#include "VegetationLandscape.h"
#include "WaterLandscape.h"
#include "Weather.h"

using namespace std;

WaterLandscape* waterLandscape;
VegetationLandscape* vegLandscape;
Weather* weather;
Parameter* parameter;
BiocrustLandscape* crustLandscape;
CrustParameters* crust_parameter;

/***************************************************************************************
 * \brief Delete the objects.
 * is called in --> controller::runSimulation
 *******************************************************************************************/
void deleteObjects() {
  delete weather;
  weather = NULL;
  delete waterLandscape;
  waterLandscape = NULL;
  delete vegLandscape;
  vegLandscape = NULL;
  delete crustLandscape;
  crustLandscape = NULL;
}

/***************************************************************************************
 * \brief Start and control the simulation.
 * is called in --> controller::main
 *******************************************************************************************/
void runSimulation() {
  // just a time stamp for the start of the simulation
  time_t t = time(0);  // get time now
  struct tm* now = localtime(&t);

  // variable to calculate duration of a simulation
  time_t t2;
  double sim_duration;

  cout << "Simulation for Sc. " << parameter->scenario + 1 << ", Climate Rep. "
       << parameter->climateRepetition + 1 << " and Model Rep. "
       << parameter->modelRepetition + 1 << " started at: " << now->tm_hour
       << ":" << now->tm_min << ":" << now->tm_sec << endl;

  ///----------------------------YEARLY
  /// SIMULATION-------------------------------------------
  for (int year = 0; year < parameter->simYears; year++) {
    cout << "Year:"
         << "\t" << year << "\n";
    vegLandscape->wetStart = weather->growStart[year];
    vegLandscape->wetEnd = weather->growEnd[year] + parameter->vegTimeStep;

    ///----------------------------DAILY
    /// SIMULATION-------------------------------------------
    for (int day = 0; day < daysPerYear; day++) {
      std::cout << day << ", ";
      //---------------------------------------------
      /// UPDATE CALCULATED WATER DATA TO VEGETATION MODEL,
      /// which is necessary to calculate vegetation processes.
      //---------------------------------------------
      for (int x = 0; x < parameter->xsize; x++) {
        for (int y = 0; y < parameter->ysize; y++) {
          /// \implements \fn VegetationCell::setMoist(double mL1, double mL2,
          /// int day,int grStart, int grEnd, int wetStart, int wetEnd)
          vegLandscape->veggrid[x][y]->setMoist(
              waterLandscape->soilgrid[x][y]->waterL1_rel,
              waterLandscape->soilgrid[x][y]->waterL2_rel, day,
              vegLandscape->growStart, vegLandscape->growEnd,
              vegLandscape->wetStart, vegLandscape->wetEnd);

          /// To show annual evaporation and transpiration in the vegetation
          /// result file \implements \fn
          /// VegetationCell::setEvaporationandTranspiration(double Evaporation,
          /// double TranspirationL1, double TranspirationL2, int day)
          vegLandscape->veggrid[x][y]->setEvaporationandTranspiration(
              waterLandscape->soilgrid[x][y]->EL1,
              waterLandscape->soilgrid[x][y]->TL1,
              waterLandscape->soilgrid[x][y]->TL2, day);

          /// To transmit the calculated transpiration of each PFT in water
          /// process to vegetation cell. \implements \fn
          /// VegetationCell::setTransFactors(double TransL1[nPFTs], double
          /// TransL2[nPFTs], int day)
          vegLandscape->veggrid[x][y]->setTransFactors(
              waterLandscape->soilgrid[x][y]->PftTL1,
              waterLandscape->soilgrid[x][y]->PftTL2, day);
        }
      }

      //---------------------------------------------
      /// CALCULATE VEGETATION PROCESSES.
      /// All concerning vegetation changes are controlled here.
      /// \implements \fn VegetationLandscape::calculateProcesses(int year, int
      /// day, WaterLandscape* waterlandscape, Weather* we, Parameter* p)
      //---------------------------------------------
      vegLandscape->calculateProcesses(year, day, waterLandscape, weather,
                                       parameter);

      //---------------------------------------------
      /// UPDATE CALCULATED VEGETATION DATA TO WATER MODEL,
      /// which is necessary to calculate water processes.
      //---------------------------------------------
      for (int x = 0; x < parameter->xsize; x++) {
        for (int y = 0; y < parameter->ysize; y++) {
          /// \implements \fn VegetationCell::setRootL2()
          vegLandscape->veggrid[x][y]->setRootL2();
          /// \implements \fn VegetationCell::setRootL1()
          vegLandscape->veggrid[x][y]->setRootL1();
          /// \implements \fn VegetationCell::setCover()
          vegLandscape->veggrid[x][y]->setCover();
          /// \implements \fn VegetationCell::setpftUse()
          vegLandscape->veggrid[x][y]->setpftUse();
          /// \implements \fn WaterCell::setRootL2(double RootL2_[])
          waterLandscape->soilgrid[x][y]->setRootL2(
              vegLandscape->veggrid[x][y]->RootL2);
          /// \implements \fn WaterCell::setRootL1(double RootL1_[])
          waterLandscape->soilgrid[x][y]->setRootL1(
              vegLandscape->veggrid[x][y]->RootL1);
          /// \implements \fn WaterCell::setCover(double shrubcover[], double
          /// grasscover[], double annualcover[])
          waterLandscape->soilgrid[x][y]->setCover(
              vegLandscape->veggrid[x][y]->Scover,
              vegLandscape->veggrid[x][y]->Pcover,
              vegLandscape->veggrid[x][y]->Acover);
          /// \implements \fn WaterCell::totalRoots()
          waterLandscape->soilgrid[x][y]->totalRoots();
          /// \implements \fn WaterCell::totalCover()
          waterLandscape->soilgrid[x][y]->totalCover();
          /// \implements \fn WaterCell::setVegUseFactors(double PftUseL1[],
          /// double PftUseL2[])
          waterLandscape->soilgrid[x][y]->setVegUseFactors(
              vegLandscape->veggrid[x][y]->UseL1,
              vegLandscape->veggrid[x][y]->UseL2);
        }
      }

      //---------------------------------------------
      /// CALCULATE WATER PROCESSES
      /// Calculate soil moisture processes.
      /// \implements \fn WaterLandscape::calculateProcesses(int year, int day,
      /// Weather* w)
      //---------------------------------------------
      waterLandscape->calculateProcesses(year, day, weather, crustLandscape);

      //   if ((year == 99) && (day == 263)) {
      //     cout << "spatial analysis\n";
      //     waterLandscape->exportCurrentMoisture(year, parameter);
      //     // vegLandscape->exportActualVegetation(year,parameter); // later
      //     when
      //     // veg file is adjusted
      //   }
    }
    ///----------------------------END OF DAILY
    /// SIMULATION--------------------------------
  }
  ///----------------------------END OF YEARLY
  /// SIMULATION--------------------------------

  //---------------------------------------------
  /// FINISH SIMULATION (WATER).
  /// Close files etc.
  /// \implements \fn WaterLandscape::finishSimulation()
  //---------------------------------------------
  waterLandscape->finishSimulation();

  //---------------------------------------------
  /// FINISH SIMULATION (VEGETATION).
  /// Close files etc.
  /// \implements \fn VegetationLandscape::finishSimulation()
  //---------------------------------------------
  vegLandscape->finishSimulation();

  //---------------------------------------------
  /// FINISH SIMULATION (VEGETATION).
  /// Close files etc.
  /// \implements \fn BiocrustLandscape::finishSimulation()
  //---------------------------------------------
  crustLandscape->finishSimulation();

  //---------------------------------------------
  /// DELETE OBJECTS.
  /// Deletes all objects.
  /// \implements \fn  deleteObjects()
  //---------------------------------------------
  deleteObjects();

  // calculate duration of simulation
  t2 = time(0);
  sim_duration = difftime(t2, t);
  t = time(0);  // get time now
  now = localtime(&t);
  cout << "Simulation for Sc. " << parameter->scenario + 1 << ", Climate Rep. "
       << parameter->climateRepetition + 1 << " and Model Rep. "
       << parameter->modelRepetition + 1 << " ended at: " << now->tm_hour << ":"
       << now->tm_min << ":" << now->tm_sec << endl
       << endl;

  cout << "Duration of simulation: " << sim_duration << " seconds " << endl;
}

/***************************************************************************************
 * \brief Create the objects.
 * New landscapes are generated for the different modules.
 * Additionally a weather object which contains the weather data
 * and a parameter object which contains general parameters which
 * should be available in other files as well.
 * is called in --> controller::main
 *******************************************************************************************/
void initialiseSimulation() {
  //---------------------------------------------
  /// NEW WEATHER.
  /// Constructor for weather.
  /// \implements \fn Weather::Weather(Parameter* p)
  //---------------------------------------------
  weather = new Weather(parameter);

  //---------------------------------------------
  /// NEW VEGLANDSCAPE.
  /// Constructor for landscape.
  /// \implements \fn VegetationLandscape::VegetationLandscape(Parameter* p,
  /// WaterLandscape* waterLandscape)
  //--------------------------------------------*
  vegLandscape = new VegetationLandscape(parameter, waterLandscape);

  //---------------------------------------------
  /// NEW WATERLANDSCAPE.
  /// Constructor for landscape.
  /// \implements \fn WaterLandscape::WaterLandscape(Parameter* p, VegLandscape*
  /// vegLandscape)
  //---------------------------------------------
  waterLandscape = new WaterLandscape(parameter, vegLandscape);

  //---------------------------------------------
  /// NEW BIOCRUST LANDSCAPE
  /// Constructor for biocrust landscape.
  /// \implements \fn BiocrustLandscape::BiocrustLandscape()
  //---------------------------------------------
  crustLandscape = new BiocrustLandscape(parameter, crust_parameter);

  //---------------------------------------------
  /// GET WEATHER.
  /// Read parameters from weather file.
  /// \implements \fn Weather::getWeather(Parameter* p)
  //---------------------------------------------
  weather->getWeather(parameter);

  //---------------------------------------------
  /// MEAN WEATHER.
  /// Calculate daily mean temperature, radiation etc.
  /// \implements \fn Weather::meanWeather(Parameter* p)
  //---------------------------------------------
  weather->meanWeather(parameter);

  //---------------------------------------------
  /// GROWING SEASON.
  /// Find first and last day of water season.
  /// Growing season read from "vegetationparameters.txt".
  /// \implements \fn Weather::growingSeason()
  //---------------------------------------------
  weather->growingSeason(parameter);
}

/***************************************************************************************
 * \brief Main model routine.
 * The routine runSimulation() starts the simulation itself.
 *******************************************************************************************/
int main() {
  //---------------------------------------------
  /// NEW PARAMETER.
  /// Constructor for parameter.
  /// \implements \fn Parameter::Parameter()
  //---------------------------------------------
  parameter = new Parameter();

  //---------------------------------------------
  /// NEW CRUST PARAMETER.
  /// Constructor for crust parameter.
  /// \implements \fn CrustParameters::CrustParameters()
  //---------------------------------------------
  crust_parameter = new CrustParameters();

  //---------------------------------------------
  /// READ IN MODEL SCENARIOS.
  /// Read in scenarios for climate, soil and vegetation
  /// and number of repetitions from file.
  /// \implements \fn Parameter::readInModelScenarios()
  //---------------------------------------------
  parameter->readInModelScenarios();

  //---------------------------------------------
  /// INITIALIZE RANDOM NUMBER GENERATOR
  /// \implements \fn RandomNumberGenerator::set_random_seed(seed)
  //---------------------------------------------
  // \todo: Add random seed as parameter or command line option, don't used
  // hard-coded seed
  RandomNumberGenerator::set_random_seed(1234);

  for (parameter->scenario = 0; parameter->scenario < parameter->scenarios;
       parameter->scenario++) {
    for (parameter->climateRepetition = 0;
         parameter->climateRepetition < parameter->climateRepetitions;
         parameter->climateRepetition++) {
      for (parameter->modelRepetition = 0;
           parameter->modelRepetition < parameter->modelRepetitions;
           parameter->modelRepetition++) {
        cout << "____ Parameters for Sc. " << parameter->scenario + 1
             << ", Climate Rep. " << parameter->climateRepetition + 1
             << " and Model Rep. " << parameter->modelRepetition + 1 << " ____"
             << endl
             << endl;

        //---------------------------------------------
        /// READ IN GENERAL MODEL PARAMETERS.
        /// Read general model parameters.
        /// \implements \fn Parameter::readInModelParameters()
        //---------------------------------------------
        parameter->readInModelParameters();

        //---------------------------------------------
        /// INITIALISE SIMULATION.
        /// Create the objects.
        /// Initialise fields, parameters etc.
        /// \implements \fn  initialiseSimulation()
        //---------------------------------------------
        initialiseSimulation();

        //---------------------------------------------
        /// RUN SIMULATION.
        /// Start, control and end the simulation.
        /// \implements \fn  runSimulation()
        //---------------------------------------------
        runSimulation();
      }
    }
  }
  cout << "_______________________________________________________________"
       << endl;
  cout << ">>> Model Run(s) for " << parameter->scenarios << " Scenario(s), "
       << parameter->climateRepetitions << " Climate Rep(s) and "
       << parameter->modelRepetitions << " Model Rep(s) successfully completed!"
       << endl;

  delete parameter;
  parameter = NULL;

  return 0;
}

#pragma once
#include <map>
#include <string>
#include <vector>

using namespace std;

class CrustParameters {
 private:
  /*************************************************************************************
   \author Selina Baldauf
   \date October 2021
   \brief Read one crust parameter input file into a map. Is called in
   CrustParameters->readAllCrustParameters
   \param filename Filename of crust parameters to read
  ****************************************************************************************/
  map<string, double> readSingleCrustParameterFile(const string& filename);

  /*************************************************************************************
    \author Selina Baldauf
    \date October 2021
    \brief Read all Crust Parameters from crust_default.yml. Is called in
    CrustParameters constructor
    \todo At the moment only one set of parameters is read in from fixed file ->
  should become more flexible to also read multiple crust parameter files
  ****************************************************************************************/
  void readAllCrustParameters();

  /*************************************************************************************
   \author Selina Baldauf
   \date October 2021
   \brief Calculate pore disconnectedness from pore distribution index b. Is
   called in CrustParameters->readAllCrustParameters.
   \param b pore distribution index; crust input parameter (-)
 ****************************************************************************************/
  double calculatePoreDisconnectedness(double b);

  /*************************************************************************************
    \author Selina Baldauf
    \date October 2021
    \brief Calculate pore size distribution from pore disconnectedness. Is
   called in CrustParameters->readAllCrustParameters.
    \param c pore disconnectedness (-)-> can be calculated using the function
    calculatePoreDisconnectedness
  ****************************************************************************************/
  double calculatePoreSizeDistribution(double c);

  /*************************************************************************************
    \author Selina Baldauf
    \date October 2021
    \brief Helper function to print crust parameters to the console (to check if
     they are correct). Can be called in CrustParameters constructor
    \param par_maps vector of parameter maps to be printed sequentially
  ****************************************************************************************/
  void printCrustParameters(vector<map<string, double>>& par_maps);

  /*************************************************************************************
    \author Selina Baldauf
    \date January 2022
    \brief Function to check if all crust parameters are in the parameter file
      or if there are typos or wrong numers. Is called in
      CrustParameters->readSingleCrustParameter
    \param par_map parameter map to check
  ****************************************************************************************/
  void checkSingleCrustParameters(map<string, double>& par_map);

 public:
  vector<map<string, double>> crust_params;

  // constructor: Reads all crust parameters
  CrustParameters();
  // destructor (default)
  ~CrustParameters();
};
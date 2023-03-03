#include "CrustParameters.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "Parameters.h"

/*************************************************************************************
Read one crust parameter input file into a map. Is called in
CrustParameters->readAllCrustParameters
****************************************************************************************/
map<string, double> CrustParameters::readSingleCrustParameterFile(
    const string& filename) {
  string line, name;
  double value;
  map<string, double> parameters;

  // open the file
  ifstream parameterFile;
  cout << "Read crust parameters from " << filename << endl;

  parameterFile.open(filename.c_str(), ios::binary | ios::in);
  // if something goes wrong...
  if (!parameterFile) {
    cerr << filename << " could not be opened!\n";
    exit(-1);
  }

  // loop through each parameter in the file

  // now loop through the file for each parameter
  while (getline(parameterFile, line)) {
    if (line.rfind('#', 0) == 0 || line.empty()) {
      // header or comment line, don't parse this line
    } else {
      // read individual parameters
      stringstream curr_par(line);
      // read parameter name
      curr_par >> name;
      // delete last ':' character< of the par_name string
      name.pop_back();
      // read parameter value
      curr_par >> value;
      // add parameter to map
      parameters.emplace(name, value);
    }
  }

  parameterFile.close();

  // check map for completeness:
  checkSingleCrustParameters(parameters);

  return parameters;
}
/*************************************************************************************
\brief Function to check if all crust parameters are in the parameter file
       or if there are typos or wrong numers. Is called in
CrustParameters->readSingleCrustParameter
****************************************************************************************/
void CrustParameters::checkSingleCrustParameters(map<string, double>& par_map) {
  vector<string> to_find{"Ks",  "nc",     "Z",         "b",        "sfc",
                         "shc", "map_id", "EP_factor", "manning_n"};
  for (auto& elem : to_find) {
    if (par_map.find(elem) == par_map.end()) {
      cerr << "Could not find parameter " << elem << " in crust parameter map"
           << endl;
      exit(-1);
    }
  }
}

/*************************************************************************************
Read all Crust Parameters from crust_default.yml. Is called in CrustParameters
constructor
****************************************************************************************/
void CrustParameters::readAllCrustParameters() {
  // loop over all files with crust_ in the parameters folder and read them in
  string dir_path = "Parameters" + delimiter + "crust";
  string file_path = "";

  // initialize a map for crust type 0 (no crust)
  map<string, double> no_crust{{"map_id", 0}};
  crust_params.push_back(no_crust);

  for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
    cout << entry.path().generic_string() << endl;
    crust_params.push_back(
        readSingleCrustParameterFile(entry.path().generic_string()));
  }
  // // add the pore disconnectedness and pore size distribution
  double m, c;
  for (size_t i = 0; i < crust_params.size(); i++) {
    if ((int)crust_params[i].find("map_id")->second != 0) {
      c = calculatePoreSizeDistribution(crust_params[i]["b"]);
      m = calculatePoreDisconnectedness(c);
      crust_params[i].emplace("c", c);
      crust_params[i].emplace("m", m);
    }
  }
}

/*************************************************************************************
Calculate pore disconnectedness from pore distribution index b. Is called in
CrustParameters->readAllCrustParameters.
****************************************************************************************/
double CrustParameters::calculatePoreDisconnectedness(double b) {
  return 2 * b + 3;
}

/*************************************************************************************
Calculate pore size distribution from pore disconnectedness. Is called in
CrustParameters->readAllCrustParameters.
****************************************************************************************/
double CrustParameters::calculatePoreSizeDistribution(double c) {
  return 2 / (c - 3);
}

/*************************************************************************************
Helper function to print crust parameters to the console (to check if they are
correct). Can be called in CrustParameters constructor
****************************************************************************************/
void CrustParameters::printCrustParameters(
    vector<map<string, double>>& par_maps) {
  for (size_t i = 0; i < par_maps.size(); i++) {
    cout << "Crust parameters for crust " << i << endl;
    for (auto& x : par_maps[i]) {
      cout << x.first << ": " << x.second << endl;
    }
  }
}

/*************************************************************************************
Constructor
****************************************************************************************/
CrustParameters::CrustParameters() {
  readAllCrustParameters();
  printCrustParameters(crust_params);
  // add an empty map to the beginning of the crust_params to represent a crust
  // cell without crust cover
  map<string, double> empty_map{{"id", 0.0}};
  crust_params.insert(crust_params.begin(), empty_map);
}

/*************************************************************************************
Destructor
****************************************************************************************/
CrustParameters::~CrustParameters() = default;

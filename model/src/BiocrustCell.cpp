#include "BiocrustCell.h"

#include <math.h>

#include <algorithm>
#include <iostream>

#include "CrustParameters.h"
#include "Weather.h"

/*************************************************************************************
Calculate crust infiltration. Is called in
BiocrustCell->calculateCellProcessHour
****************************************************************************************/
double BiocrustCell::crustInfiltration(double surface_water,
                                       double moisture_crust,
                                       bool crust_wettingFront,
                                       bool with_hydrophobicity) {
  double inf_crust = 0.0;
  // calculate Kgeom from unsaturated hydraulic conductivity Ku based on
  // empirical relationship
  // double Ku = crust_params_cell["Ks"] / (pow(1.0 + ((-30.0) /
  // (-26.0)), 3.0)); double Kgeom = sqrt(crust_params_cell["Ks"] * Ku);

  // // calculate fillable pores
  // double Na = max(
  //     crust_params_cell["nc"] - crust_params_cell["nc"] * moisture_crust,
  //     0.0);

  // if (surface_water < Kgeom) {
  //   // If there is less surface water than can infiltrate in this time step:
  //   // everything can infiltrate
  //   inf_crust = surface_water;
  //   crust_wettingFront = false;
  // } else {
  //   // if there is a wettingfront, time until ponding is 0
  //   // nothing can infiltrate
  //   double Fs = 0;
  //   double ts = 0;  // time until ponding
  //   if (!crust_wettingFront) {
  //     // if there is no wetting front yet, calculate time until ponding
  //     ts = ((Na * crust_params_cell["Psi"]) / (surface_water / Kgeom - 1) /
  //           surface_water) *
  //          veg_function;
  //     Fs = min(1.0, ts) * surface_water;
  //   }
  //   // check if ponding occurred in this timestep (1 h)
  //   if (ts >= 0.0 && ts < 1.0) {
  //     // ponding occurred in this timestep
  //     double A = Kgeom * (1 - ts);
  //     double B = max(Fs + 2 * (Na * crust_params_cell["Psi"]), 0.0);
  //     double h5 = A * A / 4.0 + A * B + Fs * Fs;
  //     inf_crust = min(
  //         (A / 2.0 + (A * A + 2 * A * B) / (4.0 * sqrt(h5))) * veg_function +
  //             Fs,
  //         surface_water);
  //     crust_wettingFront = true;
  //   } else {
  //     // ponding did not occur in this timestep
  //     inf_crust = surface_water;
  //     crust_wettingFront = false;
  //   }
  // }

  // infiltration into crust is minimum of surface water and water that can
  // infiltrate into the crust

  inf_crust =
      min(surface_water, crust_params_cell["Z"] * crust_params_cell["nc"] *
                             (1 - moisture_crust));

  // If there is hydrophobicity of the crust, infiltation is reduced
  if (with_hydrophobicity) {
    inf_crust = inf_crust * calculateHydrophobicityFactor(moisture_crust);
  }
  // check saturation excess
  // make sure that only the water that actually fits into the crust can
  // infiltrate
  // inf_crust = min(inf_crust, crust_params_cell["Z"] * crust_params_cell["nc"]
  // *
  //                                (1 - moisture_crust));

  return inf_crust;
}

/*************************************************************************************
Calculate hydrophobicity factor. This factor is multiplied with infiltration
to reduce infiltration if crust is dry. Is called in
BiocrustCell->crustInfiltration
****************************************************************************************/
double BiocrustCell::calculateHydrophobicityFactor(double moisture_crust) {
  double hyd_fct = 1.0;
  if (moisture_crust < crust_params_cell["theta_crit"]) {
    // cout << "moisture smaller theta: " << moisture_crust
    //  << " theta: " << crust_params_cell["theta_crit"];
    hyd_fct =
        crust_params_cell["h_min"] * exp((-1 / crust_params_cell["theta_crit"] *
                                          log(crust_params_cell["h_min"])) *
                                         moisture_crust);
  }
  return hyd_fct;
}

/*************************************************************************************
Update crust moisture with current (positive or negative) water input. Is
called in BiocrustCell->calculateCellProcessHour and
BiocrustCell->calculateCellProcessesDay
****************************************************************************************/
double BiocrustCell::updateCrustMoisture(double current_moisture,
                                         double water_input) {
  double new_moist = current_moisture + water_input / (crust_params_cell["nc"] *
                                                       crust_params_cell["Z"]);
  return min(max(new_moist, 0.0), 1.0);
}

/*************************************************************************************
Calculate leakage of water from biocrust to soil. Is called in
BiocrustCell->calculateCellProcessHour
****************************************************************************************/
double BiocrustCell::crustLeakageIntoSoil(double current_moisture,
                                          double surface_water) {
  double beta = 2 * crust_params_cell["b"] * 4;
  double leak =
      crust_params_cell["Ks"] *
      ((exp(beta * (current_moisture - crust_params_cell["sfc"])) - 1) /
       (exp(beta * (1 - crust_params_cell["sfc"])) - 1));
  // if more water could leak than is available only return available water
  return max(min(leak, surface_water), 0.0);
}

/*************************************************************************************
Calculate evaporation from the crust layer. Is called in
BiocrustCell->calculateCellProcessDay
****************************************************************************************/
double BiocrustCell::crustEvaporation(double current_moisture, double potEP,
                                      double openDemand) {
  // no evaporation from crust if if moisture is below hygroscopic point
  if (current_moisture < crust_params_cell["shc"]) {
    return 0.0;
  } else {
    // potential biocrust evaporation
    double crustEvap = crust_params_cell["kc"] * potEP;
    // biocrust evaporation can only be as large as demand
    crustEvap = min(crustEvap, openDemand);
    // biocrust evaporation cannot exceed available water content in the crust
    crustEvap =
        min(crustEvap, (current_moisture - crust_params_cell["shc"]) *
                           crust_params_cell["nc"] * crust_params_cell["Z"]);
    return crustEvap;
  }
}

/*************************************************************************************
Calculate evaporation from the ponded water. All surface water can freely
evaporate if it is less than the evaporative demand. Evaporative demand is
calculated as potential evaporation. Is called in
BiocrustCell->calculateCellProcessDay
****************************************************************************************/
double BiocrustCell::pondedEvaporation(double openDemand,
                                       double surface_water) {
  return min(surface_water, openDemand);
}

/*************************************************************************************
Calculate all biocrust cell processes that run on hourly timescale. Is called
in BiocrustLandscape->calculateProcessHour
****************************************************************************************/
void BiocrustCell::calculateCellProcessesHour() {
  // infiltration
  infiltration_crust = crustInfiltration(
      surface_water, moisture_crust, crust_wettingFront, with_hydrophobicity);
  // update moisture
  moisture_crust = updateCrustMoisture(moisture_crust, infiltration_crust);
  // update surface water
  surface_water -= infiltration_crust;
  // leakage into the soil
  leakage_crust = crustLeakageIntoSoil(moisture_crust, surface_water);
  // update surface water
  surface_water -= leakage_crust;
  // evaporation crust
}

/*************************************************************************************
Calculate all biocrust cell processes that run on daily timescale. Is called
in BiocrustLandscape->calculateProcessDay
****************************************************************************************/
void BiocrustCell::calculateCellProcessesDay(double potEP) {
  double openDemand = potEP;
  ponded_evaporation = pondedEvaporation(potEP, surface_water);
  surface_water -= ponded_evaporation;
  openDemand -= ponded_evaporation;
  evaporation_crust = crustEvaporation(moisture_crust, potEP, openDemand);
  // update crust moisture
  moisture_crust = updateCrustMoisture(moisture_crust, -evaporation_crust);
}

/*************************************************************************************
Constructor
****************************************************************************************/
BiocrustCell::BiocrustCell(vector<map<string, double>> crust_parameters,
                           unsigned int crust_type, double moisture_init) {
  bool found_params = false;
  // find the correct parameters based on map id
  for (auto& params : crust_parameters) {
    if ((int)params.find("map_id")->second == crust_type) {
      crust_params_cell = params;
      found_params = true;
      break;
    }
  }
  // if paramters for crust_type were not found: error
  if (!found_params) {
    cerr << "Crust parameters with map_id " << crust_type
         << " not found! Please check you crust parameter input" << endl;
    exit(-1);
  }

  moisture_crust = moisture_init;
  this->crust_type = crust_type;
  if (crust_params_cell.find("theta_crit") == crust_params_cell.end() ||
      crust_params_cell.find("h_min") == crust_params_cell.end()) {
    this->with_hydrophobicity = false;
  } else {
    this->with_hydrophobicity = true;
  }
}
/********************************************************************************************
 * WaterCell.cpp
 * This file describes the water dynamics in one grid cell.
 ********************************************************************************************/

#include "WaterCell.h"

#include <math.h>  //to use mathematical functions

#include <algorithm>  //to use min/max
#include <cstdlib>    //to use exit()
#include <fstream>    //to read and write files
#include <iostream>   //to use cin / cout
#include <string>     //to use strings

#include "Parameters.h"
#include "WaterLandscape.h"
#include "Weather.h"

using namespace std;

/*******************************************************************************************
 * calculate actual evapotranspiration in two layers
 * here, also the aspect and the inclination are taken into account
 * which normally belong to the potential evapotranspiration
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::actLayerEvapotranspiration(
    double rainTot, double potEP, bool crust,
    double ep_factor_cell) { /*Tong, based on PFT-oriented*/

  // decreasing function with increasing cover
  double vegfunctionL1 =
      Vegcoef1 - Vegcoef2 * (cGtotal + cAtotal + cStotal); /*Tong*/

  double vegfunctionL2 =
      (double)(PtotRootL2 +
               StotRootL2);  // This parameter has considered the root fraction
                             // and add the part of perennial grass /*Tong*/

  double cos_inclination = cos(inclination);

  if (potEP != nodata) {
    // aspect factor
    potEP = potEP * aspectFactor;

    ///--------------------------evapotranspiration from surface
    /// water------------------------------------------------
    // dependent on aspect factor and on vegetation cover (shading effect

    double potL0EP =
        potEP * (1 - (cStotal + cAtotal + cGtotal) / 2.0) * cos_inclination;

    // only calculate Evaporation from upper layer if there is no crust in this
    // cell
    if (!crust) {
      if (waterL0_abs > 0) {
        // if there is more surface water than there could
        // potentially evaporate some surface water is left
        if (waterL0_abs - potL0EP > 0) {
          EPL0 = potL0EP;
          waterL0_abs -= EPL0;
        }
        // otherwise all surface water evaporates
        else {
          EPL0 = waterL0_abs;
          waterL0_abs = 0;
        }
      } else {
        EPL0 = 0;
        waterL0_abs = 0;
      }
    } else {
      EPL0 = 0;
    }

    ///--------------------------evapotranspiration from soil soil
    /// layer----------------------------------------------
    // this routine is oriented according to
    // Porporato's models but also to the HBV model
    // it is flattened between 0 and wp and constant above the wp
    ///--------------------------evapotranspiration in upper soil
    /// layer-----------------------------------------------
    if (waterL1_rel > WP) {
      // constant value evaporates above wp, but a residual water content must
      // always be there
      EPL1 = min(potEP * cos_inclination * vegfunctionL1 * ep_factor_cell,
                 (waterL1_rel - resWater) * depthL1);
    } else
      EPL1 = max(0.0, min(potEP * pow((double)waterL1_rel / WP, 2) *
                              cos_inclination * vegfunctionL1 * ep_factor_cell,
                          (waterL1_rel - resWater) * depthL1));

    ///--------------------------transpiration in lower soil
    /// layer----------------------------------------------------
    if (waterL2_rel <= WP)
      EPL2 = max(0.0, min(potEP * pow((double)waterL2_rel / WP, 2) *
                              vegfunctionL2 * ep_factor_cell * cos_inclination,
                          (waterL2_rel - resWater) * depthL2));

    else
      EPL2 = min(potEP * cos_inclination * vegfunctionL2 * ep_factor_cell,
                 (waterL2_rel - resWater) * depthL2);

    // absolute moisture
    waterL1_abs = max(waterL1_abs - EPL1, 0.0);
    waterL2_abs = max(waterL2_abs - EPL2, 0.0);
    // relative moisture
    waterL1_rel = waterL1_abs / depthL1;
    waterL2_rel = waterL2_abs / depthL2;

    EPtot = EPL0 + EPL1 + EPL2;
  }  // end if (potEP != nodata)

  else {
    EPL0 = nodata;
    EPL1 = nodata;
    EPL2 = nodata;
  }
}

/*******************************************************************************************
 * transpiration and evaporation in the upper layer
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::transpirationupper() { /*Tong, separation of
                                          evapotranspiration*/

  if ((cStotal + cGtotal + cAtotal) >= 0 && EPL1 >= 0) {
    TL1 = max(0.0, Veg_ETcoef1 * (cStotal + cGtotal + cAtotal) - Veg_ETcoef2) *
          EPL1;
  } else
    TL1 = 0;
}

/*******************************************************************************************
 * transpiration in the lower layer, evaporation equals to 0
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::transpirationlower() { /*Tong, separation of
                                          evapotranspiration*/

  if ((cStotal + cGtotal + cAtotal) > 0 && EPL2 > 0) {
    TL2 = EPL2;
  } else
    TL2 = 0;
}

/*******************************************************************************************
 * evaporation in the upper layer
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::evaporationupper() { /*Tong, separation of evapotranspiration*/

  if (EPL1 > TL1) {
    EL1 = EPL1 - TL1;
  } else
    EL1 = 0;
}

/*******************************************************************************************
 * transpiration of PFTs in the upper layer
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::transpirationPFTU() { /*Tong, separation of evapotranspiration*/

  for (int i = 0; i < PFTs; i++) {
    PftTL1[i] = UseL1[i] * TL1 * coverPFT[i];
  }
}

/*******************************************************************************************
 * transpiration of PFTs in the lower layer
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::transpirationPFTL() { /*Tong, separation of evapotranspiration*/

  for (int i = 0; i < PFTs; i++) {
    PftTL2[i] = UseL2[i] * TL2 * coverPFT[i];
  }
}

/*******************************************************************************************
 * total soil water content
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::totalWater() { totWater = waterL1_abs + waterL2_abs; }

/*******************************************************************************************
 * set shrub and grass cover
 * this routine is called by the landscape
 * values for shrub cover and grass cover should be between 0 and 1
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void WaterCell::setCover(Double1D shrubcover, Double1D grasscover,
                         Double1D annualcover) { /*Tong, array transmission*/

  for (int i = 0; i < PFTs; i++) {
    if (i < Shrubs)
      coverPFT[i] = shrubcover[i];  // between 0 and 1
    else if (i >= Shrubs && i < Shrubs + Perennials)
      coverPFT[i] = grasscover[i - Shrubs];  // between 0 and 1
    else
      coverPFT[i] = annualcover[i - (Shrubs + Perennials)];  // between 0 and 1
  }
}

/*******************************************************************************************
 * calculate total cover for each PFT
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void WaterCell::totalCover() {
  cStotal = 0;
  cGtotal = 0;
  cAtotal = 0;

  for (int i = 0; i < PFTs; i++) {
    if (i < Shrubs)
      cStotal += coverPFT[i];
    else if (i >= Shrubs && i < Shrubs + Perennials)
      cGtotal += coverPFT[i];
    else
      cAtotal += coverPFT[i];
  }
}

/*******************************************************************************************
 * set vegetation use factors
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void WaterCell::setVegUseFactors(
    Double1D PftUseL1, Double1D PftUseL2) { /*Tong, array transmission*/

  for (int i = 0; i < PFTs; i++) {
    UseL1[i] = PftUseL1[i];
    UseL2[i] = PftUseL2[i];
  }
}

/*******************************************************************************************
 * the aspect factor determines how the radiation is affected by the aspect
 * of the slope
 * the values for this factor were found in
 * Shevenell L., 1999. Regional potential evapotranspiration in arid climates
 *    based on temperature, topography and calculated solar radiation.
 *    Hydrological Processes 13, 577-596
 * this aspect factor is later on multiplied to the radiation to calculate
 * the potential evapotranspiration
 * called in --> WaterLandscape::initialiseSoilParameters
 *******************************************************************************************/
void WaterCell::setAspectFactor() {
  switch (aspect) {
    case FL:
      aspectFactor = 1.0;
      break;
    case NN:
      aspectFactor = 0.90;
      break;
    case NE:
      aspectFactor = 0.95;
      break;
    case EE:
      aspectFactor = 0.98;
      break;
    case SE:
      aspectFactor = 1.03;
      break;
    case SS:
      aspectFactor = 1.10;
      break;
    case SW:
      aspectFactor = 1.05;
      break;
    case WW:
      aspectFactor = 1.02;
      break;
    case NW:
      aspectFactor = 0.97;
      break;
  }
}

/*******************************************************************************************
 * set roots layer two
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void WaterCell::setRootL2(Double1D RootL2_) { /*Tong, array transmission*/

  for (int i = 0; i < PFTs; i++) {
    if (i < Shrubs + Perennials)
      RootL2[i] = RootL2_[i];
    else
      RootL2[i] = 0;
  }
}

/*******************************************************************************************
 * set roots layer one
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void WaterCell::setRootL1(Double1D RootL1_) { /*Tong, array transmission*/

  for (int i = 0; i < PFTs; i++) {
    RootL1[i] = RootL1_[i];
  }
}

/*******************************************************************************************
 * calculate total roots in upper and lower soil layer for all PFTs per cell
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void WaterCell::totalRoots() {
  // Calculate total roots per PFT in layer one and two
  StotRootL1 = 0;
  PtotRootL1 = 0;
  AtotRootL1 = 0;
  StotRootL2 = 0;
  PtotRootL2 = 0;

  for (int i = 0; i < PFTs; i++) {
    if (i < Shrubs) {
      StotRootL1 += RootL1[i] * coverPFT[i];
      StotRootL2 += RootL2[i] * coverPFT[i];
    } else if (i >= Shrubs && i < Shrubs + Perennials) {
      PtotRootL1 += RootL1[i] * coverPFT[i];
      PtotRootL2 += RootL2[i] * coverPFT[i];
    } else {
      AtotRootL1 += RootL1[i] * coverPFT[i];
    }
  }
}

/*******************************************************************************************
 * run infiltration routine into two layers:
 *   1. layer (upper): infiltration until ponding (saturation),
 *   				   afterwards Green & Ampt infiltration -->
 *wetting front
 *   2. layer (lower): deep infiltration via root channels
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
bool WaterCell::twoLayerInfiltration(
    double meanWaterL1) { /*Tong, based on PFT-oriented*/

  double A, B;           //!< variables to calculate the cumulative infiltration
  double vegFunction;    //!< function to describe the functional response of
                         //!< infiltration to vegetation
  double shrubFunction;  //!< same, exclusively for shrubs
  double a3;  //!< parameter for FL2, determines infiltration without shrub
  double a1;  //!< vegFunc(0.0, 0.0) = 0.2; ... at least some infiltration, if
              //!< there is no vegetation
  double a2;  //!< vegFunc(0.5, 0.5) = 0.95; ... almost full infiltration at a
              //!< total cover of 1

  // infiltration without shrubs is determined by hydraulic conductivity
  if (K_geom < 1)
    a3 = 0.1;
  else if (K_geom < 10)
    a3 = 0.5;
  else
    a3 = 0.8;

  // cout << "Kgeom " << K_geom << endl;
  // cout << "Ku " << Ku << endl;
  shrubFunction =
      (a3 + (StotRootL2 + PtotRootL2)) /
      (1 + (StotRootL2 + PtotRootL2));  //_dirk01112010 to include grasses for
                                        // infiltration  /*Tong*/

  FL2 = min((waterL0_abs)*shrubFunction * maxFL2, maxAmountFL2);

  // remaining surface water after infiltration
  waterL0_abs -= FL2;

  ///----------------- slower infiltration into the upper layer
  ///-------------------

  // fillable pores
  Na = waterL1_sat - meanWaterL1;

  // Holling type II function: the more vegetation the bigger vegFunction
  // until saturation occurs, the strange value is chosen in a way that for
  // a cover of cGrass+cShrub=1 the function equals 0.95
  // and without vegetation the function equals 0.2
  a1 = 0.013;  // vegFunc(0.0, 0.0) = 0.2; ... at least some infiltration, if
               // there is no vegetation // previous 0.013
  a2 = 0.067;  // vegFunc(0.5, 0.5) = 0.95; ... almost full infitration at a
               // total cover of 1 // previous 0.067

  vegFunction =
      (a1 + (StotRootL1 + PtotRootL1 + AtotRootL1)) /
      (a2 + (StotRootL1 + PtotRootL1 + AtotRootL1));  // adjusted by /*Tong*/

  if (waterL0_abs <= K_geom) {
    Fs = waterL0_abs;
    ts = -1.0;
    FL1 = Fs;
    wettingFrontFlag = false;
  }

  // else ponding can theoretically occur
  else {
    // there was no ponding in the previous time step
    if (wettingFrontFlag == false) {
      ts = ((Na * suction) / (waterL0_abs / K_geom - 1) / waterL0_abs) *
           (vegFunction);
      // ts =
      // ((Na*suction)/(((waterL0_abs/K_geom)-1)*waterL0_abs))*(vegFunction);
      Fs = min(1.0, ts) * waterL0_abs;

    }
    // there was already ponding in the previous time step
    else {
      ts = 0.0;
      Fs = 0.0;
    }
  }

  // if ponding occurs in this time step
  if ((ts >= 0) && (ts < 1)) {  // previous ts < 1
    // this routine is taken from WaSiM-ETH
    // after Peschke (in Dyck, Peschke 1989)
    A = K_geom * (1 - ts);
    B = Fs + 2 * (Na * suction);
    if (B < 0) B = 0;
    double h5 = (double)A * A / 4.0 + A * B + Fs * Fs;

    if (h5 < 0) cout << "SQRT bei h5 " << h5 << "\n";
    FL1 = min(
        (A / 2.0 + (A * A + 2 * A * B) / (4.0 * sqrt(h5))) * (vegFunction) + Fs,
        waterL0_abs);
    wettingFrontFlag = true;
  }

  // if ponding does not occur
  else {
    FL1 = Fs;
    wettingFrontFlag = false;
  }

  waterL0_abs -= FL1;
  waterL1_abs += FL1;
  waterL2_abs += FL2;
  waterL1_rel = waterL1_abs / depthL1;
  waterL2_rel = waterL2_abs / depthL2;

  // if the field capacity of the first layer is exceeded,
  // there is water flow to the lower layer

  if (waterL1_rel > FC) {
    drainL1L2 = (waterL1_rel - FC) * depthL1;
    waterL2_abs += drainL1L2;
    waterL1_rel = FC;
    waterL1_abs = waterL1_rel * depthL1;
    waterL2_rel = waterL2_abs / depthL2;
  }

  // if the field capacity of the second layer is exceeded,
  // there is deep drainage
  if (waterL2_rel > FC) {
    draindeep = (waterL2_rel - FC) * depthL2;
    waterL2_rel = FC;
    waterL2_abs = waterL2_rel * depthL2;
  }

  if (waterL0_abs > 0)
    return true;
  else
    return false;
}

/*******************************************************************************************
 * diffusion between the upper and the lower soil layer
 * is called in --> WaterLandscape::calculateProcesses
 *******************************************************************************************/
void WaterCell::layerDiffusion() {
  // water transport between layers (darcy's law)
  if (waterL1_rel > waterL2_rel) {
    double helpbalance =
        diffConst * K_geom * ((waterL1_rel - waterL2_rel) * depthL1) /
        ((depthL1 + depthL2) / 2.0);  // Darcy's Law, Maidment 5.17
    waterL1_abs -= helpbalance;
    waterL2_abs += helpbalance;
    waterL1_rel = waterL1_abs / depthL1;
    waterL2_rel = waterL2_abs / depthL2;
    layerDiff = helpbalance;

  } else if (waterL2_rel > waterL1_rel) {
    double helpbalance = diffConst * K_geom *
                         ((waterL2_rel - waterL1_rel) * depthL2) /
                         ((depthL1 + depthL2) / 2.0);
    waterL2_abs -= helpbalance;
    waterL1_abs += helpbalance;
    waterL2_rel = waterL2_abs / depthL2;
    waterL1_rel = waterL1_abs / depthL1;
    layerDiff = -helpbalance;

  } else
    layerDiff = 0;
}

/*******************************************************************************************
 * regular Constructor for  soil parameters
 * is called in --> WaterLandscape::initialiseSoilParameters
 *******************************************************************************************/   //constructor has been changed by /*Tong*/
WaterCell::WaterCell(int x, int y, double suction_, double Ks_,
                     double waterL1_sat_, double WP_, double FC_,
                     double depthL1_, double depthL2_, double iniWaterL0abs_,
                     double iniWaterL1_, double iniWaterL2_, double EPfactor_,
                     double resWater_, double maxFL2_, double maxAmountFL2_,
                     double diffConst_, double TE_, int PFTs_, int Perennials_,
                     int Shrubs_, int Annuals_, double Vegcoef1_,
                     double Vegcoef2_, double Veg_ETcoef1_,
                     double Veg_ETcoef2_) {
  suction = suction_;
  Ks = Ks_;
  WP = WP_;
  FC = FC_;
  depthL1 = depthL1_;
  depthL2 = depthL2_;
  waterL1_rel = iniWaterL1_;
  waterL2_rel = iniWaterL2_;
  waterL1_sat = waterL1_sat_;
  waterL0_abs = iniWaterL0abs_;
  waterL1_abs = waterL1_rel * depthL1;
  waterL2_abs = waterL2_rel * depthL2;
  EPfactor = EPfactor_;
  resWater = resWater_;
  maxFL2 = maxFL2_;
  maxAmountFL2 = maxAmountFL2_;
  diffConst = diffConst_;
  TE = TE_;
  Vegcoef1 = Vegcoef1_;       /*Tong*/
  Vegcoef2 = Vegcoef2_;       /*Tong*/
  Veg_ETcoef1 = Veg_ETcoef1_; /*Tong*/
  Veg_ETcoef2 = Veg_ETcoef2_; /*Tong*/

  PFTs = PFTs_;
  Perennials = Perennials_; /*Tong*/
  Shrubs = Shrubs_;         /*Tong*/
  Annuals = Annuals_;       /*Tong*/

  xpos = x;
  ypos = y;

  // in: Kemp et al. 1997, p.81
  // following Gardner 1958
  // unsaturated hydraulic conductivity
  //  = Ks / (1 + Psi / Psi*)^p
  Ku = (double)Ks / (pow(1.0 + (double)((-30.0) / (-26.0)), 3.0));
  double h6 = Ks * Ku;
  if (h6 < 0) cout << "Problem bei h6 " << h6 << "\n";
  K_geom = sqrt(h6);

  wettingFrontFlag = false;
  EPL1 = 0;
  EPL2 = 0;
  ts = -1;
  Fs = 0;
  FL1 = 0;
  FL2 = 0;
  QD = 0;

  coverPFT.resize(PFTs, 0);
  UseL1.resize(PFTs, 0);
  UseL2.resize(PFTs, 0);
  PftTL1.resize(PFTs, 0);
  PftTL2.resize(PFTs, 0);
  RootL1.resize(PFTs, 0);
  RootL2.resize(PFTs, 0);
}

/*******************************************************************************************
 * default constructor
 *******************************************************************************************/
WaterCell::WaterCell() {}

/*******************************************************************************************
 * destructor for soil
 *******************************************************************************************/
WaterCell::~WaterCell() {}

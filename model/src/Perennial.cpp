/*******************************************************************************************
 * Perennial.cpp
 * Methods defined here are unique to perennial grasses.
 *******************************************************************************************/

#include "Perennial.h"
#include "RandomNumberGenerator.h"

/*******************************************************************************************
 *  Calculate average moisture for layer two
 *******************************************************************************************/
void Perennial::calculateAvMoistL2(VegetationCell* vegetation) {

	AvMoistL2 = max((vegetation->moistL2 - WPWaterL2) * vegetation->soilDepthL2,0.0);
}

/*******************************************************************************************
 *  Calculate growth depending on transpiration
 *  layer one
 *******************************************************************************************/
void Perennial::calculateGrowthL1(VegetationCell* vegetation, int index) {

	GrowthL1 = vegetation->TL1[index] * Trancoef * (1.0 - min(vegetation->gCover / (max(cmax - (1 - overlap) * vegetation->sCover - vegetation->resCover, 0.0001)), 1.0));
	GrowthL1 *= (vegTimeStep/(1.0*daysPerYear));
}

/*******************************************************************************************
 *  Calculate growth depending on transpiration
 *  layer two
 *******************************************************************************************/
void Perennial::calculateGrowthL2(VegetationCell* vegetation, int index) {

	GrowthL2 = vegetation->TL2[index] * Trancoef * (1.0 - min(vegetation->gCover / (max(cmax - (1 - overlap) * vegetation->sCover, 0.0001)), 1.0));
	GrowthL2 *= (vegTimeStep/(1.0*daysPerYear));
}

/*******************************************************************************************
 *  Calculate overall growth
 *******************************************************************************************/
void Perennial::calculateGrowth() {

	Growth = GrowthL1 + GrowthL2;
		if (Growth < 0) {
			Growth = 0;
		}
}

/*******************************************************************************************
 *  Basic mortality
 *******************************************************************************************/
void Perennial::basicMortality(VegetationCell* vegetation, int index) {

	basicM = 0;
	//grasses have a gBasicMortRate  Basic mortality
    double randomNumber = RandomNumberGenerator::generate_random_float(0.0,2.0);
	basicM =  Cover * randomNumber * BasicMortRate; //equally distributed between 0 and 200% of gBasicMortRate
	Cover -= basicM;
}

/*******************************************************************************************
 *  Calculate vegetation change in one grid
 *  cell during one timestep
 *******************************************************************************************/
void Perennial::vegMortality(VegetationCell* vegetation, int day, int growStart, int growEnd, int year, Weather * weather) {

	//available water for life-support
	AvMoistSuppL1 = max((vegetation->moistSeasonL1 - WPWaterL1) * vegetation->soilDepthL1,0.0);
	AvMoistSuppL2 = max((vegetation->moistSeasonL2 - WPWaterL2) * vegetation->soilDepthL2,0.0);

	//grass mortality
	DryMortL1 = MortR * Cover * (1 - min(AvMoistSuppL1 * UseL1, 1.0)) * (RootL1 / (RootL1 + RootL2));
	DryMortL2 = MortR * Cover * (1 - min(AvMoistSuppL2 * UseL2, 1.0)) * (RootL2 / (RootL1 + RootL2));
	DryMortL1 *= (growEnd - growStart) / (daysPerYear * 1.0);   // _dirk Mort pro growing season nicht aufs Jahr
	DryMortL2 *= (growEnd - growStart) / (daysPerYear * 1.0);
	DryMort = DryMortL1 + DryMortL2;
	if (DryMortL1 == 0) DryMortL2 *= 0.5;  //_dirk03112010 to reduce mortality I divide L2 mortality by 2 if L1 is  good enough to grow undisturbed

	//shrub mortality
	totalDryMort += DryMort;

	if (Cover<=0)totalDryMort = 0;
	else  totalDryMort /= Cover;
	Cover -= DryMort;
}

/******************************************************************************************
 * calculate fire disturbance on perennial grass
 ******************************************************************************************/
void Perennial::fire(double fireIntens, double fuel, int &cohortAge, int &cohortAgeFire){

	double FireLoss;

	if (fireIntens > 300){
			FireLoss = specLoss * Cover;
			} else {
			FireLoss = 0;
		}
	Cover -= FireLoss;
}

/********************************************//**
 *  Constructor
 ***********************************************/    //constructor has been changed by /*Tong*/
Perennial::Perennial(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef,
                     double _Tcoef, double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap,
                     double _conversionC_BM, double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW,
                     double _bm_c_rain, double _MAP, double _grazeParA, double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS,
                     double _distConstS2, double _reserveLast, double _reserveThis, int _vegTimeStep, int _cellsize) : Grass(_Name, _Type, _UptakeRate, _RootL1, _RootL2, _GrowthR, _Transcoef,
                                                                                                              _Kcoef,  _ecoef, _Tcoef, _LAI_C, _MortR, _BasicMortRate, _Grazed,
                                                                                                              _WPWaterL1, _WPWaterL2, _cmax, _overlap, _conversionC_BM, _grazeParam,
                                                                                                              _grazeLim, _grazeprefer, _grazepreferintra, _disFrac, _estW, _bm_c_rain,
                                                                                                              _MAP, _grazeParA, _grazeParB, _specLoss, _firecoeff, _fuelLim,  _dist0S,
                                                                                                              _distConstS,  _distConstS2, _reserveLast, _reserveThis, _vegTimeStep, _cellsize) {

}

/********************************************//**
 *  Destructor
 ***********************************************/
Perennial::~Perennial() {

}

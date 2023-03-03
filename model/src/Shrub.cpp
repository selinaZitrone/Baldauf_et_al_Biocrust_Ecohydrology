/*******************************************************************************************
 * Shrub.cpp
 * Methods defined here are unique to all shrubs.
 *******************************************************************************************/

#include "Shrub.h"
#include "RandomNumberGenerator.h"

/*******************************************************************************************
 *  Calculate average moisture for layer two
 *******************************************************************************************/
void Shrub::calculateAvMoistL2(VegetationCell* vegetation) {

	AvMoistL2 = max((vegetation->moistL2 - WPWaterL2) * vegetation->soilDepthL2,0.0);
}

/*******************************************************************************************
 *  Calculate growth depending on transpiration
 *  layer one
 *******************************************************************************************/
void Shrub::calculateGrowthL1(VegetationCell* vegetation, int index) {

	GrowthL1 = vegetation->TL1[index] * Trancoef * (1 - min(vegetation->sCover / (max(cmax - (1 - overlap) * vegetation->gCover, 0.01)), 1.0));
	GrowthL1 *= (vegTimeStep/(1.0*daysPerYear));
}

/*******************************************************************************************
 *  Calculate growth depending on transpiration
 *  layer two
 *******************************************************************************************/
void Shrub::calculateGrowthL2(VegetationCell* vegetation, int index) {

	GrowthL2 = vegetation->TL2[index] * Trancoef * (1 - min(vegetation->sCover / (max(cmax - (1 - overlap) * vegetation->gCover, 0.01)), 1.0));
	GrowthL2 *= (vegTimeStep/(1.0*daysPerYear));
}

/*******************************************************************************************
 *  Calculate overall growth
 *******************************************************************************************/
void Shrub::calculateGrowth() {

	Growth = GrowthL1 + GrowthL2;
		if (Growth < 0) {
			Growth = 0;
		}
}

/*******************************************************************************************
 *  Calculate specific uptake rate
 *******************************************************************************************/
double Shrub::calculateUptakeRateCover(VegetationCell* vegetation, Weather* weather, int year) {

    double m = (1.0 - vegetation->bm_c_rain) / vegetation->MAP;
	UptakeRateCover = UptakeRate * (conversionC_BM * ((weather->precSumYear[year] * m) + vegetation->bm_c_rain)) * 0.5;
	return UptakeRateCover;
}

/*******************************************************************************************
 *  Basic mortality
 *******************************************************************************************/
void Shrub::basicMortality(VegetationCell* vegetation, int index) {

    /// \todo Dirk: For shrubs new additional function ageMortality, basicMortality as fraction or random process

	basicM = 0;
	//shrubs have at least a 1% basic mortality when the cohort is older than 40 years   ...quatsch
	/// \todo delete the following if-statement
		if (vegetation->cohortAge[index] > 40){//_dirk_10092010
            /// \todo rn is always 1! --> to check: new random number should be between 0 and 1
            double rn = RandomNumberGenerator::generate_random_float(0.0,1.0);
			rn = rn + 0.5;
			basicM = 0.01 * Cover * rn;
		}

	// age dependent mortality in form of dieback of complete shrubcover in the cell  according to Meyer and the Satchmo people
    double rn = RandomNumberGenerator::generate_random_float(0.0,1.0);
		if (rn < vegetation-> gBasicMortRate_min && vegetation->cohortAge[index] > 40){ /// \todo Dirk; add additional age based mortality parameter for shrubs
			basicM =  Cover ;
		}

	Cover -= basicM;
}

/******************************************************************************************
 *  Calculate vegetation change in one grid
 *  cell during one timestep
 ******************************************************************************************/
void Shrub::vegMortality(VegetationCell* vegetation, int day, int growStart, int growEnd, int year, Weather * weather) {

	AvMoistSuppL1 = max((vegetation->moistSeasonL1 - WPWaterL1) * vegetation->soilDepthL1, 0.0);
	AvMoistSuppL2 = max((vegetation->moistSeasonL2 - WPWaterL2) * vegetation->soilDepthL2, 0.0);	//shrub mortality
	DryMortL1 = MortR * Cover * (1 - min(AvMoistSuppL1 * UseL1, 1.0)) * (RootL1 / (RootL1 + RootL2));
	DryMortL2 = MortR * Cover * (1 - min(AvMoistSuppL2 * UseL2, 1.0)) * (RootL2 / (RootL1 + RootL2));
	DryMortL1 *= (growEnd - growStart) / (daysPerYear * 1.0);
	DryMortL2 *= (growEnd - growStart) / (daysPerYear * 1.0);
    double randomNumber = RandomNumberGenerator::generate_random_float(0.0,2.0);
	//randomNumber *= 2.0;/// to be deleted \todo Dirk: random Number local
	DryMortL1 *= randomNumber;
	DryMortL2 *= randomNumber;
	DryMort = DryMortL1 + DryMortL2;
	totalDryMort += DryMort;
		if (Cover <= 0) {
			totalDryMort = 0;
		} else {
			totalDryMort /= Cover;
		}
	Cover -= DryMort;
}

/****************************************************************************************
 * calculate fire disturbance on shrub
 ****************************************************************************************/
void Shrub::fire(double fireIntens, double fuel, int &cohortAge, int &cohortAgeFire){

	double FireLoss = 0.0;

	if (fireIntens > 300){
        if (cohortAgeFire < 3) {
            double zuf = RandomNumberGenerator::generate_random_float(0.0,1.0);
            if (zuf < 0.5){  //according to Joubert et al 2012 97% of these seedlings die in case of a fire
                FireLoss = Cover;
                cohortAge = 0;
                cohortAgeFire = 0;
            }
        }
        else if (cohortAgeFire < 11) {
                FireLoss = Cover * 0.3;
                cohortAge = 0;
                cohortAgeFire = 0;
        }
        else if (cohortAgeFire > 10) {
            FireLoss = min(specLoss * ((double) fuel / (cellsize * cellsize)) *firecoeff, 0.15) * Cover;
        }
        else {
            FireLoss = 0;
        }
    }
    else {
        FireLoss = 0;
    }
    Cover -= FireLoss;
}

/*******************************************************************************************
 *  Constructor
 *******************************************************************************************/   //constructor has been changed by /*Tong*/
Shrub::Shrub(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2,double _GrowthR,double _Transcoef, double _Kcoef, double _ecoef, double _Tcoef,
             double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap, double _conversionC_BM,
             double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain, double _MAP,
             double _grazeParA, double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2, double _reserveLast,
             double _reserveThis, int _vegTimeStep, int _cellsize) : Woody(_Name, _Type, _UptakeRate, _RootL1, _RootL2, _GrowthR, _Transcoef,  _Kcoef,  _ecoef, _Tcoef, _LAI_C,
                                                                           _MortR, _BasicMortRate, _Grazed, _WPWaterL1, _WPWaterL2, _cmax, _overlap, _conversionC_BM, _grazeParam,
                                                                           _grazeLim,  _grazeprefer,  _grazepreferintra, _disFrac, _estW, _bm_c_rain, _MAP, _grazeParA, _grazeParB,
                                                                           _specLoss, _firecoeff, _fuelLim,  _dist0S,  _distConstS,  _distConstS2, _reserveLast, _reserveThis,
                                                                           _vegTimeStep, _cellsize) {

}

/*******************************************************************************************
 *  Destructor
 *******************************************************************************************/
Shrub::~Shrub() {

}

/*******************************************************************************************
 * Annual.cpp
 * Methods defined here are unique for annual grasses.
 *******************************************************************************************/

#include "Annual.h"
#include "RandomNumberGenerator.h"

/*******************************************************************************************
 *  Calculate growth depending on soil
 *  layer one
 *******************************************************************************************/    //Growth has been changed by /*Tong/
void Annual::calculateGrowthL1(VegetationCell* vegetation, int index) {

    // j = random number in [0,1)
    double j = RandomNumberGenerator::generate_random_float(0.0,1.0);
    // [0,1) x number_of_annuals
    j = j * vegetation->Annuals;
    // if number_of_annuals 1, then annual_num is always 0
    int annual_num = floor(j);

    /// \todo Why this strange if condition?
	if (index - (vegetation->PFTs - vegetation->Annuals) == annual_num){
        GrowthL1 = min(AvMoistL1, 1.0) * GrowthR * (max(1.0 - vegetation->gCover - (1-overlap) * vegetation->sCover /*overlap*/ - vegetation->aCover - vegetation->resCover, 0.01) );
        GrowthL1 *= (vegTimeStep/(1.0*daysPerYear));
	}

    else {
        GrowthL1 = 0;
    }
}

/*******************************************************************************************
 *  Calculate overall growth
 *******************************************************************************************/
void Annual::calculateGrowth() {

	Growth = GrowthL1;
		if (Growth < 0) {
			Growth = 0;
		}
}

/*******************************************************************************************
 *  Calculate vegetation change in one grid
 *  cell during one timestep
 *******************************************************************************************/
void Annual::vegMortality(VegetationCell* vegetation, int day, int growStart, int growEnd, int year, Weather * weather) {

	AvMoistSuppL1 = max((vegetation->moistSeasonL1 - WPWaterL1) * vegetation->soilDepthL1, 0.0);
	DryMort = MortR * Cover * (1-min(AvMoistSuppL1 * 1, 1.0));
	DryMort *= (growEnd - growStart) / (daysPerYear * 1.0);
	Cover -= DryMort;
	if (day == daysPerYear-1) Cover = 0;
}

/*******************************************************************************************
 *  Constructor
 *******************************************************************************************/   // Constructor has been changed by /*Tong/
Annual::Annual(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef,
               double _Tcoef, double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap,
               double _conversionC_BM, double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain,
               double _MAP, double _grazeParA, double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2,
               double _reserveLast, double _reserveThis, int _vegTimeStep, int _cellsize) : Grass(_Name, _Type, _UptakeRate, _RootL1, _RootL2, _GrowthR,  _Transcoef,  _Kcoef,  _ecoef,
                                                                                                   _Tcoef, _LAI_C, _MortR, _BasicMortRate, _Grazed, _WPWaterL1, _WPWaterL2, _cmax,
                                                                                                   _overlap, _conversionC_BM, _grazeParam, _grazeLim, _grazeprefer, _grazepreferintra,
                                                                                                   _disFrac, _estW, _bm_c_rain, _MAP, _grazeParA, _grazeParB, _specLoss, _firecoeff,
                                                                                                   _fuelLim,  _dist0S, _distConstS,  _distConstS2, _reserveLast, _reserveThis,
                                                                                                   _vegTimeStep, _cellsize) {
}

/*******************************************************************************************
 *  Destructor
 *******************************************************************************************/
Annual::~Annual() {

}


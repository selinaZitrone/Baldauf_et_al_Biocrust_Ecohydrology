/*******************************************************************************************
 * Grass.cpp
 * Methods defined here are unique for all grasses.
 *******************************************************************************************/

#include "Grass.h"

/*******************************************************************************************
 *  Constructor
 *******************************************************************************************/
Grass::Grass(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef,
             double _Tcoef, double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap,
             double _conversionC_BM, double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain,
             double _MAP, double _grazeParA, double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2,
             double _reserveLast, double _reserveThis, int _vegTimeStep, int _cellsize) : PlantFunctionalType(_Name, _Type, _UptakeRate, _RootL1, _RootL2, _GrowthR,  _Transcoef,  _Kcoef,
                                                                                               _ecoef,  _Tcoef, _LAI_C, _MortR, _BasicMortRate, _Grazed, _WPWaterL1, _WPWaterL2,
                                                                                               _cmax, _overlap, _conversionC_BM, _grazeParam, _grazeLim,  _grazeprefer, _grazepreferintra,
                                                                                               _disFrac, _estW, _bm_c_rain, _MAP, _grazeParA, _grazeParB, _specLoss, _firecoeff, _fuelLim,
                                                                                               _dist0S,  _distConstS,  _distConstS2, _reserveLast, _reserveThis, _vegTimeStep, _cellsize) {
}

/*******************************************************************************************
 *  Destructor
 *******************************************************************************************/
Grass::~Grass() {

}

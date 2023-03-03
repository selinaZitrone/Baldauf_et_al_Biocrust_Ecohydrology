/****************************************************************************************//**
 * \class   Grass Grass.h
 * \brief   Methods defined here are unique for all grasses.
 * \author  Tong Guo
 *******************************************************************************************/

#ifndef GRASS_H_
#define GRASS_H_

#include "PlantFunctionalType.h"

//Forward declaration
class VegetationCell;

class Grass : public PlantFunctionalType {
    public:

        //Constructor and Destructor
        /***************************************************************************************//**
         * \brief Constructor
         * \todo include params here
         *******************************************************************************************/
        Grass(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef,
              double _Tcoef, double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap,
              double _conversionC_BM, double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain,
              double _MAP, double _grazeParA, double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2,
              double _reserveLast, double _reserveThis, int _vegTimeStep, int _cellsize);

        /***************************************************************************************//**
         * \brief Destructor
         *******************************************************************************************/
        virtual ~Grass();
};

#endif /* GRASS_H_ */

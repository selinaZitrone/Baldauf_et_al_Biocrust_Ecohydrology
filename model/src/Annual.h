/****************************************************************************************//**
 * \class   Annual Annual.h
 * \brief   Methods defined here are unique for annual grasses.
 * \author  Tong Guo
 *******************************************************************************************/

#ifndef ANNUAL_H_
#define ANNUAL_H_

#include "Grass.h"

class Annual : public Grass {
    public:
        //Constructor and Destructor
        /***************************************************************************************//**
         * \brief Constructor
         * \todo include params here or better not?
         *******************************************************************************************/
        Annual(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef, double _Tcoef,
               double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap, double _conversionC_BM,
               double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain, double _MAP,
               double _grazeParA, double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2,
               double reserveLast, double reserveThis, int _vegTimeStep, int _cellsize);

        /***************************************************************************************//**
         * \brief Destructor
         *******************************************************************************************/
        virtual ~Annual();

        //Virtual member functions
        /****************************************************************************************//**
         * \brief Calculate growth depending on soil layer one.
         * \param vegetation pointer to the VegetationCell object;
         * \param index index of the PFT array;
         *******************************************************************************************/
        virtual void calculateGrowthL1(VegetationCell* vegetation, int index);

        /****************************************************************************************//**
         * \brief Calculate overall growth.
         *******************************************************************************************/
        virtual void calculateGrowth();

        /****************************************************************************************//**
         * \brief Calculate vegetation change in one grid cell during one timestep.
         * \param vegetation pointer to the VegetationCell object;
         * \param growStart day of growing season start;
         * \param growEnd day of growing season end;
         * \param year current year;
         * \param weather pointer to the weather object;
         *******************************************************************************************/
        virtual void vegMortality(VegetationCell* vegetation, int day, int growStart, int growEnd, int year, Weather* weather);

};

#endif /* ANNUAL_H_ */

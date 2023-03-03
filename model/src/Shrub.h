/****************************************************************************************//**
 * \class   Shrub Shrub.h
 * \brief   Methods defined here are unique for shrubs.
 * \author  Tong Guo
 *******************************************************************************************/

#ifndef SHRUB_H_
#define SHRUB_H_

#include <iostream>

#include "Woody.h"

class Shrub : public Woody {
    public:

        //Constructor and Destructor
        /***************************************************************************************//**
         * \brief Constructor
         * \todo include params here
         *******************************************************************************************/
        //constructor has been changed by /*Tong*/
        Shrub(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef, double _Tcoef,
              double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap, double _conversionC_BM,
              double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain, double _MAP, double _grazeParA,
              double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2, double _reserveLast,
              double _reserveThis, int _vegTimeStep, int _cellsize);

        /***************************************************************************************//**
         * \brief Destructor
         *******************************************************************************************/
        virtual ~Shrub();

        //Virtual member functions
        /****************************************************************************************//**
         * \brief Calculate average moisture for layer two.
         * \param vegetation pointer to the VegetationCell object;
         *******************************************************************************************/
		virtual void calculateAvMoistL2(VegetationCell* vegetation);

        /****************************************************************************************//**
         * \brief Calculate specific uptake rate.
         * \param vegetation pointer to the VegetationCell object;
         * \param weather pointer to the weather object;
         * \param year current year;
         *******************************************************************************************/
		virtual double calculateUptakeRateCover(VegetationCell* vegetation, Weather* weather, int year);

        /****************************************************************************************//**
         * \brief Calculate growth depending on soil layer one.
         * \param vegetation pointer to the VegetationCell object;
         * \param index index of the PFT array;
         *******************************************************************************************/
        virtual void calculateGrowthL1(VegetationCell* vegetation, int index);

		/****************************************************************************************//**
         * \brief Calculate growth depending on soil layer two.
         * \param vegetation pointer to the VegetationCell object;
         * \param index index of the PFT array;
         *******************************************************************************************/
        virtual void calculateGrowthL2(VegetationCell* vegetation, int index);

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

        /****************************************************************************************//**
         * \brief Basic mortality.
         * \todo Drik: basic Mortality for shrubs and grasses differs, what is going on there, what about the parameters check for all types
         * \param vegetation pointer to the VegetationCell object;
         * \param index index of the PFT array;
         *******************************************************************************************/
		virtual void basicMortality(VegetationCell* vegetation, int index);

        /****************************************************************************************//**
         * \brief Calculate fire disturbance on shrubs.
         * \param fireIntens ?;
         * \param fuel ?;
         * \param cohortAge ?;
         * \param cohortAgeFire ?;
         *******************************************************************************************/
        virtual void fire(double fireIntens, double fuel, int &cohortAge, int &cohortAgeFire);
};

#endif /* SHRUB_H_ */

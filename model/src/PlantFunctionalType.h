/****************************************************************************************//**
 * \class   PlantFunctionalType PlantFunctionalType.cpp
 * \brief   Methods defined here are universal for all plant functional types.
 * \author  Tong Guo
 *******************************************************************************************/

#ifndef PLANTFUNCTIONALTYPE_H_
#define PLANTFUNCTIONALTYPE_H_

#include "Parameters.h"
#include "VegetationCell.h"
#include <iostream>  //to use cin / cout
#include <math.h>    //to use mathematical functions
#include <fstream>   //to read and write files
#include <cstdlib>   //to use exit() and for srand() and rand()
#include <algorithm> //to use min/max
#include <string.h>

/// \todo does a PFT know its own type (e.g. by string?), this would be useful!
/// \todo handing over all these parameters is really nasty, eventually we could manage to build a parameter object, that is created at the very beginning and handed over to all other classes if necessary. By this we also needed to make ONE change in case of a new Parameter

class PlantFunctionalType {
	public:

        //Constructor and Destructor
        /***************************************************************************************//**
         * \brief Constructor
         * \todo include params here
         * \todo we should get rid of this parameter mess!
         *******************************************************************************************/
        PlantFunctionalType(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef, double _Tcoef, double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap, double _conversionC_BM, double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain, double _MAP, double _grazeParA,
	    		            double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2, double _reserveLast, double _reserveThis, int _vegTimeStep, int _cellsize);

        /***************************************************************************************//**
         * \brief Destructor
         *******************************************************************************************/
	    virtual ~PlantFunctionalType();

        //Member functions
        /****************************************************************************************//**
         * \brief Output parameters!
         *******************************************************************************************/
		void output();

		//Virtual member functions
        /****************************************************************************************//**
         * \brief Calculate average moisture for layer one.
         * \param vegetation pointer to the VegetationCell object;
         *******************************************************************************************/
		virtual void calculateAvMoistL1(VegetationCell* vegetation);

        /****************************************************************************************//**
         * \brief Calculate average moisture for layer two.
         * \param vegetation pointer to the VegetationCell object;
         *******************************************************************************************/
		virtual void calculateAvMoistL2(VegetationCell* vegetation) {};

        /****************************************************************************************//**
         * \brief Calculate growth depending on soil layer one.
         * \param vegetation pointer to the VegetationCell object;
         * \param index index of the PFT array;
         *******************************************************************************************/
        virtual void calculateGrowthL1(VegetationCell* vegetation, int index) {};

        /****************************************************************************************//**
         * \brief Calculate growth depending on soil layer two.
         * \param vegetation pointer to the VegetationCell object;
         * \param index index of the PFT array;
         *******************************************************************************************/
        virtual void calculateGrowthL2(VegetationCell* vegetation, int index) {};

        /****************************************************************************************//**
         * \brief Calculate overall growth.
         *******************************************************************************************/
		virtual void calculateGrowth() {};

        /****************************************************************************************//**
         * \brief Calculate total growth.
         *******************************************************************************************/
		virtual double calculateTotalGrowth();

        /****************************************************************************************//**
         * \brief Add growth to cover.
         *******************************************************************************************/
		virtual void addGrowthToCover();

        /****************************************************************************************//**
         * \brief Calculate the corresponding biomass on a cell by a linear regression.
         * \param vegetation pointer to the VegetationCell object;
         * \param we2 pointer to the weather object;
         * \param yr2 current year;
         *******************************************************************************************/
		virtual void calculateBiomass(VegetationCell* vegetation, Weather* we2, int yr2);

        /****************************************************************************************//**
         * \brief Calculate the corresponding cover on a cell by a linear regression.
         * \param vegetation pointer to the VegetationCell object;
         * \param we2 pointer to the weather object;
         * \param yr2 current year;
         * \param factor ?;
         *******************************************************************************************/
		virtual void calculateCover(VegetationCell* vegetation, Weather* we2, int yr2, double factor);

        /****************************************************************************************//**
         * \brief Calculate specific uptake rate.
         * \param vegetation pointer to the VegetationCell object;
         * \param weather pointer to the weather object;
         * \param year current year;
         *******************************************************************************************/
		virtual double calculateUptakeRateCover(VegetationCell* vegetation, Weather* weather, int year);

        /****************************************************************************************//**
         * \brief Calculate vegetation change in one grid cell during one timestep.
         * \param vegetation pointer to the VegetationCell object;
         * \param growStart day of growing season start;
         * \param growEnd day of growing season end;
         * \param year current year;
         * \param weather pointer to the weather object;
         *******************************************************************************************/
        virtual void vegMortality(VegetationCell* vegetation, int day, int growStart, int growEnd, int year, Weather * weather) {};

        /****************************************************************************************//**
         * \brief Basic mortality.
         * \param vegetation pointer to the VegetationCell object;
         * \param index index of the PFT array;
         *******************************************************************************************/
		virtual void basicMortality(VegetationCell* vegetation, int index) {};

        /****************************************************************************************//**
         * \brief Calculate fire disturbance on perennial grass.
         * \todo add params here and their description!
         * \param fireIntens ?;
         *******************************************************************************************/
        virtual void fire(double fireIntens, double fuel, int &cohortAge, int &cohortAgeFire) {};

		//Parameters
		int vegTimeStep;
		int cellsize;
		string Name;                //!< Name of PFTs;                                                                     /*Tong*/
		string Type;                //!< Type of PFTs;                                                                     /*Tong*/
		double UptakeRate; 	        //!< water uptake rate per unit biomass [mm/(g*a)]
		double UptakeRateCover; 	//!< water uptake rate per unit cover [mm/(cover*a)]
		double URL1_scaled;		    //!< scaled such that the fraction in calculateWaterUptake is 1
		double URL2_scaled;		    //!< scaled such that the fraction in calculateWaterUptake is 1
		double URL1_help;		    //!< scaled such that the fraction in calculateWaterUptake is 1
		double URL2_help;		    //!< scaled such that the fraction in calculateWaterUptake is 1
		double RootL1;		        //!< root distribution: proportion of roots in the upper layer
		double RootL2;		        //!< root distribution: proportion of roots in the lower layer
		double GrowthR;		        //!< theoretical growth rate
		double Trancoef;            //!< Transpiration rate of plant functional types                                        /*Tong*/
		double Kcoef;               //!< coefficient of light attenuation                                                    /*Tong*/
		double ecoef;               //!< exponential coefficient                                                             /*Tong*/
		double temcoef;             //!< temperature coefficient                                                             /*Tong*/
		double LAI_C;               //!< transformed coeffcient of leaf area to vegetation cover                             /*Tong*/

		double MortR; 		        //!< mortality rate [1/a]
		double WPWaterL1;	        //!< minimal amount of water needed for growth [vol]
		double WPWaterL2;	        //!< minimal amount of water needed for growth [vol]
		double conversionC_BM;      //!<  conversion factor c (-) in biomass_g (g/ha)
		double cmax;		        //!< capacity
		double overlap;			    //!< possible overlapping of grasses and shrubs                                     /*Tong*/
		double grazeParam;          //!<  _dirk slope of function for alpha . alpha determines percentage (of expected coverdecrease) of cover decrease that is realized
		double grazeParA;
		double grazeParB;
		double grazeLim;            //!<  _dirk percentage of biomass that is grazed
		double grazeprefer;         //!< grazed preference for each PFT (inter-competition)                                /*Tong*/
		double grazepreferintra;    //!< grazed preference for each PFT (intra-competition)                              /*Tong*/
		double disFrac;             //!< Plant fraction can disperse
		double dist0S;              //!< determine the seed decline with distance
		double distConstS;          //!< determine the seed decline of the exponential function with distance
		double distConstS2;         //!< determine the seed decline of the exponential function with distance, related with grazing intensity
		double reserveLast;         //!< fraction of previous year that survives                                           /*Tong*/
		double reserveThis;         //!< fraction of this year's biomass that goes into reserve biomass                     /*Tong*/
		double estW;                //!< _dirk  establishment if  gestW* gWPWaterL1 was available on average during the respective growing season    /*Tong*/
		double bm_c_rain;           //!< rain influencing coefficient on plant biomass
		double MAP;                 //!< mean anuual precipitation

		//parameters van langevelde
		double specLoss ;	        //!< specific (fixed) loss (fraction)of grass biomass per year y^-1  (per year...
		double firecoeff;           //!<  coeff for increase in fire intensity with grass biomass W* m^-2 * g^-1
		int fuelLim;                //!<  minimum grass biomass that is needed for a fire to kill anything

		//this parameter should not be set in the vegetation model!
		double resWater;			//!< minimal water content in the soil

		//these parameters are updated according to the timestep
		double moistL1;		        //!< moisture in upper layer
		double moistL2;		        //!< moisture in lower layer
		double AvMoistL1;           //!< how much water is there for grass growth in L1
		double AvMoistL2;           //!< how much water is there for grass growth in L2
		double AvMoistSuppL1;       //!< how much water is there for grass support in L1
		double AvMoistSuppL2;       //!< how much water is there for grass support in L2
		double AvMoistSeasonL1;	    //!< how much water is there for grass growth in L1 (in whole season)
		double AvMoistSeasonL2;	    //!< how much water is there for grass growth in L2 (in whole season)

		double UseL1;	            //!< theoretical water use because of competitive effects (layer1)
		double UseL2;	            //!< theoretical water use because of competitive effects (layer2)

		double Rcover;              //!< total respiration of vegetation cover
		double MRcover;             //!< maintenance respiration of vegetation cover
		double GrowthL1;	        //!< growth because of water in L1
		double GrowthL2;	        //!< growth because of water in L2
		double BasicMortRate;       //!< basic mortality rate, independent of water \todo what is the difference between BasicMort and BasicMortRate?

		double Growth;		        //!< total growth in one time step (increase of veg cover)
		double BasicMort;	        //!< total basic mortality in one time step \todo what is the difference between BasicMort and BasicMortRate?
		double Grazed;
		bool GrazedValidation;
		double GrazedBM;
		double DryMortL1;
		double DryMortL2;
		double DryMort;		        //!< total mortality due to water shortage

		double totalDryMort;        //!< _dirk to sum up mortality
		double totalGrowth;         //!< _dirk to sum up all year growth
		double basicM;              //!< _dirk to sum up all year growth

		//changing vegetation parameters
		double Cover;	            //!< total cover of this PFT in one cell
		double Biomass;	            //!< biomass per cell
		double Biomassold;	        //!< biomass before grazing per cell
		double resBM;	            //!< _dirk_22112011 reserve biomass per cell

		//actual dispersal in each cell
		double Disp;		        //!< received
		double DispSend;	        //!< sent

		//changing grazing parameters
		double safeBM;
		double BMavailable;
		double BMres;

		int est;	                //!< establishment this season (1 true or 0 false)
		int ger;                    //!< _dirk can seeds germinate (after one good year) (important for fire management)
		int est2;                   //!< _dirk for fire initialization, memorizes the establishment conditions of one year before
		int est3;                   //!< _dirk for fire initialization, memorizes the establishment conditions of two year before

		double SpecLoss;
};

#endif /* PLANTFUNCTIONALTYPE_H_ */

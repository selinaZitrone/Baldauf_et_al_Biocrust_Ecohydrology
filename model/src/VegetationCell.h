/***************************************************************************************//**
 * \class   VegetationCell VegetationCell.h
 * \brief   Describes the vegetation dynamics in one grid cell.
 * \author  Britta Tietjen
 *******************************************************************************************/

#ifndef VEGETATIONCELL_H_
#define VEGETATIONCELL_H_

#include <iostream>  //to use cin / cout
#include <math.h>    //to use mathematical functions
#include <fstream>   //to read and write files
#include <cstdlib>   //to use exit() and for srand() and rand()
#include <algorithm> //to use min/max
#include <string>    //to use strings
#include <vector>

#include "Weather.h"

using namespace std;

//Forward declaration
class PlantFunctionalType;

class VegetationCell{
	private:

        //Private parameters
        int vegTimeStep;
        int cellsize;

	public:

        //Constructor and Destructor
        /***************************************************************************************//**
         * \brief Constructor
         * \todo include param here
         *******************************************************************************************/
		VegetationCell(int PFTs_, int Annuals_, int Perennials_, int Shrubs_, double resWater_, double moistL1_, double moistL2_, double overlap_, double grazeParam_,
             double aGrowthConst_,double bm_c_rain_, double grazeparA_, double grazeparB_, double MAP_,  double firecoeff_, int fuelLim_, double EnScover_,
             double initialshare_, int vegTimeStep_, int cellsize_);

        /***************************************************************************************//**
         * \brief Destructor
         *******************************************************************************************/
		~VegetationCell(); //Destructor


        //Public parameters
        int PFTs;
        int Perennials;                                         /*Tong*/
        int Annuals;                                            /*Tong*/
        int Shrubs;                                             /*Tong*/

        double totalCover;                                      //!< total vegetation cover
		double shrubCover;                                      //!< total shrub cover
		double perennialCover;                                  //!< total perennial cover
		double annualCover;                                     //!< total annual cover

		Double1D safeBM_reserve;                                //!< a certain fraction of reserve biomass is non removable safeBM      /*Tong*/
		Double1D safeBM_alive;                                  //!< a certain fraction of alive biomass is non removable safeBM          /*Tong*/
		Double1D edibleBM_reserve;                              //!< a certain fraction of reserve biomass is non removable safeBM    /*Tong*/
		Double1D edibleBM_alive;                                //!< a certain fraction of alive biomass is non removable safeBM           /*Tong*/
		Double1D edibleBM_total;                                //!< total edible biomass for each PFT                                /*Tong*/
		double edibleShrub;                                     //!< total edible biomass of shrubs in a cell                                    /*Tong*/
		double edibleAnnual;                                    //!< total edible biomass of annuals in a cell                                  /*Tong*/
		double ediblePerennial;                                 //!< total edible biomass of perennials in a cell                              /*Tong*/
		double reserveBM;                                       //!< reserve biomass in a cell                                                  /*Tong*/
		double aliveBM;                                         //!< alive biomass in a cell                                                    /*Tong*/
		double totalBM;                                         //!< total biomass in a cell                                                    /*Tong*/
		double AnnualresBM;                                     //!< reserve annual biomass in a cell                                                       /*Tong*/
        double PerennialresBM;                                  //!< reserve perennial biomass in a cell                                                       /*Tong*/
	    double ShrubaliveBM;                                    //!< alive shrub biomass in a cell                                                         /*Tong*/
        double AnnualaliveBM;                                   //!< alive annual biomass in a cell                                                         /*Tong*/
        double PerennialaliveBM;                                //!< alive perennial biomass in a cell                                                         /*Tong*/

		Double1D relprefer;                                     //!< relative grazing preference for each PFT                                         /*Tong*/

		//constant vegetation parameters
		double overlap;			                                //!< possible overlapping of grasses and shrubs
		double grazeParam;                                      //!< _dirk slope of function for alpha . alpha determines percentage (of expected coverdecrease) of cover decrease that is realized
		double MAP;                                             //!< _dirk mean annual precipitation
		double aGrowthConst;
		double bm_c_rain;                                       //!< _dirk
		double grazeparA;                                       //!< _dirk
		double grazeparB;                                       //!< _dirk
		double resWater;			                            //!< minimal water content in the soil
		double firecoeff;                                       //!< coeff for increase in fire intensity with grass biomass W* m^-2 * g^-1
		int fuelLim;                                            //!< minimum grass biomass that is needed for a fire to kill anything
		double EncroachScover;                                  //!< maximum encroachment first shrub cover
		double initialshare;                                    //!< initial share of grazing for each cell
		double moistL1;		                                    //!< moisture in upper layer
		double moistL2;		                                    //!< moisture in lower layer

		//these parameters are updated according to the timestep //moisture parameter
		double moistSeasonL1;	                                //!< mean moisture during growing season
		double moistSeasonL2;
		double moistWetSeasonL1;                                //!< mean moisture during wet season
		double moistWetSeasonL2;
		double resBM; 	                                        //!< _dirk_22112011 total reserve biomass per cell
		double moistL1mem[3];
		double soilDepthL1; 	                                //!< this is not nice, this should not be in here!
		double soilDepthL2;	                                    //!< since it is part of the soil model

		//changing vegetation parameters
		double totalPFTcover;                                   //!< total PFT cover in one cell                                                 /*Tong*/
		double sCover;	                                        //!< total shrub cover in one cell                                                      /*Tong*/
		double gCover;	                                        //!< total grass cover in one cell                                                       /*Tong*/
		double aCover;                                          //!< total annual grass cover in one cell                                                /*Tong*/
		double resCover;                                        //!< _dirk cover of "reserve" biomass
		double resAnnualcover;                                  //!< cover of reserve part for annual grass                                         /*Tong*/
		double resPerennialcover;                               //!< cover of reserve part for perennial grass                                   /*Tong*/
		double gBasicMortRate_min;                              //!< choose the minimum growth basic mortality rate for multi-perennials        /*Tong*/

		//shrub establishment parameter
		int SestS2;                                             //!< _dirk for fire initialization, memorizes the establishment conditions of one year before
		int SestS3;                                             //!< _dirk for fire initialization, memorizes the establishment conditions of two year before
		Int1D cohortAge;                                        //!< counts years from LAST establishment event
		Int1D cohortAgeFire;                                    //!< counts age of oldest plant in cell, so that fire can be simulated accordingly

        Double1D TL1;                                           //!< Transpiration in the upper layer                /*Tong*/
		Double1D TL2;                                           //!< Transpiration in the lower layer        /*Tong*/
	    Double1D yearlyPftTL1;                                  //!< transpiration in the upper layer;           /*Tong*/
	    Double1D yearlyPftTL2;                                  //!< transpiration in the lower layer;          /*Tong*/

        double yearlyTL1;                                       /*Tong*/
	    double yearlyTL2;                                       /*Tong*/
	    double yearlyE;                                         /*Tong*/

        Double1D Pcover;                                        /*Tong*/
        Double1D Acover;                                        /*Tong*/
        Double1D Scover;                                        /*Tong*/

        Double1D UseL1;                                         /*Tong*/
        Double1D UseL2;                                         /*Tong*/

        //Root parameters
        Double1D RootL1;                                        /*Tong*/
        Double1D RootL2;                                        /*Tong*/

        //Fields
        vector<PlantFunctionalType*> pftList;		            //!< grid of plant functional types --> these are objects

		//Public member functions
        /***************************************************************************************//**
         * \brief Is called by the controller to set the soil moisture.
         * Calculate mean soil moisture of every vegTimeStep.
         * \param mL1 waterL1_rel;
         * \param mL2 waterL2_rel;
         * \param day current day;
         * \param grStart first day of growing season;
         * \param wetStart first day of water season;
         * \param wetEnd last day of water season;
         *******************************************************************************************/
		void setMoist(double mL1, double mL2, int day,int grStart, int grEnd, int wetStart, int wetEnd);

        /***************************************************************************************//**
         * \brief Return the transpiration of PFTs in the waterlandscape process.
         * \todo explain ??? param here
         * \param Double1D TransL1 ???
         * \param Double1D TransL2 ???
         * \param day current day;
         *******************************************************************************************/
		void setTransFactors(Double1D TransL1, Double1D TransL2, int day);

        /***************************************************************************************//**
         * \brief Calculate water uptake per unit cover.
         *******************************************************************************************/
		void setCover();          /*Tong*/

        /***************************************************************************************//**
         * \brief ???
         * \todo description here!
         *******************************************************************************************/
        void setpftUse();                          /*Tong*/

         /***************************************************************************************//**
         * \brief set roots layer two
         *******************************************************************************************/
        void setRootL2();                           /*Tong*/

        /***************************************************************************************//**
         * \brief set roots layer one
         *******************************************************************************************/
        void setRootL1();                          /*Tong*/

        /***************************************************************************************//**
         * \brief Calculate cover reduction in each cell due to fire.
         * Only called, if a fire actually occurs.
         * Revised the fire process based on PFTs.
         *******************************************************************************************/
		void calculatePFTfire();  // called, if a fire occurs, revised based on former calculateFire()

        /***************************************************************************************//**
         * \brief Calculate mean soil moisture during water season.
         * \todo Explain parameters here!
         * \param mL1 ???;
         * \param mL2 ???;
         * \param day current day;
         * \param grStart growing season start;
         * \param grEnd growing season end;
         *******************************************************************************************/
		void meanGrowingSeasonMoist(double mL1, double mL2, int day, int grStart, int grEnd);

        /***************************************************************************************//**
         * \brief Calculate mean soil moisture during wet season.
         * \todo Explain parameters here!
         * \param mL1 ???;
         * \param mL2 ???;
         * \param day current day;
         * \param wStart wet season start;
         * \param wEnd wet season end;
         *******************************************************************************************/
		void meanWetSeasonMoist(double mL1, double mL2, int day,int wStart, int wEnd);

        /***************************************************************************************//**
         * \brief Calculate mean soil moisture during wet season.
         * \todo Explain parameters and include brief description here
         * \param Evaporation ???;
         * \param TranspirationL1 ???;
         * \param TranspirationL2 current day;
         * \param day current day;
         *******************************************************************************************/
		void setEvaporationandTranspiration(double Evaporation, double TranspirationL1, double TranspirationL2, int day);   /*Tong*/

        /***************************************************************************************//**
         * \brief Calculate the water uptake per unit cover.
         *******************************************************************************************/
		void calculateWaterUptake();

		/***************************************************************************************//**
         * \brief Calculate total cover of shrubs in one cell.
         *******************************************************************************************/
        double calculateShrubCover();

        /***************************************************************************************//**
         * \brief Calculate total cover of annuals in one cell.
         *******************************************************************************************/
        double calculateAnnualCover();

		/***************************************************************************************//**
         * \brief Calculate total cover of perennials in one cell.
         *******************************************************************************************/
        double calculatePerennialCover();

        /***************************************************************************************//**
         * \brief Calculate total cover of all PFTs in one cell.
         *******************************************************************************************/
        double calculateTotalCover();

        /***************************************************************************************//**
         * \brief ???
         * \todo include brief description here
         *******************************************************************************************/
        double calculateResCover();

        /***************************************************************************************//**
         * \brief ???
         * \todo include brief description here
         *******************************************************************************************/
        double calculateResAnnualcover();         /*Tong*/

        /***************************************************************************************//**
         * \brief ???
         * \todo include brief description here
         *******************************************************************************************/
        double calculateResPerennialcover();                    /*Tong*/
};

#endif //VEGETATIONCELL_H



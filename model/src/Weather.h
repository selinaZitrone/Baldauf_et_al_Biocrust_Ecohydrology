/****************************************************************************************//**
 * \class   Weather Weather.h
 * \brief   In this file, all functions are described that are related to climate
 *          climate input is read;
 *          climate functions are calculated;
 *          growing seasons are calculated;
 * \author  Britta Tietjen
 *******************************************************************************************/

#ifndef WEATHER_H_
#define WEATHER_H_

#include "Parameters.h"

class Weather{
	private:
	
		const double PI = 3.14159265358979323846264338327950288;
	
	public:

        //Constructor and Destructor
        /***************************************************************************************//**
         * \brief Constructor
         *******************************************************************************************/
		Weather(Parameter* p);

        /***************************************************************************************//**
         * \brief Destructor
         *******************************************************************************************/
		~Weather(); //Destructor

        //Public parameters
		
        Double1D tempAvg;				                    //!< mean temperature during one day [�C]
		Double1D tempMax;				                    //!< maximal temperature during one day [�C]
		Double1D tempMin;				                    //!< minimal temperature during one day [�C]
		Double1D precSum;				                    //!< total precipitation during one day [mm/d]
        Double1D precSumYear;                               //!< total precipitation during one year [mm/a]
		Double1D moistBBAvg; 			                    //!< mean soil moisture between bushes [vol%]
		Double1D moistUBAvg;			                    //!< mean soil moisture under bushes [vol%]
		Double1D relHumAvg; 				                //!< mean relative humidity during one day [�C]
		Double1D windSpeedAvg;			                    //!< mean windSpeed during one day [m/s]
		Double1D solRadSum; 				                //!< total solar radiation during one day [mm/d]
		Double1D extraterrRadiation; 		                //!< evaporation equivalent of extraterrestrial evaporation [mm/d]
        Int1D growStart;		                            //!< when is the first day of that year with rainfall above xxx mm?
		Int1D growEnd; 			                            //!< when is the last day of that year with rainfall above xxx mm?
		Double2D prec;		                                //!< hourly precipitation [mm/h]
		Double2D temp;		                                //!< hourly temperature [�C]
		Double2D relHum;		                            //!< hourly relative air humidity [%]
		Double2D solRad;		                            //!< hourly solar radiation [W/m^2]
		Double2D windSpeed;	                                //!< hourly wind speed [m/s]
		Double2D tempBB;		                            //!< hourly temperature between bushes[�C]
		Double2D tempUB;		                            //!< hourly temperature under bushes[�C]
		Double2D moistBB;	                                //!< hourly soil moisture between bushes[vol%]
		Double2D moistUB; 	                                //!< hourly soil moisture under bushes[vol%]

		//Public member functions
        /****************************************************************************************//**
         * \brief read in weather data from file
         * \param p pointer to parameter object;
         *******************************************************************************************/
		void getWeather(Parameter* p);

        /****************************************************************************************//**
         * \brief calculated daily sums and means
         * \param p pointer to parameter object;
         *******************************************************************************************/
		void meanWeather(Parameter* p);

        /****************************************************************************************//**
         * \brief determine first and last day of growing season
         * depending on effective rain (for each year)
         * \param p pointer to parameter object;
         *******************************************************************************************/
		void growingSeason(Parameter* p);
};

#endif /*WEATHER_H_*/

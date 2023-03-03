/*******************************************************************************************
 * Weather.cpp
 * In this file, all functions are described that are related to climate
 * climate input is read
 * climate functions are calculated
 *******************************************************************************************/

#include <iostream>  //to use cin / cout
#include <string>    //to use strings
#include <math.h>    //to use mathematical functions
#include <fstream>   //to read and write files
#include <cstdlib>   //to use exit()

#include "Weather.h"

using namespace std;

/*******************************************************************************************
 * determine first and last day of growing season depending on effective rain (for each year)
 * is called in --> controller::initialiseSimulation
 *******************************************************************************************/
void Weather::growingSeason(Parameter* p) {

	bool foundFirst;
	double effectiveRain = 5.0; /// \todo effectiveRain should be a parameter

	for (int year = 0; year < p->simYears; year++) {
		foundFirst = false;
		/// \todo Dtodo: question to Dirk->should this be discussed why 17?
		for (int day = 0; day < daysPerYear-17; day++){  //_dirk eingeführt weil sonst bei seltenen regenereignissen ausserhalb des üblichen ranges Probleme auftreten
			//effective rain
			if (precSum[year*daysPerYear + day] >= effectiveRain){
				//first day of the year with rain > effective rain
				if (!foundFirst){
					growStart[year] = day;
					foundFirst = true;
				}
				growEnd[year] = day;
			}
		}
	}
}

/*******************************************************************************************
 * read parameters from weather file
 * this time the much simpler one, just with the format
 * Year | Month | Day | Hour | Prec | Temp | SolRad
 * if is called in --> controller::initialiseSimulation
 *******************************************************************************************/
void Weather::getWeather(Parameter* p){

	//open precipitation file
  	ifstream parameterFile;
  	string parameterFileName;
  	string dummy1; //to read in values that are not needed
	string param; //to read in a "real" parameter
  	int day, hour;

    parameterFileName = "Weather" + delimiter + p->site + "_" + std::to_string(p->simYears) + "years_" + p->weatherFileID[p->scenario] + "_climrep-" + std::to_string(p->climateRepetition + 1) + ".txt";

    cout << "Read weather from: " << parameterFileName << endl;

    parameterFile.open(parameterFileName.c_str(), ios::binary|ios::in);

  	//if something goes wrong...
  	if(!parameterFile){
  		cerr << parameterFileName
  		     << " could not be opened!\n";
  		exit(-1);
  	}

	cout << "________________________________________________________________" << endl;
	cout << endl;

	//read header
	for (int i = 0; i < 6; i++) {
		if (!parameterFile.eof()) parameterFile >> dummy1;
        else {
            cout << "No Data in " + parameterFileName;
            exit(-1);
        }
	}

	day = hour = 0;

	for (int i=0; i<(p->simYears*daysPerYear*hoursPerDay); i++) {

		//date
		parameterFile >> dummy1 >> dummy1 >> dummy1 >> dummy1;

		//precipitation
		parameterFile >> param;
		if (param != ".") prec[day][hour] = atof((param).c_str());
		else prec[day][hour] = nodata;

		//temperature
		parameterFile >> param;
		if (param != ".") temp[day][hour] = atof((param).c_str());
		else temp[day][hour] = nodata;

		hour = (hour+1)%hoursPerDay;
		if (hour == 0) day += 1;
	}

  	parameterFile.close();
}


/*******************************************************************************************
 * calculate daily mean temperature, radiation etc
 * if there is any nodata value for one variable during one day, the whole
 * day is flagged with nodata;
 * this should be taken into account in the other routines
 * called in --> controller::initialiseSimulation
 *******************************************************************************************/
void Weather::meanWeather(Parameter* p){

	double temp_help, prec_help, relHum_help, windSpeed_help, solRad_help;
	double moistBB_help, moistUB_help;
	bool nodataflag_temp, nodataflag_prec, nodataflag_relHum, nodataflag_windSpeed, nodataflag_solRad;
	bool nodataflag_moistBB, nodataflag_moistUB;
	double distEarthSun, solDeclination, sunSetHour;
    int yrcount = 0; //_dirk to count years

    for (int i = 0; i < p->simYears; i++) {
        precSumYear[i]=0;
    }

	for (int day = 0; day < daysPerYear*p->simYears; day++) {

        temp_help = prec_help = relHum_help = windSpeed_help = solRad_help = 0;
		moistBB_help = moistUB_help = 0;
		//flags for nodata values: true if there is a nodata value
		nodataflag_temp = nodataflag_moistBB = nodataflag_moistUB = nodataflag_prec = false;
		nodataflag_relHum = nodataflag_windSpeed = nodataflag_solRad = false;

		//calculation of evaporation equivalent of extraterrestrial evaporation
		//routine found in Maidment "Handbook of Hydrology", page 4.31
        //northern hemispere: 1st of January   (day 0)
        //southern hemisphere: 1st of July
		distEarthSun = 1.0 + 0.033*cos(2.0*PI/(daysPerYear*1.0)*(day%daysPerYear+1));
		solDeclination = 0.4093*sin(2.0*PI/(daysPerYear*1.0)*(day%daysPerYear+1) - 1.405);
		sunSetHour = acos(-1.0*tan(p->latitude*PI/180.0)*tan(solDeclination));
		extraterrRadiation[day] = 15.392*distEarthSun*(sunSetHour*sin(p->latitude*PI/180.0)*sin(solDeclination)+cos(p->latitude*PI/180.0)*cos(solDeclination)*sin(sunSetHour));
		tempMin[day] = temp[day][0];	 //just set the first value
		tempMax[day] = temp[day][0];

		for (int hour = 0; hour < hoursPerDay; hour++) {  //calculated daily sum...
			if (temp[day][hour]!=nodata)		//if everything is right
				temp_help += temp[day][hour];
			else								//if there is a nodata value
				nodataflag_temp = true;
			if (prec[day][hour]!=nodata)
				prec_help += prec[day][hour];
			else
				nodataflag_prec = true;
			if (moistBB[day][hour]!=nodata)
				moistBB_help += moistBB[day][hour];
			else
				nodataflag_moistBB = true;
			if (moistUB[day][hour]!=nodata)
				moistUB_help += moistUB[day][hour];
			else
				nodataflag_moistUB = true;
			if (relHum[day][hour]!=nodata)
				relHum_help += relHum[day][hour];
			else
				nodataflag_relHum = true;
			if (windSpeed[day][hour]!=nodata)
				windSpeed_help += windSpeed[day][hour];
			else
				nodataflag_windSpeed = true;
			if (solRad[day][hour]!=nodata)
				solRad_help += solRad[day][hour];
			else
				nodataflag_solRad = true;
			if (temp[day][hour] < tempMin[day])
				tempMin[day] = temp[day][hour];
			if (temp[day][hour] > tempMax[day])
				tempMax[day] = temp[day][hour];
		}

		if (nodataflag_temp == false)	//and divide it by 24 if necessary
			tempAvg[day] = temp_help/24.0;
		else
			tempAvg[day] = nodata;
		if (nodataflag_prec == false) {
			precSum[day] = prec_help;
            precSumYear[yrcount] += prec_help; //_dirk annual rainfall
        }
		else {
            precSum[day] = nodata;
            precSumYear[yrcount] += 0;
        }

		if (nodataflag_moistBB == false)
			moistBBAvg[day] = moistBB_help/24.0;
		else
			moistBBAvg[day] = nodata;
		if (nodataflag_moistUB == false)
			moistUBAvg[day] = moistUB_help/24.0;
		else
			moistUBAvg[day] = nodata;
		if (nodataflag_relHum == false)
			relHumAvg[day] = relHum_help/24.0;
		else
			relHumAvg[day] = nodata;
		if (nodataflag_windSpeed == false)
			windSpeedAvg[day] = windSpeed_help/24.0;
		else
			windSpeedAvg[day] = nodata;
		if (nodataflag_solRad == false)
			solRadSum[day] = solRad_help;
		else
			solRadSum[day] = nodata;

        if (day%daysPerYear == daysPerYear-1) {
            yrcount++;
        }
	}
}

/*******************************************************************************************
 * constructor for weather
 * is called in --> controller::initialise Simulation
 *******************************************************************************************/
Weather::Weather(Parameter* p){

    int totalSimDays = daysPerYear * p->simYears;

    //initialize vectors
    prec.resize(totalSimDays, Double1D(hoursPerDay, 0));
    temp.resize(totalSimDays, Double1D(hoursPerDay, 0));
    relHum.resize(totalSimDays, Double1D(hoursPerDay, 0));
    solRad.resize(totalSimDays, Double1D(hoursPerDay, 0));
    windSpeed.resize(totalSimDays, Double1D(hoursPerDay, 0));
    tempBB.resize(totalSimDays, Double1D(hoursPerDay, 0));
    tempUB.resize(totalSimDays, Double1D(hoursPerDay, 0));
    moistBB.resize(totalSimDays, Double1D(hoursPerDay, 0));
    moistUB.resize(totalSimDays, Double1D(hoursPerDay, 0));
    tempAvg.resize(totalSimDays, 0);
    tempMax.resize(totalSimDays, 0);
    tempMin.resize(totalSimDays, 0);
    precSum.resize(totalSimDays, 0);
    precSumYear.resize(p->simYears, 0);
    moistBBAvg.resize(totalSimDays, 0);
    moistUBAvg.resize(totalSimDays, 0);
    relHumAvg.resize(totalSimDays, 0);
    windSpeedAvg.resize(totalSimDays, 0);
    solRadSum.resize(totalSimDays, 0);
    extraterrRadiation.resize(totalSimDays, 0);
    growStart.resize(p->simYears, 0);
    growEnd.resize(p->simYears, 0);
}

/*******************************************************************************************
 * destructor for weather
 *******************************************************************************************/
Weather::~Weather(){

}

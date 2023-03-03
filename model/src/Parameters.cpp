/*******************************************************************************************
 * Parameters.cpp
 * The parameter object is responsible for communication.
 * Here, general model parameters are stored (in parameter.h)
 *******************************************************************************************/

#include <string>    //to use strings
#include <iostream>	 //to read and write files
#include <sstream>	 //to convert integers to strings...
#include <cstdlib>   //to use exit()
#include <fstream>   //to read and write files

#include "Parameters.h"

using namespace std;

/*******************************************************************************************
 * read in general model parameters from file
 * is called in --> controller::initialiseSimulation
 *******************************************************************************************/
void Parameter::readInModelParameters() {

	ifstream parameterFile;    //input file
	string modelParametersName;

    modelParametersName = "Parameters" + delimiter + "modelparameters_" + modelFileID[scenario] + ".txt";

    cout << "Read general model parameters from: " << modelParametersName << endl;
    cout << endl;

	parameterFile.open(modelParametersName.c_str(), ios::binary|ios::in);
	parameterFile.seekg(0); //position 0

    //if something goes wrong...
	if(!parameterFile){
        cerr << modelParametersName
		 << " could not be opened!\n";
		exit(-1);
    }

	int fsize = 10000;
	parameterFile.ignore(fsize, '#');
	parameterFile.ignore(fsize, ':');               //parses until end of file is reached (INT_MAX) OR until next ":" is found
	parameterFile >> site;
	parameterFile.ignore(fsize, ':');
	parameterFile >> latitude;
	parameterFile.ignore(fsize, ':');
	parameterFile >> simYears;
	parameterFile.ignore(fsize, ':');
	parameterFile >> xsize;
    parameterFile.ignore(fsize, ':');
	parameterFile >> ysize;
    parameterFile.ignore(fsize, ':');
	parameterFile >> cellsize;
    parameterFile.ignore(fsize, ':');
	parameterFile >> vegTimeStep;

	parameterFile.close();

	cout << "Site:\t\t\t\t" << site << endl;
	cout << "Latitude:\t\t\t" << latitude << " deg" << endl;
	cout << "Simulation time:\t\t" << simYears << " years" << endl;
	cout << "Grid scale:\t\t\t" << xsize << " x " << ysize << endl;
	cout << "Cell size:\t\t\t" << cellsize << " m" << endl;
	cout << "Vegetation time steps:\t\t" << vegTimeStep << " days" << endl;
	cout << endl;
}

/*******************************************************************************************
 * read in scenarios for climate, soil and vegetation and number of repetitions from file
 * is called in --> controller::initialiseModel
 *******************************************************************************************/
void Parameter::readInModelScenarios() {

	ifstream parameterFile;    //input file

    cout << "__________________________ Scenarios __________________________" << endl;
    cout << "Read model scenarios from: " << modelScenariosName << endl;
    cout << endl;

	parameterFile.open(modelScenariosName.c_str(), ios::binary|ios::in);
	parameterFile.seekg(0); //position 0

    //if something goes wrong...
	if(!parameterFile){
        cerr << modelScenariosName
		<< " could not be opened!\n";
		exit(-1);
    }

	int fsize = 10000;
	parameterFile.ignore(fsize, '#');
	parameterFile.ignore(fsize, ':');               //parses until end of file is reached (INT_MAX) OR until next ":" is found
	parameterFile >> scenarios;
    parameterFile.ignore(fsize, ':');
	parameterFile >> climateRepetitions;
	parameterFile.ignore(fsize, ':');
	parameterFile >> modelRepetitions;

	if (scenarios > 0) {
        outputFileID.resize(scenarios, "NA");
        modelFileID.resize(scenarios, "NA");
        weatherFileID.resize(scenarios, "NA");
        elevationFileID.resize(scenarios, "NA");
        soilFileID.resize(scenarios, "NA");
        vegetationFileID.resize(scenarios, "NA");
	}

    for (int i = 0; i < scenarios; i++) {
        parameterFile.ignore(fsize, '*');
        parameterFile.ignore(fsize, ':');
        parameterFile >> outputFileID[i];
        parameterFile.ignore(fsize, ':');
        parameterFile >> modelFileID[i];
        parameterFile.ignore(fsize, ':');
        parameterFile >> weatherFileID[i];
        parameterFile.ignore(fsize, ':');
        parameterFile >> elevationFileID[i];
        parameterFile.ignore(fsize, ':');
        parameterFile >> soilFileID[i];
        parameterFile.ignore(fsize, ':');
        parameterFile >> vegetationFileID[i];
    }

	parameterFile.close();

	cout << "Total scenarios:\t\t" << scenarios << endl;
	if (scenarios <= 0) {
        cerr << endl << ">>> Error: There must be at least one scenario!\n";
		exit(-1);
	}

    cout << "Climate repetitions:\t\t" << climateRepetitions << endl;
	if (climateRepetitions <= 0) {
        cerr << endl << ">>> Error: There must be at least one repetition!\n";
		exit(-1);
	}

    cout << "Model repetitions:\t\t" << modelRepetitions << endl;
	if (modelRepetitions <= 0) {
        cerr << endl << ">>> Error: There must be at least one repetition!\n";
		exit(-1);
	}
    cout << "Total model runs:\t\t" << scenarios * climateRepetitions * modelRepetitions << endl << endl;

    for (int i = 0; i < scenarios; i++) {
        cout << "--------------------" << endl;
        cout << "Scenario " << i + 1 << endl;
        cout << "--------------------" << endl;
        cout << "OutputFileID:\t\t\t" << outputFileID[i] << endl;
        cout << "ModelParametersID:\t\t" << modelFileID[i] << endl;
        cout << "WeatherFileID:\t\t\t" << weatherFileID[i] << endl;
        cout << "ElevationFileID:\t\t" << elevationFileID[i] << endl;
        cout << "SoilParametersID:\t\t" << soilFileID[i] << endl;
        cout << "VegetationParametersID:\t\t" << vegetationFileID[i] << endl;
    }
    cout << endl;
}

/*******************************************************************************************
 * constructor for parameters
 * is called in --> controller::initialiseSimulation
 *******************************************************************************************/
Parameter::Parameter() {

}

/*******************************************************************************************
 * destructor for parameters
 *******************************************************************************************/
Parameter::~Parameter() {

}

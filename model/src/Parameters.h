/****************************************************************************************//**
 * \class   Parameter Parameters.h
 * \brief   The parameter object is responsible for communication.
 *          Here, general model parameters are stored (in parameter.h)
 *          The only functions are string conversion functions.
 * \author  Britta Tietjen
 *******************************************************************************************/

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>    //to use strings
#include <fstream>	 //to read and write files
#include <vector>

using namespace std;

//General and constant variables, which should be visible everywhere
const int daysPerYear = 365;
const int hoursPerDay = 24;
const int nodata = -99999;

//To define vector types
typedef vector<double> Double1D;
typedef vector< vector<double> > Double2D;
typedef vector< vector <vector<double> > > Double3D;
typedef vector< vector <vector <vector<double> > > > Double4D;
typedef vector<int> Int1D;
typedef vector< vector<int> > Int2D;
typedef vector<bool> Bool1D;
typedef vector< vector<bool> > Bool2D;
typedef vector<string> String1D;

//To define delimiter depending on operating system
const string delimiter =
    #ifdef _WIN32
        "\\";        //for windows
    #else
        "/";         //for linux
    #endif

//Parameter paths which are not dependent on general model parameters
const string modelScenariosName = "Parameters" + delimiter + "modelscenarios.txt";                          //!< Path name of scenarios and repetitions

class Parameter{
	private:
	public:
	    //Constructor and Destructor
        /***************************************************************************************//**
         * \brief Constructor
         *******************************************************************************************/
		Parameter();

        /***************************************************************************************//**
         * \brief Destructor
         *******************************************************************************************/
		~Parameter();

		//Public model parameters
		string site;                        //!< site name
		double latitude;                    //!< latitude
		int simYears;                       //!< total simulation time
		int xsize;                          //!< x cells of the grid
        int ysize;                          //!< y cells of the grid
        int cellsize;                       //!< size of a single cell (length of a side) [m]
        int vegTimeStep;                    //!< number of days after which vegetation is updated

        //Public scenario parameters
        int modelRepetitions;               //!< number of total repetitions for each scenario
        int modelRepetition;                //!< current repetition of scenario
        int climateRepetitions;             //!< number of total repetitions for each scenario
        int climateRepetition;              //!< current repetition of scenario
        int scenarios;                      //!< total number of scenarios
        int scenario;                       //!< current number of scenario
        String1D outputFileID;		        //!< ID of output files for each scenario
        String1D modelFileID;		        //!< ID of model parameters file for each scenario
        String1D weatherFileID;		        //!< ID of weather file for each scenario
        String1D elevationFileID;		    //!< ID of elevation file for each scenario
        String1D soilFileID;		        //!< ID of soil parameters file for each scenario
        String1D vegetationFileID;		    //!< ID of vegetation parameters file for each scenario

        //Public member functions
        /***************************************************************************************//**
         * \brief read in model specific parameters from file
         *******************************************************************************************/
		void readInModelParameters();

        /***************************************************************************************//**
         * \brief read in scenarios for climate, soil and vegetation and number of repetitions from file
         *******************************************************************************************/
		void readInModelScenarios();
};

#endif /*PARAMETERS_H_*/


/*******************************************************************************************
 * PlantFunctionalType.cpp
 * Methods defined here are universal for all plant functional types.
 *******************************************************************************************/

#include "PlantFunctionalType.h"

using namespace std;

/*******************************************************************************************
 *  Output parameters
 *******************************************************************************************/
void PlantFunctionalType::output() {

	cout << endl;
	cout << "____________Plant Functional Type Parameters____________" << endl;
	cout << "Name: " << Name << endl;
	cout << "Type: " << Type << endl;
	cout << "UptakeRate: " << UptakeRate << endl;
	cout << "RootL1:" << RootL1 << endl;
	cout << "RootL2:" << RootL2 << endl;
	cout << "GrowthR: " << GrowthR << endl;
	cout << "Trancoef: " << Trancoef << endl;
	cout << "Kcoef: " << Kcoef << endl;
	cout << "ecoef: " << ecoef << endl;
	cout << "temcoef: " << temcoef << endl;
	cout << "LAI_C: " << LAI_C << endl;
	cout << "MortR: " << MortR << endl;
	cout << "BasicMortRate:" << BasicMortRate << endl;
	cout << "Grazed: " << Grazed << endl;
	cout << "WPWaterL1: " << WPWaterL1 << endl;
	cout << "WPWaterL2: " << WPWaterL2 << endl;
	cout << "cmax: " << cmax << endl;
	cout << "overlap: " << overlap << endl;
	cout << "conversionC_BM: " << conversionC_BM << endl;
	cout << "grazeParam: " << grazeParam << endl;
	cout << "grazeLim: " << grazeLim << endl;
	cout << "grazeprefer: " << grazeprefer << endl;
	cout << "disFrac: " << disFrac << endl;
	cout << "estW: " << estW << endl;
	cout << "bm_c_rain: " << bm_c_rain << endl;
	cout << "MAP: " << MAP << endl;
	cout << "grazeParA: " << grazeParA << endl;
	cout << "grazeParB: " << grazeParB << endl;
	cout << "specLoss: " << specLoss << endl;
	cout << "firecoeff: " << firecoeff << endl;
	cout << "fuelLim: " << fuelLim << endl;
	cout << "dist0S: " << dist0S << endl;
	cout << "distConstS: " << distConstS << endl;
	cout << "distConstS2: " << distConstS2 << endl;
	cout << "________________________________________________________" << endl;
	system("pause");
}

/*******************************************************************************************
 *  Calculate average moisture for layer one
 *******************************************************************************************/
void PlantFunctionalType::calculateAvMoistL1(VegetationCell* vegetation) {

	AvMoistL1 = max((vegetation->moistL1 - WPWaterL1) * vegetation->soilDepthL1,0.0);
}

/*******************************************************************************************
 *  Add growth to cover
 *******************************************************************************************/
void PlantFunctionalType::addGrowthToCover() {

	Cover += Growth;
}

/*******************************************************************************************
 *  calculate the corresponding biomass on a cell by a linear regression
 *******************************************************************************************/
void PlantFunctionalType::calculateBiomass(VegetationCell* vegetation, Weather *we2, int yr2) {

	double m = (1.0 - vegetation->bm_c_rain) / vegetation->MAP;
	Biomass = (conversionC_BM * Cover * ((we2->precSumYear[yr2]* m) + vegetation->bm_c_rain)) * (cellsize*cellsize)/(100.0*100.0); //biomass in g/cell
}

/*******************************************************************************************
 *  Calculate total growth
 *******************************************************************************************/
double PlantFunctionalType::calculateTotalGrowth() {

	totalGrowth += Growth;
	return totalGrowth;
}

/*******************************************************************************************
 *  Calculate specific uptake rate
 *******************************************************************************************/
double PlantFunctionalType::calculateUptakeRateCover(VegetationCell* vegetation, Weather* weather, int year) {

	double m = (1.0 - vegetation->bm_c_rain) / vegetation->MAP;
	UptakeRateCover = UptakeRate * (conversionC_BM * ((weather->precSumYear[year] * m) + vegetation->bm_c_rain));
	return UptakeRateCover;
}

/******************************************************************************************
 *  calculate the corresponding cover on
 *  a cell by a linear regression
 ******************************************************************************************/
void PlantFunctionalType::calculateCover(VegetationCell* vegetation, Weather* we2, int yr2, double factor) {

    // ggG: relative damage caused to plant cover due to grazing, which increases linearly with an increasing fraction of biomass that was removed by the grazers
	// the parameter grazeParA and grazeParB define the shape of this quadratic function of cover reduction
	//new version of _dirk03112010   to reduce impact of heavy grazing a little (reduced max damage)
	double gg = factor * grazeParA + grazeParB;
	double gg1 = Biomassold - (Biomassold * factor *  gg); // Biomassold * factor = Biomass eaten ;
	double m = (1.0 - vegetation->bm_c_rain) / vegetation->MAP;
	Cover = ((double) gg1 /(conversionC_BM* ((we2->precSumYear[yr2]* m) + vegetation->bm_c_rain))) * (100.0*100.0 / (cellsize*cellsize*1.0));
}

/*******************************************************************************************
 *  Constructor
 *******************************************************************************************/
PlantFunctionalType::PlantFunctionalType(string _Name, string _Type, double _UptakeRate, double _RootL1, double _RootL2, double _GrowthR, double _Transcoef, double _Kcoef, double _ecoef, double _Tcoef, double _LAI_C, double _MortR, double _BasicMortRate, double _Grazed, double _WPWaterL1, double _WPWaterL2, double _cmax, double _overlap, double _conversionC_BM, double _grazeParam, double _grazeLim, double _grazeprefer, double _grazepreferintra, double _disFrac, double _estW, double _bm_c_rain, double _MAP, double _grazeParA,
		                                 double _grazeParB, double _specLoss, double _firecoeff, int _fuelLim, double _dist0S, double _distConstS, double _distConstS2, double _reserveLast, double _reserveThis, int _vegTimeStep, int _cellsize) {

    /// \todo we should get rid of this parameter mess!
	vegTimeStep = _vegTimeStep;
	cellsize = _cellsize;
	Name = _Name;
	Type = _Type;
	UptakeRate = _UptakeRate;
	RootL1 = _RootL1;
	RootL2 = _RootL2;
	GrowthR =_GrowthR;
	Trancoef = _Transcoef;
	Kcoef = _Kcoef;
	ecoef = _ecoef;
	temcoef = _Tcoef;
	LAI_C   =_LAI_C;
	MortR = _MortR;
	BasicMortRate = _BasicMortRate;
	Grazed = _Grazed;
	WPWaterL1 = _WPWaterL1;
	WPWaterL2 = _WPWaterL2;
	cmax = _cmax;
	overlap = _overlap;
	conversionC_BM = _conversionC_BM;
	grazeParam = _grazeParam;
	grazeLim = _grazeLim;
	grazeprefer = _grazeprefer;
	grazepreferintra = _grazepreferintra;
	disFrac = _disFrac;
	estW = _estW;
	bm_c_rain = _bm_c_rain;
	MAP = _MAP;
	grazeParA = _grazeParA;
	grazeParB = _grazeParB;
	specLoss = _specLoss;
	firecoeff = _firecoeff;
	fuelLim = _fuelLim;
	dist0S = _dist0S;
	distConstS = _distConstS;
	distConstS2 = _distConstS2;
	reserveLast = _reserveLast;
	reserveThis = _reserveThis;

	Growth = 0;
	BasicMort = 0;
	DryMort = 0;
	DryMortL1 = 0;
	DryMortL2 = 0;
	Disp = 0;
	GrowthL1 = 0;
	GrowthL2 = 0;
	AvMoistL1 = 0;
	AvMoistL2 = 0;
	AvMoistSuppL1 = 0;
	AvMoistSuppL2 = 0;
	Biomass = 0;
	Biomassold = 0;
	resBM = 0;
	Cover = 0;
	basicM = 0;
	totalGrowth = 0;
}

/*******************************************************************************************
 *  Destructor
 *******************************************************************************************/
PlantFunctionalType::~PlantFunctionalType() {

}

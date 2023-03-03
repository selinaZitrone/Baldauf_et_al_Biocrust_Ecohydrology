/*******************************************************************************************
 * VegetationCell.cpp
 * Describes the vegetation dynamics in one grid cell.
 *******************************************************************************************/

#include "VegetationCell.h"
#include "PlantFunctionalType.h"

using namespace std;

/*******************************************************************************************
 * calculate total cover of shrubs in one cell
 * is called in -->
 *******************************************************************************************/
double VegetationCell::calculateShrubCover(){

    shrubCover = 0;
	for (int i = 0; i < Shrubs; i++) {
		shrubCover += pftList[i]->Cover;
	}
	return shrubCover;
}

/*******************************************************************************************
 * calculate total cover of perennials in one cell
 * is called in -->
 *******************************************************************************************/
double VegetationCell::calculatePerennialCover(){

	perennialCover = 0;
	for (int i = Shrubs; i < Shrubs + Perennials; i++) {
		perennialCover += pftList[i]->Cover;
	}
	return perennialCover;
}

/*******************************************************************************************
 * calculate total cover of annuals in one cell
 * is called in -->
 *******************************************************************************************/
double VegetationCell::calculateAnnualCover(){

	annualCover = 0;
	for (int i = Shrubs + Perennials; i < Shrubs + Perennials + Annuals; i++) {
		annualCover += pftList[i]->Cover;
	}
	return annualCover;
}

/*******************************************************************************************
 * calculate total cover of all PFTs in one cell
 * is called in -->
 *******************************************************************************************/
double VegetationCell::calculateTotalCover(){

    //-------------------------------------------------------
    /// CALCULATE SHRUB COVER.
    /// \implements \fn VegetationCell::calculateShrubCover()
    //-------------------------------------------------------
	calculateShrubCover();

    //-------------------------------------------------------
    /// CALCULATE PERENNIAL COVER.
    /// \implements \fn VegetationCell::calculatePerennialCover()
    //-------------------------------------------------------
	calculatePerennialCover();

    //-------------------------------------------------------
    /// CALCULATE ANNUAL COVER.
    /// \implements \fn VegetationCell::calculateAnnualCover()
    //-------------------------------------------------------
	calculateAnnualCover();

	totalCover =  shrubCover + perennialCover + annualCover;
	return totalCover;
}

/*******************************************************************************************
 * calculate mean soil moisture during water season
 * is called in --> VegetationCell::setMoist
 *******************************************************************************************/
void VegetationCell::meanGrowingSeasonMoist(double mL1, double mL2, int day, int grStart, int grEnd){/*Tong, based on PFT-oriented*/

    if (day == grStart) {
        moistSeasonL1 = 0;
        moistSeasonL2 = 0;
        moistL1mem[0] = moistL1mem[1] ; //moisture memory for shrub establishment
        moistL1mem[1] = moistL1mem[2] ;
    }

    moistSeasonL1 += mL1;
    moistSeasonL2 += mL2;

    if (day == grEnd){
    /// \todo Btodo: this code should go into PFT establishment functions!
        moistSeasonL1 /= ((grEnd-grStart+1)*1.0);
        moistSeasonL2 /= ((grEnd-grStart+1)*1.0);
        moistL1mem[2] = moistSeasonL1;    //_rule

        double Scover = 0;
        for (int i = 0; i < Shrubs; i++){
             Scover +=  pftList[i]->Cover;
        }

        //Here we determine whether Seeds of Acacia could germinate.
        //This doesn't mean that they establish a year later, but in reality in this season you would see seedlings emerge, so for management modelling we include this variable
        for (int i = 0; i < Shrubs; i++){
            if (moistSeasonL1 > pftList[i]->estW *  pftList[i]->WPWaterL1) {
                pftList[i]->ger = 1;
            } else pftList[i]->ger = 0;

            if (moistL1mem[1] > pftList[i]->estW * pftList[i]->WPWaterL1 && moistL1mem[2] >  pftList[i]->estW * pftList[i]->WPWaterL1 ){
                if (Scover < EncroachScover) {
                    cohortAge[i] = 0; // if cell is not already full of shrubs establishment can take place
                    pftList[i]->est = 1;
                }
            } else pftList[i]->est = 0;

			pftList[i]->est3 = pftList[i]->est2;
            //S1estS2 for history of estS (fire)
            if (moistL1mem[1] > pftList[i]->estW * pftList[i]->WPWaterL1 && moistL1mem[0] > pftList[i]->estW * pftList[i]->WPWaterL1 ) {
                pftList[i]->est2 = 1;
            } else SestS2=0;
        }

        //estG and estS: is there enough moisture for establishment of plants?
        for (int i = Shrubs; i < Shrubs + Perennials; i++){
            if (moistSeasonL1 > pftList[i]->estW * pftList[i]->WPWaterL1) {
                pftList[i]->est = 1;
            } else pftList[i]->est = 0;
        }
    }
}

/*******************************************************************************************
 * \todo include brief description here
 * is called in --> controller::runSimulation
 *******************************************************************************************/
void VegetationCell::setEvaporationandTranspiration(double Evaporation, double TranspirationL1, double TranspirationL2, int day){/*Tong, transmission of evaporation and transpiration*/

    if (day == 0){
		yearlyE = 0;
		yearlyTL1 = 0;
		yearlyTL2 = 0;
	}

	yearlyE += Evaporation;
	yearlyTL1 += TranspirationL1;
	yearlyTL2 += TranspirationL2;
}

/*******************************************************************************************
 * calculate mean soil moisture during wet season
 * is called in --> VegetationCell::setMoist
 *******************************************************************************************/
void VegetationCell::meanWetSeasonMoist(double mL1, double mL2, int day, int wStart, int wEnd) {

	if (day == wStart) {
		moistWetSeasonL1 = 0;
		moistWetSeasonL2 = 0;
	}

	moistWetSeasonL1 += mL1;
	moistWetSeasonL2 += mL2;

	if (day == wEnd) {
		moistWetSeasonL1 /= ((wEnd-wStart+1)*1.0);
		moistWetSeasonL2 /= ((wEnd-wStart+1)*1.0);
	}
}

/*******************************************************************************************
 * set soil moisture of both layers
 * is called in --> controller::runSimulation
 *******************************************************************************************/
void VegetationCell::setMoist(double mL1, double mL2, int day, int grStart, int grEnd, int wetStart, int wetEnd) {

    //new counting of meanMoist starts at day 1
    //(evaluation is at day 0)
    if (day%vegTimeStep == 0){
        moistL1 = 0;
        moistL2 = 0;
    }

    //in each year, moisture is added
	moistL1 += mL1;
	moistL2 += mL2;

    //at day 0, the mean moisture is calculated
    if (day%vegTimeStep == vegTimeStep-1){
        moistL1 /= (vegTimeStep*1.0);
        moistL2 /= (vegTimeStep*1.0);
    }
    //or at the end of one simulation year
    else if (day == (daysPerYear-1)){
        moistL1 /= (((daysPerYear%vegTimeStep)+1)*1.0);
        moistL2 /= (((daysPerYear%vegTimeStep)+1)*1.0);
    }

    //---------------------------------------------
    /// MEAN GROWING SEASON MOIST.
    //---------------------------------------------
    if ((day >= grStart) && (day <= grEnd)) {
        /// \implements \fn VegetationCell::meanGrowingSeasonMoist(double mL1, double mL2, int day, int grStart, int grEnd)
        meanGrowingSeasonMoist(mL1, mL2, day, grStart, grEnd);
    }

    //---------------------------------------------
    /// MEAN WATER SEASON MOIST.
    //---------------------------------------------
    if ((day >= wetStart) && (day <= wetEnd)) {
        /// \implements \fn VegetationCell::meanWaterSeasonMoist(double mL1, double mL2, int day, int wetStart, int wetEnd)
        meanWetSeasonMoist(mL1, mL2, day, wetStart, wetEnd);
    }
}

/*******************************************************************************************
 * \todo include brief description here
 * is called in --> VegetationLandscape::Grazing
 *******************************************************************************************/
double VegetationCell::calculateResAnnualcover(){/*Tong, calculate the total reserve cover of annual grass*/

	double conversionC_BM_Annual = 0;

	for (int i = Shrubs + Perennials; i < PFTs; i++){
		conversionC_BM_Annual += pftList[i]->conversionC_BM;
	}

    conversionC_BM_Annual /= Annuals;
	resAnnualcover = (double) (AnnualresBM / conversionC_BM_Annual )  * (100.0*100.0 / (cellsize*cellsize*1.0));
	return resAnnualcover;
}

/*******************************************************************************************
 * \todo include brief description here
 * is called in --> VegetationLandscape::Grazing
 *******************************************************************************************/
double VegetationCell::calculateResPerennialcover(){/*Tong, calculate the total reserve cover of perennial grass*/

	double conversionC_BM_Perennial = 0;

	for (int i = Shrubs; i < Shrubs + Perennials; i++){
		conversionC_BM_Perennial += pftList[i]->conversionC_BM;
	}

	conversionC_BM_Perennial /= Perennials;
	resPerennialcover = (double) (PerennialresBM / conversionC_BM_Perennial )  * (100.0*100.0 / (cellsize*cellsize*1.0));
	return resPerennialcover;
}

/*******************************************************************************************
 * \todo include brief description here
 * is called in --> VegetationLandscape::Grazing
 *******************************************************************************************/
double VegetationCell::calculateResCover() {/*Tong, calculate the total reserve cover of annual grass and perennial grass*/

	resCover = resAnnualcover + resPerennialcover;

	return resCover;
}

/*******************************************************************************************
 * Calculate water uptake per unit cover
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void VegetationCell::setCover(){ /*Tong*/

    for (int i = 0; i < PFTs; i++){
        if (i < Shrubs) Scover[i] = pftList[i]->Cover;
        else if (i >= Shrubs && i < Shrubs + Perennials) Pcover[i-Shrubs] = pftList[i]->Cover;
        else Acover[i-(Shrubs+Perennials)] = pftList[i]->Cover;
	}
}

/*******************************************************************************************
 * \todo description here!
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void VegetationCell::setpftUse(){   /*Tong*/

	for (int i = 0; i < PFTs; i++){
        UseL1[i] = pftList[i]->UseL1;
        UseL2[i] = pftList[i]->UseL2;
	}
}

/*******************************************************************************************
 * set roots layer two
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void VegetationCell::setRootL2(){   /*Tong*/

	for (int i = 0; i < PFTs; i++){
		RootL2[i] = pftList[i]->RootL2;
	}
}

/*******************************************************************************************
 * set roots layer one
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void VegetationCell::setRootL1(){  /*Tong*/

	for (int i = 0; i < PFTs; i++){
		RootL1[i] = pftList[i]->RootL1;
	}
}

/*******************************************************************************************
 * calculate the water uptake per unit cover
 * \todo needs to be reviewed!!!
 * is called in --> VegetationLandscape::calculateProcesses
 *******************************************************************************************/
void VegetationCell::calculateWaterUptake() {  /*Tong, based on PFT-oriented*/

	double Perennial_URL1_scaled = 0;
	double Perennial_URL2_scaled = 0;
	double Annual_URL1_scaled = 0;
	double Shrub_URL1_scaled = 0;
	double Shrub_URL2_scaled = 0;

    // calculate the fraction of the available water for each plant functional type, it depends on the potential uptake rate per cover, root fraction and vegetation cover for each PFT
	for (int i = 0; i < PFTs; i++){
		pftList[i]->URL1_help = pftList[i]->UptakeRateCover * pftList[i]->RootL1;
		pftList[i]->URL1_scaled = pftList[i]->URL1_help * pftList[i]->Cover;
		pftList[i]->URL2_help = pftList[i]->UptakeRateCover * pftList[i]->RootL2;
        pftList[i]->URL2_scaled = pftList[i]->URL2_help * pftList[i]->Cover;

		if (i < Shrubs) {
            Shrub_URL1_scaled += pftList[i]->URL1_scaled;
            Shrub_URL2_scaled += pftList[i]->URL2_scaled;
		}
		else if (i >= Shrubs && i < Shrubs + Perennials) {
            Perennial_URL1_scaled += pftList[i]->URL1_scaled;
            Perennial_URL2_scaled += pftList[i]->URL2_scaled;
		}
		else {
            Annual_URL1_scaled += pftList[i]->URL1_scaled;
		}
	}

	for (int i = 0; i < PFTs; i++){
		if (pftList[i]->Cover > 0) {
            if (i < Shrubs + Perennials) {
                pftList[i]->UseL1 = pftList[i]->URL1_help / (Perennial_URL1_scaled + Annual_URL1_scaled + Shrub_URL1_scaled);
                pftList[i]->UseL2 = pftList[i]->URL2_help / (Perennial_URL2_scaled + Shrub_URL2_scaled);
            }
            else {
                pftList[i]->UseL1 = pftList[i]->URL1_help / (Perennial_URL1_scaled + Annual_URL1_scaled + Shrub_URL1_scaled);
            }
		}
		else {
            if (i < Shrubs + Perennials) {
                pftList[i]->UseL1 = 0;
                pftList[i]->UseL2 = 0;
            }
            else {
                pftList[i]->UseL1 = 0;
            }
		}
	}
}

/*******************************************************************************************
 * return the transpiration of PFTs in the waterlandscape process
 * is called in --> Controller::runSimulation
 *******************************************************************************************/
void VegetationCell::setTransFactors(Double1D TransL1, Double1D TransL2, int day){ /*Tong, transmission of transpiration to vegetation cell*/

	if (day == 0){
        for (int i = 0; i < PFTs; i++){
            yearlyPftTL1[i] = 0;
            yearlyPftTL2[i] = 0;
        }
	}

    // 	transparent the transpiration calculated in soil cell to vegetation cell
	for (int i = 0; i < PFTs; i++){
		TL1[i] = TransL1[i];
		TL2[i] = TransL2[i];
	}

    // sum up the transpiration of each plant functional type
	for (int i = 0; i < PFTs; i++){
		 yearlyPftTL1[i] += TL1[i];
		 yearlyPftTL2[i] += TL2[i];
	}
}

/*******************************************************************************************
 * calculate cover reduction in each cell due to fire
 * only called, if a fire actually occurs
 * revised the fire process based on PFTs
 * is called in --> VegetationLandscape::Fire
 *******************************************************************************************/
void VegetationCell::calculatePFTfire(){

	double fuel = 0;
	for (int i = Shrubs; i < Perennials; i++) {
		fuel = resBM + pftList[i]->Biomass ;
	}

    double windspeed = 1.4;     //windspeed at Sandfeld in October: 1.43 m/s
    double relhum = 35.6;       //relhum at sandveld  35.6 %
    double FM = 19.6;           //MAX Fuelmoisture from Dave 19.6%
    double fuel2 = (100.0*100)* fuel /1000;
    fuel2 = fuel2 / (cellsize*cellsize) ;
    double fireIntens = 2729 + 0.8684 * fuel2 - (530*sqrt(FM)) - (0.1907*relhum*relhum) - (596/windspeed);   //fire intensity is:  fuel in kg/ha

    if (fireIntens > 300) {
        resBM = 0;
    }

    for (int i = 0; i < Shrubs + Perennials; i++) {
        /// \implements \fn PlantFunctionalType::fire(double fireIntens, double fuel, Int1D &cohortAge, Int1D &cohortAgeFire)
        pftList[i]->fire(fireIntens, fuel, cohortAge[i], cohortAgeFire[i]);
    }
}

/*******************************************************************************************
 * constructor for class vegetation
 * is called in --> VegetationLandscape::initialiseVegParams
 *******************************************************************************************/
VegetationCell::VegetationCell(int PFTs_, int Annuals_, int Perennials_, int Shrubs_, double resWater_, double moistL1_, double moistL2_, double overlap_, double grazeParam_,
                       double aGrowthConst_, double bm_c_rain_, double grazeparA_, double grazeparB_, double MAP_, double firecoeff_, int fuelLim_, double EnScover_,
                       double initialshare_, int vegTimeStep_, int cellsize_){//constructor has been changed by Tong

	vegTimeStep = vegTimeStep_;
	cellsize = cellsize_;
	PFTs = PFTs_;
	Annuals = Annuals_;
	Perennials = Perennials_;
	Shrubs = Shrubs_;

	resWater = resWater_;
	moistL1 = moistL1_;
	moistL2 = moistL2_;
	overlap = overlap_;
	aGrowthConst = aGrowthConst_;
	grazeParam = grazeParam_;
	bm_c_rain = bm_c_rain_;
	grazeparA = grazeparA_;
	grazeparB = grazeparB_;
	MAP = MAP_;
	firecoeff = firecoeff_ ;
	fuelLim = fuelLim_;
	EncroachScover = EnScover_;
	initialshare = initialshare_;

	moistL1mem [0] = 0;
	moistL1mem [1] = 0;
	moistL1mem [2] = 0;
	moistSeasonL1 = 0;
	moistSeasonL2 = 0;
	moistWetSeasonL1 = 0;
	moistWetSeasonL2 = 0;
	reserveBM = 0;

	safeBM_reserve.resize(PFTs, 0);
	safeBM_alive.resize(PFTs, 0);
	edibleBM_reserve.resize(PFTs, 0);
	edibleBM_alive.resize(PFTs, 0);
	edibleBM_total.resize(PFTs, 0);
	relprefer.resize(PFTs, 0);
	TL1.resize(PFTs, 0);
	TL2.resize(PFTs, 0);
    cohortAge.resize(Shrubs, 0);
    cohortAgeFire.resize(Shrubs, 0);
	yearlyPftTL1.resize(PFTs, 0);
	yearlyPftTL2.resize(PFTs, 0);
	Pcover.resize(Perennials, 0);
	Acover.resize(Annuals, 0);
	Scover.resize(Shrubs, 0);
	UseL1.resize(PFTs, 0);
	UseL2.resize(PFTs, 0);
	RootL1.resize(PFTs, 0);
	RootL2.resize(PFTs, 0);
	pftList.resize(PFTs);
}

/*******************************************************************************************
 * destructor for class vegetation
 *******************************************************************************************/
VegetationCell::~VegetationCell(){

    for(unsigned int i = 0; i < pftList.size(); i++)
    {
        delete pftList[i];
    }
}

/****************************************************************************************
 * \class   WaterCell WaterCell.h
 * \brief   Describes the water dynamics in one grid cell
 * \author  Britta Tietjen
 *******************************************************************************************/

#ifndef WATERCELL_H_
#define WATERCELL_H_

#include <fstream>  //to read and write files
#include <string>   //to use strings

#include "Parameters.h"
#include "Weather.h"

class WaterCell {
 private:
  // Private parameters
  double StotRootL1;  //!< total roots of shrubs in layer 1 in one cell
  double PtotRootL1;  //!< total roots of perennials in layer 1 in one cell
  double AtotRootL1;  //!< total roots of annuals in layer 1 in one cell
  double StotRootL2;  //!< total roots of shrubs in layer 2 in one cell
  double PtotRootL2;  //!< total roots of perennials in layer 2 in one cell
  double cStotal;     //!< total cover of shrubs
  double cGtotal;     //!< total cover of perennials
  double cAtotal;     //!< total cover of annuals

 public:
  // Constructor and Destructor
  /***************************************************************************************
   * \brief Constructor
   * \todo include params here
   *******************************************************************************************/
  WaterCell(int x, int y, double suction, double Ks, double waterL1_sat,
            double WP, double FC, double depthL1, double depthL2,
            double iniWaterL0abs, double iniWaterL1, double iniWaterL2,
            double EPfactor, double resWater, double maxFL2,
            double maxAmountFL2, double diffConst, double TE_, int PFTs_,
            int Perennials_, int Shrubs_, int Annuals_, double Vegcoef1_,
            double Vegcoef2_, double Veg_ETcoef1_,
            double Veg_ETcoef2_);  // Constructor has been changed by /*Tong*/

  /***************************************************************************************
   * \brief Default constructor
   *******************************************************************************************/
  WaterCell();

  /***************************************************************************************
   * \brief Destructor
   *******************************************************************************************/
  ~WaterCell();  // Destructor

  // Public parameters
  typedef vector<double> Double1D;

  int xpos;
  int ypos;

  double suction;    //!< matrix potential (suction?) [mm]
  double Ks;         //!< saturated hydraulic conductivity [mm/h]
  double Ku;         //!< unsaturated hydraulic conductivity [mm/h]
  double K_geom;     //!< geometric mean of saturated and unsaturated hydraulic
                     //!< conductivity [mm/h]
  double porosity;   //!< fraction of pores [mm³/mm³]
  double Na;         //!< space in pores, which can be filled [mm]
  double depthL1;    //!< depth of the upper soil layer [mm]
  double depthL2;    //!< depth of the lower soil layer [mm]
  double drainL1L2;  //!< drainage from layer1 to layer2 due to saturation [mm]
  double draindeep;  //!< deep drainage from layer2 due to saturation [mm]
  double WP;         //!< wilting point [vol%]
  double FC;         //!< field capacity [vol%]
  double waterL1_sat;   //!< saturated water soil water content [-]
  double inclination;   //!< inclination of cell to lowest neighbour (radians)
  int aspect;           //!< direction of slope inclination
  double aspectFactor;  //!< impact of aspect on radiation
  double EPfactor;  //!< additional factor to account for e.g. ET reduction by
                    //!< soil crusts etc.
  double resWater;  //!< residual soil water content
  double maxFL2;    //!< maximal fraction of surface water that infiltrates into
                    //!< lower layer
  double maxAmountFL2;  //!< maximal total amount of infiltrated water into the
                        //!< lower layer
  double diffConst;     //!< diffusion constant that influences the speed of
                        //!< vertical water movement
  double Vegcoef1;      //!< total vegetation cover coefficient in vegetation
                    //!< function of upper layer for evapotranspiration /*Tong*/
  double Vegcoef2;  //!< total vegetation cover coefficient in vegetation
                    //!< function of upper layer for evapotranspiration /*Tong*/
  double Veg_ETcoef1;  //!< the transforming coefficient of vegetation cover to
                       //!< evapotranspiration /*Tong*/
  double Veg_ETcoef2;  //!< the transforming coefficient of vegetation cover to
                       //!< evapotranspiration /*Tong*/

  int PFTs;
  int Perennials; /*Tong*/
  int Shrubs;     /*Tong*/
  int Annuals;    /*Tong*/

  Double1D RootL2;  //!< root fraction of each PFT in the lower layer /*Tong*/
  Double1D RootL1;  //!< root fraction of each PFT in the upper layer /*Tong*/

  double TE;  //!< TE constant for decreased evapotranspiration at increased CO2
              //!< levels

  // calculated soil parameters dependent on precipitation
  double ts;     //!< time until saturation/ponding [h]
  double Fs;     //!< infiltrated water until saturation/ponding [mm/h]
  double FL1;    //!< total infiltrated water in layer1 [mm/h]
  double FL2;    //!< total infiltrated water in layer2 [mm/h]
  double F;      //!< infiltrated water via wetting front [mm/h]
  double QD;     //!< total run-off [mm/h]
  double EPL0;   //!< evaporated water (surface) [mm/d]
  double EPL1;   //!< evaporated water (layer1) [mm/d]
  double EPL2;   //!< evaporated water (layer2) [mm/d]
  double EPtot;  //!< total evaporated water from all layers [mm/d]
  double TL1;    //!< transpiration in the upper layer [mm/d] /*Tong*/
  double TL2;    //!< transpiration in the lower layer [mm/d] /*Tong*/
  double EL1;    //!< evaporation in the upper layer [mm/d] /*Tong*/
  Double1D
      PftTL1;  //!< transpiration of each PFT in the upper layer [mm/d] /*Tong*/
  Double1D
      PftTL2;  //!< transpiration of each PFT in the lower layer [mm/d] /*Tong*/

  // water content in soil layers
  double waterL0_abs;  //!< actual water content above the soil [mm]
  double waterL1_abs;  //!< actual water content in layer 1 [mm]
  double waterL2_abs;  //!< actual water content in layer 2 [mm]
  double waterL1_rel;  //!< actual water content in layer 1 [%]
  double waterL2_rel;  //!< actual water content in layer 2 [%]

  // monthly and daily means for each cell
  double waterL1_rel_day, waterL2_rel_day = 0.0;
  double waterL1_rel_month, waterL2_rel_month = 0.0;

  double totalInfiltration;  //!< infiltration in all layers
  double totWater;
  Double1D coverPFT;
  Double1D UseL1;  //!< Water using rate of each PFT in the upper layer /*Tong*/
  Double1D UseL2;  //!< Water using rate of each PFT in the lower layer /*Tong*/
  double layerDiff;  //!< water that flows from the upper to the lower layer
                     //!< (other way round if negative)

  bool wettingFrontFlag;

  // Public member functions
  /***************************************************************************************
   * \brief Infiltration routine into two layers:
   *   1. layer (upper): infiltration until ponding (saturation),
   *   				   afterwards Green & Ampt infiltration -->
   *wetting front
   *   2. layer (lower): deep infiltration via root channels
   * \param meanWaterL1 mean water in the upper layer in one timestep;
   *******************************************************************************************/
  bool twoLayerInfiltration(double meanWaterL1);

  /***************************************************************************************
   * \brief Calculate actual evapotranspiration in two layers.
   * Here, also the aspect and the inclination are taken into account
   * which normally belong to the potential evapotranspiration.
   * \param rainTot ?;
   * \param potEvapo potential Evapotranspiration;
   * \param crust is this cell covered by crust? If yes, evaporation from L0 is
   *skipped because this is managed by crustCell
    \param ep_factor_cell EP factor of the cell from which the evaporation is
   calculated
   *******************************************************************************************/
  void actLayerEvapotranspiration(double rainTot, double potEvapo, bool crust,
                                  double ep_factor_cell);

  /***************************************************************************************
   * \brief Total soil water content
   *******************************************************************************************/
  void totalWater();

  /***************************************************************************************
   * \brief Set PFT cover.
   * \todo Add params here!
   *******************************************************************************************/
  void setCover(Double1D shrubcover, Double1D grasscover,
                Double1D annualcover); /*Tong*/

  /*******************************************************************************************
   * \brief Calculate total cover for each PFT.
   *******************************************************************************************/
  void totalCover();

  /***************************************************************************************
   * \brief Set vegetation use factors.
   * \todo Add params here!
   * \todo why these parameters again?
   *******************************************************************************************/
  void setVegUseFactors(Double1D PftUseL1, Double1D PftUseL2);

  /***************************************************************************************
   * \brief The aspect factor determines how the radiation is affected by the
   *aspect of the slope. The values for this factor were found in Shevenell L.,
   *1999. Regional potental evapotranspiration in arid climates based on
   *temperature, topography and calculated solar radiation. Hydrological
   *Processes 13, 577-596 This aspect factor is later on multiplied to the
   *radiation to calculate the potential evapotranspiration.
   *******************************************************************************************/
  void setAspectFactor();

  /***************************************************************************************
   * \brief Set roots layer two.
   * \author Tong Guo
   * \todo Add params here!
   * \todo why these parameters again?
   *******************************************************************************************/
  void setRootL2(Double1D RootL2_);

  /***************************************************************************************
   * \brief Set roots layer one.
   * \author Tong Guo
   * \todo Add params here!
   * \todo why these parameters again?
   *******************************************************************************************/
  void setRootL1(Double1D RootL1_);

  /***************************************************************************************
   * \brief calculate total roots in upper and lower soil layer for all PFTs per
   *cell
   *******************************************************************************************/
  void totalRoots();

  /***************************************************************************************
   * \brief diffusion between the upper and the lower soil layer
   *******************************************************************************************/
  void layerDiffusion();

  /***************************************************************************************
   * \brief transpiration and evaporation in the upper layer
   * \author Tong Guo
   *******************************************************************************************/
  void transpirationupper();

  /***************************************************************************************
   * \brief evaporation in the upper layer
   * \author Tong Guo
   *******************************************************************************************/
  void evaporationupper();

  /***************************************************************************************
   * \brief transpiration and evaporation in the lower layer
   * \author Tong Guo
   *******************************************************************************************/
  void transpirationlower();

  /***************************************************************************************
   * \brief transpiration of PFTs in the upper layer
   * \author Tong Guo
   *******************************************************************************************/
  void transpirationPFTU();

  /***************************************************************************************
   * \brief transpiration of PFTs in the lower layer
   * \author Tong Guo
   *******************************************************************************************/
  void transpirationPFTL();
};

#endif  // WATERCELL_H

//The following stuff is for the doxygen documentation webpage
/// \todo Needs to be updated by someone!
 //-----------------------------------------------------------------------------
/*! \mainpage EcoHyD by Britta Tietjen

 \section pub_sec Publications

 \subsection hydro_sub 1. Hydrological submodel of EcoHyD
 Tietjen B, Zehe E, Jeltsch F (2009)
 Simulating plant water availability in dry lands under climate change: A generic model of two soil layers.
 Water Resources Research 45:W01418.

 \subsection veg_sub 2. VegetationCell submodel of EcoHyD
 Tietjen B, Jeltsch F, Zehe E, Classen N, Groengroeft A, Schiffers K, Oldeland J (2010)
 Effects of climate change on the coupled dynamics of water and vegetation in drylands.
 Ecohydrol. 3:226-237.

 \subsection mod_changes 3. EcoHyD model changes
 Lohmann, D., Tietjen, B., Blaum, N., Joubert, D., Jeltsch, F. (2012)
 Shifting thresholds and changing degradation patterns: Climate change effects on the simulated long-term response of a semi-arid savanna to grazing.
 Journal of Applied Ecology 49, 814-823

 Lohmann D, Tietjen B, Blaum N, Joubert D, Jeltsch F. (2014)
 Prescribed fires as a tool for sustainable management of semi-arid savanna rangelands.
 Journal of Arid Environments 107, 49-56

 <BR>
 <BR>

 \section process_sec Process overview
 \image html ProcessOverview.jpg \todo include JPG in folder structure

 see also \ref page_1 "some equations"  you need to adapt, if you want to include new PFTs and \ref page_2 "the order of function calls"

 <BR>
 <BR>

 \subsection gen_sub General Information:
 - cell size: 5x5 m
 - each cell is characterized by a topographic height
 - open boundary conditions -> water loss due to runoff
 - unique soil texture for the whole grid
 - each cell covered with grasses, shrubs, both or without vegetation
 - no individuals, calculations depend on cover or biomass

 \subsection hydro2_sub  Hydrological Processes:
 - hourly \ref WaterLandscape::precipitation(double actualPrec) "precipitation" accumulates as surface water
 - surface water infiltrates into both layers
 - not infiltrated water is lost due to runoff to lower cells and/or evaporation
 - runoff is calculated hourly
 - evapotranspiration is calculated every 24 hours
 - \ref WaterCell::layerDiffusion() "diffusive flux" between both layers and deep drainage from lower layer(hourly)

 \subsection veg2_sub  Vegetation Processes:
 - two life forms: shrubs and perennial grasses
 - vegetation cover determined by growth, mortality, and dispersal/establishment
 - vegetation growth calculated biweekly in the growing season
 - three types of mortality: \ref Vegetation::winterMortality(int day)   ????
    "winter", water-independent mortality calculated biweekly,
    mortality due to water shortage calculated once a year
 - dispersal and establishment calculated once a year in a joined function
 - grazing can be simulated by an increased ??? mortality (see Table 4.2 of Britta's phd)

 \subsection feed_sub  Feebacks between water and vegetation:
 - \ref WaterCell::twoLayerInfiltration(double meanWaterL1) "infiltration" ->
    promoted by vegetation; further dependent on soil texture
 - \ref WaterLandscape::calculateRunoff() "runoff" -> reduced by vegetation;
    further dependent on inclination
 - \ref WaterCell::actLayerEvapotranspiration(double rainTot,int ep_type,double potEP)
   "evapotranspiration" -> reduced by vegetation due to shading (mainly upper layer)
    and increased due to transpiration (mainly lower layer);
    further dependent on \ref Weather::meanWeather(Parameter* p) "radiation",
    \ref Weather::getWeather(Parameter* p) "temperature",
    \ref WaterLandscape::setAspectAndInclination() "slope, aspect factor"
 - \ref Vegetation::vegGrowth(int day) "growth" ->  follows logistic behaviour with growth rate r
    and is reduced by limited water availability; further depended on own cover,
    cover of other plant type and overlap (20 %)
 - \ref Vegetation::vegMortality(int day, int waterStart, int waterEnd) "mortality" ->
    calculated as a result of unfavourable soil moisture conditions during the whole season
 - \ref VegetationLandscape::dispersal() "dispersal/establishment" -> dependent on water conditions in target cell;
    grasses dispers homogenous, shrubs have an exponential decrease with distance

 \subsection par_sub  Parameter Input:
 -  Parameters/<EM>Site</EM>\ref page_3 "soilparameters.txt"
 - Parameters/<EM>Site</EM>\ref page_4 "vegetationparameters.txt"
 - Parameters/grasscover_<EM>site_size</EM>.txt -> initial grass cover of each cell
 - Parameters/shrubcover_<EM>site_size</EM>.txt -> initial shrub cover of each cell
 - Parameters/\ref page_5 "modelparameters.txt"
 - Parameters/elevation_<EM>site_size</EM>.txt -> elevation of each cell
 - Weather/<EM> \ref page_6 "Site_simYears"</EM>.txt   (not real simYears, but simulation years from "modelparameters.txt")


 */

 //-----------------------------------------------------------------------------
 /*! \page page_2 Soil parameters

<h1>Parameters/<EM>Site</EM>soilparameters.txt</h1>

<table border="0">
  <tr>
    <td><h2>sf</h2></td>
    <td>suction [mm]</td>
    <td>Rawls et al. (1992)</td>
    <td>see also Table 4.5 of Britta's phd for values of different soils</td>
  </tr>


  <tr>
    <td><h2>Ks</h2></td>
    <td>(saturated) hydraulic conductivity [mm/h]</td>
    <td>Rawls et al. (1992)</td>
    <td>see also Table 4.5 of Britta's phd for values of different soils</td>
  </tr>

  <tr>
    <td><h2>WP</h2></td>
    <td>wilting point [vol%] / water content at beginning of stomatal closure</td>
    <td>Rawls et al. (1992)</td>
    <td>value for the calculation of Evapotranspiration; see also Table 4.5 of Britta's phd for values of different soils</td>
  </tr>

  <tr>
    <td><h2>FC</h2></td>
    <td>field capacity [vol%]</td>
    <td>Rawls et al. (1992)</td>
    <td> see also Table 4.5 of Britta's phd for values of different soils</td>
  </tr>

  <tr>
    <td><h2>sat</h2></td>
    <td>water content at saturation [vol%]</td>
    <td>called FC in her thesis, so take same value as for FC if you don't know the right one</td>
    <td></td>
  </tr>

  <tr>
    <td><h2>depL1</h2></td>
    <td>depth of the upper soil layer [mm]</td>
    <td>default: 200 mm </td>
    <td></td>
  </tr>

  <tr>
    <td><h2>depL2</h2></td>
    <td>depth of the lower soil layer [mm]</td>
    <td>default: 800 mm</td>
    <td> </td>
  </tr>

  <tr>
    <td><h2>iniSur</h2></td>
    <td>initial surface Water [mm]</td>
    <td>default: 0 mm</td>
    <td> </td>
  </tr>

  <tr>
    <td><h2>iniL1</h2></td>
    <td>initial moisture in upper soil layer [vol%]</td>
    <td>default: value of residual water content </td>
    <td> </td>
  </tr>

  <tr>
    <td><h2>iniL2</h2></td>
    <td>initial moisture in lowersoil layer [vol%]</td>
    <td>default: value of residual water content </td>
    <td> </td>
  </tr>

  <tr>
    <td><h2>EP</h2></td>
    <td>evaporation factor to reduce / increase EP from the upper layer [-]</td>
    <td>"Constant aET accounts for evaporation reduction e.g. caused by soil crusts. Normally, aET should not influence evaporation and should be set to 1."</td>
    <td>see p. 47 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>rw</h2></td>
    <td>residual water content, the soil cannot dry out below this content [vol%]</td>
    <td>Rawls et al. (1992)</td>
    <td> see also Table 4.5 of Britta's phd for values of different soils</td>
  </tr>

  <tr>
    <td><h2>maxFL2</h2></td>
    <td>rate of maximal infiltration [-]</td>
    <td>"As the parameters FL2,frac and FL2,bare are difficult to observe in the field, they have to be estimated
with great care." Had in the end been set to 0</td>
    <td> see p. 43 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>maxAmountFL2</h2></td>
    <td>maximal total infiltration [mm/h]</td>
    <td>estimated, had in the end been set to 0</td>
    <td> see p. 49 of Britta's phd</td>
  </tr>


  <tr>
    <td><h2>diffConst</h2></td>
    <td>diffusion speed between layers [-] / water balancing constant between layers</td>
    <td>estimated according to Bai et al. (2007)</td>
    <td> see also Table 4.5 of Britta's phd </td>
  </tr>
</table>

 */


 //-----------------------------------------------------------------------------
 /*! \page page_3 Vegetation parameters

<h1>Parameters/<EM>Site</EM>vegetationparameters.txt</h1>

<table border="0">
  <tr>
    <td><h2>gUptakeRate</h2></td>
    <td>relative uptake rate per gram grass biomass [mm_Water/(g_Biomass*Year)]</td>
    <td>van Langevelde et al., 2003</td>
    <td>0.9; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>sUptakeRate</h2></td>
    <td>relative uptake rate per gram shrub biomass [mm_Water/(g_Biomass*Year)]</td>
    <td>van Langevelde et al., 2003</td>
    <td>0.5; see also Table 2 of Tietjen et. al 2010</td>
  </tr>

  <tr>
    <td><h2>gRoot1</h2></td>
    <td>relative fraction of grass roots in the upper layer [-]</td>
    <td>root_g,L1 = 1-(b_g)^d; b_g = 0.952, d = 20 cm - Jackson et al. (1996)</td>
    <td>0.63; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>gRoot2</h2></td>
    <td>relative fraction of grass roots in the lower layer [-]</td>
    <td>gRoot2 = 1-gRoot1</td>
    <td>0.37; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>sRoot1</h2></td>
    <td>relative fraction of grass roots in the upper layer [-]</td>
    <td>root_s,L1 = 1-(b_g)^d; b_g = 0.978, d = 20 cm - Jackson et al. (1996)</td>
    <td>0.36; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>sRoot2</h2></td>
    <td>relative fraction of grass roots in the lower layer [-]</td>
    <td>sRoot2 = 1-sRoot1</td>
    <td>0.64; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>gGrowthR</h2></td>
    <td>potential growth rate of grass [mm^-1 yr^-1]</td>
    <td>calibrated to gain long-term of grass values with a mean grass cover of 40-80% and a mean cover of 15-25%</td>
    <td>0.8; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>sGrowthR</h2></td>
    <td>potential growth rate of shrubs [mm^-1 yr^-1]</td>
    <td>calibrated to gain long-term of grass values with a mean grass cover of 40-80% and a mean cover of 15-25%</td>
    <td>0.2; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>gMortR</h2></td>
    <td>mortality rate of grass due to water stress [mm^-1 yr^-1]</td>
    <td>van Langevelde et al. (2003)</td>
    <td>0.9; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>sMortR</h2></td>
    <td>mortality rate of shrub due to water stress [mm^-1 yr^-1]</td>
    <td>van Langevelde et al. (2003)</td>
    <td>0.4; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>gMortWinter</h2></td>
    <td>winter mortality of grass, independent of water (fraction of removed cover) [-]</td>
    <td></td>
    <td>default: 0.5; can be increased to simulate grazing (see Table 4.2 of Britta's phd)</td>
  </tr>

    <tr>
    <td><h2>gMortGrazing</h2></td>
    <td>grazing mortality of grass, independent of water (fraction of removed cover) [-]</td>
    <td></td>
    <td>default: 0.5; can be increased to simulate grazing (see Table 4.2 of Britta's phd)</td>
  </tr>

  <tr>
    <td><h2>sMortWinter</h2></td>
    <td>winter mortality of shrubs, independent of water (fraction of removed cover) [-]</td>
    <td></td>
    <td>default: 0</td>
  </tr>

    <tr>
    <td><h2>sMortGrazing</h2></td>
    <td>grazing mortality of shrubs, independent of water (fraction of removed cover) [-]</td>
    <td></td>
    <td>default: 0</td>
  </tr>

  <tr>
    <td><h2>gMinWaterL1</h2></td>
    <td>Wilting point, herbaceous vegetation [vol%]</td>
    <td>estimated based on Sala et al. (1989) and Neilson (1995)</td>
    <td>value for growth, mortality, establishment; see also Table 4.5 of Britta's phd for values of different soils
    (as it's only estimated, you should lower it, if nothing's growing, and increase it, if water content is above most of the time;
    see also http://www.terragis.bees.unsw.edu.au/terraGIS_soil/sp_water-soil_moisture_classification.html for ranges)</td>
  </tr>

  <tr>
    <td><h2>gMinWaterL2</h2></td>
    <td>Wilting point, herbaceous vegetation [vol%]</td>
    <td>estimated based on Sala et al. (1989) and Neilson (1995)</td>
    <td>value for growth, mortality; see also Table 4.5 of Britta's phd for values of different soils
    (as it's only estimated, you should lower it, if nothing's growing, and increase it, if water content is above most of the time;
    see also http://www.terragis.bees.unsw.edu.au/terraGIS_soil/sp_water-soil_moisture_classification.html for ranges)</td>
  </tr>

  <tr>
    <td><h2>sMinWaterL1</h2></td>
    <td>Wilting point, woody vegetation [vol%]</td>
    <td>estimated based on Sala et al. (1989) and Neilson (1995)</td>
    <td>value for growth, mortality, establishment; see also Table 4.5 of Britta's phd for values of different soils
    (as it's only estimated, you should lower it, if nothing's growing, and increase it, if water content is above most of the time;
    see also http://www.terragis.bees.unsw.edu.au/terraGIS_soil/sp_water-soil_moisture_classification.html for ranges)</td>
  </tr>

  <tr>
    <td><h2>sMinWaterL2</h2></td>
    <td>Wilting point, woody vegetation [vol%]</td>
    <td>estimated based on Sala et al. (1989) and Neilson (1995)</td>
    <td>value for growth, mortality; see also Table 4.5 of Britta's phd for values of different soils
    (as it's only estimated, you should lower it, if nothing's growing, and increase it, if water content is above most of the time;
    see also http://www.terragis.bees.unsw.edu.au/terraGIS_soil/sp_water-soil_moisture_classification.html for ranges)</td>
  </tr>

  <tr>
    <td><h2>resWater</h2></td>
    <td>residual water content, the soil cannot dry out below this content [vol%]</td>
    <td>Rawls et al. (1992)</td>
    <td> see also Table 4.5 of Britta's phd for values of different soils (not sure why it's needed here AND in soilparameters.txt!)</td>
  </tr>

  <tr>
    <td><h2>moistL1</h2></td>
    <td>initial moisture in upper soil layer [vol%]</td>
    <td>default: value of residual water content </td>
    <td>(not sure why it's needed here AND in soilparameters.txt!)</td>
  </tr>

  <tr>
    <td><h2>moistL2</h2></td>
    <td>initial moisture in upper soil layer [vol%]</td>
    <td>default: value of residual water content </td>
    <td>(not sure why it's needed here AND in soilparameters.txt!)</td>
  </tr>

  <tr>
    <td><h2>dispFracG</h2></td>
    <td>rate of successful establishment of grasses [yr^-1]</td>
    <td>estimated value, low sensitivity</td>
    <td>0.05; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>dispFracS</h2></td>
    <td>rate of successful establishment of shrubs [yr^-1]</td>
    <td>estimated value, low sensitivity</td>
    <td>0.005; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>dist0S</h2></td>
    <td>constant for exponential dispersal decline with distance [-]</td>
    <td>values are chosen in a way that the cumulative dispersal is 75% at dist = 2.5m (border and 99% at distmax</td>
    <td>0.5; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>distConstS</h2></td>
    <td>constant for exponential dispersal decline with distance [-]</td>
    <td>values are chosen in a way that the cumulative dispersal is 75% at dist = 2.5m (border and 99% at distmax</td>
    <td>0.53; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>maxDistS</h2></td>
    <td>maximal distance of shrub seed dispersal [m]</td>
    <td>short distance dispersal only, e.g. Jeltsch et al. (1996)</td>
    <td>20; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>growStart</h2></td>
    <td>starting day of the growing season [-]</td>
    <td>long-term soil moisture observations in the model</td>
    <td>150; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>growEnd</h2></td>
    <td>last day of the growing season [-]</td>
    <td>long-term soil moisture observations in the model</td>
    <td>330; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>cmaxg</h2></td>
    <td>maximum value for grass cover [-]</td>
    <td></td>
    <td>default: 1</td>
  </tr>

  <tr>
    <td><h2>cmaxs</h2></td>
    <td>maximal extent to which grass and shrubs can overlap [%]</td>
    <td>observations (J. Oldeland, unpublished)</td>
    <td>20; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>conversionC_BM_G</h2></td>
    <td>conversion: grass biomass to relative cover [g]</td>
    <td>J. Oldeland, unpublished</td>
    <td>8000000; see also Table 3.2 of Britta's phd</td>
  </tr>

  <tr>
    <td><h2>conversionC_BM_S</h2></td>
    <td>conversion: shrub biomass to relative cover [g]</td>
    <td>J. Oldeland, unpublished</td>
    <td>50000000; see also Table 2 of Tietjen et. al 2010</td>
  </tr>

</table>
  */

//------------------------------------------------------------------------------
  /*! \page page_4 Model parameters

<h1>Parameters/modelparameters.txt</h1>

<table border="0">
  <tr>
    <td><h2>site</h2></td>
    <td>name of the site (or soil type or whatever you want, important for names of the other input files)</td>
  </tr>

  <tr>
    <td><h2>latitude</h2></td>
    <td>latitude of the region the model is applied to</td>
  </tr>

  <tr>
    <td><h2>simulation years</h2></td>
    <td>just for weather input file, log file and output files, real simulation years are changed in parameters.h</td>
  </tr>

  <tr>
    <td><h2>latitude</h2></td>
    <td>calculation type of potential evapotranspiration (model just validated for "hargreaves")</td>
  </tr>

  <tr>
    <td><h2>weather data</h2></td>
    <td>weather data, 0: original data (with lots of entries) 1: default; simple data (yy mm dd hh prec temp, tab-separated)</td>
  </tr>

</table>
  */


//------------------------------------------------------------------------------
  /*! \page page_5 Weather

<h1>Weather/<EM>Site_simYears</EM>.txt </h1>

Format is Year | Month | Day | Hour | Precipitation [mm] | Temperature [°C]

in Pictures/TemperatureCurveExample.xls is an example of how to create hourly
temperature data with a cosinus function in case you don't have real values  (see p. 55 of Britta's phd)

Filename: Site and simYears from modelparameters.txt


  */


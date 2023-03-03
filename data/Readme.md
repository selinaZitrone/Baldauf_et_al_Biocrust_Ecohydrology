The `model_results` folder contains all model output for the scenarios 
analysed in the paper. To reproduce the results, please refer to `/model`.
All `transect.*` files are needed to extract the transect from the model results.

The result folders are all named `Results_*` where `*` represents the sceanario name.
Scenarios are:

- cyano: landscape only covered by cyanobacteria biocrust
- lichen: landscape only covered by lichen biocrust
- physical: landscape only covered by incipient cyanobacteria biocrust
- tabernas: current surfacecover and weather from 2010 for El Cautivo
- tabernas_halfrain: current surfacecover and 50% reduced rainfall for El Cautivo
- zeroCrust: no biocrust cover (bare soil) and weather from 2010 for El Cautivo
- zeroCrust_halfrain: no biocrust cover (bare soil) and 50% reduced rainfall for El Cautivo

In every result folder, the following files are contained:

- `modelscenarios.txt`: Scenario file that shows which scenario was run
- `/Spatial`: Folder with spatial output for the months
- `*_1years_CS_*`: Hourly model output for the crust processes
- `*_1years_hourly_*`: Hourly model output of hydrology

The folder `spatial` contains shape files and tifs to extract the hillslope and 
the transect from the results for the plots presented in the paper.
# Reproduce EcohyD model results

This folder contains all scripts and the pipeline necessary to run the model and
reproduce the results from the script. Below you find a description of all files
and in the end a guide on how to run the model on your machine.

## Contents of this folder

### Model folders and files

#### Folder `/src` with model

Contains the Ecohyd model source files. The biocrust functions are contained in 
`BiocrustCell.cpp` for cell functions and `BiocrustLandscape.cpp` for landscape
functions.

#### Model input folders

##### `Parameters_all/`

Folder that contains all model input parameter files. When running the different
scenarios, the currently needed parameter files are picked from here and moved
to the `Parameters/` directory to run the current simulation.

- `/crust`: Crust parameters for incipient cyanobacteria (`crust_1.yml`), cyanobacteria (`crust_2.yml`) and
lichen (`crust_3.yml`) biocrusts
- `crustcover_74_*`: Map of biocrust cover for the different scenarios
  - `_cyano`/`_incipient`/`_lichen`/`_nocrust`: Scenarios where the entire landscape is covered by only one biocrust type or no biocrust at all
  - `_tabernas`: current biocrust cover at the study site
- `vegperennialcover_74_tabernas`: location of annual herbaceous vegetation
- `vegshrubcover_74_tabernas`: location of shrub vegetation

##### `Parameters/`

Temporary folder that is used by the model to access the current parameter files.
Parameter files come from `Parameters_all/`

##### `Weather`

- Weather files used to run the model for standard and 50% reduced rainfall (file
ending `_half`)

#### Folders and files to run the model

##### `R/`

- `write_scenarios.R`: R script that writes the scenario definition file that is 
used by the model. This script is automatically called by `run_scenario.sh`

##### `compile.sh`

Compiles the model using g++ (called automatically in `run_all_scenarios.sh`).

##### `make_output_dirs.sh`

Shell file that creates output folders (called automatically in `run_all_scenarios.sh`).

##### `run_scenario.sh`

Shell script to run a single model scenario (called automatical3ly in `run_all_scenarios.sh`).

#### `run_all_scenarios.sh`

Shell script that runs all model scenarios analysed in the paper.

## How to build and run the model

A prerequisite to build the model is the installation and setup of the C++ compiler
g++. It can be downloaded from [here](https://gcc.gnu.org/). Installation guides can 
be found online.

The script `run_all_scenarios.sh` contains the workflow to compile the model and then
subsequently run all scenarios that are analyed in the paper:
- standard scenario with current conditions in El Cautivo
- standard scenario with current conditions in El Cautivo with 50% reduced rainfall
- scenarios only covered with one of the 3 biocrust types
- scenario with only bare soil (standard climate and 50% reduced climate)

You can start the script in the terminal with:

```sh
./run_all_scenarios.sh
```

The results are also already included in the `data/model_results/` folder of this repository.
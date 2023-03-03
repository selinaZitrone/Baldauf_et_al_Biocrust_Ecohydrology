#!/bin/bash -e

# read the command line parameters
# s = soil
# e = elevation
# p = parameter
# v = vegetation

while getopts ":d:n:s:e:p:v:x:y:z:" option; do
    case ${option} in

    d) dim=${OPTARG} ;;
    n) name=${OPTARG} ;;
    e) elevation=${OPTARG} ;;
    s) soil=${OPTARG} ;;
    p) parameter=${OPTARG} ;;
    v) vegetation=${OPTARG} ;;
    x) mapgrass=${OPTARG} ;;
    y) mapshrub=${OPTARG} ;;
    z) mapcrust=${OPTARG} ;;
    esac
done

# set default values if not given in command line
if [ "$name" == "" ]; then
    echo "No name given for scenario. Use [-n] option"
    exit 1
fi

if [ "$dim" == "" ]; then
    dim="30"
fi

if [ "$elevation" == "" ]; then
    elevation="steep"
fi

if [ "$soil" == "" ]; then
    soil="tabernas"
fi

if [ "$parameter" == "" ]; then
    parameter="tabernas"
fi

if [ "$vegetation" == "" ]; then
    vegetation="control"
fi

echo "Scenario   : $name"
echo "Dimension  : $dim"
echo "Elevation  : $elevation"
echo "Soil       : $soil"
echo "Vegetation : $vegetation"
echo "Parameter  : $parameter"
echo "Grass_map  : $mapgrass"
echo "Shrub_map  : $mapshrub"
echo "Crust_map  : $mapcrust"

# make the output directories

# echo "Create output directories ------------------"

./make_output_dirs.sh

# write correct model scenarios file using R

echo "Run Rscript to create modelscenarios ------------------"

Rscript ./R/write_scenarios.R "$elevation" "$soil" "$vegetation" "$parameter" "$name"

# copy the model scenarios that will be run into the ouput dir

echo "Copy model scenarios to output ------------------"
cp ./Parameters/modelscenarios.txt ./Results/

# Move the vegetation maps if necessary
echo "Move vegetation maps ------------------"

if [ "$mapgrass" != "" ]; then
    echo "Move ./Parameters_all/vegperennialcover_${dim}_${mapgrass}.txt to ./Parameters/vegperennialcover_$dim.txt"
    cp "./Parameters_all/vegperennialcover_${dim}_${mapgrass}.txt" "./Parameters/vegperennialcover_$dim.txt"
fi

if [ "$mapshrub" != "" ]; then
    echo "Move ./Parameters_all/vegshrubcover_${dim}_${mapshrub}.txt to ./Parameters/vegshrubcover_$dim.txt"
    cp "./Parameters_all/vegshrubcover_${dim}_${mapshrub}.txt" "./Parameters/vegshrubcover_$dim.txt"
fi

if [ "$mapcrust" != "" ]; then
    echo "Move ./Parameters_all/crustcover_${dim}_${mapcrust}.txt to ./Parameters/crustcover_$dim.txt"
    cp "./Parameters_all/crustcover_${dim}_${mapcrust}.txt" "./Parameters/crustcover_$dim.txt"
fi

# run ecohyd
echo "Run model ------------------"
./Ecohyd.exe

# rename the results folder
echo "Rename results ------------------"
mv Results/ "Results_$name"

# move back the vegetation maps

echo "Move back vegetation maps------------------"

if [ "$mapgrass" != "" ]; then
    echo "Remove ./Parameters/vegperennialcover_$dim.txt"
    rm "./Parameters/vegperennialcover_$dim.txt"
fi

if [ "$mapshrub" != "" ]; then
    echo "Remove ./Parameters/vegshrubcover_$dim.txt"
    rm "./Parameters/vegshrubcover_$dim.txt"
fi

if [ "$mapcrust" != "" ]; then
    echo "Remove ./Parameters/crustcover_$dim.txt"
    rm "./Parameters/crustcover_$dim.txt"
fi

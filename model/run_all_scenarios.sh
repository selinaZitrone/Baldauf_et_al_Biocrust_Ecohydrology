#!/bin/bash -e

# --------------------------------------------------------------------------
# Compile the model
# --------------------------------------------------------------------------

./compile.sh

# --------------------------------------------------------------------------
# run all scenarios with standard rainfall
# --------------------------------------------------------------------------

./run_scenario.sh -n cyano -s tabernas -x tabernas -y tabernas -z cyano -e tabernas -d 74
./run_scenario.sh -n lichen -s tabernas -x tabernas -y tabernas -z lichen -e tabernas -d 74
./run_scenario.sh -n incipient -s tabernas -x tabernas -y tabernas -z incipient -e tabernas -d 74
./run_scenario.sh -n tabernas -s tabernas -x tabernas -y tabernas -z tabernas -e tabernas -d 74
./run_scenario.sh -n zeroCrust -s tabernas -x tabernas -y tabernas -z nocrust -e tabernas -d 74

mv ./Weather/tabernas_1years_tabernas_climrep-1.txt ./Weather/tabernas_1years_tabernas_climrep-1_def.txt

mv ./Weather/tabernas_1years_tabernas_climrep-1_half.txt ./Weather/tabernas_1years_tabernas_climrep-1.txt

# --------------------------------------------------------------------------
# run all scenarios with half rainfall
# --------------------------------------------------------------------------

./run_scenario.sh -n tabernas_halfrain -s tabernas -x tabernas -y tabernas -z tabernas -e tabernas -d 74
./run_scenario.sh -n zeroCrust_halfrain -s tabernas -x tabernas -y tabernas -z nocrust -e tabernas -d 74
./run_scenario.sh -n cyano_halfrain -s tabernas -x tabernas -y tabernas -z cyano -e tabernas -d 74
./run_scenario.sh -n lichen_halfrain -s tabernas -x tabernas -y tabernas -z lichen -e tabernas -d 74
./run_scenario.sh -n incipient_halfrain -s tabernas -x tabernas -y tabernas -z incipient -e tabernas -d 74

# move weather files back
mv ./Weather/tabernas_1years_tabernas_climrep-1.txt ./Weather/tabernas_1years_tabernas_climrep-1_half.txt

mv ./Weather/tabernas_1years_tabernas_climrep-1_def.txt ./Weather/tabernas_1years_tabernas_climrep-1.txt

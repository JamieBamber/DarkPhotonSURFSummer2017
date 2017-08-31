#!/bin/bash

########### Full Input -> output simulation process (from the .pythia file provided by Maurizio) ########

# © Jamie Bamber 2017

## un-comment the lines required

export NAME=HiggsZH
export PYTHIA_CARD=HiggsZH.pythia 

### .pythia to .lhe file

cd [path]/BSMatLHC/BSMGen/
source setup.sh
./GenPythiaToLHE data/pythiaCards/HiggsInclusive/$PYTHIA_CARD ../../lhe_hepmc_files/$NAME.lhe

### modify the .lhe file

cd /Users/Jamie/Documents/University/Caltech/Project
python3 ./Scripts/LHE_photon_neutralino_v2.py "./lhe_hepmc_files/${NAME}.lhe" "./lhe_hepmc_files/${NAME}_modified.lhe"

### make the hepmc file

cd /Users/Jamie/Documents/University/Caltech/Project/BSMatLHC/BSMGen
./LHEGenToHepMC data/pythiaCards/HiggsInclusive/$PYTHIA_CARD ../../lhe_hepmc_files/${NAME}_modified.lhe ../../lhe_hepmc_files/HiggsZH_modified.hepmc

### make the new directory (starting from Scripts directory)

cd [path]
mkdir $NAME

### make root file

rm [path]/ROOT_outputs/ZH/${NAME}_modified.root
#
cd [path]
./BSMatLHC/delphes/DelphesHepMC BSMatLHC/delphes/cards/delphes_card_CMS_Jamie.tcl ROOT_outputs/ZH/${NAME}_modified.root lhe_hepmc_files/${NAME}_modified.hepmc
cd ROOT_outputs/ZH

### open root file

root -l << EOF
TFile::Open("HiggsZH_modified.root");
Delphes->MakeClass();
.q
EOF

cd [scripts directory]




#!/bin/bash

#cd /Users/Jamie/Documents/University/Caltech/Project/ROOT_outputs/VBF

### /// This script produces a csv file containing the Trigger Table data

### run root scripts 

### clear file 

export OutFile=Trigger_Table_data_VBF.csv

>$OutFile

### create headings

printf "\nTrigger Table1: Background\n" >> $OutFile
printf "Cut no., Selected events (raw), Selected events (weighted, lum. = 100fb^-1)" >> $OutFile

###

root -l DarkPhotonAnalyzer_GJets_HT-40ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root << EOF
.q
EOF

for NUM in 0 1 2 3 4 5 7
do
root -l  << EOF >> $OutFile
.L Trigger_Table_VBF_bkg.C
BkgTree t;
t.Loop($NUM);
.q
EOF
done

echo "Background done"

### create headings

printf "\n" >> $OutFile
printf "Trigger Table 2: Signal\n" >> $OutFile
printf "Cut no., Selected events (total events = 50,000), Selected events (weighted, lum = 100fb^-1)" >> $OutFile

###

root -l HiggsVBF_modified.root << EOF
.q
EOF

for NUM in 0 1 2 3 4 5 7
do
root -l  << EOF >> $OutFile
.L Trigger_Table_VBF_sig.C
Delphes t;
t.Loop($NUM);
.q
EOF
done

#cd /Users/Jamie/Documents/University/Caltech/Project/ROOT_outputs/VBF



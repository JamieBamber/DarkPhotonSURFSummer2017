#!/bin/bash

### /// This script produces a csv file containing the Trigger Table data for the signal and background, for the all inclusive Higgs production
# 
# © Jamie Bamber 2017
# 

## clear file 

export OutFile=Trigger_Table_data_v2.csv

> $OutFile

### create headings

printf "Trigger Table 1: Background\n" >> $OutFile
printf "Cut no., Selected events (raw), Selected events (weighted lum. = 100fb^-1)\n" >> $OutFile

###

root -l <Background_file> << EOF
.q
EOF

for NUM in 0 1 2 3 4 5
do
root -l << EOF >> $OutFile
.L Trigger_Table_bkg.C
BkgTree t;
t.Loop($NUM);
.q
EOF
done

echo "Background done"

### create headings

printf "\n" >> $OutFile
printf "Trigger Table 2: Signal\n" >> $OutFile
printf "Cut no., Selected events (total events = 50000), Selected events (weighted lum = 100fb^-1)\n" >> $OutFile

###

root -l <Signal_file> << EOF
.q
EOF

for NUM in 0 1 2 3 4 5
do
root -l << EOF >> $OutFile
.L Trigger_Table_sig.C
Delphes t;
t.Loop($NUM);
.q
EOF
done





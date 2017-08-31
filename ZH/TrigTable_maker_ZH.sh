#!/bin/bash

### /// This script produces a csv file containing the Trigger Table data

### run root scripts 

### clear file 

export OutFile=Trigger_Table_data_ZH.csv

> $OutFile

### create headings

#printf "Trigger Table1: Background\n" >> $OutFile
#printf "Cut no., Selected events (raw), Selected events (weighted, lum. = 100fb^-1)\n" >> $OutFile

###

root -l ZGamma_Bkg_1pb_weighted.root << EOF
.q
EOF

for NUM in 7 8 9 10 11
do
root -l << EOF >> $OutFile
.L Trigger_Table_bkg_ZH.C
BkgTree t;
t.Loop($NUM);
.q
EOF
COUNT=$(($NUM+1))
echo "done cut no." $COUNT
done

echo "Background done"

### create headings

#printf "\n" >> $OutFile
#printf "Trigger Table 2: Signal\n" >> $OutFile
#printf "Cut no., Selected OutTrees (total OutTrees = 50,000), Selected OutTrees (weighted, lum = 100fb^-1), Efficiency / %%, Cuts\n" >> $OutFile

###

#root -l HiggsZH_modified.root << EOF
#.q
#EOF

#for NUM in 0 1 2 3 4 5 6 7 8 9 10
#do
#root -l << EOF >> $OutFile
#.L Trigger_Table_sig_ZH_v2.C
#Delphes t;
#t.Loop($NUM);
#.q
#EOF
#let COUNT=$NUM+1 
#echo "done " $COUNT "of 12"
#done

#cd /Users/Jamie/Documents/University/Caltech/Project/ROOT_outputs/ZH

#!/bin/bash

###################

### /// This script produces my charts from the relevent .C ROOT scripts

### run root scripts 

for FILE in PT_ZH MT_ZH MET_ZH 2D_PT_vs_MT_ZH 2D_PT_vs_MET_ZH 2D_MT_vs_MET_ZH  
do
root -l << EOF
.L $FILE.C
Delphes t;
t.Loop();
c1->SaveAs("$FILE.pdf");
.q
EOF
done 





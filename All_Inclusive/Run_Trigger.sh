#!/bin/bash

###

### 

export b_file="Trigger_bkg.C"
export s_file="Trigger_sig.C"

root -l << EOF
TFile::Open("GJets_Bkg_1pb_weighted.root")
.L $b_file
BkgTree t1;
t1.Loop();
// ############# Signal file next
TFile::Open("HiggsInclusive_modified_CMS_Jamie_50k.root")
.L $s_file
Delphes t2;
t2.Loop();
// 
EOF
 





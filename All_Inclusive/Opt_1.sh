#!/bin/bash

###########  ########

## This script produces a 2D histogram for the all inclusive signal / gamma+jet background ratio vs the pT and MET cuts.
## additional cuts are applied detailed in the relevent signal and background scripts
##
## © Jamie Bamber 2017

### run root scripts 

export b_file="Opt1_bkg.C" 	# choose background script
export s_file="Opt1_sig.C"	# choose background script
export data_file="Opt_1.root"

root -l << EOF
TFile::Open("GJets_Bkg_1pb_weighted.root");
.L $b_file
BkgTree t1;
// ############# Signal file next
TFile::Open("HiggsInclusive_modified_CMS_Jamie_50k.root");
.L $s_file
Delphes t2;
//
// ### parameters
Float_t PT0;		// photon pT cut threshold
Float_t E0;			// MET cut threshold
Float_t MT_UP;		// upper MT (transverse mass) bound
Float_t MT_LW;		// lower MT (transverse mass) bound
Bool_t Iso_leptons;
Float_t Eta0;
Float_t HPT0;
Float_t sig_value;
Float_t bkg;
// #### parameter array Param_arr
Float_t Best_Param[10];	// [ PT0, E0, MT_UP, MT_LW, Iso_leptons, Eta0, HPT0, signal (weighted), background (weighted), ratio ]  
Best_Param[9] = 0;		// set initial best ratio to 0
//
TH2F* h_ratio = new TH2F("h_ratio", "SB Ratio vs pT & MET cuts", 40, 0, 200, 40, 0, 200);
Double_t SBratio;
Double_t max_ratio;
// ### Baseline parameter settings
Iso_leptons = true;
Eta0 = 1.442;
HPT0 = 0;
// #####

//have a requirement for minimum signal efficiency ??



// ##### Loop over variable parameters
for (int i1=0; i1<=16; i1++) {
	MT_LW = 70 + i1*5;
	for (int i2=i1+4; i2<=36; i2++) {
		MT_UP = 70 + i2*5;
		t1.Loop(50, 50, MT_UP, MT_LW, Iso_leptons, Eta0, HPT0);
		t2.Loop(50, 50, MT_UP, MT_LW, Iso_leptons, Eta0, HPT0);
		TFile::Open("$data_file");
		for (Int_t i=1; i<=41; i++) {
			for (Int_t j=1; j<=41; j++) {
				//SBratio = static_cast<Double_t>(h_sig->GetBinContent(i,j))/TMath::Sqrt(h_bkg->GetBinContent(i,j) + pow(0.1*h_bkg->GetBinContent(i,j),2));
				//h_ratio->SetBinContent(i, j, SBratio);
			}
		}
		Int_t bin = h_ratio->GetMaximumBin();
		Int_t binx, biny, binz;
		max_ratio = h_ratio->GetBinContent(binx,biny);
		if (max_ratio>Best_Param[9]) { // if the signal / background is the best value so far
 			Best_Param[9] = max_ratio;
 			Best_Param[0] = (biny-1)*(200.0/40);
 			Best_Param[1] = (binx-1)*(200.0/40);
			Best_Param[2] = MT_UP;
			Best_Param[3] = MT_LW;
			Best_Param[7] = h_sig->GetBinContent(binx,biny);
			Best_Param[8] =	h_bkg->GetBinContent(binx,biny);
		}
	}
}
// ### Ratio histogram
cout << "P0 = " << Best_Param[0] << endl;
cout << "E0 = " << Best_Param[1] << endl;
cout << "MT_LW = " << Best_Param[3] << endl;
cout << "MT_UP = " << Best_Param[2] << endl;
cout << "Background = " << Best_Param[8] << endl;
cout << "Signal = " << Best_Param[7] << endl;
cout << "ratio = " << Best_Param[9] << endl;
EOF
 





#!/bin/bash

###########  ########

## This script produces a 2D histogram for the all inclusive signal / gamma+jet background ratio vs the pT and MET cuts.
## additional cuts are applied detailed in the relevent signal and background scripts
##
## © Jamie Bamber 2017

### run root scripts 

export b_file="2D_E0_P0_bkg_extra.C" 	# choose background script
export s_file="2D_E0_P0_sig_extra.C"	# choose background script
export data_file="2D_E0_P0.root"
export o_file="2D_E0_P0_extra_cuts.pdf"	# output filename

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
// ############# Draw Histograms
// ### Ratio histogram
TCanvas *c1 = new TCanvas();
c1->SetCanvasSize(1000,700);
TFile::Open("$data_file");
TH2F* h_ratio = new TH2F("h_ratio", "SB Ratio vs pT & MET cuts, with added cuts", 40, 0, 200, 40, 0, 200);
for (Int_t i=1; i<=41; i++) {
	for (Int_t j=1; j<=41; j++) {
		Double_t ratio = static_cast<Double_t>(h_sig->GetBinContent(i,j))/TMath::Sqrt(h_bkg->GetBinContent(i,j) + pow(0.1*h_bkg->GetBinContent(i,j),2));
		h_ratio->SetBinContent(i, j, ratio);
		//cout << h_bkg->GetBinContent(i,j) << ", " << h_sig->GetBinContent(i,j) << ", " << ratio << endl;		// compute S/B ratio and output
	}
}
h_ratio->GetXaxis()->SetTitle("MET cut / GeV");
h_ratio->GetXaxis()->CenterTitle();
h_ratio->GetYaxis()->SetTitle("PT cut / GeV");
h_ratio->GetYaxis()->CenterTitle();
h_ratio->GetZaxis()->SetTitle("R = S/#sqrt{B + (0.1B)^{2}}");
h_ratio->GetZaxis()->CenterTitle();
//c1->SetLogz();
h_ratio->SetTitleOffset(1.25,"X");
h_ratio->SetTitleOffset(1.25,"Y");
h_ratio->SetTitleOffset(1.25,"Z");
c1->SetLeftMargin(0.14);
c1->SetRightMargin(0.14);
h_ratio->Draw("COLZ");
TPaveStats *st1 = (TPaveStats*)h_ratio->FindObject("stats");
st1->SetX1NDC(0.6);
st1->SetY1NDC(0.72);
st1->SetX2NDC(0.8);
st1->SetY2NDC(0.88);
// ############ Save File
c1->SaveAs("$o_file");
EOF
 





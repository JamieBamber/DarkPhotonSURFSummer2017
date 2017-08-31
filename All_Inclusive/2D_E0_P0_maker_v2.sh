#!/bin/bash

###########  ########

## This script produces a plot with multiple histograms for the all inclusive signal / gamma+jet background ratio vs the pT cut for different MET cuts.
## cuts of |eta|<1.44 and 100GeV < MT < 130GeV are also used (although this can be modified)
##
## MET cuts of 50, 60, 70, 80, ... 150, 160 GeV are shown
##
## © Jamie Bamber 2017

### run root scripts 

export b_file="2D_E0_P0_bkg.C"
export s_file="2D_E0_P0_sig.C"
export data_file="2D_E0_P0.root"
export o_file="2D_E0_P0_v2.pdf"

root -l << EOF
/*TFile::Open("GJets_Bkg_1pb_weighted.root")
.L $b_file
BkgTree t1;
t1.Loop2();
// ############# Signal file next
TFile::Open("HiggsInclusive_modified_CMS_Jamie_50k.root")
.L $s_file
Delphes t2;
t2.Loop2();*/
// ############# Draw Histograms
// ### Ratio histogram
TFile::Open("$data_file");
TH2F* h_ratio = new TH2F("h_ratio", "Signal to Background Ratio vs pT & MET threshold", 40, 0, 200, 40, 0, 200);
for (Int_t i=1; i<=41; i++) {
	for (Int_t j=1; j<=41; j++) {
		Double_t ratio = static_cast<Double_t>(h_sig->GetBinContent(i,j))/TMath::Sqrt(h_bkg->GetBinContent(i,j) + pow(0.1*h_bkg->GetBinContent(i,j),2));
		h_ratio->SetBinContent(i, j, ratio);
		//cout << h_bkg->GetBinContent(i,j) << ", " << h_sig->GetBinContent(i,j) << ", " << ratio << endl;		// compute S/B ratio and output
	}
}
// ######## make multi-line plot of many histograms
TCanvas *c1 = new TCanvas();
c1->SetCanvasSize(1000,700);
TH1F* h = new TH1F("h", "Signal to Backgroud ratio vs pT cut for different MET cuts", 40, 0, 200);
auto legend = new TLegend(0.8,0.45,0.95,0.95);
legend->SetHeader("Key","C");
legend->SetTextSize(0.02);
Double_t MAX=1;
Double_t MIN=0;
h->SetAxisRange(MIN,MAX,"Y");
h->GetXaxis()->SetTitle("PT cut / GeV");
h->GetXaxis()->CenterTitle();
h->GetYaxis()->SetTitle("R = S/#sqrt{B + (0.1B)^{2}}");
h->GetYaxis()->CenterTitle();
h->SetTitleOffset(1.25,"X");
h->SetTitleOffset(1.25,"Y");
gPad->Update();
h->Draw("HIST");
TH1F *hist_arr[12];
for (int n=0; n<12; n++) {
	TString h_name = TString(Form("MET>%dGeV",(10*n+50)));
	hist_arr[n] = new TH1F(h_name,"", 40, 0, 200);
	for (int j=1; j<=41; j++) {
		hist_arr[n]->SetBinContent(j,h_ratio->GetBinContent((2*n+10),j));
	}
	hist_arr[n]->SetLineColor(1+(11-n)%9);
	if( n<=1 ) hist_arr[n]->SetLineStyle(8);
	hist_arr[n]->Draw("HIST SAME");
	gPad->Update();
	legend->AddEntry(hist_arr[n],h_name);
}
legend->Draw();
c1->SetLeftMargin(0.14);
c1->SetRightMargin(0.14);
TPaveStats *st1 = (TPaveStats*)h->FindObject("stats");
st1->SetX1NDC(1);
st1->SetY1NDC(1);
st1->SetX2NDC(1.2);
st1->SetY2NDC(1.2);
// ############ Save File
TFile *hfile = new TFile("2D_E0_P0.root", "UPDATE");
c1->Write("RECREATE");
hfile->Close();
c1->SaveAs("$o_file");

 





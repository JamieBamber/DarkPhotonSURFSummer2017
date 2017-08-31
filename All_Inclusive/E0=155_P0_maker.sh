#!/bin/bash

###########  ########

### /// This script produces histograms for the signal / gamma+jet background ratio vs different pT cuts for MET>155

### run root scripts 

export b_file="E0=155_P0_bkg.C"
export s_file="E0=155_P0_sig.C"
export data_file="E0=155_P0.root"
export o_file1="E0=155_P0_1.pdf"
export o_file2="E0=155_P0_2.pdf"

root -l << EOF
TFile::Open("GJets_Bkg_1pb_weighted_old3.root")
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
c1->SetLeftMargin(0.07);
c1->SetRightMargin(0.03);
TFile::Open("$data_file");
TH1F* h_ratio = new TH1F("h_ratio", "Signal to Background Ratio vs pT cut, MET>155GeV", 40, 0, 200);
for (Int_t i=1; i<=41; i++) {
	Double_t ratio = static_cast<Double_t>(h_sig->GetBinContent(i))/TMath::Sqrt(h_bkg->GetBinContent(i) + pow(0.1*h_bkg->GetBinContent(i),2));
	h_ratio->SetBinContent(i, ratio);
	cout << h_bkg->GetBinContent(i) << ", " << h_sig->GetBinContent(i) << ", " << ratio << endl;		// compute S/B ratio and output
}
h_ratio->GetXaxis()->SetTitle("p_{t} & MET cut threshold / GeV");
h_ratio->GetXaxis()->CenterTitle();
h_ratio->GetYaxis()->SetTitle("R = S/#sqrt{B + (0.1B)^{2}}");
h_ratio->GetYaxis()->CenterTitle();
h_ratio->SetTitleOffset(1.25,"X");
h_ratio->SetTitleOffset(1.35,"Y");
h_ratio->Draw("HIST");
TPaveStats *st1 = (TPaveStats*)h_ratio->FindObject("stats");
st1->SetX1NDC(0.6);
st1->SetY1NDC(0.72);
st1->SetX2NDC(0.8);
st1->SetY2NDC(0.88);
// ############ Save File
c1->SaveAs("$o_file1");
// ### Higgs PT histogram
TCanvas *c2 = new TCanvas();
c2->SetCanvasSize(1000,700);
// normalize
auto norm1 = h_bkg->GetBinContent(1);
h_bkg->Scale(100.0/norm1);
auto norm2 = h_sig->GetBinContent(1);
h_sig->Scale(100.0/norm2);
//
h_bkg->SetLineColor(2);
h_sig->SetLineColor(4);
//c2->SetLogy(1);
//gPad->SetLogy(1);
Double_t MAX=100;
Double_t MIN=0;
//MAX=TMath::Max(h_bkg->GetBinContent(h_bkg->GetMaximumBin()),h_sig->GetBinContent(h_sig->GetMaximumBin()));
//MIN=0.1*(h_sig->GetBinContent(h_sig->GetMinimumBin()));
h_bkg->SetAxisRange(MIN,MAX,"Y");
h_bkg->GetXaxis()->SetTitle("p_{t} & MET cut threshold / GeV");
h_bkg->GetYaxis()->SetTitle("efficiency / %");
h_bkg->SetTitleOffset(1.25,"X");
h_bkg->SetTitleOffset(1.25,"Y");
h_bkg->Draw("HIST");
gPad->Update();
//gPad->SetLogy(1);
h_sig->Draw("HIST SAME");
gPad->Update();
TPaveStats *st2 = (TPaveStats*)h_bkg->FindObject("stats");
st2->SetX1NDC(0.6);
st2->SetY1NDC(0.72);
st2->SetX2NDC(0.8);
st2->SetY2NDC(0.88);
auto legend = new TLegend(0.5,0.5,0.9,0.67);
legend->SetHeader("Key","C");
legend->AddEntry(h_sig,"Higgs signal");
legend->SetTextSize(0.035);
legend->AddEntry(h_bkg,"Jet+photon background");
legend->Draw();
// ############ Save File
c2->SaveAs("$o_file2");
EOF
 





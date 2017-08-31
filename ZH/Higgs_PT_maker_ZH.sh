#!/bin/bash

###########  ########

### /// This script produces histograms for the signal and gamma+jet background vs different Higgs pT cuts

### run root scripts 

export b_file="Higher_Higgs_pT_bkg_ZH.C"
export s_file="Higher_Higgs_pT_sig_ZH.C"
export o_file="Higgs_PT_ZH.pdf"

root -l << EOF
TFile::Open("ZGamma_Bkg_1pb_weighted.root")
.L $b_file
BkgTree t1;
t1.Loop();
// ############# Signal file next
TFile::Open("HiggsZH_modified.root")
.L $s_file
Delphes t2;
t2.Loop();
// ############# Draw Histograms
// ### Ratio histogram
TCanvas *c1 = new TCanvas();
c1->SetCanvasSize(1000,700);
TFile::Open("Higgs_PT_ZH.root");
TH1F* h_ratio = new TH1F("h_ratio", "Signal to Background Ratio (ZH) vs Higgs p_{T}", 50, 0, 500);
for (Int_t i=1; i<=51; i++) {
	Double_t ratio;
	if(h_bkg->GetBinContent(i)==0){
		ratio = 0;
	} else {
		ratio = static_cast<Double_t>(h_sig->GetBinContent(i))/TMath::Sqrt(h_bkg->GetBinContent(i) + pow(0.1*h_bkg->GetBinContent(i),2));
	}
	h_ratio->SetBinContent(i, ratio);
	cout << h_bkg->GetBinContent(i) << ", " << h_sig->GetBinContent(i) << ", " << ratio << endl;		// compute S/B ratio and output
}
h_ratio->GetXaxis()->SetTitle("Higgs p_{t} cut threshold / GeV");
h_ratio->GetXaxis()->CenterTitle();
h_ratio->GetYaxis()->SetTitle("R = S/#sqrt{B + (0.1B)^{2}}");
h_ratio->GetYaxis()->CenterTitle();
h_ratio->SetTitleOffset(1.25,"X");
h_ratio->SetTitleOffset(1.25,"Y");
h_ratio->Draw("HIST");
TPaveStats *st1 = (TPaveStats*)h_ratio->FindObject("stats");
st1->SetX1NDC(0.2);
st1->SetY1NDC(0.72);
st1->SetX2NDC(0.4);
st1->SetY2NDC(0.88);
// ############ Save File
c1->SaveAs("SB_ratio_Higgs_PT.pdf");
// ### Higgs PT histogram
TCanvas *c2 = new TCanvas();
c2->SetCanvasSize(1000,700);
// normalize
auto norm1 = h_bkg->Integral();
h_bkg->Scale(1.0/norm1);
norm1 = h_sig->Integral();
h_sig->Scale(1.0/norm1);
//
h_bkg->SetLineColor(2);
h_sig->SetLineColor(4);
//c2->SetLogy(1);
//gPad->SetLogy(1);
Double_t MAX;
Double_t MIN;
MAX=TMath::Max(h_bkg->GetBinContent(h_bkg->GetMaximumBin()),h_sig->GetBinContent(h_sig->GetMaximumBin()));
//MIN=0.1*(h_sig->GetBinContent(h_sig->GetMinimumBin()));
MIN=0
h_bkg->SetAxisRange(MIN,MAX,"Y");
h_bkg->GetXaxis()->SetTitle("Higgs p_{t} cut threshold / GeV");
h_bkg->GetYaxis()->SetTitle("normalized frequency");
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
legend->AddEntry(h_sig,"Higgs ZH signal");
legend->SetTextSize(0.035);
legend->AddEntry(h_bkg,"Z+#gamma background");
legend->Draw();
// ############ Save File
c2->SaveAs("$o_file");
EOF
 





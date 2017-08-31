#!/bin/bash

###################

### /// This script produces histograms showing signal and Gamma+Jet background from the r .C ROOT scripts

# © Jamie Bamber 2017

for TYPE in MT PT MET
do
export b_file="$TYPE.bkg.C"
export s_file="$TYPE.sig.C"
export o_file="$TYPE.both.pdf"

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
TCanvas *c1 = new TCanvas();
TFile::Open("$TYPE.root");
c1->SetCanvasSize(1000,700);
h_bkg->SetLineColor(2);
h_bkg->GetYaxis()->SetTitle("normalized freq. density");
//h_bkg->GetXaxis()->SetTitle("MET / GeV");
h_sig->SetLineColor(4);
// ##### Normalize histograms
Float_t binw1 = 250.0/200;
Float_t binw2 = 250.0/50;
auto norm1 = binw1*(h_bkg->Integral());
//if (norm1!=0) {
	h_bkg->Scale(1.0/norm1);
//}
auto norm2 = binw2*(h_sig->Integral());
//if (norm2!=0) {
	h_sig->Scale(1.0/norm2);
//}
// ########## 
Double_t MAX;
Double_t MIN;
MAX=TMath::Max( (h_bkg->GetBinContent(h_bkg->GetMaximumBin())),(h_sig->GetBinContent(h_sig->GetMaximumBin())) );
//MAX=0.02;
//MIN=h_sig->GetBinContent(h_sig->GetMinimumBin());
MIN=0;
h_bkg->SetAxisRange(MIN,MAX,"Y");
gPad->Update();
h_bkg->Draw("HIST");
gPad->Update();
h_sig->Draw("HIST SAME");
gPad->Update();
TPaveStats *st = (TPaveStats*)h_bkg->FindObject("stats");
st->SetX1NDC(0.65);
st->SetY1NDC(0.7);
st->SetX2NDC(0.85);
st->SetY2NDC(0.85);
auto legend = new TLegend(0.63,0.5,0.88,0.6);
legend->SetHeader("Key","C");
legend->AddEntry(h_sig,"Higgs signal (|#eta|<1.44)");
legend->SetTextSize(0.025);
legend->AddEntry(h_bkg,"Jet+photon background");
legend->Draw();
gPad->Update();
// ############ Save File
c1->SaveAs("$o_file")
EOF
done
 





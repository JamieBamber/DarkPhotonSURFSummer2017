#define Delphes_cxx
#include "Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

void Delphes::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   // Create Variables
   Float_t PT;
   Float_t E;
   Float_t Eta;
   Float_t Phi;
   Float_t M_T;
   Float_t Pi = TMath::Pi();
   //
   //Define Histograms
   TH1F* h2 = new TH1F("h2", "M_{#gamma#bar{#gamma}}^{T} (ZH)", 50, 0, 250);
   TH1F* h = new TH1F("h", "", 50, 0, 250);
   // 
   Float_t w = 5; // bin width
   //
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     
     // Fill Histogram 
     PT = TMath::MaxElement(Photon_,Photon_PT);
     E = MissingET_MET[0];
     Phi = Photon_Phi[TMath::LocMax(Photon_,Photon_PT)] - MissingET_Phi[0];
     Eta = Photon_Eta[TMath::LocMax(Photon_,Photon_PT)];
     if (Phi > Pi) {
     	 Phi = Phi - 2*Pi;
     } else if (Phi < -Pi) {
     	 Phi = Phi + 2*Pi;
     }
     M_T = sqrt(2*PT*E*(1-cos(Phi)));
     h->Fill(M_T);
     if ( Eta < 1.44 ) { // if eta < 1.44 fill second histogram
     	 h2->Fill(M_T);
     }
     // progress indicator
     if (((jentry+1)*100) % nentries == 0) {
     	 Float_t percent;
     	 percent = static_cast<Float_t>((jentry+1)*100)/nentries;
     	 printf("\r%.0f%% complete",percent);
     	 cout.flush();
     }
   }
   cout << endl;
   // Normalise histogram
   auto norm1 = h->Integral();
   h->Scale(1.0/norm1);
   auto norm2 = h2->Integral();
   TCanvas *c1 = new TCanvas();
   h2->Scale(1.0/norm2);
   // Draw Histogram
   h->Draw("HIST");
   h2->Draw("HIST");
   h2->GetXaxis()->SetTitle("M_{#gamma#bar{#gamma}}^{T} / GeV");
   h2->GetXaxis()->CenterTitle();
   h2->GetYaxis()->SetTitle("normalised frequency");
   h2->GetYaxis()->CenterTitle();
   h2->SetLineColor(2);
   c1->Update();
   gPad->Update();
   TPaveStats *st = (TPaveStats*)h2->FindObject("stats");
   st->SetX1NDC(0.6);
   st->SetY1NDC(0.72);
   st->SetX2NDC(0.8);
   st->SetY2NDC(0.88);
   h2->Draw("HIST");
   gPad->Update();
   h->Draw("HIST SAME");
   // 
   auto legend = new TLegend(0.6,0.5,0.8,0.67);
   legend->SetHeader("Key","C");
   legend->AddEntry(h,"no |#eta| cut");
   legend->AddEntry(h2,"|#eta| < 1.44");
   legend->Draw();
}


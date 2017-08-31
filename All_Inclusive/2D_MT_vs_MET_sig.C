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
   Float_t Phi;
   Float_t M_T;
   Float_t Pi = TMath::Pi();
   Int_t N = 1050; // number of events 
   //
   //Define Histogram
   TH2F* h = new TH2F("h", "M_{#gamma#bar{#gamma}}^{T} vs MET", 50, 0, 250, 50, 0, 250);
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
     if (Phi > Pi) {
     	 Phi = Phi - 2*Pi;
     } else if (Phi < -Pi) {
     	 Phi = Phi + 2*Pi;
     }
     M_T = sqrt(2*PT*E*(1-cos(Phi)));
     h->Fill(E, M_T);
     //
     // progress indicator
     if (((jentry+1)*100) % nentries == 0) {
     	 Float_t percent;
     	 percent = static_cast<Float_t>((jentry+1)*100)/nentries;
     	 printf("\r%.0f%% complete",percent);
     	 cout.flush();
     }
   }
   // Draw Histogram
   TCanvas *c1 = new TCanvas();
   h->Draw();
   h->GetXaxis()->SetTitle("MET / GeV");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->SetTitle("M_{#gamma#bar{#gamma}}^{T} / GeV");
   h->GetYaxis()->CenterTitle();
   h->GetZaxis()->SetTitle("frequency");
   h->GetZaxis()->CenterTitle();
   c1->SetLogz(1);
   c1->Update();
   gPad->SetRightMargin(0.14);
   gPad->Update();
   TPaveStats *st = (TPaveStats*)h->FindObject("stats");
   st->SetX1NDC(0.67);
   st->SetY1NDC(0.81);
   st->SetX2NDC(0.85);
   st->SetY2NDC(0.98);
   //gPad->Update();
   //TPaletteAxis *palette = (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
   h->Draw("COlZ");
}

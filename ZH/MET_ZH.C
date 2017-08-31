#define Delphes_cxx
#include "Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

void Delphes::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Delphes.C
//      root> Delphes t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   // Create Variables
   Float_t PT;
   Float_t E;
   Float_t Phi;
   Float_t Eta;
   Float_t M_T;
   Float_t Pi = TMath::Pi();
   //
   //Define Histogram
   TH1F* h = new TH1F("h", "MET (ZH)", 50, 0, 250);
   TH1F* h2 = new TH1F("h2", "", 50, 0, 250);
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
     h->Fill(E);
     if ( Eta < 1.44 ) { // if eta < 1.44 fill second histogram
     	 h2->Fill(E);
     }
     //
     // progress indicator
     if (((jentry+1)*100) % nentries == 0) {
     	 Float_t percent;
     	 percent = static_cast<Float_t>((jentry+1)*100)/nentries;
     	 printf("\r%.0f%% complete",percent);
     	 cout.flush();
     }
   }
   // Normalise histogram
   auto norm = h->Integral();
   h->Scale(1.0/norm);
   norm = h2->Integral();
   h2->Scale(1.0/norm);
   // Draw Histogram
   TCanvas *c1 = new TCanvas();
   h->Draw("HIST");
   h->GetXaxis()->SetTitle("MET / GeV");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->SetTitle("normalised frequency");
   h->GetYaxis()->CenterTitle();
   h2->SetLineColor(2);
   c1->Update();
   gPad->Update();
   TPaveStats *st = (TPaveStats*)h->FindObject("stats");
   st->SetX1NDC(0.6);
   st->SetY1NDC(0.72);
   st->SetX2NDC(0.8);
   st->SetY2NDC(0.88);
   h->Draw("HIST");
   h2->Draw("HIST SAME");
   // 
   auto legend = new TLegend(0.6,0.5,0.8,0.67);
   legend->SetHeader("Key","C");
   legend->AddEntry(h,"no |#eta| cut");
   legend->AddEntry(h2,"|#eta| < 1.44");
   legend->Draw();
}

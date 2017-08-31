#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"
/*
#########
Script to create a histogram of the MT variable & save it to a root file 

© Jamie Bamber 2017
#########
*/

void BkgTree::Loop(Int_t index=0)
{
//   In a ROOT session, you can do:
//      root> .L event.C
//      root> event t
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
   Float_t M_T;
   Float_t Pi = TMath::Pi();
   Float_t Lumin = 100; // luminosity in fb^-1
   //
   //Define Histogram for MT solely
   TH1F* h_bkg = new TH1F("h_bkg", "M_{#gamma#bar{#gamma}}^{T}, signal and #gamma+jet background", 200, 0, 250);
   // 
   Float_t w = 5; // bin width

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
     
      // Fill Histogram 
     PT = PhotonPt;
     E = MET;
     Phi = PhotonPhi - METPhi;
     if (Phi > Pi) {
     	 Phi = Phi - 2*Pi;
     } else if (Phi < -Pi) {
     	 Phi = Phi + 2*Pi;
     }
     M_T = sqrt(2*PT*E*(1-cos(Phi)));
     h_bkg->Fill(M_T, weight*1000*Lumin);
     // progress indicator
     Float_t percent;
     Int_t number;
     percent = static_cast<Float_t>((jentry+1)*100)/nentries;
     number = int(10*percent);
     if ((number % 10) == 0) {
     	 printf("\r%.0f%% complete",percent);
     	 cout.flush();
     }
   }
   cout << endl;
   // Normalise histogram
   //auto norm = h_bkg->Integral();
   //h_bkg->Scale(1.0/norm);
  
   // Make Histogram
   h_bkg->GetXaxis()->SetTitle("M_{#gamma#bar{#gamma}}^{T} / GeV");
   h_bkg->GetXaxis()->CenterTitle();
   h_bkg->GetYaxis()->SetTitle("frequency");
   h_bkg->GetYaxis()->CenterTitle();
   /*gPad->Update();
   TPaveStats *st = (TPaveStats*)h_bkg->FindObject("stats");
   st->SetX1NDC(0.6);
   st->SetY1NDC(0.72);
   st->SetX2NDC(0.8);
   st->SetY2NDC(0.88);*/
   TFile *hfile = new TFile("MT.root", "RECREATE");
   h_bkg->Write();
   hfile->Close();
}

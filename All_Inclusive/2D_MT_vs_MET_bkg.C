#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void BkgTree::Loop(Int_t index=0)
{
//   In a ROOT session, you can do:
//      root> .L BkgTree.C
//      root> BkgTree t
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
     PT = PhotonPt;
     E = MET;
     Phi = PhotonPhi - METPhi;
     if (Phi > Pi) {
     	 Phi = Phi - 2*Pi;
     } else if (Phi < -Pi) {
     	 Phi = Phi + 2*Pi;
     }
     M_T = sqrt(2*PT*E*(1-cos(Phi)));
     h->Fill(E, M_T, weight);
     //
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


#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

/*
##########
Script to produce a weighted output of events that pass the cuts for the background

© Jamie Bamber 2017
##########
*/

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
   Float_t PT;				// photon transverse momentum 
   Float_t E;				// Missing transverse energy (MET)
   Float_t Lumin=100;			// luminosity in fb^-1 (standard is 1pb^-1)
   Float_t E_0=150;
   Float_t PT_0=70;
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   Float_t MET_sgnf;		// MET significance
   Float_t DR;				// Delta R displacement variable 
   Float_t Pi = TMath::Pi();
   //
   Long64_t nbytes = 0, nb = 0;
   Int_t Ns = 0;
   Double_t Nw = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		// ## Trigger cuts
		
		if (no_jets >= 2) continue;			// no more than 1 jet (with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5)
		
		//if (no_jets==1) {
		//	if (D_phi_g_jet>2.5) continue; 	// Δφ(γ,jet) < 2.5 
		//}
		
		//if (HT>100) continue;
		
		//if (MET_sgnf < 3) continue; // MET significance

		//if ( ((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) continue; // least stringent eta cut
		
		if (Eta > 1.44) continue; // more stringent eta cut
		
		if (PT < PT_0) continue;	// PT criteria 
	
		if (E < E_0) continue;	// MET criteria 	
		
		Phi = PhotonPhi - METPhi;
		 if (Phi > Pi) {
			 Phi = Phi - 2*Pi;
		 } else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		 }
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ((MT < 100) || (MT > 130)) continue;	// Transverse Mass cut 
		
		if (Iso_lepton2==1) continue; 		    // Isolated lepton criterion
			 
		Ns++;
		Nw = Nw + weight;
		
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
   //
   //
   // print out results
   cout << "Background" << endl;
   cout << "Selected events (raw) , Selected events (weighted)" << endl;
   cout << Ns << ", " << Nw*1000*Lumin << endl;
}

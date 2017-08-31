#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

/* 
###########
This script outputs the data on selected events for successive numbers of cuts, acting on the Gamma+Jet background data produced from DarkPhotonAnalyzer.cc RazorAnalyzer. BkgTree->MakeClass() should be 
applied to generate BkgTree.h and BkgTree.C files. The number of cuts applied is given as an integer input "index". The script can be used in conjuction with "TrigTable_maker.sh"

This version my optimised cuts (or the most optimal cuts I found for my project) for use on the γ+jet background 

© Jamie Bamber 2017
###########
*/ 

void BkgTree::Loop(Int_t index = 0)
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
   Float_t E_0 = 150;		// energy threshold
   Float_t Lumin=100;		// luminosity in fb^-1 (standard is 1pb^-1)
   Float_t PT_0 = 70;
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   //
   Float_t D_phi_g_jet;		// Δφ(γ,jet)
   Int_t no_jets;			// no jets with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5
   Float_t HT;				// scalar sum of transverse momentum of jets that satisfy the conditions
   //
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of BkgTrees 
   Double_t N_w[6];			// number of selected BkgTrees (weighted)
   Int_t N_s[6];			// number of selected BkgTrees (raw)
   Float_t Eff;				// efficiency
   //
   N_s[index] = 0;
   N_w[index] = 0;
   PT_0 = E_0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		// ###### My calculations 
		// Variable assignment
		E = MET;
		PT = PhotonPt; 
		Eta = TMath::Abs(PhotonEta);
		
		Phi = PhotonPhi - METPhi;
		if (Phi > Pi) {
			Phi = Phi - 2*Pi;
		} else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		}
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		// Trigger cuts
		
		//if ( (((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) && (index==1) ) continue; // least stringent eta cut
		
		if ( (Eta > 1.44) && (index>=1) ) continue; 				// more stringent eta cut
		
		if ( ((MT < 100) || (MT > 130)) && (index>=2) ) continue;	// Transverse Mass cut
		
		if ( (Iso_lepton2==1) && (index>=3) ) continue; 			// Isolated lepton cut
		
		if ( (no_jets >= 2) && (index>=4) ) continue;				// no more than 1 jet (with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5)
		
		if ( (no_jets==1) && (index>=5) ) {
			if (D_phi_g_jet>2.5) continue; 							// Δφ(γ,jet) < 2.5 
		}
		
		if ( (HT>100) && (index>=6) ) continue;						// require HT (sum of Jet pT) to be less than 100GeV
		
		if ( (PT < PT_0) && (index>=7) ) continue;	// Photon PT cut
	
		if ( (E < E_0) && (index>=8) ) continue;	// MET cut	
			 
		N_w[index] = N_w[index] + weight;							// increment element by weight
		N_s[index] = N_s[index] + 1;
		//fail:;													// goto this point if the event fails the test
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
   // Print out data 
   // outputs: "no. of selected events (raw)" , "no. of selected events (weighted)" ,
   cout << N_s[index] << ", " << N_w[index]*1000*Lumin << ", " << endl;
}

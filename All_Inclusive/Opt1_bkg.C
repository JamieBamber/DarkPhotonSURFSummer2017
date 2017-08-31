#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

/*
##########
Optimisation script for the γ+jet background

used with Opt_1.sh

© Jamie Bamber 2017
##########
*/


void BkgTree::Loop(Float_t PT0=50, Float_t E0=50, Float_t MT_UP=130, Float_t MT_LW=100, Bool_t Iso_Leptons=true, Float_t Eta0=1.44,Float_t HPT0=0)
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
   Float_t PT;				// photon transverse momentum 
   Float_t E;				// Missing transverse energy (MET)
   //
   Float_t Lumin=100;		// Luminosity
   Float_t sigma=50000;		// cross section in fb
   Float_t BR=0.01;			// Branching ratio
   Float_t Event_weight = Lumin*sigma*BR/50000;
   //

   //
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of events 
   Float_t bin_width = 5;
   Int_t nbins = 40;	
   //
   //Define Histogram 
   TH2F* h_bkg = new TH2F("h_sig", "", nbins, 0, nbins*bin_width, nbins, 0, nbins*bin_width);
   //
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		// My calculations 
		E = MET;
		PT = PhotonPt;
		Eta = PhotonEta;
		
		/*TVector3 photon_v;
		photon_v.SetPtEtaPhi(PT,Eta,PhotonPhi);
		TVector3 MET_v;
		MET_v.SetPtEtaPhi(E,0,METPhi);
		
		TVector3 Higgs_v = photon_v + MET_v;
		Float_t HPT = TMath::Abs(Higgs_v.Perp());
		
		sum_ET = ScalarHT_HT[0];*/
		
		// ### Trigger Cuts
		
		// ### Baseline cuts
		if (Eta > Eta0) continue; // more stringent eta cut
		
		if (no_jets >= 2) continue;			// no more than 1 jet (with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5)
		
		if (no_jets==1) {
			if (D_phi_g_jet>2.5) continue; 	// Δφ(γ,jet) < 2.5 
		}
		
		if ((Iso_Leptons) && (Iso_lepton2==1)) continue; // Isolated lepton cut
		
		// #### Cuts to be optimised 
		//if (HT>100) continue;				// require HT (sum of Jet pT) to be less than 100GeV
		
		//if (PT < PT_0) continue;	// PT criteria 
	
		//if (E < E_0) continue;	// MET criteria 	
		
		Phi = PhotonPhi - METPhi;
		 if (Phi > Pi) {
			 Phi = Phi - 2*Pi;
		 } else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		 }
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ((MT < MT_LW) || (MT > MT_UP)) continue;	// Transverse Mass cut 
		
		// Fill Histogram
		Int_t PT_min_bin;
		if((PT >= 0) && (PT <= nbins*bin_width)) {
			PT_min_bin = TMath::FloorNint(PT/bin_width)+1;
		} else {
			PT_min_bin = nbins + 1;
		}
		Int_t E_min_bin;
		if((E >= 0) && (E <= nbins*bin_width)) {
			E_min_bin = TMath::FloorNint(E/bin_width)+1;
		} else {
			E_min_bin = nbins + 1;
		}
		for(Int_t i=1; i<=nbins+1; i++) {
			if (i <= E_min_bin) {
				for(Int_t j=1; j<=nbins+1; j++) {
					if (j <= PT_min_bin) {
						h_bkg->SetBinContent(i,j,h_bkg->GetBinContent(i,j)+1000*Lumin*weight);
					}
				}
			}
		}
		fail:;
		//
		// progress indicator
		 Float_t percent;
		 Int_t number;
		 percent = static_cast<Float_t>((jentry+1)*100)/nentries;
		 number = int(10*percent);
		 if ((number % 10) == 0) {
			 printf("\r%.0f%% complete Bkg, MT_UP = %.0f, MT_LW = %.0f",percent, MT_UP, MT_LW);
			 cout.flush();
		}
   }
   cout << endl;
   //
   //
   /*cout << "Signal:" << endl;
   for(Int_t i=1; i<=(nbins+1); i++) {
   	   cout << h_sig->GetBinContent(i) << endl;
   }*/
   //
   TFile *hfile = new TFile("Opt_1.root", "RECREATE");
   h_bkg->Write();
   hfile->Close();
}

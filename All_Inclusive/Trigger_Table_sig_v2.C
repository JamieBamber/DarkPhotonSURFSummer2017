#define Delphes_cxx
#include "Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

/* 
###########
This script outputs the data on selected events for successive numbers of cuts, acting on the signal data produced from Delphes. Delphes->MakeClass() should be 
applied to generate Delphes.h and Delphes.C files. The number of cuts applied is given as an integer input "index". The script can be used in conjuction with "TrigTable_maker.sh"

This version my optimised cuts (or the most optimal cuts I found for my project) for use on the all inclusive Higgs signal

© Jamie Bamber 2017
###########
*/ 
void Delphes::Loop(Int_t index)
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
   Float_t Lumin=100;		 // luminosity in fb^-1
   Float_t sigma=50000; 	// cross section in fb for the all inclusive channel
   Float_t BR=0.01;			 // Branching ratio
   Float_t Event_weight = Lumin*sigma*BR/50000;
   //
   Float_t E_0 = 150;		// limit for MET
   Float_t PT_0 = 70;		// limit for PT
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   Float_t Phi_l;			// lepton phi
   Float_t Eta_l; 			// lepton eta
   Float_t DR;				// Delta R displacement variable 
   //
   Float_t D_phi_g_jet;		// Δφ(γ,jet)
   Int_t no_jets;			// no jets with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5
   Float_t HT;				// scalar sum of transverse momentum of jets that satisfy the conditions
   //
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of events 
   Float_t N_s[7];			// number of selected events (array for different cuts) 
   //
   N_s[index] = 0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (ientry < 0) break;
		// if (Cut(ientry) < 0) continue;
		
		// ###### My calculations 
		// Variable assignment
		E = MissingET_MET[0];
		PT = TMath::MaxElement(Photon_,Photon_PT);
		Eta = TMath::Abs(Photon_Eta[TMath::LocMax(Photon_,Photon_PT)]);
		
		Phi = Photon_Phi[TMath::LocMax(Photon_,Photon_PT)] - MissingET_Phi[0];
		if (Phi > Pi) {
			Phi = Phi - 2*Pi;
		} else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		}
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		HT = 0;
		no_jets=0;
		for (int i=0; i<Jet_; i++) {
			D_phi_g_jet = TMath::Abs(Jet_Phi[i] - PhotonPhi);
			if (D_phi_g_jet>Pi) D_phi_g_jet = 2*Pi - D_phi_g_jet;
			if( (Jet_PT[i]>30) && (TMath::Abs(Jet_Eta[i])<2.4) && (D_phi_g_jet>0.5) ) {
				HT = HT + Jet_PT[i];	
				no_jets++;
			}
		}
		
		// Trigger cuts
		
		//if ( (((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) && (index==1) ) continue; // least stringent eta cut
		
		if ( (Eta > 1.44) && (index>=1) ) continue; 				// more stringent eta cut
		
		if ( ((MT < 100) || (MT > 130)) && (index>=2) ) continue;	// Transverse Mass cut
		
		float e_Eta;
		if (index >= 3) {	
			for (Int_t j=0; j<Electron_; j++) {	
				Phi_l = TMath::Abs(Electron_Phi[j] - Photon_Phi[TMath::LocMax(Photon_,Photon_PT)]);
				if (Phi_l > Pi) {
					Phi_l = 2*Pi - Phi_l;
				}
				Eta_l = Electron_Eta[j] - Photon_Eta[TMath::LocMax(Photon_,Photon_PT)];
				DR = sqrt(pow(Phi_l,2) + pow((Eta_l - Eta),2));
				e_Eta = Electron_Eta[j];
				if (DR > 0.3) {																							// Isolated electron cut
					if ((Electron_PT[j]>10) && ( ((e_Eta < 1.4442) || (e_Eta > 1.566)) && (e_Eta<2.5) ) ) goto fail;
				}
			}
			for (Int_t j=0; j<Muon_; j++) {	
				Phi_l = Muon_Phi[j] - Photon_Phi[TMath::LocMax(Photon_,Photon_PT)];
				if (Phi_l > Pi) {
					Phi_l = 2*Pi - Phi_l;
				}
				Eta_l = Muon_Eta[j] - Photon_Eta[TMath::LocMax(Photon_,Photon_PT)];
				DR = sqrt(pow(Phi_l,2) + pow(Eta_l,2));
				if (DR > 0.3) {																							// Isolated muon cut
					if ((Muon_PT[j]>10) && (TMath::Abs(Muon_Eta[j])<2.1)) goto fail;
				}
			}	
		}
		
		if ( (no_jets >= 2) && (index>=4) ) continue;				// no more than 1 jet (with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5)
		
		if ( (no_jets==1) && (index>=5) ) {
			if (D_phi_g_jet>2.5) continue; 							// Δφ(γ,jet) < 2.5 
		}
		
		if ( (HT>100) && (index>=6) ) continue;						// require HT (sum of Jet pT) to be less than 100GeV
		
		if ( (PT < PT_0) && (index>=7) ) continue;	// Photon PT cut
	
		if ( (E < E_0) && (index>=8) ) continue;	// MET cut	
		
		N_s[index] = N_s[index] + 1;								// increment element 
		fail:;													// goto this point if the event fails the test
		// progress indicator
		if (((jentry+1)*100) % nentries == 0) {
			Float_t percent;
			percent = static_cast<Float_t>((jentry+1)*100)/nentries;
			printf("\r%.0f%% complete",percent);
			cout.flush();
		}
   }
   //cout << endl;
   // Print out data 
   // outputs: "no. of selected events (raw)" , "no. of selected events (weighted)"
   cout << N_s[index] << ", " << N_s[index]*Event_weight << ", " << 100*N_s[index]/N << ", " << endl;
}

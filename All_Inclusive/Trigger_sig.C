#define Delphes_cxx
#include "Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"
/*
##########
Script to produce a 2D Histogram of the signal weighted selected events, after applying cuts of |eta|<1.44, 100<MT<130, plotting the PT>... cut vs 
the MET>... cut. (i.e. PT > [y value] & MET > [x value])

used with 2D_E0_P0_maker.sh

© Jamie Bamber 2017
##########
*/


void Delphes::Loop(Int_t index=0)
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
   Float_t Phi_l;			// lepton phi
   Float_t Eta_l;			// lepton Eta
   Float_t DR;				// DeltaR = Sqrt( DeltaEta^2 + DeltaPhi^2 ) 
   //
   Float_t E_0 = 150;
   Float_t PT_0 = 70;
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   Float_t PhotonPhi;
   Float_t METPhi;
   Float_t MET_sgnf;		// MET significance 
   Float_t sum_ET;			// sum of transverse energy
   Float_t D_phi_g_jet;		// Δφ(γ,jet)
   Int_t no_jets;			// no jets with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5
   Float_t HT;				// scalar sum of transverse momentum of jets that satisfy the conditions
   Float_t Pi = TMath::Pi();	
   //
   PT_0 = E_0;
   //
   Int_t Ns = 0; 
   //
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
		
		// My calculations 
		E = MissingET_MET[0];
		PT = TMath::MaxElement(Photon_,Photon_PT);
		Eta = TMath::Abs(Photon_Eta[TMath::LocMax(Photon_,Photon_PT)]);
		PhotonPhi = Photon_Phi[TMath::LocMax(Photon_,Photon_PT)];
		METPhi = MissingET_Phi[0]; 
		//
		//
		TVector3 photon_v;
		photon_v.SetPtEtaPhi(PT,Eta,PhotonPhi);
		TVector3 MET_v;
		MET_v.SetPtEtaPhi(E,0,METPhi);
		//
		TVector3 Higgs_v = photon_v + MET_v;
		Float_t HPT = TMath::Abs(Higgs_v.Perp());
		
		sum_ET = ScalarHT_HT[0];
     
		MET_sgnf = E/sqrt(TMath::Abs(sum_ET));
		
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
		
		// ### Trigger Cuts
		
		//if (MET_sgnf < 3) continue; // MET significance cut
		
		//if ( ((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) continue; // least stringent eta cut
		
		if (Eta > 1.44) continue; // more stringent eta cut
		
		if (no_jets >= 2) continue;			// no more than 1 jet (with pT>30GeV |η|<2.4 Δφ(γ,jet)>0.5)
		
		if (no_jets==1) {
			if (D_phi_g_jet>2.5) continue; 	// Δφ(γ,jet) < 2.5 
		}
		
		if (HT>100) continue;				// require HT (sum of Jet pT) to be less than 100GeV
		
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
			 	
		float e_Eta;	
		for (Int_t j=0; j<Electron_; j++) {	
			Phi_l = TMath::Abs(Electron_Phi[j] - Photon_Phi[TMath::LocMax(Photon_,Photon_PT)]);
			if (Phi_l > Pi) {
				Phi_l = 2*Pi - Phi_l;
			}
			Eta_l = Electron_Eta[j] - Photon_Eta[TMath::LocMax(Photon_,Photon_PT)];
			DR = sqrt(pow(Phi_l,2) + pow((Eta_l - Eta),2));
			e_Eta = Electron_Eta[j];
			if (DR > 0.3) {																							// Isolated electron criterion
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
			if (DR > 0.3) {																							// Isolated muon criterion
				if ((Muon_PT[j]>10) && (TMath::Abs(Muon_Eta[j])<2.1)) goto fail;
			}
		}		
		
		Ns++;					
		 
		fail:;
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
   //
   //
   cout << "Signal" << endl;
   cout << "Selected Events (raw) , Selected events (weighted)" << endl;
   cout << Ns << ", " << Ns*Event_weight << endl;
}

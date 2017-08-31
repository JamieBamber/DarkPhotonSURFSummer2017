#define Delphes_cxx
#include "Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"
/*
##########
Script to produce a Histogram of the signal weighted selected events, with photon |η|<1.44, 100<MT<130, no isolated leptons
MET>155 vs photon pT cut

used with Higgs_PT_maker.sh

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
   Float_t E_0 = 155;
   Float_t PT_0;
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   Float_t PhotonPhi;
   Float_t METPhi;
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of events 
   Float_t bin_width = 5;
   Int_t nbins = 40;	
   //
   //Define Histogram 
   TH1F* h_sig = new TH1F("h_sig", "", nbins, 0, nbins*bin_width);
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
		
		if ( ((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) continue; // least stringent eta cut
		
		if (Eta > 1.44) continue; // more stringent eta cut
		
		//if (PT < PT_0) continue;	// PT criteria 
	
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
		
		//h_sig->Fill(HPT, Event_weight);						
		
		// Fill Histogram
		Int_t min_bin;
		if((PT >= 0) && (PT <= nbins*bin_width)) {
			min_bin = TMath::FloorNint(PT/bin_width)+1;
		} else {
			min_bin = nbins + 1;
		}
		for(Int_t i=1; i<=nbins+1; i++) {
			if (i <= min_bin) {	
				h_sig->SetBinContent(i,h_sig->GetBinContent(i)+Event_weight);
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
			 printf("\r%.0f%% complete",percent);
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
   TFile *hfile = new TFile("E0=155_P0.root", "UPDATE");
   h_sig->Write();
   hfile->Close();
}

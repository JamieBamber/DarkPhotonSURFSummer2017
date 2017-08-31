#define Delphes_cxx
#include "Delphes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

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
   Float_t Lumin=100;
   Float_t sigma=3.779*1000; // cross section in fb for the VBF channel
   Float_t BR=0.01;
   Float_t Event_weight = Lumin*sigma*BR/50000;
   //
   Float_t E_0 = 30;
   Float_t PT_0;
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable
   Float_t Phi_l;			// lepton phi
   Float_t Eta_l; 			// lepton eta
   Float_t DR;				// Delta R displacement variable 
   //
   Float_t Copy_Jet_PT[kMaxJet];	// copy of Jet_PT array
   Int_t j1_index;					// j1 index
   Int_t j2_index;					// j2 index
   Float_t j1_Eta;					// j1 eta
   Float_t j2_Eta;					// j2 eta
   Float_t j1_PT;					// j1 PT
   Float_t j2_PT;					// j2 PT
   //
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of events 
   Float_t N_s[7];			// number of selected events (array for different cuts) 
   //
   N_s[index] = 0;
   PT_0 = E_0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (ientry < 0) break;
		// if (Cut(ientry) < 0) continue;
		
		// My calculations 
		E = MissingET_MET[0];
		PT = TMath::MaxElement(Photon_,Photon_PT);
		Eta = TMath::Abs(Photon_Eta[TMath::LocMax(Photon_,Photon_PT)]);
		if ( (((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) && (index==1) ) continue; // least stringent eta cut
		
		if ( (Eta > 1.44) && (index>=2) ) continue; // least stringent eta cut
		
		if ( (PT < PT_0) && (index>=3) ) continue;	// Photon criteria 
	
		if ( (E < E_0) && (index>=4) ) continue;	// Photon criteria 	
		
		Phi = Photon_Phi[TMath::LocMax(Photon_,Photon_PT)] - MissingET_Phi[0];
		if (Phi > Pi) {
			Phi = Phi - 2*Pi;
		} else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		}
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ( ((MT < 100) || (MT > 130)) && (index>=5) ) continue;	// Transverse momentum criterion 
			 
		float e_Eta;
		if (index >= 6) {	
			for (Int_t j=0; j<Electron_; j++) {	
				Phi_l = TMath::Abs(Electron_Phi[j] - Photon_Phi[TMath::LocMax(Photon_,Photon_PT)]);
				if (Phi_l > Pi) {
					Phi_l = 2*Pi - Phi_l;
				}
				Eta_l = Electron_Eta[j] - Photon_Eta[TMath::LocMax(Photon_,Photon_PT)];
				DR = sqrt(pow(Phi_l,2) + pow((Eta_l - Eta),2));
				e_Eta = Electron_Eta[j];
				if (DR > 0.3) {									// Isolated electron criterion
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
				if (DR > 0.3) {									// Isolated muon criterion
					if ((Muon_PT[j]>10) && (TMath::Abs(Muon_Eta[j])<2.1)) goto fail;
				}
			}	
		}
		
		// ## Find two jets with highest pT: j1, j2 with j1_pT > j2_pT
		// make array copy of Jet
		if (index>=7) {
			if (Jet_>=2) {
				for(int i=0; i<Jet_; i++) {
					Copy_Jet_PT[i] = Jet_PT[i];	
				}
				// find indices
				j1_index = TMath::LocMax(Jet_,Jet_PT);	// find j1 index
				Copy_Jet_PT[j1_index] = 0;				// set j1 value in Copy to zero
				j2_index = TMath::LocMax(Jet_,Copy_Jet_PT);	// find j2 index*/
				// assign variables
				j1_Eta = Jet_Eta[j1_index];
				j2_Eta = Jet_Eta[j2_index];
			
				if ( ((j1_Eta*j2_Eta)>0) || (TMath::Abs(j1_Eta-j2_Eta)<4.0) ) continue; // apply rapidity cuts*/
					
			} else if (Jet_<2) {
				continue;
			}
		}
		
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
   cout << endl;
   // Print out data in csv format
   // using eff = no. selected / no. events 
   cout << N_s[index] << ", " << N_s[index]*Event_weight << ", " << 100*N_s[index]/N << ", " << endl;
}

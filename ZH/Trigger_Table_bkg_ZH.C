#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

void BkgTree::Loop(const Int_t index)
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
   //
   Float_t Eta;				// Psuedo-rapidity of the photon
   Float_t Phi;				// difference in azimuthal angle from the photon to the MET
   Float_t MT;				// Transverse momentum variable 
   //
   Float_t Me=0.5109989461*0.001;	// electron invarient mass in GeV
   Float_t Mm=105.6583745*0.001;	// muon invarient mass in GeV
   Float_t Mz=91.1876;				// Z boson invarient mass in GeV
   TVector3 l1;
   TVector3 l2;
   TVector3 ZCandidate;
   TVector3 photon_vec;
   TVector3 MET_vec;
   Float_t dPhi_ll_EE;				// azimuthal angle between the sum of the lepton Pt and the sum of the MET and photon Pt
   Float_t dPhi_ll;					// azimuthal angle between leptons
   Float_t p_variable;
   Float_t DeltaM;
   //
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of events 
   Int_t N_s[12];			// number of selected events (raw)
   Double_t N_w[12];			// number of selected events (weighted)
   //
   for(int n=0; n<12; n++) {
   	   N_s[n] = 0;
   	   N_w[n] = 0;
   }
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if (ientry < 0) break;
		// if (Cut(ientry) < 0) continue;
		
		// My calculations 
		E = MET;
		PT = PhotonPt; 
		Eta = TMath::Abs(PhotonEta);
		
		// #### Trigger cuts
		//if ( (((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) && (index==1) ) continue; // least stringent eta cut
		
		//if ( (Eta > 1.44) && (index>=2) ) continue; // barrel region eta cut
		
		if ( (PT < 20) && (index>=3) ) continue;	// Photon PT cut ### 20 value taken from CMS paper
	
		if ( (E < 60) && (index>=4) ) continue;	// MET cut	### 60 value taken from CMS paper
		
		Phi = PhotonPhi - METPhi;
		if (Phi > Pi) {
			Phi = Phi - 2*Pi;
		} else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		}
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ( ((MT < 100) || (MT > 130)) && (index>=5) ) continue;	// Transverse momentum criterion ### kept over from Gabrielli et al. paper
		
		//
		
		// ## further cuts also from CMS paper
		
		if( index>=6 ) {
			l1.SetPtEtaPhi(lep1Pt,lep1Eta,lep1Phi);
			l2.SetPtEtaPhi(lep2Pt,lep2Eta,lep2Phi);
			ZCandidate = l1 + l2;
			if(Zdecay == 0) continue;
		}
		
		if ( (index>=7) && (Jet_PT_max > 30) ) continue; // no jets with PT greater than 30GeV cut
		
		if (index>=8) {														// first Delta Phi cut
			photon_vec.SetPtEtaPhi(PT,PhotonEta,PhotonPhi);
			MET_vec.SetPtEtaPhi(E,0,METPhi);
			dPhi_ll_EE = TMath::Abs( (photon_vec + MET_vec).Phi() - ZCandidate.Phi() );
			if ( dPhi_ll_EE > Pi) dPhi_ll_EE = 2*Pi - dPhi_ll_EE;
			if ( dPhi_ll_EE < 2.7) continue;
		}	
		
		if (index>=9) {		// pt variable cut
			p_variable = TMath::Abs( (photon_vec+MET_vec).Pt() - ZCandidate.Pt() )/ZCandidate.Pt();
			if (p_variable > 0.5) continue;
		}
		
		if (index>=10) {		// azimuthal angle between leptons cut
			dPhi_ll = TMath::Abs(lep1Phi-lep2Phi);
			if ( dPhi_ll > Pi ) dPhi_ll = 2*Pi - dPhi_ll;
			if ( dPhi_ll > 2.25) continue;
		}
		
		if (index>=11) {		// lepton combined pT cut
			if ( ZCandidate.Pt() < 60) continue;
		}
		
		N_s[index] = N_s[index] + 1;								// increment element
		N_w[index] = N_w[index] + weight;							// increment element
		
		fail:;													// goto this point if the event fails the test
		// progress indicator
		/*if (((jentry+1)*100) % nentries == 0) {
			Float_t percent;
			percent = static_cast<Float_t>((jentry+1)*100)/nentries;
			printf("\r%.0f%% complete",percent);
			cout.flush();
		}*/
   }
   //cout << endl;
   // Print out data in csv format
   // using eff = no. selected / no. events
   for(int n=index; n<=index; n++) {
   	   cout << N_s[n] << ", " << N_w[n]*Lumin*1000 << ", " << 100*N_s[n]/N << ", " << endl;
   }
}

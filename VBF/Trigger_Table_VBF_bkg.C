#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

void BkgTree::Loop(Int_t index = 0)
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
   //cout << nentries << endl;

   // Create Variables
   Float_t PT;				// photon transverse momentum 
   Float_t E;				// Missing transverse energy (MET)
   //
   Float_t Lumin=100;
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
   Float_t Pi = TMath::Pi();
   Int_t N = nentries; 		// number of events 
   Int_t N_s[7];			// number of selected events (array for different cuts) - raw
   Long_t N_l[7];
   Double_t N_w[7]; 			// number of selected events (array for different cuts) - weighted
   Float_t Eff;				// efficiency
   //
   N_s[index] = 0;
   N_w[index] = 0;
   PT_0 = E_0;
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
		if ( ( ((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) && (index==1) ) continue; // least stringent eta cut
		
		if ( (Eta > 1.44) && (index>=2) ) {
			continue; // least stringent eta cut
		}
		
		if ( (PT < PT_0) && (index>=3) ) continue;	// Photon criteria 
	
		if ( (E < E_0) && (index>=4) ) continue;	// Photon criteria 	
		
		Phi = PhotonPhi - METPhi;
		if (Phi > Pi) {
			Phi = Phi - 2*Pi;
		} else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		}
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ( ((MT < 100) || (MT > 130)) && (index>=5) ) continue;	// Transverse momentum criterion
		
		if ( (Iso_lepton2==1) && (index>=6) ) continue; 				// Isolated lepton criterion
		
		if ( ( ((j1_Eta*j2_Eta)>0) || (TMath::Abs(j1_Eta-j2_Eta)>4.0) ) && (index>=7) ) continue;
		
		N_s[index]++;
		N_w[index] = N_w[index] + weight;								// increment element by weight
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
   // Print out data in csv format
   // using eff = no. selected / no. events 
   cout << N_s[index] << ", " << N_w[index]*1000*Lumin << endl;
}

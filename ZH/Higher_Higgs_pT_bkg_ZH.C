#define BkgTree_cxx
#include "BkgTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "Riostream.h"

void BkgTree::Loop(const Int_t index=0)
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
   Int_t nbins = 50; // number of bins
   Float_t bin_width = 10;
   //
   //Define Histogram 
   TH1F* h_bkg = new TH1F("h_bkg", "Higgs p_{T}, signal and Z+#gamma background", nbins, 0, nbins*bin_width);
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
		Eta = TMath::Abs(PhotonEta);
		
		// #### Trigger cuts
		//if (((Eta > 1.4442) && (Eta < 1.566)) || (Eta > 2.5)) continue; // least stringent eta cut
		
		//if (Eta > 1.44) continue; // barrel region eta cut
		
		if (PT < 20) continue;	// Photon PT cut ### 20 value taken from CMS paper
	
		if (E < 60) continue;	// MET cut	### 60 value taken from CMS paper
		
		Phi = PhotonPhi - METPhi;
		if (Phi > Pi) {
			Phi = Phi - 2*Pi;
		} else if (Phi < -Pi) {
			 Phi = Phi + 2*Pi;
		}
		MT = sqrt(2*PT*E*(1-cos(Phi)));	
		
		if ((MT < 100) || (MT > 130)) continue;	// Transverse momentum criterion ### kept over from Gabrielli et al. paper
		
		//
		TVector3 photon_vec;
		photon_vec.SetPtEtaPhi(PT,Eta,PhotonPhi);
		TVector3 MET_vec;
		MET_vec.SetPtEtaPhi(E,0,METPhi);
		//
		TVector3 Higgs_v = photon_vec + MET_vec;
		Float_t HPT = TMath::Abs(Higgs_v.Perp());
		
		// ## further cuts also from CMS paper
		
		l1.SetPtEtaPhi(lep1Pt,lep1Eta,lep1Phi);
		l2.SetPtEtaPhi(lep2Pt,lep2Eta,lep2Phi);
		ZCandidate = l1 + l2;
		if(Zdecay == 0) continue;
		
		if (Jet_PT_max > 30) continue; // no jets with PT greater than 30GeV cut
															
		dPhi_ll_EE = TMath::Abs( (photon_vec + MET_vec).Phi() - ZCandidate.Phi() ); // first Delta Phi cut
		if ( dPhi_ll_EE > Pi) dPhi_ll_EE = 2*Pi - dPhi_ll_EE;
		if ( dPhi_ll_EE < 2.7) continue;
		
		p_variable = TMath::Abs( (photon_vec+MET_vec).Pt() - ZCandidate.Pt() )/ZCandidate.Pt();
		if (p_variable > 0.5) continue;
		
				
		dPhi_ll = TMath::Abs(lep1Phi-lep2Phi);			// azimuthal angle between leptons cut
		if ( dPhi_ll > Pi ) dPhi_ll = 2*Pi - dPhi_ll;
		if ( dPhi_ll > 2.25) continue;
			
		if ( ZCandidate.Pt() < 60) continue; // lepton combined pT cut
			 
		//h_bkg->Fill(HPT, weight*1000*Lumin);						// Fill Histogram
		Int_t min_bin;
		if ((HPT >= 0) && (HPT <= nbins*bin_width)) {
			min_bin = TMath::FloorNint(HPT/bin_width)+1;
		} else if (HPT > nbins*bin_width) {
			min_bin = nbins + 1;
		}
		for(Int_t i=1; i<=nbins+1; i++) {
			if (i <= min_bin) {	
				h_bkg->SetBinContent(i,h_bkg->GetBinContent(i)+1000*Lumin*weight);
			}
		}
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
   h_bkg->GetXaxis()->SetTitle("Higgs p_{t} / GeV");
   h_bkg->GetXaxis()->CenterTitle();
   h_bkg->GetYaxis()->SetTitle("frequency");
   h_bkg->GetYaxis()->CenterTitle();
   //
   /*cout << "Background:" << endl;
   for(Int_t i=1; i<=(nbins+1); i++) {
   	   cout << h_bkg->GetBinContent(i) << endl;
   }*/
   //
   TFile *hfile = new TFile("Higgs_PT_ZH.root", "RECREATE");
   h_bkg->Write();
   hfile->Close();
}
